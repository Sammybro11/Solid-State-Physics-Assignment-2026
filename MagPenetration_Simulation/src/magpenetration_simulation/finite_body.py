"""
Two-dimensional scalar London solver for a rectangular cross-section.

Solves d2B/dx2 + d2B/dy2 = B/lam**2 on the rectangle
-width/2 <= x <= width/2 and -height/2 <= y <= height/2, with B = b0 held
on all four edges. The field Bz(x, y) is one Cartesian component in a field
parallel to z. The grid is uniform, the stencil is the five-point form, and
red-black successive over-relaxation (SOR) does the sweeps.

What the model does not include:

- This is the scalar London equation on the sample only. There is no
  outside field, no demagnetization, and no turning of the screening
  current near an edge. That turning is exactly what the Dirichlet edge
  value b0 leaves out, so treat every edge value as an imposed datum.
- The valid-zone mask, from rect_valid_mask, marks where the 1D slab limit
  is a useful approximation: at least margin * lam from any edge.
- Box and disk geometries are later slices. This module has no box solver
  and no disk solver, and it draws no field-line plot.

Note on the relaxation factor (this departs from docs/q6_geometry_models.md,
section "Model 4"): the interior update divides by 2/hx**2 + 2/hy**2
+ 1/lam**2. The report and the task sheet write the last term with a minus
sign. The minus sign turns the discrete operator into a Helmholtz form whose
solution grows instead of decaying, and the sweep then diverges: for a
20 x 5 rectangle in lam units the centre value runs away to 1e21. The plus
sign is the one that reproduces the 1D slab limit.

All inputs and outputs are in SI units: lengths in m, fields in T,
curvature terms in T/m**2.
"""

from dataclasses import dataclass

import numpy as np

__all__ = [
    "RectSolution",
    "rect_1d_limit",
    "rect_center_field",
    "rect_edge_distance",
    "rect_field_finite_difference",
    "rect_valid_mask",
]


@dataclass
class RectSolution:
    """
    Result of the 2D rectangle solve.

    Attributes
    ----------
    x : numpy.ndarray
        Grid positions along x in m, shape (nx,), from -width/2 to width/2.
    y : numpy.ndarray
        Grid positions along y in m, shape (ny,), from -height/2 to height/2.
    b : numpy.ndarray
        Field in T, shape (ny, nx), indexed b[j, i] at (x[i], y[j]).
    residual : float
        Final maximum absolute finite-difference residual in the interior,
        in T/m**2. This is max |(B_E - 2B + B_W)/hx**2 +
        (B_N - 2B + B_S)/hy**2 - B/lam**2| over the interior nodes.
    sweeps : int
        Number of red-black sweep pairs performed.
    hx : float
        Grid step along x in m.
    hy : float
        Grid step along y in m.
    converged : bool
        True if the residual fell to the tolerance, False if the sweeps
        stopped at max_sweeps.
    """

    x: np.ndarray
    y: np.ndarray
    b: np.ndarray
    residual: float
    sweeps: int
    hx: float
    hy: float
    converged: bool


def _check_rect(width, height, lam, b0):
    """
    Shared input checks for the rectangle model.

    Parameters
    ----------
    width : float
        Rectangle width in m, positive and finite.
    height : float
        Rectangle height in m, positive and finite.
    lam : float
        London penetration depth in m, positive and finite.
    b0 : float
        Field on the four edges in T, non-negative and finite.

    Raises
    ------
    ValueError
        If an input is not finite, width or height <= 0, lam <= 0, or b0 < 0.
    """
    if not np.isfinite(width) or width <= 0:
        raise ValueError(f"width must be positive and finite, got {width}")
    if not np.isfinite(height) or height <= 0:
        raise ValueError(f"height must be positive and finite, got {height}")
    if not np.isfinite(lam) or lam <= 0:
        raise ValueError(f"lam must be positive and finite, got {lam}")
    if not np.isfinite(b0) or b0 < 0:
        raise ValueError(f"b0 must be non-negative and finite, got {b0}")


def rect_edge_distance(x, y, width, height):
    """
    Return the distance from each grid point to the nearest rectangle edge.

    The interior distance to an edge is min(x + width/2, width/2 - x,
    y + height/2, height/2 - y). The result is 0 on an edge and negative
    outside the rectangle.

    Parameters
    ----------
    x : float or numpy.ndarray
        Positions along x in m, measured from the rectangle centre.
    y : float or numpy.ndarray
        Positions along y in m. A meshgrid pair from np.meshgrid(x, y)
        gives one distance per grid point.
    width : float
        Rectangle width in m, positive.
    height : float
        Rectangle height in m, positive.

    Returns
    -------
    numpy.ndarray
        Distance to the nearest edge in m, broadcast over x and y.

    Raises
    ------
    ValueError
        If width <= 0, height <= 0, or an input is not finite.
    """
    _check_rect(width, height, 1.0, 0.0)

    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    return np.minimum(
        np.minimum(x_arr + width / 2, width / 2 - x_arr),
        np.minimum(y_arr + height / 2, height / 2 - y_arr),
    )


def rect_valid_mask(x, y, width, height, lam, margin=3.0):
    """
    Return the mask where the 1D slab limit is a useful approximation.

    The mask keeps the points at least margin * lam from every edge, so the
    edge zone where the screening current turns is excluded. The centre
    survives unless the rectangle is thinner than 2 * margin * lam.

    Parameters
    ----------
    x : float or numpy.ndarray
        Positions along x in m, measured from the rectangle centre.
    y : float or numpy.ndarray
        Positions along y in m, broadcast against x.
    width : float
        Rectangle width in m, positive.
    height : float
        Rectangle height in m, positive.
    lam : float
        London penetration depth in m, positive.
    margin : float, optional
        Edge band width in units of lam. Non-negative.

    Returns
    -------
    numpy.ndarray
        Boolean mask, broadcast over x and y, True where the point is valid.

    Raises
    ------
    ValueError
        If width <= 0, height <= 0, lam <= 0, margin < 0, or an input is
        not finite.
    """
    _check_rect(width, height, lam, 0.0)
    if not np.isfinite(margin) or margin < 0:
        raise ValueError(f"margin must be non-negative and finite, got {margin}")

    return rect_edge_distance(x, y, width, height) >= margin * lam


def rect_field_finite_difference(width, height, lam, b0, nx=81, ny=81,
                                 omega=1.7, tol=1e-9, max_sweeps=20000):
    """
    Solve the scalar London equation on a rectangle by red-black SOR.

    Each interior node is updated from its four neighbours:

        B = (B_E/hx**2 + B_W/hx**2 + B_N/hy**2 + B_S/hy**2)
            / (2/hx**2 + 2/hy**2 + 1/lam**2)

    with the SOR relaxation B <- (1 - omega) * B + omega * B_new. The
    denominator is positive at every grid step, so the discrete operator is
    symmetric and positive definite and the sweep converges for
    0 < omega < 2. The four edges stay at b0. Points are updated in
    red-black order (even index sum, then odd), so each half sweep reads
    only values the other half own.

    tol is relative: the sweeps stop when
    residual <= tol * b0 / lam**2, i.e. when the largest residual is tol of
    the natural curvature scale of the problem. With b0 = 0 the solution is
    zero everywhere and the residual is exactly 0.

    Parameters
    ----------
    width : float
        Rectangle width in m, positive.
    height : float
        Rectangle height in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Field held on the four edges in T, non-negative.
    nx : int, optional
        Number of grid points along x, at least 21. Both edges are
        included, so hx = width / (nx - 1). Odd and even are both allowed.
    ny : int, optional
        Number of grid points along y, at least 21, with
        hy = height / (ny - 1).
    omega : float, optional
        SOR relaxation factor, strictly between 0 and 2. 1 is plain
        Gauss-Seidel.
    tol : float, optional
        Relative residual tolerance, positive.
    max_sweeps : int, optional
        Maximum number of red-black sweep pairs, at least 1.

    Returns
    -------
    RectSolution
        Grid, field, final residual, sweep count, grid steps, and the
        convergence flag.

    Raises
    ------
    ValueError
        If width or height <= 0, lam <= 0, b0 < 0, nx < 21, ny < 21,
        omega outside (0, 2), tol <= 0, max_sweeps < 1, or an input is
        not finite.
    """
    _check_rect(width, height, lam, b0)
    if not isinstance(nx, (int, np.integer)) or nx < 21:
        raise ValueError(f"nx must be an integer >= 21, got {nx}")
    if not isinstance(ny, (int, np.integer)) or ny < 21:
        raise ValueError(f"ny must be an integer >= 21, got {ny}")
    if not np.isfinite(omega) or not 0.0 < omega < 2.0:
        raise ValueError(f"omega must satisfy 0 < omega < 2, got {omega}")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError(f"tol must be positive and finite, got {tol}")
    if not isinstance(max_sweeps, (int, np.integer)) or max_sweeps < 1:
        raise ValueError(f"max_sweeps must be an integer >= 1, got {max_sweeps}")

    x = np.linspace(-width / 2, width / 2, nx)
    y = np.linspace(-height / 2, height / 2, ny)
    hx = width / (nx - 1)
    hy = height / (ny - 1)

    ax = 1.0 / hx**2
    ay = 1.0 / hy**2
    denom = 2.0 * ax + 2.0 * ay + 1.0 / lam**2

    b = np.full((ny, nx), float(b0))
    inner = b[1:-1, 1:-1]
    # Start from the 1D slab profile in y: it already fits the top and bottom
    # edges, so the sweeps do not open with an SOR overshoot.
    ht = height / lam
    yt = y[1:-1] / lam
    inner[:] = b0 * (np.exp(yt - ht / 2) + np.exp(-yt - ht / 2))[:, None] / (1.0 + np.exp(-ht))
    ir, jr = np.indices(inner.shape)
    red = (ir + jr) % 2 == 0

    def residual():
        """Largest interior finite-difference residual in T/m**2."""
        lap = ((b[1:-1, 2:] - 2.0 * inner + b[1:-1, :-2]) * ax
               + (b[2:, 1:-1] - 2.0 * inner + b[:-2, 1:-1]) * ay)
        return float(np.abs(lap - inner / lam**2).max())

    target = tol * b0 / lam**2
    res = residual()

    for sweep in range(1, max_sweeps + 1):
        for mask in (red, ~red):
            num = ((b[1:-1, 2:] + b[1:-1, :-2]) * ax
                   + (b[2:, 1:-1] + b[:-2, 1:-1]) * ay)
            inner[mask] = ((1.0 - omega) * inner[mask]
                           + omega * num[mask] / denom)
        res = residual()
        if res <= target:
            return RectSolution(x, y, b, res, sweep, hx, hy, True)

    return RectSolution(x, y, b, res, max_sweeps, hx, hy, False)


def rect_center_field(solution):
    """
    Return the field at the grid node nearest the rectangle centre.

    With odd nx and ny one node sits exactly at the centre. With even counts
    the nearest of the four middle nodes is used, so the value carries the
    grid error of that node.

    Parameters
    ----------
    solution : RectSolution
        Result of rect_field_finite_difference.

    Returns
    -------
    float
        Field at the centre node in T.
    """
    return float(solution.b[solution.b.shape[0] // 2, solution.b.shape[1] // 2])


def rect_1d_limit(width, height, lam, b0, n_points=401):
    """
    Return the 1D slab profile across y, the height direction.

    The slab form B(y) = b0 * cosh(y/lam) / cosh(height/(2*lam)) is the
    limit of this model at the centre line x = 0 when width >> lam. It is
    the comparison curve for the rectangle centre value, and it ignores the
    width band where the 1D picture fails.

    Parameters
    ----------
    width : float
        Rectangle width in m, positive. Checked only; the slab form does
        not depend on it.
    height : float
        Rectangle height in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Field on the edges in T, non-negative.
    n_points : int, optional
        Number of grid points, odd and at least 3, so one node sits at 0.

    Returns
    -------
    y : numpy.ndarray
        Positions in m, shape (n_points,), from -height/2 to height/2.
    b : numpy.ndarray
        Slab field in T, shape (n_points,).

    Raises
    ------
    ValueError
        If an input is invalid, or n_points is even or less than 3.
    """
    _check_rect(width, height, lam, b0)
    if n_points % 2 == 0 or n_points < 3:
        raise ValueError(f"n_points must be odd and >= 3, got {n_points}")

    y = np.linspace(-height / 2, height / 2, n_points)
    yt = y / lam
    ht = height / lam
    # Multiply top and bottom by exp(-ht/2) so no exponent is positive.
    return y, b0 * (np.exp(yt - ht / 2) + np.exp(-yt - ht / 2)) / (1 + np.exp(-ht))
