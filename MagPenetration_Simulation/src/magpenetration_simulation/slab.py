"""
One-dimensional slab model for London-equation field penetration.

Solves lam**2 * d2B/dx2 = B on -d/2 <= x <= d/2 with B(+/-d/2) = b0,
using four methods that share one grid: analytic, finite difference,
relaxation (SOR), and shooting.

All inputs and outputs are in SI units: lengths in m, fields in T.
"""

import numpy as np
from scipy.linalg import solve_banded
from scipy.integrate import solve_ivp

def _grid(d, n_points):
    """
    Return the shared grid for every slab method.

    Parameters
    ----------
    d : float
        Slab thickness in m.
    n_points : int
        Number of grid points. Must be odd and at least 3, so that
        one node sits exactly at x = 0.

    Returns
    -------
    numpy.ndarray
        Positions x in m, shape (n_points,), evenly spaced from
        -d/2 to +d/2.

    Raises
    ------
    ValueError
        If n_points is even or less than 3.
    """
    if n_points % 2 == 0 or n_points < 3:
        raise ValueError(f"n_points must be odd and >= 3, got {n_points}")
    return np.linspace(-d/2, d/2, n_points)

def slab_field_analytic(x, d, lam, b0):
    """
    Return the exact London field B(x) inside a slab.

    Evaluates B(x) = b0 * cosh(x/lam) / cosh(d/(2*lam)) in a form
    that cannot overflow: the numerator and denominator are both
    multiplied by exp(-d/(2*lam)), so every exponent is <= 0.

    Parameters
    ----------
    x : float or numpy.ndarray
        Position from the slab centre in m, with -d/2 <= x <= d/2.
    d : float
        Slab thickness in m.
    lam : float
        London penetration depth in m.
    b0 : float
        Applied field at both faces in T.

    Returns
    -------
    float or numpy.ndarray
        Field B in T, same shape as x.
    """
    xt = x / lam
    dt = d / lam
    
    return b0 * (np.exp(xt - dt / 2) + np.exp(-xt - dt / 2)) / (1 + np.exp(-dt))

def slab_field_finite_difference(d, lam, b0, n_points):
    """
    Solve the slab field with second-order central differences.

    Discretizes lam**2 * d2B/dx2 = B on the shared grid. Each interior
    node gives -B[i-1] + (2 + h**2) * B[i] - B[i+1] = 0, with
    h = dx/lam. The face values B(+/-d/2) = b0 move to the right-hand
    side, and scipy.linalg.solve_banded solves the tridiagonal system
    in O(n_points) time. The error falls as dx**2.

    Parameters
    ----------
    d : float
        Slab thickness in m.
    lam : float
        London penetration depth in m.
    b0 : float
        Applied field at both faces in T.
    n_points : int
        Number of grid points, odd and at least 3.

    Returns
    -------
    x : numpy.ndarray
        Positions in m, shape (n_points,).
    b : numpy.ndarray
        Field in T, shape (n_points,).
    """
    x = _grid(d, n_points)
    h = (x[1] - x[0]) / lam

    ab = np.zeros((3, n_points - 2))
    ab[0, 1:] = -1.0
    ab[1, :] = 2.0 + h**2
    ab[2, :-1] = -1.0

    rhs = np.zeros(n_points - 2)
    rhs[0] += 1.0
    rhs[-1] += 1.0

    b_inner = solve_banded((1, 1), ab, rhs)
    b = np.concatenate(([1.0], b_inner, [1.0]))
    return x, b0 * b

def slab_field_relaxation(d, lam, b0, n_points, tol=1e-10, max_sweeps=10000):
    """
    Solve the slab field by red-black successive over-relaxation (SOR).

    Solves the same discrete equations as slab_field_finite_difference,
    but by repeated sweeps instead of a matrix solve, so it is an
    independent check on the matrix assembly. Each sweep updates the
    odd interior nodes, then the even ones, with
    B[i] <- (1 - omega) * B[i] + omega * c * (B[i-1] + B[i+1]),
    where c = 1 / (2 + h**2) and h = dx/lam. omega is the optimal SOR
    factor for this tridiagonal system. Plain Jacobi sweeps, as in the
    techniques document, need about n_points**2 / 2 sweeps; SOR needs
    about 2 * n_points.

    Parameters
    ----------
    d : float
        Slab thickness in m.
    lam : float
        London penetration depth in m.
    b0 : float
        Applied field at both faces in T.
    n_points : int
        Number of grid points, odd and at least 3.
    tol : float, optional
        Stop when the largest change in one sweep falls below tol,
        measured in units of b0 (dimensionless).
    max_sweeps : int, optional
        Maximum number of sweeps before giving up.

    Returns
    -------
    x : numpy.ndarray
        Positions in m, shape (n_points,).
    b : numpy.ndarray
        Field in T, shape (n_points,).

    Raises
    ------
    RuntimeError
        If the sweeps do not converge within max_sweeps.
    """
    x = _grid(d, n_points)
    h = (x[1] - x[0]) / lam
    c = 1.0 / (2.0 + h**2)

    rho = 2.0 * c * np.cos(np.pi / (n_points - 1))
    omega = 2.0 / (1.0 + np.sqrt(1.0 - rho**2))

    b = np.ones(n_points)
    red = np.arange(1, n_points - 1, 2) # odd indices
    black = np.arange(2, n_points - 1, 2) # even indices

    for sweep in range(1, max_sweeps + 1):
        old = b.copy()
        for idx in (red, black):
            b[idx] = (1.0 - omega) * b[idx] + omega * c * (b[idx - 1] + b[idx + 1])
        change = np.abs(b - old).max()
        if change < tol:
            return x, b0 * b

    raise RuntimeError(
        f"relaxation did not converge in {max_sweeps} sweeps "
        f"(last change {change:.2e}, tol {tol:.0e})")


def _london_rhs(xt, y):
    """Right-hand side of the dimensionless slab equation.

    With y = [B/b0, d(B/b0)/d(x/lam)], the equation d2y/dxt2 = y
    becomes the first-order pair dy0/dxt = y1, dy1/dxt = y0.
    """
    return [y[1], y[0]]


def slab_field_shooting(d, lam, b0, n_points=401, tol=1e-10):
    """Solve the slab field by linear shooting from the slab centre.

    By symmetry dB/dx = 0 at x = 0, so the only unknown starting value
    is B(0). The equation is linear, so one integration from
    B(0) = 1, dB/dx(0) = 0 out to the face, followed by division by
    the face value, gives the exact correction in a single shot. The
    right half is then mirrored onto the left half of the shared grid.
    Starting from the centre keeps the growing exp(x/lam) part of the
    solution from swamping the answer; the techniques document shot
    from the left face, which fails for d/lam above about 20. This
    version holds to d/lam of about 1400, beyond that the unnormalized
    solution, which grows like cosh(x/lam), overflows before the
    division and solve_ivp fails with a RuntimeError.

    Parameters
    ----------
    d : float
        Slab thickness in m.
    lam : float
        London penetration depth in m.
    b0 : float
        Applied field at both faces in T.
    n_points : int, optional
        Number of grid points, odd and at least 3.
    tol : float, optional
        Relative tolerance passed to scipy.integrate.solve_ivp.

    Returns
    -------
    x : numpy.ndarray
        Positions in m, shape (n_points,).
    b : numpy.ndarray
        Field in T, shape (n_points,).

    Raises
    ------
    RuntimeError
        If the integrator reports a failure.
    """
    x = _grid(d, n_points)
    mid = n_points // 2
    xt_half = x[mid:] / lam

    sol = solve_ivp(_london_rhs, (xt_half[0], xt_half[-1]), [1.0, 0.0],
                    t_eval=xt_half, method="DOP853", rtol=tol, atol=tol * 1e-3) # DOP853 is scipy's high-order RK method, a better fit than the default for a tolerance as tight as 1e-10
    if not sol.success:
        raise RuntimeError(f"shooting integration failed: {sol.message}")

    half = sol.y[0] / sol.y[0, -1]
    b = np.concatenate((half[:0:-1], half))
    return x, b0 * b
