"""
Axial infinite cylinder for London-equation field penetration.

Solves lam**2 * (1/r) * d/dr (r * dB/dr) = B on 0 <= r <= R with B(R) = b0.
The bounded solution is B(r) = b0 * I0(r/lam) / I0(R/lam), where I0 is the
modified Bessel function of the first kind of order zero.

The Bessel functions are evaluated in the exponentially scaled form
i0e(x) = exp(-x) * I0(x), so large R/lam cannot overflow.

All inputs and outputs are in SI units: lengths in m, fields in T.
"""

import numpy as np
from scipy.special import i0e, i1e

__all__ = [
    "cylinder_center_field",
    "cylinder_field_analytic",
    "cylinder_field_profile",
    "cylinder_mean_field",
]


def _check(radius, lam, b0):
    """
    Shared input checks for the cylinder model.

    Parameters
    ----------
    radius : float
        Cylinder radius in m, positive and finite.
    lam : float
        London penetration depth in m, positive and finite.
    b0 : float
        Applied field at the surface in T, non-negative and finite.

    Raises
    ------
    ValueError
        If an input is not finite, radius <= 0, lam <= 0, or b0 < 0.
    """
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError(f"radius must be positive and finite, got {radius}")
    if not np.isfinite(lam) or lam <= 0:
        raise ValueError(f"lam must be positive and finite, got {lam}")
    if not np.isfinite(b0) or b0 < 0:
        raise ValueError(f"b0 must be non-negative and finite, got {b0}")


def cylinder_field_analytic(r, radius, lam, b0):
    """
    Return the exact London field B(r) inside an axial infinite cylinder.

    Evaluates B(r) = b0 * I0(r/lam) / I0(R/lam) as
    b0 * i0e(r/lam) / i0e(R/lam) * exp(r/lam - R/lam). The exponent is never
    positive, so nothing overflows, and at r = R the ratio and the exponent
    are exactly 1, so the surface value is exactly b0.

    Parameters
    ----------
    r : float or numpy.ndarray
        Radius from the axis in m, with 0 <= r <= radius.
    radius : float
        Cylinder radius R in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field at the surface in T, non-negative.

    Returns
    -------
    float or numpy.ndarray
        Field B in T, same shape as r.

    Raises
    ------
    ValueError
        If radius <= 0, lam <= 0, b0 < 0, or an input is not finite.
    """
    _check(radius, lam, b0)

    r_arr = np.asarray(r, dtype=float)
    x = r_arr / lam
    x_r = radius / lam

    out = b0 * (i0e(x) / i0e(x_r)) * np.exp(x - x_r)
    return float(out) if out.ndim == 0 else out


def cylinder_field_profile(radius, lam, b0, n_points=401):
    """
    Return the radial field profile of the cylinder on a uniform grid.

    Parameters
    ----------
    radius : float
        Cylinder radius R in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field at the surface in T, non-negative.
    n_points : int, optional
        Number of grid points. Must be odd and at least 3.

    Returns
    -------
    r : numpy.ndarray
        Radii in m, shape (n_points,), from 0 to radius.
    b : numpy.ndarray
        Field in T, shape (n_points,).

    Raises
    ------
    ValueError
        If n_points is even or less than 3, or another input is invalid.
    """
    _check(radius, lam, b0)
    if n_points % 2 == 0 or n_points < 3:
        raise ValueError(f"n_points must be odd and >= 3, got {n_points}")

    r = np.linspace(0.0, radius, n_points)
    return r, cylinder_field_analytic(r, radius, lam, b0)


def cylinder_center_field(radius, lam, b0):
    """
    Return the field on the cylinder axis, B(0) = b0 / I0(R/lam).

    Evaluated as b0 * exp(-R/lam) / i0e(R/lam), so the exponent is never
    positive.

    Parameters
    ----------
    radius : float
        Cylinder radius R in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field at the surface in T, non-negative.

    Returns
    -------
    float
        Field at the axis in T.

    Raises
    ------
    ValueError
        If radius <= 0, lam <= 0, b0 < 0, or an input is not finite.
    """
    _check(radius, lam, b0)

    x_r = radius / lam
    return float(b0 * np.exp(-x_r) / i0e(x_r))


def cylinder_mean_field(radius, lam, b0):
    """
    Return the volume-averaged field, <B> = b0 * (2*lam/R) * I1(R/lam)/I0(R/lam).

    Evaluated as b0 * (2 / (R/lam)) * i1e(R/lam) / i0e(R/lam). The scaled
    ratio i1e/i0e is bounded by 1, so large R/lam cannot overflow.

    Parameters
    ----------
    radius : float
        Cylinder radius R in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field at the surface in T, non-negative.

    Returns
    -------
    float
        Volume-averaged field in T.

    Raises
    ------
    ValueError
        If radius <= 0, lam <= 0, b0 < 0, or an input is not finite.
    """
    _check(radius, lam, b0)

    x_r = radius / lam
    return float(b0 * (2.0 / x_r) * i1e(x_r) / i0e(x_r))
