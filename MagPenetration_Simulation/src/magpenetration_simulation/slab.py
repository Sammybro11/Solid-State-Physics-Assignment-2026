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
    """Return the shared grid for every slab method.

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
    """Return the exact London field B(x) inside a slab.

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




