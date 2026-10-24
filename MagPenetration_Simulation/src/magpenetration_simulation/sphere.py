"""
Exact local London vector solution for a sphere in a uniform applied field.

A sphere of radius R sits in a uniform field B0 along the polar axis. The
London equation separates with a toroidal potential psi(r) sin(theta) in the
phi direction, and its regular radial solution is the modified spherical
Bessel function i1(r/lam). Continuity of both field components at the surface
fixes the prefactor, so demagnetization is included through that vector
boundary matching rather than through a scalar demagnetizing factor.

All inputs and outputs are in SI units: lengths in m, fields in T, the
potential psi in T m, and the dipole moment in A m^2.

i1 and its derivative are evaluated in the scaled form exp(-x) * i1(x) and
exp(-x) * i1'(x), which cannot overflow for large R/lam.
"""

import numpy as np

__all__ = [
    "sphere_center_field",
    "sphere_dipole_moment",
    "sphere_equator_field",
    "sphere_field",
    "sphere_mean_field",
    "sphere_potential",
    "sphere_screening_factor",
]

# Vacuum permeability in N/A^2, same value as the other geometry modules.
MU0 = 4.0 * np.pi * 1e-7


def _check(radius, lam, b0):
    """
    Shared input checks for the sphere model.

    Parameters
    ----------
    radius : float
        Sphere radius R in m, positive and finite.
    lam : float
        London penetration depth in m, positive and finite.
    b0 : float
        Applied field in T, non-negative and finite.

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


def _i1_scaled(x):
    """
    Return exp(-x) * i1(x) for x >= 0, stable at both small and large x.

    Uses the power series for x < 0.5 and the closed form
    exp(-x) * i1(x) = ((x - 1) + (x + 1) exp(-2x)) / (2 x**2) for x >= 0.5.
    """
    x = np.asarray(x, dtype=float)
    small = x < 0.5

    xs = np.where(small, x, 0.0)
    series = np.exp(-xs) * (xs / 3.0 + xs**3 / 30.0 + xs**5 / 840.0
                            + xs**7 / 45360.0)

    xl = np.where(small, 1.0, x)
    large = ((xl - 1.0) + (xl + 1.0) * np.exp(-2.0 * xl)) / (2.0 * xl**2)

    return np.where(small, series, large)


def _i1p_scaled(x):
    """
    Return exp(-x) * i1'(x) for x >= 0, stable at both small and large x.

    Uses i1'(x) = i0(x) - 2 i1(x)/x with i0(x) = sinh(x)/x, evaluated in the
    scaled form (1 - exp(-2x)) / (2x) - 2 exp(-x) i1(x) / x, and the power
    series i1'(x) = 1/3 + x**2/10 + x**4/168 + ... for x < 0.5.
    """
    x = np.asarray(x, dtype=float)
    small = x < 0.5

    xs = np.where(small, x, 0.0)
    series = np.exp(-xs) * (1.0 / 3.0 + xs**2 / 10.0 + xs**4 / 168.0
                            + xs**6 / 6480.0)

    xl = np.where(small, 1.0, x)
    large = (1.0 - np.exp(-2.0 * xl)) / (2.0 * xl) - 2.0 * _i1_scaled(xl) / xl

    return np.where(small, series, large)


def _x_over_sinh(x):
    """Return x / sinh(x) for x >= 0 without overflowing sinh."""
    if x < 1.0:
        return x / np.sinh(x)
    e = np.exp(-x)
    return 2.0 * x * e / (1.0 - e * e)


def _screening_factor(x):
    """
    Return S(x) = 1 - (3/x) coth(x) + 3/x**2 for x = R/lam > 0.

    Uses the power series S = x**2/15 - 2 x**4/315 + ... for x < 0.2, which
    avoids the catastrophic cancellation of the direct form at small x.
    """
    if x < 0.2:
        x2 = x * x
        return x2 / 15.0 - 2.0 * x2**2 / 315.0 + x2**3 / 1575.0 \
            - 2.0 * x2**4 / 31185.0
    coth = 1.0 / np.tanh(x)
    return 1.0 - 3.0 * coth / x + 3.0 / x**2


def _interior_psi(r_arr, radius, lam, b0):
    """
    Return psi, dpsi/dr, and psi/r inside the sphere.

    psi(r) = 3 B0 R / (2 sinh(R/lam)) * i1(r/lam) is evaluated as
    3 B0 R / (1 - exp(-2 R/lam)) * exp(-(R/lam - r/lam)) * exp(-r/lam) i1,
    so the decaying exponent never overflows. At r = 0 the ratio psi/r takes
    its limit B0 R / (2 lam sinh(R/lam)).
    """
    x_r = radius / lam
    x = r_arr / lam

    factor = 3.0 * b0 * radius / (1.0 - np.exp(-2.0 * x_r))
    decay = np.exp(-(x_r - x))
    psi = factor * decay * _i1_scaled(x)
    psi_prime = factor * decay * _i1p_scaled(x) / lam

    q0 = 0.5 * b0 * _x_over_sinh(x_r)
    q = np.empty_like(psi)
    np.divide(psi, r_arr, out=q, where=r_arr > 0)
    q = np.where(r_arr > 0, q, q0)
    return psi, psi_prime, q


def sphere_potential(r, radius, lam, b0):
    """
    Return the toroidal London potential psi(r) inside the sphere.

    psi(r) = 3 B0 R / (2 sinh(R/lam)) * i1(r/lam), in T m.

    Parameters
    ----------
    r : float or numpy.ndarray
        Radius from the sphere center in m, with 0 <= r <= radius.
    radius : float
        Sphere radius R in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field in T, non-negative.

    Returns
    -------
    float or numpy.ndarray
        Potential psi in T m, same shape as r.

    Raises
    ------
    ValueError
        If radius <= 0, lam <= 0, b0 < 0, an input is not finite, or r is
        outside [0, radius].
    """
    _check(radius, lam, b0)

    r_arr = np.asarray(r, dtype=float)
    if not np.all(np.isfinite(r_arr)):
        raise ValueError("r must be finite")
    if np.any(r_arr < 0) or np.any(r_arr > radius * (1.0 + 1e-12)):
        raise ValueError(f"r must satisfy 0 <= r <= radius, got {r}")

    return float(_interior_psi(r_arr, radius, lam, b0)[0]) \
        if r_arr.ndim == 0 else _interior_psi(r_arr, radius, lam, b0)[0]


def sphere_field(r, theta, radius, lam, b0):
    """
    Return the interior London field (B_r, B_theta) in the sphere.

    B_r = 2 psi(r)/r * cos(theta)
    B_theta = -(psi(r) + r psi'(r))/r * sin(theta)

    Both formulas are finite at r = 0 and are evaluated with the limiting
    value of psi/r there.

    Parameters
    ----------
    r : float or numpy.ndarray
        Radius from the sphere center in m, with 0 <= r <= radius.
    theta : float or numpy.ndarray
        Polar angle from the applied-field axis in rad. Broadcastable with r.
    radius : float
        Sphere radius R in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field in T, non-negative.

    Returns
    -------
    B_r : float or numpy.ndarray
        Radial field component in T.
    B_theta : float or numpy.ndarray
        Polar field component in T.

    Raises
    ------
    ValueError
        If radius <= 0, lam <= 0, b0 < 0, an input is not finite, or r is
        outside [0, radius].
    """
    _check(radius, lam, b0)

    r_arr, th_arr = np.broadcast_arrays(np.asarray(r, dtype=float),
                                        np.asarray(theta, dtype=float))
    if not np.all(np.isfinite(r_arr)) or not np.all(np.isfinite(th_arr)):
        raise ValueError("r and theta must be finite")
    if np.any(r_arr < 0) or np.any(r_arr > radius * (1.0 + 1e-12)):
        raise ValueError(f"r must satisfy 0 <= r <= radius, got {r}")

    _, psi_prime, q = _interior_psi(r_arr, radius, lam, b0)
    b_r = 2.0 * q * np.cos(th_arr)
    b_theta = -(q + psi_prime) * np.sin(th_arr)

    if b_r.ndim == 0:
        return float(b_r), float(b_theta)
    return b_r, b_theta


def sphere_center_field(radius, lam, b0):
    """
    Return the field at the sphere center, B(0) = B0 * (R/lam)/sinh(R/lam).

    The form 2 (R/lam) exp(-R/lam) / (1 - exp(-2 R/lam)) is used for large
    R/lam, so sinh cannot overflow.

    Parameters
    ----------
    radius : float
        Sphere radius R in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field in T, non-negative.

    Returns
    -------
    float
        Field at the center in T.

    Raises
    ------
    ValueError
        If radius <= 0, lam <= 0, b0 < 0, or an input is not finite.
    """
    _check(radius, lam, b0)
    return float(b0 * _x_over_sinh(radius / lam))


def sphere_mean_field(radius, lam, b0):
    """
    Return the axial mean field, <B_z> = B0 (1 - S).

    This equals B0 * ((3 lam/R) coth(R/lam) - 3 lam**2/R**2) and reuses the
    stable screening factor S.

    Parameters
    ----------
    radius : float
        Sphere radius R in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field in T, non-negative.

    Returns
    -------
    float
        Mean field along the polar axis in T.

    Raises
    ------
    ValueError
        If radius <= 0, lam <= 0, b0 < 0, or an input is not finite.
    """
    _check(radius, lam, b0)
    return float(b0 * (1.0 - _screening_factor(radius / lam)))


def sphere_screening_factor(radius, lam, b0):
    """
    Return the dimensionless screening factor S of the exterior dipole.

    S = 1 - (3 lam/R) coth(R/lam) + 3 lam**2/R**2, which tends to R**2/15 at
    small R/lam and to 1 in the ideal limit. It does not depend on b0; the
    argument is kept for a uniform signature with the other functions.

    Parameters
    ----------
    radius : float
        Sphere radius R in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field in T, non-negative. Unused.

    Returns
    -------
    float
        Screening factor S, dimensionless.

    Raises
    ------
    ValueError
        If radius <= 0, lam <= 0, b0 < 0, or an input is not finite.
    """
    _check(radius, lam, b0)
    return float(_screening_factor(radius / lam))


def sphere_dipole_moment(radius, lam, b0):
    """
    Return the exterior dipole moment of the sphere.

    m = -(2 pi R**3 B0 / mu0) * S, in A m^2.

    Parameters
    ----------
    radius : float
        Sphere radius R in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field in T, non-negative.

    Returns
    -------
    float
        Dipole moment in A m^2, negative for b0 > 0 (diamagnetic).

    Raises
    ------
    ValueError
        If radius <= 0, lam <= 0, b0 < 0, or an input is not finite.
    """
    _check(radius, lam, b0)
    s = _screening_factor(radius / lam)
    return float(-(2.0 * np.pi * radius**3 * b0 / MU0) * s)


def sphere_equator_field(radius, lam, b0):
    """
    Return the interior field magnitude at the equator, r = R, theta = pi/2.

    There B_r = 0 and |B_theta| = (3/2) B0 (1 - i1(R/lam)/sinh(R/lam)),
    evaluated in the scaled form 2 exp(-x) i1(x) / (1 - exp(-2x)).

    Parameters
    ----------
    radius : float
        Sphere radius R in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field in T, non-negative.

    Returns
    -------
    float
        Equator field magnitude in T. Approaches 1.5 B0 in the ideal limit.

    Raises
    ------
    ValueError
        If radius <= 0, lam <= 0, b0 < 0, or an input is not finite.
    """
    _check(radius, lam, b0)
    x_r = radius / lam
    ratio = 2.0 * float(_i1_scaled(x_r)) / (1.0 - np.exp(-2.0 * x_r))
    return float(1.5 * b0 * (1.0 - ratio))
