"""
Thin films and the Pearl regime for London-equation field penetration.

The Pearl length Lambda = 2 * lam**2 / d replaces the London depth lam as the
in-plane screening length of a film of thickness d. The crossover is
d = 2 * lam, where Lambda = lam:

- d < 2 * lam ("pearl"): Lambda > lam, so in-plane fields and currents
  spread over Lambda instead of lam. lam alone is not enough there.
- d = 2 * lam ("crossover"): Lambda = lam.
- d > 2 * lam ("thick"): the film is thick and lam sets the scale.

For a parallel field the 1D slab profile stays exact at every thickness, so
film_parallel_profile delegates to slab.slab_field_analytic.

The perpendicular-field film is not solved here. Its demagnetizing factor
approaches 1, the field concentrates at the edge, and the edge currents need
a separate calculation.

All inputs and outputs are in SI units: lengths in m, fields in T.
"""

import numpy as np

from .slab import slab_field_analytic

__all__ = [
    "film_crossover",
    "film_parallel_center_field",
    "film_parallel_mean_field",
    "film_parallel_profile",
    "film_regime",
    "inplane_screening_length",
    "pearl_length",
    "pearl_ratio",
]

# Relative tolerance for d == 2 * lam in film_regime.
_REGIME_RTOL = 1e-12


def _check(lam, d, b0=None):
    """
    Shared input checks for the film model.

    Parameters
    ----------
    lam : float
        London penetration depth in m, positive and finite.
    d : float
        Film thickness in m, positive and finite.
    b0 : float, optional
        Applied field in T, non-negative and finite. Checked only when
        given.

    Raises
    ------
    ValueError
        If an input is not finite, lam <= 0, d <= 0, or b0 < 0.
    """
    if not np.isfinite(lam) or lam <= 0:
        raise ValueError(f"lam must be positive and finite, got {lam}")
    if not np.isfinite(d) or d <= 0:
        raise ValueError(f"d must be positive and finite, got {d}")
    if b0 is not None and (not np.isfinite(b0) or b0 < 0):
        raise ValueError(f"b0 must be non-negative and finite, got {b0}")


def pearl_length(lam, d):
    """
    Return the Pearl length Lambda = 2 * lam**2 / d.

    Parameters
    ----------
    lam : float
        London penetration depth in m, positive.
    d : float
        Film thickness in m, positive.

    Returns
    -------
    float
        Pearl length in m. Equals lam at d = 2 * lam.

    Raises
    ------
    ValueError
        If lam <= 0, d <= 0, or an input is not finite.
    """
    _check(lam, d)
    return float(2.0 * lam**2 / d)


def film_crossover(lam):
    """
    Return the film thickness at which Lambda = lam.

    Parameters
    ----------
    lam : float
        London penetration depth in m, positive.

    Returns
    -------
    float
        Crossover thickness 2 * lam in m.

    Raises
    ------
    ValueError
        If lam <= 0 or lam is not finite.
    """
    if not np.isfinite(lam) or lam <= 0:
        raise ValueError(f"lam must be positive and finite, got {lam}")
    return float(2.0 * lam)


def pearl_ratio(lam, d):
    """
    Return the dimensionless Pearl ratio Lambda/lam = 2 * lam / d.

    Parameters
    ----------
    lam : float
        London penetration depth in m, positive.
    d : float
        Film thickness in m, positive.

    Returns
    -------
    float
        Lambda/lam, dimensionless. Greater than 1 for d < 2 * lam.

    Raises
    ------
    ValueError
        If lam <= 0, d <= 0, or an input is not finite.
    """
    _check(lam, d)
    return float(2.0 * lam / d)


def inplane_screening_length(lam, d):
    """
    Return the length that sets in-plane screening, max(Lambda, lam).

    For d > 2 * lam the Pearl length falls below lam, and lam sets the
    scale, so the two agree with the thick-film limit.

    Parameters
    ----------
    lam : float
        London penetration depth in m, positive.
    d : float
        Film thickness in m, positive.

    Returns
    -------
    float
        Screening length in m.

    Raises
    ------
    ValueError
        If lam <= 0, d <= 0, or an input is not finite.
    """
    _check(lam, d)
    return float(max(pearl_length(lam, d), lam))


def film_regime(lam, d):
    """
    Return the film regime name for a given thickness.

    Parameters
    ----------
    lam : float
        London penetration depth in m, positive.
    d : float
        Film thickness in m, positive.

    Returns
    -------
    str
        "pearl" if d < 2 * lam, "crossover" if d == 2 * lam within a
        relative tolerance of 1e-12, otherwise "thick".

    Raises
    ------
    ValueError
        If lam <= 0, d <= 0, or an input is not finite.
    """
    _check(lam, d)
    dc = film_crossover(lam)
    if abs(d - dc) <= _REGIME_RTOL * dc:
        return "crossover"
    return "pearl" if d < dc else "thick"


def film_parallel_profile(x, d, lam, b0):
    """
    Return the parallel-field profile of a film of any thickness.

    The 1D slab solution is exact at every thickness, so this delegates to
    slab.slab_field_analytic.

    Parameters
    ----------
    x : float or numpy.ndarray
        Position from the film centre in m, with -d/2 <= x <= d/2.
    d : float
        Film thickness in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field at both faces in T, non-negative.

    Returns
    -------
    float or numpy.ndarray
        Field B in T, same shape as x.

    Raises
    ------
    ValueError
        If lam <= 0, d <= 0, b0 < 0, or an input is not finite.
    """
    _check(lam, d, b0)
    return slab_field_analytic(x, d, lam, b0)


def film_parallel_center_field(d, lam, b0):
    """
    Return the field at the film centre, B(0) = b0 / cosh(d/(2*lam)).

    Evaluated as b0 * 2 * m / (1 + m**2) with m = exp(-d/(2*lam)), so large
    d/lam cannot overflow and m stays non-zero.

    Parameters
    ----------
    d : float
        Film thickness in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field at both faces in T, non-negative.

    Returns
    -------
    float
        Centre field in T.

    Raises
    ------
    ValueError
        If lam <= 0, d <= 0, b0 < 0, or an input is not finite.
    """
    _check(lam, d, b0)
    m = np.exp(-d / (2.0 * lam))
    return float(b0 * 2.0 * m / (1.0 + m * m))


def film_parallel_mean_field(d, lam, b0):
    """
    Return the mean parallel field, <B> = b0 * (2*lam/d) * tanh(d/(2*lam)).

    Evaluated with m = exp(-d/(2*lam)) and tanh = (1 - m**2)/(1 + m**2),
    so large d/lam cannot overflow.

    Parameters
    ----------
    d : float
        Film thickness in m, positive.
    lam : float
        London penetration depth in m, positive.
    b0 : float
        Applied field at both faces in T, non-negative.

    Returns
    -------
    float
        Thickness-averaged field in T.

    Raises
    ------
    ValueError
        If lam <= 0, d <= 0, b0 < 0, or an input is not finite.
    """
    _check(lam, d, b0)
    u = d / (2.0 * lam)
    m = np.exp(-u)
    return float(b0 * ((1.0 - m * m) / (1.0 + m * m)) / u)
