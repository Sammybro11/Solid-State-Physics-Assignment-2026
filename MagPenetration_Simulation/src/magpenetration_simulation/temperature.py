"""
Two-fluid temperature dependence of the London penetration depth.

The superfluid density follows the two-fluid law n_s(T)/n = 1 - (T/Tc)**4,
so the London depth grows as lambda0 / sqrt(1 - (T/Tc)**4). The slab profile
of slab.py then takes that depth as its only temperature input.

All inputs and outputs are in SI units: lengths in m, fields in T,
temperatures in K.
"""

import numpy as np

from .slab import slab_field_analytic

__all__ = [
    "clamp_temperature",
    "guard_metrics",
    "lambda_curve",
    "lambda_of_temperature",
    "temperature_profile",
]


def lambda_of_temperature(t, lambda_0, t_c):
    """
    Return the London penetration depth from the two-fluid law.

    lambda(T) = lambda0 / sqrt(1 - (T/Tc)**4), which has no value at
    T = Tc: the superfluid fraction vanishes and screening fails there.
    Points at or above Tc return inf, per the guard contract.

    Parameters
    ----------
    t : float or numpy.ndarray
        Temperature in K, with t >= 0.
    lambda_0 : float
        Zero-temperature penetration depth in m, positive.
    t_c : float
        Critical temperature in K, positive.

    Returns
    -------
    float or numpy.ndarray
        Penetration depth in m. A Python float for scalar input, an array
        of the same shape for array input. Entries with t >= t_c are inf.

    Raises
    ------
    ValueError
        If t < 0, lambda_0 <= 0, or t_c <= 0.
    """
    if lambda_0 <= 0:
        raise ValueError(f"lambda_0 must be positive, got {lambda_0}")
    if t_c <= 0:
        raise ValueError(f"t_c must be positive, got {t_c}")

    t_arr = np.asarray(t, dtype=float)
    if np.any(t_arr < 0):
        raise ValueError("t must be non-negative")

    scalar = t_arr.ndim == 0
    t_arr = np.atleast_1d(t_arr)

    out = np.full(t_arr.shape, np.inf)
    valid = t_arr < t_c
    # Only the valid branch is evaluated, so t == t_c cannot divide by zero.
    out[valid] = lambda_0 / np.sqrt(1.0 - (t_arr[valid] / t_c) ** 4)

    return float(out[0]) if scalar else out


def lambda_curve(t_values, lambda_0, t_c):
    """
    Return penetration depths for a temperature array, for the lambda(T) plot.

    Parameters
    ----------
    t_values : array_like
        Temperatures in K, all >= 0.
    lambda_0 : float
        Zero-temperature penetration depth in m, positive.
    t_c : float
        Critical temperature in K, positive.

    Returns
    -------
    numpy.ndarray
        Penetration depths in m, same shape as t_values. Entries at or
        above t_c are inf.

    Raises
    ------
    ValueError
        If any t < 0, lambda_0 <= 0, or t_c <= 0.
    """
    return np.asarray(lambda_of_temperature(np.asarray(t_values, dtype=float),
                                            lambda_0, t_c), dtype=float)


def clamp_temperature(t_c, factor=0.995):
    """
    Return the clamp value of the temperature slider, factor * t_c.

    Parameters
    ----------
    t_c : float
        Critical temperature in K, positive.
    factor : float
        Fraction of t_c that the slider may reach, default 0.995.

    Returns
    -------
    float
        Clamped temperature in K, equal to factor * t_c.

    Raises
    ------
    ValueError
        If t_c <= 0.
    """
    if t_c <= 0:
        raise ValueError(f"t_c must be positive, got {t_c}")
    return factor * t_c


def temperature_profile(x, d, lambda_0, t_c, t, b_0):
    """
    Return the slab field B(x) at temperature t.

    Below t_c the field is the exact London profile with lam = lambda(T).
    At or above t_c screening fails, so the field is flat at b0.

    Parameters
    ----------
    x : float or numpy.ndarray
        Position from the slab centre in m, with -d/2 <= x <= d/2.
    d : float
        Slab thickness in m, positive.
    lambda_0 : float
        Zero-temperature penetration depth in m, positive.
    t_c : float
        Critical temperature in K, positive.
    t : float
        Temperature in K, with t >= 0.
    b_0 : float
        Applied field at both faces in T, non-negative.

    Returns
    -------
    numpy.ndarray
        Field B in T, same shape as x.

    Raises
    ------
    ValueError
        If t < 0, d <= 0, or b_0 < 0.
    """
    if t < 0:
        raise ValueError(f"t must be non-negative, got {t}")
    if d <= 0:
        raise ValueError(f"d must be positive, got {d}")
    if b_0 < 0:
        raise ValueError(f"b_0 must be non-negative, got {b_0}")

    x_arr = np.asarray(x, dtype=float)
    lam = lambda_of_temperature(t, lambda_0, t_c)

    if t >= t_c:
        return np.full(x_arr.shape, float(b_0))
    return np.asarray(slab_field_analytic(x_arr, d, lam, b_0), dtype=float)


def guard_metrics(lambda_0, t_c, d):
    """
    Return the warm-end numbers the app reports for the slab.

    Parameters
    ----------
    lambda_0 : float
        Zero-temperature penetration depth in m, positive.
    t_c : float
        Critical temperature in K, positive.
    d : float
        Slab thickness in m, positive.

    Returns
    -------
    dict
        clamp_t : float
            Clamped temperature in K.
        lambda_clamp : float
            Penetration depth at the clamp in m.
        center_ratio : float
            B(0)/B0 at the clamp, dimensionless.

    Raises
    ------
    ValueError
        If lambda_0 <= 0, t_c <= 0, or d <= 0.
    """
    if lambda_0 <= 0:
        raise ValueError(f"lambda_0 must be positive, got {lambda_0}")
    if t_c <= 0:
        raise ValueError(f"t_c must be positive, got {t_c}")
    if d <= 0:
        raise ValueError(f"d must be positive, got {d}")

    clamp_t = clamp_temperature(t_c)
    lambda_clamp = lambda_of_temperature(clamp_t, lambda_0, t_c)
    # The app keeps d / (2 * lambda) below 50 for the slider ranges, so
    # plain cosh stays far from overflow.
    center_ratio = 1.0 / np.cosh(d / (2.0 * lambda_clamp))

    return {
        "clamp_t": float(clamp_t),
        "lambda_clamp": float(lambda_clamp),
        "center_ratio": float(center_ratio),
    }
