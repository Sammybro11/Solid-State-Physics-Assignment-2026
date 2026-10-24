"""Smoke check and figures for the 1D slab model.

Runs the four slab methods at the agreed defaults, prints the maximum error of
each numerical method against the analytic answer, writes the slab figures to
docs/figures/, writes the temperature lambda and profile figures, checks the
temperature guard numbers, writes the cylinder geometry figures, checks the
cylinder and slab acceptance values, and exits with status 1 if any error
exceeds the 1e-3 target or a guard leaves its bound.

Run from MagPenetration_Simulation/ with:
    uv run python scripts/make_figures.py
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # write files only; never open a window
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import i0, i1

from magpenetration_simulation.cylinder import (
    cylinder_center_field,
    cylinder_field_analytic,
    cylinder_field_profile,
    cylinder_mean_field,
)
from magpenetration_simulation.slab import (
    _grid,
    slab_field_analytic,
    slab_field_finite_difference,
    slab_field_relaxation,
    slab_field_shooting,
)
from magpenetration_simulation.temperature import (
    clamp_temperature,
    guard_metrics,
    lambda_curve,
    lambda_of_temperature,
    temperature_profile,
)

# Agreed defaults (SI units)
LAM = 50e-9        # penetration depth lambda_0, m
D = 5 * LAM        # slab thickness, m
B0 = 50e-3         # applied field, T
N_POINTS = 401     # grid points
TARGET = 1e-3      # maximum allowed relative error

# Geometry defaults (SI units). T = 0 K, so lambda = LAM exactly.
SIZE = 5                      # slab d/lambda and cylinder R/lambda
FAMILY_R = (1.0, 2.0, 5.0, 10.0, 20.0)   # cylinder R/lambda family

# Acceptance values from the report, as B/B0 ratios
CENTER_5 = 0.0367108923
MEAN_5 = 0.3573532548
CENTER_20 = 2.2957746294e-8
MEAN_20 = 0.0974670508
SLAB_CENTER_5 = 0.1630712319
SLAB_MEAN_5 = 0.3946457193
VALUE_TOL = 1e-9              # decimal places printed in the report
WALL_TOL = 1e-12              # surface field must equal B0
SCALED_TOL = 1e-12            # scaled Bessel forms vs direct ratios
WALL_FORM_TOL = 1e-3          # near-wall approximation
WALL_FORM_S = (0.5, 1.0, 2.0)

# Temperature defaults (SI units), separate from the slab defaults above
TC = 9.2                 # K
D_TEMP = 200e-9          # m
T_TEMP = 4.2             # K
FAMILY_T = (0.0, 0.217, 0.435, 0.652, 0.870, 0.957)
TEMP_POINTS = 401

# Guard bounds from the report
LAMBDA_CLAMP_TARGET_NM = 355.0
LAMBDA_CLAMP_TOL_NM = 3.55      # 1 %
CENTER_TARGET = 0.96
CENTER_TOL = 0.005

METHODS = {
    "Finite difference": slab_field_finite_difference,
    "Relaxation (SOR)": slab_field_relaxation,
    "Shooting": slab_field_shooting,
}
COLORS = {"Finite difference": "#2563eb", "Relaxation (SOR)": "#d97706", "Shooting": "#059669"}
MARKERS = {"Finite difference": "o", "Relaxation (SOR)": "s", "Shooting": "D"}

FIG_DIR = Path(__file__).resolve().parent.parent / "docs" / "figures"


def run_methods():
    """Return the shared grid, the analytic field and each method's field, all in SI units."""
    x = _grid(D, N_POINTS)
    exact = slab_field_analytic(x, D, LAM, B0)
    fields = {}
    for name, method in METHODS.items():
        x_method, b = method(D, LAM, B0, N_POINTS)
        assert np.array_equal(x_method, x), f"{name} did not use the shared grid"
        fields[name] = b
    return x, exact, fields


def style_axes(ax, xt):
    """Shade the slab, mark one lambda in from each face, and tidy the axes."""
    half = D / LAM / 2
    ax.axvspan(-half, half, color="#3b82f6", alpha=0.06, lw=0)
    for s in (-1, 1):
        ax.axvline(s * (half - 1), color="#94a3b8", ls=":", lw=1)
    ax.set_xlim(xt[0], xt[-1])
    ax.set_xlabel(r"$x/\lambda$")
    ax.grid(color="#e5e7eb", lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def profile_figure(x, exact, fields, path):
    """B/B0 against x/lambda on linear axes: analytic line, numerical methods as markers."""
    xt = x / LAM
    fig, ax = plt.subplots(figsize=(7, 4.2), layout="constrained")
    style_axes(ax, xt)
    ax.plot(xt, exact / B0, color="#0f172a", lw=2, label="Analytic", zorder=3)
    step = N_POINTS // 20
    for k, (name, b) in enumerate(fields.items()):
        s = slice(step // 2 + (k - 1) * step // 3, None, step)  # stagger so markers do not stack
        ax.plot(xt[s], b[s] / B0, MARKERS[name], ms=6, mfc="white", mec=COLORS[name], mew=1.6,
                label=name, zorder=4)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel(r"$B/B_0$")
    ax.set_title(rf"Slab field profile, $d = {D / LAM:g}\lambda$, $N = {N_POINTS}$", fontsize=11)
    ax.legend(frameon=False, loc="upper center", ncol=4, fontsize=9)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def error_figure(x, exact, fields, path):
    """Pointwise |B - B_analytic| / B0 on a log axis, with the 1e-3 target."""
    xt = x / LAM
    fig, ax = plt.subplots(figsize=(7, 4.2), layout="constrained")
    style_axes(ax, xt)
    for name, b in fields.items():
        err = np.abs(b - exact) / B0
        ax.semilogy(xt, np.where(err > 0, err, np.nan), color=COLORS[name], lw=1.6,
                    ls="--" if name == "Relaxation (SOR)" else "-", label=name)
    ax.axhline(TARGET, color="#ef4444", ls="--", lw=1, label=r"Target $10^{-3}$")
    ax.set_ylim(1e-16, 1)
    ax.set_ylabel(r"$|B - B_\mathrm{analytic}|\,/\,B_0$")
    ax.set_title(rf"Pointwise error vs analytic, $d = {D / LAM:g}\lambda$, $N = {N_POINTS}$", fontsize=11)
    ax.legend(frameon=False, loc="lower center", ncol=4, fontsize=9)
    fig.savefig(path, dpi=200)
    plt.close(fig)
def temperature_lambda_figure(path):
    """lambda(T)/lambda0 against T/Tc, up to the slider clamp, with the current T marked."""
    clamp = clamp_temperature(TC) / TC
    t_ratio = np.linspace(0.0, clamp, TEMP_POINTS)
    lam = lambda_curve(t_ratio * TC, LAM, TC) / LAM
    fig, ax = plt.subplots(figsize=(7, 4.2), layout="constrained")
    ax.plot(t_ratio, lam, color="#2563eb", lw=2, label=r"$\lambda(T)/\lambda_0$")
    ax.plot([T_TEMP / TC], [lambda_of_temperature(T_TEMP, LAM, TC) / LAM], "o",
            color="#b91c1c", ms=7, label=rf"$T = {T_TEMP:g}$ K", zorder=3)
    ax.axvline(clamp, color="#64748b", ls=":", lw=1.2, label=rf"clamp $T/T_c = {clamp:g}$")
    ax.set_xlim(0, clamp * 1.02)
    ax.set_xlabel(r"$T / T_c$")
    ax.set_ylabel(r"$\lambda(T) / \lambda_0$")
    ax.grid(color="#e5e7eb", lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_title(rf"Two-fluid penetration depth, $T_c = {TC:g}$ K, "
                 rf"$\lambda_0 = {LAM * 1e9:g}$ nm", fontsize=11)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def temperature_profiles_figure(path):
    """B/B0 across the slab for a family of T/Tc values, plus the current T as a bold curve."""
    x_nm = np.linspace(-D_TEMP / 2, D_TEMP / 2, TEMP_POINTS) * 1e9
    cmap = matplotlib.colormaps["viridis"]
    fig, ax = plt.subplots(figsize=(7, 4.2), layout="constrained")
    for k, f in enumerate(FAMILY_T):
        b = temperature_profile(x_nm * 1e-9, D_TEMP, LAM, TC, f * TC, B0) / B0
        ax.plot(x_nm, b, color=cmap(k / (len(FAMILY_T) - 1)), lw=1.6,
                label=rf"$T/T_c = {f:g}$")
    b_now = temperature_profile(x_nm * 1e-9, D_TEMP, LAM, TC, T_TEMP, B0) / B0
    ax.plot(x_nm, b_now, color="#0f172a", lw=3, label=rf"$T = {T_TEMP:g}$ K", zorder=4)
    ax.set_xlim(x_nm[0], x_nm[-1])
    ax.set_xlabel("x (nm)")
    ax.set_ylabel(r"$B / B_0$")
    ax.grid(color="#e5e7eb", lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_title(rf"Slab field vs temperature, $d = {D_TEMP * 1e9:g}$ nm, "
                 rf"$\lambda_0 = {LAM * 1e9:g}$ nm", fontsize=11)
    ax.legend(frameon=False, loc="upper center", ncol=4, fontsize=8)
    fig.savefig(path, dpi=200)
    plt.close(fig)

def tidy(ax):
    """Grid below the data, top and right spines off."""
    ax.grid(color="#e5e7eb", lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def cylinder_family_figure(path):
    """B/B0 against r/R for the R/lambda family, at T = 0 K so lambda = LAM."""
    cmap = matplotlib.colormaps["viridis"]
    fig, ax = plt.subplots(figsize=(7, 4.2), layout="constrained")
    for k, f in enumerate(FAMILY_R):
        r, b = cylinder_field_profile(f * LAM, LAM, B0, N_POINTS)
        ax.plot(r / (f * LAM), b / B0, color=cmap(k / (len(FAMILY_R) - 1)), lw=1.8,
                label=rf"$R/\lambda = {f:g}$")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel(r"$r/R$")
    ax.set_ylabel(r"$B/B_0$")
    tidy(ax)
    ax.set_title(rf"Cylinder field family, $\lambda = {LAM * 1e9:g}$ nm", fontsize=11)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def geometry_comparison_figure(path):
    """Slab and cylinder B/B0 against depth from the nearest face, d = R = SIZE * lambda."""
    d = SIZE * LAM
    s = np.linspace(0.0, SIZE / 2, N_POINTS)
    slab = slab_field_analytic(d / 2 - s * LAM, d, LAM, B0) / B0
    r, b = cylinder_field_profile(d, LAM, B0, N_POINTS)
    fig, ax = plt.subplots(figsize=(7, 4.2), layout="constrained")
    ax.plot(s, slab, color="#2563eb", lw=2, label=r"Slab, $d = 5\lambda$")
    ax.plot(SIZE - r / LAM, b / B0, color="#059669", lw=2, ls="--", label=r"Cylinder, $R = 5\lambda$")
    ax.set_xlim(0, SIZE)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel(r"$s/\lambda$ from the nearest face")
    ax.set_ylabel(r"$B/B_0$")
    tidy(ax)
    ax.set_title(rf"Slab vs cylinder, $\lambda = {LAM * 1e9:g}$ nm, size $= {SIZE:g}\lambda$",
                 fontsize=11)
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def geometry_checks():
    """Return the labels of every geometry check that leaves its bound."""
    x_5, x_20, x_50 = 5.0, 20.0, 50.0
    failed = []

    def compare(label, value, target, tol):
        ok = bool(abs(value - target) <= tol)
        print(f"{label:<40}{value:>22.10g}{target:>22.10g}   {'pass' if ok else 'FAIL'}")
        if not ok:
            failed.append(label)

    wall = max(abs(cylinder_field_analytic(f * LAM, f * LAM, LAM, B0) / B0 - 1.0) for f in FAMILY_R)
    print(f"Geometry defaults: lambda = {LAM * 1e9:g} nm, B0 = {B0 * 1e3:g} mT, "
          f"size = {SIZE:g} lambda, N = {N_POINTS}")
    print(f"{'Check':<40}{'value':>22}{'expected':>22}   result")
    compare("wall B(R)/B0 - 1 (max, R/lambda family)", wall, 0.0, WALL_TOL)
    compare(f"cylinder centre, R = {x_5:g} lambda", cylinder_center_field(x_5 * LAM, LAM, B0) / B0,
            CENTER_5, VALUE_TOL)
    compare(f"cylinder mean, R = {x_5:g} lambda", cylinder_mean_field(x_5 * LAM, LAM, B0) / B0,
            MEAN_5, VALUE_TOL)
    compare(f"cylinder centre, R = {x_20:g} lambda", cylinder_center_field(x_20 * LAM, LAM, B0) / B0,
            CENTER_20, VALUE_TOL)
    compare(f"cylinder mean, R = {x_20:g} lambda", cylinder_mean_field(x_20 * LAM, LAM, B0) / B0,
            MEAN_20, VALUE_TOL)

    d_slab = x_5 * LAM
    compare(f"slab centre, d = {x_5:g} lambda", slab_field_analytic(0.0, d_slab, LAM, B0) / B0,
            SLAB_CENTER_5, VALUE_TOL)
    compare(f"slab mean, d = {x_5:g} lambda",
            2.0 / x_5 * np.tanh(x_5 / 2.0), SLAB_MEAN_5, VALUE_TOL)

    radius = x_20 * LAM
    for s in WALL_FORM_S:
        bessel = cylinder_field_analytic(radius - s * LAM, radius, LAM, B0) / B0
        approx = np.sqrt(x_20 / (x_20 - s)) * np.exp(-s)
        compare(f"near-wall form, R = {x_20:g} lambda, s/lambda = {s:g}", approx, bessel,
                WALL_FORM_TOL)

    # The scaled forms must reproduce the direct Bessel ratios where those still fit.
    x = x_50
    bessel_mean = (2.0 / x) * i1(x) / i0(x)
    compare(f"scaled mean vs direct ratios, R/lambda = {x:g}",
            cylinder_mean_field(x * LAM, LAM, B0) / B0, bessel_mean, SCALED_TOL)
    compare(f"scaled field vs direct ratios, R/lambda = {x:g}",
            cylinder_field_analytic(0.5 * x * LAM, x * LAM, LAM, B0) / B0,
            i0(0.5 * x) / i0(x), SCALED_TOL)

    return failed


def main():
    x, exact, fields = run_methods()

    print(f"Slab defaults: lambda = {LAM * 1e9:g} nm, d = {D / LAM:g} lambda, "
          f"B0 = {B0 * 1e3:g} mT, N = {N_POINTS}")
    print(f"{'Method':<20}{'max |dB| (T)':>15}{'max |dB| / B0':>16}   result")
    failed = []
    for name, b in fields.items():
        abs_err = np.abs(b - exact).max()
        rel_err = abs_err / B0
        ok = rel_err < TARGET
        if not ok:
            failed.append(name)
        print(f"{name:<20}{abs_err:>15.3e}{rel_err:>16.3e}   {'pass' if ok else 'FAIL'}")

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    profile_figure(x, exact, fields, FIG_DIR / "slab_profile.png")
    error_figure(x, exact, fields, FIG_DIR / "slab_error.png")
    print(f"Wrote {FIG_DIR / 'slab_profile.png'}")
    print(f"Wrote {FIG_DIR / 'slab_error.png'}")

    print()
    print(f"Temperature defaults: lambda0 = {LAM * 1e9:g} nm, Tc = {TC:g} K, "
          f"d = {D_TEMP * 1e9:g} nm, B0 = {B0 * 1e3:g} mT, T = {T_TEMP:g} K")
    g = guard_metrics(LAM, TC, D_TEMP)
    lambda_clamp_nm = g["lambda_clamp"] * 1e9
    center_ratio = g["center_ratio"]
    print(f"lambda at clamp: {lambda_clamp_nm:.1f} nm "
          f"(expected {LAMBDA_CLAMP_TARGET_NM:g} +/- {LAMBDA_CLAMP_TOL_NM:g} nm)")
    print(f"B(0)/B0 at clamp: {center_ratio:.4f} "
          f"(expected {CENTER_TARGET:g} +/- {CENTER_TOL:g})")

    guard_failed = []
    if abs(lambda_clamp_nm - LAMBDA_CLAMP_TARGET_NM) > LAMBDA_CLAMP_TOL_NM:
        guard_failed.append("lambda at clamp")
    if abs(center_ratio - CENTER_TARGET) > CENTER_TOL:
        guard_failed.append("B(0)/B0 at clamp")
    if abs(lambda_of_temperature(0.0, LAM, TC) - LAM) > 1e-15:
        guard_failed.append("lambda(T = 0) != lambda_0")
    x_temp = np.linspace(-D_TEMP / 2, D_TEMP / 2, TEMP_POINTS)
    flat = temperature_profile(x_temp, D_TEMP, LAM, TC, TC, B0)
    if np.abs(flat - B0).max() > 1e-15:
        guard_failed.append("field at T = Tc not flat at B0")

    temperature_lambda_figure(FIG_DIR / "temperature_lambda.png")
    temperature_profiles_figure(FIG_DIR / "temperature_profiles.png")
    print(f"Wrote {FIG_DIR / 'temperature_lambda.png'}")
    print(f"Wrote {FIG_DIR / 'temperature_profiles.png'}")

    print()
    geom_failed = geometry_checks()
    cylinder_family_figure(FIG_DIR / "geometry_cylinder.png")
    geometry_comparison_figure(FIG_DIR / "geometry_comparison.png")
    print(f"Wrote {FIG_DIR / 'geometry_cylinder.png'}")
    print(f"Wrote {FIG_DIR / 'geometry_comparison.png'}")

    if failed:
        print(f"FAIL: {', '.join(failed)} above the {TARGET:g} target", file=sys.stderr)
    if guard_failed:
        print(f"FAIL: {', '.join(guard_failed)} outside bound", file=sys.stderr)
    if geom_failed:
        print(f"FAIL: {', '.join(geom_failed)} outside bound", file=sys.stderr)
    if failed or guard_failed or geom_failed:
        sys.exit(1)
    print(f"All three methods within the {TARGET:g} target.")
    print("Temperature guards within bound.")
    print("Geometry checks within bound.")

if __name__ == "__main__":
    main()
