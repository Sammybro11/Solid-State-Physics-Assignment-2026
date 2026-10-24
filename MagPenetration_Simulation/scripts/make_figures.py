"""Smoke check and figures for the 1D slab model.

Runs the four slab methods at the agreed defaults, prints the maximum error of
each numerical method against the analytic answer, writes two figures to
docs/figures/, and exits with status 1 if any error exceeds the 1e-3 target.

Run from MagPenetration_Simulation/ with:
    uv run python scripts/make_figures.py
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # write files only; never open a window
import matplotlib.pyplot as plt
import numpy as np

from magpenetration_simulation.slab import (
    _grid,
    slab_field_analytic,
    slab_field_finite_difference,
    slab_field_relaxation,
    slab_field_shooting,
)

# Agreed defaults (SI units)
LAM = 50e-9        # penetration depth lambda_0, m
D = 5 * LAM        # slab thickness, m
B0 = 50e-3         # applied field, T
N_POINTS = 401     # grid points
TARGET = 1e-3      # maximum allowed relative error

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

    if failed:
        print(f"FAIL: {', '.join(failed)} above the {TARGET:g} target", file=sys.stderr)
        sys.exit(1)
    print(f"All three methods within the {TARGET:g} target.")


if __name__ == "__main__":
    main()
