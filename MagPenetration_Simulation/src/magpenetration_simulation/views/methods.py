"""Slab numerical methods page: how the three solvers compare with the analytic field."""

import math

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from magpenetration_simulation import ui
from magpenetration_simulation.slab import (
    _grid,
    slab_field_analytic,
    slab_field_finite_difference,
    slab_field_relaxation,
    slab_field_shooting,
)

METHODS = {
    "Finite difference": slab_field_finite_difference,
    "Relaxation (SOR)": slab_field_relaxation,
    "Shooting": slab_field_shooting,
}
SHORT = {"Finite difference": "FD", "Relaxation (SOR)": "SOR", "Shooting": "Shooting"}
LABEL = {"Analytic": "Analytic", **SHORT}
TARGET = 1e-3
MU0 = 1.25663706212e-6  # vacuum permeability in H/m


@st.cache_data
def solve(d_over_lam, b_0, lam, n_points):
    """Run all four methods on the shared grid.

    Returns positions in units of lam and fields and errors in units of b0,
    so the charts do not depend on the absolute scales.
    """
    d = d_over_lam * lam
    x = _grid(d, n_points)
    exact = slab_field_analytic(x, d, lam, b_0)
    res = {}
    for name, method in METHODS.items():
        try:
            _, b = method(d, lam, b_0, n_points)
            res[name] = {"b": b / b_0, "signed": (b - exact) / b_0,
                         "rel": np.abs(b - exact).max() / b_0, "abs": np.abs(b - exact).max()}
        except RuntimeError as err:
            res[name] = {"error": str(err)}
    return {"xt": x / lam, "dt": d_over_lam, "n": n_points, "exact": exact / b_0, "m": res,
            "lam": lam, "b0": b_0, "d": d}


def x_axis(r, units):
    """Positions and axis title: x/lam when scaled, x in nm when physical."""
    factor = 1.0 if units == "Scaled" else r["lam"] * 1e9
    return r["xt"] * factor, factor, ("x / λ" if units == "Scaled" else "x (nm)")


def slab_backdrop(fig, r, t, xf):
    """Shade the slab and mark one penetration depth in from each face.

    xf converts lengths in units of lam to the chart's x units.
    """
    half = r["dt"] / 2
    fig.add_vrect(x0=-half * xf, x1=half * xf, fillcolor=t["slab"], line_width=0, layer="below")
    if 2 <= r["dt"] <= 40:
        for s in (-1, 1):
            fig.add_vline(x=s * (half - 1) * xf, line={"color": t["axis"], "dash": "dot", "width": 1})
        fig.add_annotation(x=-(half - 1) * xf, y=1, yref="y domain", text="λ", showarrow=False,
                           xshift=-7, yshift=-2, font={"size": 12, "color": t["muted"]})


def profile_chart(r, shown, units, t):
    """Field profile, as B/b0 against x/lam (scaled) or B in mT against x in nm (physical).

    Markers are staggered so methods do not stack; an invisible full-resolution
    line per method drives the hover readout.
    """
    x, xf, x_title = x_axis(r, units)
    yf = 1.0 if units == "Scaled" else r["b0"] * 1e3
    fmt = "%{y:.6f}" if units == "Scaled" else "%{y:.4f} mT"
    half = r["dt"] / 2
    pad = max(0.6, 0.08 * r["dt"])
    fig = go.Figure()
    slab_backdrop(fig, r, t, xf)
    for xs in ([-(half + pad) * xf, -half * xf], [half * xf, (half + pad) * xf]):
        fig.add_scatter(x=xs, y=[yf, yf], mode="lines", hoverinfo="skip",
                        line={"color": t["analytic"], "width": 1.5, "dash": "dot"})
    if "Analytic" in shown:
        fig.add_scatter(x=x, y=r["exact"] * yf, name="Analytic", mode="lines",
                        line={"color": t["analytic"], "width": 2.2}, hovertemplate=fmt)
    methods = [name for name in shown if name in METHODS]
    step = max(4, r["n"] // 20)
    for k, name in enumerate(methods):
        m = r["m"][name]
        if "error" in m:
            continue
        offset = (k - (len(methods) - 1) / 2) * step / max(len(methods), 1)
        s = slice(int(step / 2 + offset), None, step)
        fig.add_scatter(x=x[s], y=m["b"][s] * yf, name=SHORT[name], mode="markers",
                        marker={"size": 8, "color": t["card"], "line": {"color": t[name], "width": 2}},
                        hoverinfo="skip")
        fig.add_scatter(x=x, y=m["b"] * yf, name=SHORT[name], mode="lines",
                        line={"color": t[name], "width": 0}, hovertemplate=fmt)
    fig.update_xaxes(title=x_title, range=[-(half + pad) * xf, (half + pad) * xf])
    fig.update_yaxes(title="B / B₀" if units == "Scaled" else "B (mT)", range=[0, 1.05 * yf])
    return ui.style(fig, t, ui.CHART_H, uirevision=f"profile-{units}")


def error_chart(r, shown, scale, units, t):
    """Pointwise error of each chosen method, signed (linear) or absolute (log),
    in units of b0 (scaled) or in µT (physical)."""
    x, xf, x_title = x_axis(r, units)
    yf = 1.0 if units == "Scaled" else r["b0"] * 1e6
    fig = go.Figure()
    slab_backdrop(fig, r, t, xf)
    for name in shown:
        m = r["m"][name]
        if "error" in m:
            continue
        y = m["signed"] * yf if scale == "Signed" else np.abs(m["signed"]) * yf
        if scale == "Log":
            y = np.where(y > 0, y, np.nan)
        fill = "tozeroy" if scale == "Signed" and len(shown) == 1 else None
        fig.add_scatter(x=x, y=y, name=SHORT[name], mode="lines", fill=fill,
                        fillcolor=ui.alpha(t[name], 0.12), line={"color": t[name], "width": 2},
                        hovertemplate="%{y:.2e}" + ("" if units == "Scaled" else " µT"))
    if scale == "Log":
        top = np.log10(yf)
        fig.add_hline(y=TARGET * yf, line={"color": "#ef4444", "dash": "dash", "width": 1})
        fig.update_yaxes(type="log", range=[top - 16, top], dtick=2, exponentformat="power",
                         title="|B − B_analytic| / B₀" if units == "Scaled" else "|B − B_analytic| (µT)")
    else:
        fig.add_hline(y=0, line={"color": t["axis"], "width": 1})
        fig.update_yaxes(type="linear", autorange=True, range=None, dtick=None,
                         exponentformat="power",
                         title="(B − B_analytic) / B₀" if units == "Scaled" else "B − B_analytic (µT)")
    fig.update_xaxes(title=x_title)
    return ui.style(fig, t, ui.CHART_H, uirevision=f"error-{units}-{scale}")


def render(material):
    """Compare the slab solvers at the shared lambda(T) and B_0."""
    t = ui.theme()
    lam, b_0 = material["lambda"], material["b_0"]
    ui.page_title("Methods", "Three slab solvers on one shared grid.")
    st.markdown(
        "**Finite difference** solves the tridiagonal London system in one pass. "
        "**Red-black SOR** relaxes the same grid in two interleaved colours. "
        "**Shooting** integrates the ODE from one face and matches B = B₀ at the other. "
        "The analytic cosh profile, B = B₀ cosh(x/λ) / cosh(d/2λ), is the reference for "
        "every comparison. The target is a maximum relative error below 1e-3, drawn as the "
        "dashed line on the log scale."
    )

    with ui.card("controls"):
        c = st.columns(2, gap="large")
        with c[0]:
            d_over_lam = ui.live_slider("Slab thickness d/λ", 0.5, 100.0, 5.0, 0.5, "d_over_lam", t, decimals=1)
        with c[1]:
            n_points = int(ui.live_slider("Grid points N (odd)", 41, 2001, 401, 20, "n_points", t))

    if not math.isfinite(lam):
        st.warning("Normal state: λ is infinite, so there is no screening profile to solve.")
        return

    r = solve(d_over_lam, b_0, lam, n_points)

    centre = r["exact"][n_points // 2]
    per_lam = (n_points - 1) / d_over_lam
    dx_nm = r["d"] / (n_points - 1) * 1e9
    surface_current = r["b0"] / (MU0 * r["lam"]) * np.tanh(d_over_lam / 2)

    k = st.columns(7)
    k[0].metric("Centre field B(0)", ui.field_text(centre * r["b0"]),
                f"{centre:.3f} B₀" if centre > 1e-3 else f"{centre:.1e} B₀",
                delta_color="off", delta_arrow="off")
    k[1].metric("Slab thickness d", f"{r['d'] * 1e9:,.0f} nm", f"λ = {r['lam'] * 1e9:.0f} nm",
                delta_color="off", delta_arrow="off")
    k[2].metric("Surface current J", f"{ui.sci(surface_current)} A/m²", "at each face",
                delta_color="off", delta_arrow="off",
                help="Screening current density at each face: J = (1/μ₀) dB/dx = B₀ tanh(d/2λ) / (μ₀λ).")
    k[3].metric("Resolution", f"{per_lam:.0f} pts/λ" if per_lam >= 10 else f"{per_lam:.1f} pts/λ",
                f"Δx = {dx_nm:.3g} nm", delta_color="off", delta_arrow="off")
    for col, name in zip(k[4:], METHODS):
        m = r["m"][name]
        if "error" in m:
            col.metric(f"{SHORT[name]} max error", "failed", help=m["error"])
        else:
            col.metric(f"{SHORT[name]} max error", ui.field_text(m["abs"]), f"{m['rel']:.1e} B₀",
                       delta_color="off", delta_arrow="off",
                       help="Maximum |B − B_analytic| over the grid, absolute and as a fraction of B₀.")

    UNITS = ["Scaled", "Physical"]
    left, right = st.columns(2, gap="medium")
    with left, ui.card("card_profile"):
        a, b, c1 = st.columns([0.9, 2, 1.1], vertical_alignment="center")
        ui.card_title("Field profile", a)
        shown_p = b.pills("Shown on profile", ["Analytic", *METHODS], format_func=LABEL.__getitem__,
                          selection_mode="multi", default=["Analytic", "Finite difference"],
                          label_visibility="collapsed", key="pills_profile")
        units_p = c1.segmented_control("Profile units", UNITS, default="Scaled",
                                       label_visibility="collapsed", key="units_profile") or "Scaled"
        ui.live_chart(profile_chart(r, shown_p, units_p, t), key="profile_chart")
    with right, ui.card("card_error"):
        a, b, c2, c3 = st.columns([0.5, 1.45, 0.85, 1.1], vertical_alignment="center")
        ui.card_title("Error", a)
        shown_e = b.pills("Shown on error", list(METHODS), format_func=LABEL.__getitem__, selection_mode="multi",
                          default=["Finite difference"], label_visibility="collapsed", key="pills_error")
        scale = c2.segmented_control("Scale", ["Signed", "Log"], default="Signed",
                                     label_visibility="collapsed", key="scale") or "Signed"
        units_e = c3.segmented_control("Error units", UNITS, default="Scaled",
                                       label_visibility="collapsed", key="units_error") or "Scaled"
        ui.live_chart(error_chart(r, shown_e, scale, units_e, t), key="error_chart")

    bad = [SHORT[n] for n, m in r["m"].items() if m.get("rel", 1) >= TARGET]
    if bad:
        st.warning(f"{', '.join(bad)} above the 10⁻³ target. Increase N or reduce d/λ.")
