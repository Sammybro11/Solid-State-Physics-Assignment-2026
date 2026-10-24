"""Slab page: the exact London field through an infinite plate."""

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from magpenetration_simulation import ui, visualization
from magpenetration_simulation import ui_config as cfg
from magpenetration_simulation.slab import slab_field_analytic


def ratio_text(value):
    """A field ratio with three decimals, or scientific notation when tiny."""
    return f"{value:.3f}" if value >= 1e-2 else f"{value:.2e}"


def profile_chart(x_nm, ratio, b_0, units, t):
    """B against depth through the full slab, from -d/2 to d/2."""
    yf = b_0 * 1e3 if units == "B (mT)" else 1.0
    fmt = "%{y:.4f} mT" if units == "B (mT)" else "B/B₀ = %{y:.4f}"
    fig = go.Figure()
    fig.add_vrect(x0=x_nm[0], x1=x_nm[-1], fillcolor=t["slab"], line_width=0, layer="below")
    fig.add_scatter(x=x_nm, y=ratio * yf, mode="lines", name="B",
                    line={"color": t["analytic"], "width": 2.4}, hovertemplate=fmt)
    fig.update_xaxes(title="x (nm)", range=[x_nm[0], x_nm[-1]])
    fig.update_yaxes(title="B (mT)" if units == "B (mT)" else "B / B₀", range=[0, cfg.PLOT["y_headroom"] * yf])
    return ui.style(fig, t, ui.CHART_H, uirevision=f"slab-profile-{units}-{x_nm[0]:.6e}-{x_nm[-1]:.6e}")


def render():
    """Show the slab field map, its profile, its mean, and the Pearl crossover."""
    t = ui.theme()
    ui.page_title("Slab", "Exact London field through an infinite plate: the profile, its mean, and the thin-film Pearl crossover.")

    left, right = st.columns([0.34, 0.66], gap="medium")
    with left:
        material = ui.material_controls("slab", compact=True)
        lam, lambda_0 = material["lambda"], material["lambda_0"]
        t_c, t_k, b_0 = material["t_c"], material["t"], material["b_0"]
        with ui.card("controls"):
            size = ui.size_input("slab", cfg.SLAB["d_over_lam"], cfg.SLAB["d_nm"],
                                 "Slab thickness d/λ", "Slab thickness d (nm)", lam, t)
            n_points = int(ui.live_slider("Grid points N (odd)", cfg.GRID["min"], cfg.GRID["max"],
                                          cfg.GRID["default"], cfg.GRID["step"], "slab_points", t))
            units = st.segmented_control("Units", ui.GEOM_VIEWS, default=ui.GEOM_VIEWS[0],
                                         key="slab_units") or ui.GEOM_VIEWS[0]

    normal = material["normal_state"]
    # λ is infinite in the normal state, so a ratio-sized slab takes its length from λ₀.
    d = size["size_m"] if np.isfinite(size["size_m"]) else size["ratio"] * lambda_0
    d_over_lam = size["ratio"]

    with right:
        with ui.card("card_visual"):
            ui.card_title("Field map")
            ui.live_chart(visualization.slab_section(d, lam, b_0, t), key="slab_visual")
            st.caption("Longitudinal section: the field points along z, the heatmap is |B|/B₀, "
                       "and the field lines are straight.")
        k = st.columns(4, gap="small")
        if normal:
            k[0].metric("λ(T)", "∞", "normal state", delta_color="off", delta_arrow="off")
            k[1].metric("Centre field B(0)/B₀", "1.000", "no screening", delta_color="off", delta_arrow="off")
            k[2].metric("Mean field ⟨B⟩/B₀", "1.000", "volume average", delta_color="off", delta_arrow="off")
            k[3].metric("Screened fraction", "0.000", "1 − B(0)/B₀", delta_color="off", delta_arrow="off")
        else:
            centre = float(slab_field_analytic(0.0, d, lam, b_0)) / b_0
            mean = (2.0 / d_over_lam) * np.tanh(d_over_lam / 2.0)
            k[0].metric("λ(T)", f"{lam * 1e9:.1f} nm", f"λ₀ = {lambda_0 * 1e9:.0f} nm",
                        delta_color="off", delta_arrow="off")
            k[1].metric("Centre field B(0)/B₀", ratio_text(centre), ui.field_text(centre * b_0),
                        delta_color="off", delta_arrow="off",
                        help="Field at the slab midplane, as a fraction of the applied field.")
            k[2].metric("Mean field ⟨B⟩/B₀", ratio_text(mean), "volume average",
                        delta_color="off", delta_arrow="off",
                        help="⟨B⟩/B₀ = (2λ/d) tanh(d/2λ).")
            k[3].metric("Screened fraction", f"{1.0 - centre:.3f}", "1 − B(0)/B₀",
                        delta_color="off", delta_arrow="off")
            if d_over_lam < 2.0:
                pearl_nm = 2.0 * lam**2 / d * 1e9
                st.warning(f"Pearl regime: d/λ = {d_over_lam:g} < 2. Λ = 2λ²/d = {pearl_nm:.3g} nm "
                           "replaces λ for in-plane fields in a perpendicular field.")

    c1, c2 = st.columns(2, gap="medium")
    if normal:
        x_nm = np.linspace(-d * 1e9 / 2.0, d * 1e9 / 2.0, n_points)
        with c1, ui.card("card_profile"):
            ui.card_title("Field profile")
            ui.live_chart(profile_chart(x_nm, np.ones_like(x_nm), b_0, units, t), key="slab_profile")
            st.caption("Normal state: λ → ∞, so the size sets the window. The field is flat at B₀.")
        with c2:
            ui.normal_state_infographic(t)
        return

    with c1, ui.card("card_profile"):
        ui.card_title("Field profile")
        ui.live_chart(ui.temperature_profile_chart(t_c, lambda_0, d, b_0, t_k, units, t),
                      key="slab_profile")
        st.caption("All temperatures meet at B₀ at both faces. Faint: the T/T_c family. "
                   "Bold: the current T.")
    with c2, ui.card("card_lambda"):
        ui.card_title("λ(T)")
        ui.live_chart(ui.lambda_chart(t_c, lambda_0, t_k, t), key="slab_lambda")
