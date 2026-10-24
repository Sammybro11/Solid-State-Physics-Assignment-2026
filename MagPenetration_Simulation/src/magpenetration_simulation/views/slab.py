"""Slab page: the exact London field through an infinite plate."""

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from magpenetration_simulation import ui
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
    fig.update_yaxes(title="B (mT)" if units == "B (mT)" else "B / B₀", range=[0, 1.05 * yf])
    return ui.style(fig, t, ui.CHART_H, uirevision=f"slab-profile-{units}")


def normal_view(t, b_0, lambda_0, d_over_lam, n_points, units):
    """Flat field at B0 and the why, when screening has vanished."""
    k = st.columns(4)
    k[0].metric("λ(T)", "∞", "normal state", delta_color="off", delta_arrow="off")
    k[1].metric("Centre field B(0)/B₀", "1.000", "no screening", delta_color="off", delta_arrow="off")
    k[2].metric("Mean field ⟨B⟩/B₀", "1.000", "volume average", delta_color="off", delta_arrow="off")
    k[3].metric("Screened fraction", "0.000", "1 − B(0)/B₀", delta_color="off", delta_arrow="off")

    # λ is infinite, so the window is sized in units of λ0 instead of λ.
    span_nm = d_over_lam * lambda_0 * 1e9
    x_nm = np.linspace(-span_nm / 2.0, span_nm / 2.0, n_points)
    with ui.card("card_profile"):
        ui.card_title("Field profile")
        ui.live_chart(profile_chart(x_nm, np.ones_like(x_nm), b_0, units, t), key="slab_profile")
    st.caption("Normal state: λ → ∞, so the d/λ slider sets the window in units of λ₀. The field is flat at B₀.")
    ui.normal_state_infographic(t)


def render(material):
    """Show the slab field, its mean, and the Pearl crossover."""
    t = ui.theme()
    lam, lambda_0 = material["lambda"], material["lambda_0"]
    t_c, t_k, b_0 = material["t_c"], material["t"], material["b_0"]

    ui.page_title("Slab", "Exact London field through an infinite plate: the profile, its mean, and the thin-film Pearl crossover.")

    with ui.card("controls"):
        c = st.columns([1, 1, 0.7], gap="large", vertical_alignment="bottom")
        with c[0]:
            d_over_lam = ui.live_slider("Slab thickness d/λ", 0.5, 100.0, 5.0, 0.5, "d_over_lam", t, decimals=1)
        with c[1]:
            n_points = int(ui.live_slider("Grid points N (odd)", 41, 2001, 401, 20, "n_points", t))
        units = c[2].segmented_control("Units", ui.GEOM_VIEWS, default=ui.GEOM_VIEWS[0],
                                       key="slab_units") or ui.GEOM_VIEWS[0]

    if material["normal_state"]:
        normal_view(t, b_0, lambda_0, d_over_lam, n_points, units)
        return

    d = d_over_lam * lam
    x = np.linspace(-d / 2.0, d / 2.0, n_points)
    ratio = slab_field_analytic(x, d, lam, b_0) / b_0
    centre = float(slab_field_analytic(0.0, d, lam, b_0)) / b_0
    mean = (2.0 / d_over_lam) * np.tanh(d_over_lam / 2.0)

    left, right = st.columns([0.9, 1.1], gap="medium")
    with left:
        st.metric("λ(T)", f"{lam * 1e9:.1f} nm", f"λ₀ = {lambda_0 * 1e9:.0f} nm",
                  delta_color="off", delta_arrow="off")
        st.metric("Centre field B(0)/B₀", ratio_text(centre), ui.field_text(centre * b_0),
                  delta_color="off", delta_arrow="off",
                  help="Field at the slab midplane, as a fraction of the applied field.")
        st.metric("Mean field ⟨B⟩/B₀", ratio_text(mean), "volume average",
                  delta_color="off", delta_arrow="off",
                  help="⟨B⟩/B₀ = (2λ/d) tanh(d/2λ).")
        st.metric("Screened fraction", f"{1.0 - centre:.3f}", "1 − B(0)/B₀",
                  delta_color="off", delta_arrow="off")
        if d_over_lam < 2.0:
            pearl_nm = 2.0 * lam**2 / d * 1e9
            st.warning(f"Pearl regime: d/λ = {d_over_lam:g} < 2. Λ = 2λ²/d = {pearl_nm:.3g} nm "
                       "replaces λ for in-plane fields in a perpendicular field.")
    with right, ui.card("card_lambda"):
        ui.card_title("λ(T)")
        ui.live_chart(ui.lambda_chart(t_c, lambda_0, t_k, t), key="slab_lambda")

    profile_col, family_col = st.columns([1.35, 1], gap="medium")
    with profile_col, ui.card("card_profile"):
        ui.card_title("Field profile")
        ui.live_chart(profile_chart(x * 1e9, ratio, b_0, units, t), key="slab_profile")
    with family_col, ui.card("card_family"):
        ui.card_title("Temperature family")
        ui.live_chart(ui.temperature_profile_chart(t_c, lambda_0, d, b_0, t_k, units, t),
                      key="slab_family")