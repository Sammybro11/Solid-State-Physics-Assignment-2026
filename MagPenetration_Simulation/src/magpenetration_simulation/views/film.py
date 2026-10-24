"""Film page: a thin film, the Pearl regime, and the exact parallel profile."""

import math

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from magpenetration_simulation import ui, visualization
from magpenetration_simulation import ui_config as cfg
from magpenetration_simulation.film import (
    film_crossover,
    film_parallel_center_field,
    film_parallel_mean_field,
    film_regime,
    inplane_screening_length,
    pearl_length,
    pearl_ratio,
)


def ratio_text(value):
    """A field ratio with three decimals, or scientific notation when tiny."""
    return f"{value:.3f}" if value >= 1e-2 else f"{value:.2e}"


def pearl_chart(d_over_lam, n_points, t):
    """Lambda/lambda = 2/(d/lambda), with the crossover and the current thickness marked."""
    x = np.linspace(cfg.FILM["thickness"]["min"], cfg.FILM["thickness"]["max"], n_points)
    fig = go.Figure()
    fig.add_scatter(x=x, y=2.0 / x, mode="lines", name="Λ/λ",
                    line={"color": t["analytic"], "width": 2.2}, hovertemplate="%{y:.3f}")
    fig.add_vline(x=2.0, line={"color": t["axis"], "dash": "dot", "width": 1})
    fig.add_scatter(x=[d_over_lam], y=[2.0 / d_over_lam], mode="markers", hoverinfo="skip",
                    marker={"size": 11, "color": t["card"],
                            "line": {"color": t["Finite difference"], "width": 2}})
    fig.update_xaxes(title="d / λ")
    fig.update_yaxes(title="Λ / λ")
    return ui.style(fig, t, ui.CHART_H, uirevision="film-pearl")


def profile_chart(x_nm, ratio, b_0, units, t):
    """Parallel field against x through the full film, from -d/2 to d/2."""
    yf = b_0 * 1e3 if units == "B (mT)" else 1.0
    fmt = "%{y:.4f} mT" if units == "B (mT)" else "B/B₀ = %{y:.4f}"
    fig = go.Figure()
    fig.add_vrect(x0=x_nm[0], x1=x_nm[-1], fillcolor=t["slab"], line_width=0, layer="below")
    fig.add_scatter(x=x_nm, y=ratio * yf, mode="lines", name="B",
                    line={"color": t["analytic"], "width": 2.4}, hovertemplate=fmt)
    fig.update_xaxes(title="x (nm)", range=[x_nm[0], x_nm[-1]])
    fig.update_yaxes(title="B (mT)" if units == "B (mT)" else "B / B₀", range=[0, cfg.PLOT["y_headroom"] * yf])
    return ui.style(fig, t, ui.CHART_H, uirevision=f"film-profile-{units}-{x_nm[0]:.6e}-{x_nm[-1]:.6e}")


def map_caption():
    """The section caption: what the picture shows and what it leaves out."""
    st.caption("Longitudinal section: the parallel field runs along z, and the colour is |B|/B₀. "
               "There is no substrate because the film model has none.")


def normal_view(col, t, b_0, d, n_points, units):
    """Flat field at B0 and the why, when screening has vanished."""
    with col, ui.card("card_map"):
        ui.card_title("Field map")
        ui.live_chart(visualization.film_section(d, math.inf, b_0, t), key="film_visual")
        map_caption()

    k = col.columns(4)
    k[0].metric("λ(T)", "∞", "normal state", delta_color="off", delta_arrow="off")
    k[1].metric("Pearl length Λ", "∞", "2λ²/d", delta_color="off", delta_arrow="off")
    k[2].metric("Centre field B(0)/B₀", "1.000", "no screening", delta_color="off", delta_arrow="off")
    k[3].metric("Mean field ⟨B⟩/B₀", "1.000", "volume average", delta_color="off", delta_arrow="off")

    d_nm = d * 1e9
    x_nm = np.linspace(-d_nm / 2.0, d_nm / 2.0, n_points)
    profile_col, info_col = st.columns(2, gap="medium")
    with profile_col, ui.card("card_profile"):
        ui.card_title("Parallel profile")
        ui.live_chart(profile_chart(x_nm, np.ones_like(x_nm), b_0, units, t), key="film_profile")
        st.caption("Normal state: λ → ∞ and Λ → ∞, so the field is flat at B₀ and the film screens nothing.")
    with info_col:
        ui.normal_state_infographic(t)


def render():
    """Show the Pearl length, the exact parallel profile, and the regime."""
    t = ui.theme()
    ui.page_title("Film", "Thin film of thickness d: the Pearl length, the regime, and the exact parallel profile.")

    left, right = st.columns([0.34, 0.66], gap="medium")
    with left:
        material = ui.material_controls("film", compact=True)
        with ui.card("controls"):
            size = ui.size_input("film", cfg.FILM["thickness"], cfg.FILM["d_nm"],
                                 "Film thickness d/λ", "Film thickness d (nm)", material["lambda"], t)
            n_points = int(ui.live_slider("Grid points N (odd)", cfg.GRID["min"], cfg.GRID["max"],
                                          cfg.GRID["default"], cfg.GRID["step"], "film_points", t))
            units = st.segmented_control("Units", ui.GEOM_VIEWS, default=ui.GEOM_VIEWS[0],
                                         key="film_units") or ui.GEOM_VIEWS[0]

    lam, lambda_0 = material["lambda"], material["lambda_0"]
    t_c, t_k, b_0 = material["t_c"], material["t"], material["b_0"]
    # λ is infinite in the normal state, so a ratio-sized film is taken in units of λ0.
    d = size["size_m"] if math.isfinite(size["size_m"]) else size["ratio"] * lambda_0
    thickness_ratio = size["ratio"]

    if material["normal_state"]:
        normal_view(right, t, b_0, d, n_points, units)
        st.caption("The perpendicular-field film, its edge currents, and demagnetization "
                   "are not solved here.")
        return

    pearl = pearl_length(lam, d)
    crossing = film_crossover(lam)
    regime = film_regime(lam, d)
    centre = film_parallel_center_field(d, lam, b_0) / b_0
    mean = film_parallel_mean_field(d, lam, b_0) / b_0

    with right:
        with ui.card("card_map"):
            ui.card_title("Field map")
            ui.live_chart(visualization.film_section(d, lam, b_0, t), key="film_visual")
            map_caption()
        m = st.columns(4)
        m[0].metric("λ(T)", f"{lam * 1e9:.1f} nm", f"λ₀ = {lambda_0 * 1e9:.0f} nm",
                    delta_color="off", delta_arrow="off")
        m[1].metric("Pearl length Λ", f"{pearl * 1e9:.3g} nm", "2λ²/d",
                    delta_color="off", delta_arrow="off",
                    help="Λ sets in-plane screening in a perpendicular field.")
        m[2].metric("Λ / λ", ratio_text(pearl_ratio(lam, d)), "2λ/d",
                    delta_color="off", delta_arrow="off",
                    help="Above 1 for d < 2λ: the film screens in plane less than λ suggests.")
        m[3].metric("In-plane screening length", f"{inplane_screening_length(lam, d) * 1e9:.3g} nm",
                    "max(Λ, λ)", delta_color="off", delta_arrow="off")
        m = st.columns(4)
        m[0].metric("Regime", regime.capitalize(), f"d/λ = {thickness_ratio:.1f}",
                    delta_color="off", delta_arrow="off")
        m[1].metric("Centre field B(0)/B₀", ratio_text(centre), ui.field_text(centre * b_0),
                    delta_color="off", delta_arrow="off",
                    help="Parallel field at the film midplane.")
        m[2].metric("Mean field ⟨B⟩/B₀", ratio_text(mean), "volume average",
                    delta_color="off", delta_arrow="off",
                    help="⟨B⟩/B₀ = (2λ/d) tanh(d/2λ).")
        m[3].metric("Crossover thickness 2λ", f"{crossing * 1e9:.1f} nm", "Λ = λ",
                    delta_color="off", delta_arrow="off")
        if regime == "pearl":
            st.warning(f"Pearl regime: d/λ = {thickness_ratio:g} < 2. Λ = 2λ²/d = {pearl * 1e9:.3g} nm "
                       "replaces λ for in-plane fields in a perpendicular field.")

    profile_col, pearl_col = st.columns(2, gap="medium")
    with profile_col, ui.card("card_profile"):
        ui.card_title("Parallel profile")
        ui.live_chart(ui.temperature_profile_chart(t_c, lambda_0, d, b_0, t_k, units, t),
                      key="film_profile")
        st.caption("All temperatures meet at B₀ at both faces; the current T is bold. "
                   "The parallel-field profile is exact at every thickness.")
    with pearl_col, ui.card("card_pearl"):
        ui.card_title("Pearl length")
        ui.live_chart(pearl_chart(thickness_ratio, n_points, t), key="film_pearl")
        st.caption("Λ/λ = 2/(d/λ). The dot is the current thickness and the dashed line is the crossover at d = 2λ.")

    st.caption("The perpendicular-field film, its edge currents, and demagnetization "
               "are not solved here.")
