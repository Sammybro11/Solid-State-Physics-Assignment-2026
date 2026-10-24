"""Film page: a thin film, the Pearl regime, and the exact parallel profile."""

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from magpenetration_simulation import ui
from magpenetration_simulation import ui_config as cfg
from magpenetration_simulation.film import (
    film_crossover,
    film_parallel_center_field,
    film_parallel_mean_field,
    film_parallel_profile,
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


def normal_view(t, b_0, lambda_0, d_over_lam, n_points, units):
    """Flat field at B0 and the why, when screening has vanished."""
    k = st.columns(4)
    k[0].metric("λ(T)", "∞", "normal state", delta_color="off", delta_arrow="off")
    k[1].metric("Pearl length Λ", "∞", "2λ²/d", delta_color="off", delta_arrow="off")
    k[2].metric("Centre field B(0)/B₀", "1.000", "no screening", delta_color="off", delta_arrow="off")
    k[3].metric("Mean field ⟨B⟩/B₀", "1.000", "volume average", delta_color="off", delta_arrow="off")

    # λ is infinite, so the window is sized in units of λ0 instead of λ.
    span_nm = d_over_lam * lambda_0 * 1e9
    x_nm = np.linspace(-span_nm / 2.0, span_nm / 2.0, n_points)
    with ui.card("card_profile"):
        ui.card_title("Parallel field profile")
        ui.live_chart(profile_chart(x_nm, np.ones_like(x_nm), b_0, units, t), key="film_profile")
    st.caption("Normal state: λ → ∞ and Λ → ∞, so the field is flat at B₀ and the film screens nothing. "
               "The d/λ slider sets the window in units of λ₀.")
    ui.normal_state_infographic(t)


def render(material):
    """Show the Pearl length, the exact parallel profile, and the regime."""
    t = ui.theme()
    lam, lambda_0 = material["lambda"], material["lambda_0"]
    t_c, t_k, b_0 = material["t_c"], material["t"], material["b_0"]

    ui.page_title("Film", "Thin film of thickness d: the Pearl length, the regime, and the exact parallel profile.")

    with ui.card("controls"):
        c = st.columns([1, 1, 0.7], gap="large", vertical_alignment="bottom")
        with c[0]:
            f_cfg = cfg.FILM["thickness"]
            d_over_lam = ui.live_slider("Film thickness d/λ", f_cfg["min"], f_cfg["max"],
                                        f_cfg["default"], f_cfg["step"], "d_over_lam", t,
                                        decimals=1)
        with c[1]:
            n_points = int(ui.live_slider("Grid points N (odd)", cfg.GRID["min"], cfg.GRID["max"],
                                          cfg.GRID["default"], cfg.GRID["step"], "n_points", t))
        units = c[2].segmented_control("Units", ui.GEOM_VIEWS, default=ui.GEOM_VIEWS[0],
                                       key="film_units") or ui.GEOM_VIEWS[0]

    if material["normal_state"]:
        normal_view(t, b_0, lambda_0, d_over_lam, n_points, units)
        return

    d = d_over_lam * lam
    pearl = pearl_length(lam, d)
    crossing = film_crossover(lam)
    regime = film_regime(lam, d)
    centre = film_parallel_center_field(d, lam, b_0) / b_0
    mean = film_parallel_mean_field(d, lam, b_0) / b_0

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
    m[0].metric("Regime", regime.capitalize(), f"d/λ = {d_over_lam:.1f}",
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
        st.warning(f"Pearl regime: d/λ = {d_over_lam:g} < 2. Λ = 2λ²/d = {pearl * 1e9:.3g} nm "
                   "replaces λ for in-plane fields in a perpendicular field.")

    pearl_col, profile_col = st.columns([1, 1.35], gap="medium")
    with pearl_col, ui.card("card_pearl"):
        ui.card_title("Pearl length")
        ui.live_chart(pearl_chart(d_over_lam, n_points, t), key="film_pearl")
        st.caption("Λ/λ = 2/(d/λ). The dot is the current thickness and the dashed line is the crossover at d = 2λ.")
    with profile_col, ui.card("card_profile"):
        ui.card_title("Parallel field profile")
        x = np.linspace(-d / 2.0, d / 2.0, n_points)
        ui.live_chart(profile_chart(x * 1e9, film_parallel_profile(x, d, lam, b_0) / b_0,
                                    b_0, units, t), key="film_profile")
        st.caption("This parallel-field profile is exact at every thickness.")

    with ui.card("card_family"):
        ui.card_title("Temperature family")
        ui.live_chart(ui.temperature_profile_chart(t_c, lambda_0, d, b_0, t_k, units, t),
                      key="film_family")
        st.caption("All temperatures meet at B₀ at both faces. The current T is bold.")

    st.caption("The perpendicular-field film, its edge currents, and demagnetization are not solved here. "
               "The 2D field-line view is a later slice.")
