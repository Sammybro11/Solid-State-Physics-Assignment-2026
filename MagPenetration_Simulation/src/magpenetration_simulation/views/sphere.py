"""Sphere page: the exact London sphere solution, shown without a solver."""

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from magpenetration_simulation import ui, visualization
from magpenetration_simulation import ui_config as cfg
from magpenetration_simulation.sphere import (
    sphere_center_field,
    sphere_dipole_moment,
    sphere_equator_field,
    sphere_field,
    sphere_mean_field,
    sphere_screening_factor,
)
from magpenetration_simulation.temperature import lambda_of_temperature

NOTE = ("Exact vector London solution for a sphere: local theory, demagnetization "
        "included through the vector boundary match. Vortices "
        "are not modeled.")


def ratio_text(value):
    """A field ratio with three decimals, or scientific notation when tiny."""
    return f"{value:.3f}" if value >= 1e-2 else f"{value:.2e}"


def y_top(ratio):
    """Y top that never clips the peak: the headroom, or 5% above the data."""
    return max(cfg.PLOT["y_headroom"], 1.05 * float(np.max(ratio)))

def surface_chart(theta, ratio, view, b_0, t):
    """|B(R, theta)| against polar angle, with the equator and both poles marked."""
    yf = b_0 * 1e3 if view == ui.GEOM_VIEWS[1] else 1.0
    fig = go.Figure()
    fig.add_scatter(x=theta, y=ratio * yf, mode="lines", name="B",
                    line={"color": t["analytic"], "width": 2.4}, hovertemplate="%{y:.4f}")
    for x, text, anchor in ((0.0, "pole", "left"), (np.pi / 2.0, "equator", "center"),
                            (np.pi, "pole", "right")):
        fig.add_vline(x=x, line={"color": t["axis"], "dash": "dot", "width": 1})
        fig.add_annotation(x=x, y=0.04, yref="paper", text=text, showarrow=False,
                           xanchor=anchor, font={"size": 11, "color": t["muted"]})
    fig.update_xaxes(title="θ (rad, 0 = north pole)", range=[0, np.pi])
    fig.update_yaxes(title="B (mT)" if view == ui.GEOM_VIEWS[1] else "B / B₀",
                     range=[0, y_top(ratio) * yf])
    return ui.style(fig, t, ui.CHART_H, uirevision=f"sphere-surface-{view}")


def temperature_chart(radius, lambda_0, t_c, b_0, t_k, view, t):
    """Equatorial radial profiles at fixed R over the T/Tc family, current T in bold."""
    yf = b_0 * 1e3 if view == ui.GEOM_VIEWS[1] else 1.0
    last = len(ui.FAMILY_T) - 1
    r = np.linspace(0.0, radius, ui.TEMP_POINTS)
    x = r / radius
    fig = go.Figure()
    peak = 0.0
    for i, f in enumerate(ui.FAMILY_T):
        lam_f = lambda_of_temperature(f * t_c, lambda_0, t_c)
        _, b_theta = sphere_field(r, np.pi / 2.0, radius, lam_f, b_0)
        ratio = np.abs(b_theta) / b_0
        peak = max(peak, float(np.max(ratio)))
        fig.add_scatter(x=x, y=ratio * yf, mode="lines", hoverinfo="skip",
                        line={"color": ui.alpha(t["Finite difference"], 0.22 + 0.5 * i / last),
                              "width": 1.4})
    lam = lambda_of_temperature(t_k, lambda_0, t_c)
    _, b_theta = sphere_field(r, np.pi / 2.0, radius, lam, b_0)
    ratio = np.abs(b_theta) / b_0
    fig.add_scatter(x=x, y=ratio * yf, name="T", mode="lines",
                    line={"color": t["analytic"], "width": 3}, hovertemplate="%{y:.4f}")
    fig.update_xaxes(title="r / R")
    fig.update_yaxes(title="B (mT)" if view == ui.GEOM_VIEWS[1] else "B / B₀",
                     range=[0, max(cfg.PLOT["y_headroom"], 1.05 * peak) * yf])
    return ui.style(fig, t, ui.CHART_H,
                    uirevision=f"sphere-temp-{view}-{radius}-{t_k}")


def normal_view(t):
    """The exact sphere has no normal-state limit here, so say so and stop."""
    with ui.card("card_sphere_normal"):
        ui.card_title("Normal state")
        st.info("Above T_c the screening length is infinite, so this exact London sphere "
                "has no solution to show. The normal-state demagnetized sphere is not "
                "modeled in this slice.")


def render():
    """Show the field map and metrics, the equator family, and the surface field."""
    t = ui.theme()
    ui.page_title("Sphere", "A sphere of radius R in a uniform field B₀: the equator, the surface, and the demagnetization.")
    left, right = st.columns([0.34, 0.66], gap="medium")

    with left:
        material = ui.material_controls("sphere", compact=True)
        lam = material["lambda"]
        lambda_0, t_c = material["lambda_0"], material["t_c"]
        t_k, b_0 = material["t"], material["b_0"]
        with ui.card("controls"):
            size = ui.size_input("sphere", cfg.SPHERE["size"], cfg.SPHERE["r_nm"],
                                 "Radius R/λ", "Radius R (nm)", lam, t)
            n_points = int(ui.live_slider("Grid points N (odd)", cfg.GRID["min"], cfg.GRID["max"],
                                          cfg.GRID["default"], cfg.GRID["step"], "sphere_points", t))
            view = st.segmented_control("Units", ui.GEOM_VIEWS, default=ui.GEOM_VIEWS[0],
                                        key="sphere_units") or ui.GEOM_VIEWS[0]

    if material["normal_state"]:
        with right:
            normal_view(t)
        st.caption(NOTE)
        return

    radius, ratio = size["size_m"], size["ratio"]
    centre_ratio = sphere_center_field(radius, lam, b_0) / b_0
    mean_ratio = sphere_mean_field(radius, lam, b_0) / b_0
    screening = sphere_screening_factor(radius, lam, b_0)
    equator_ratio = sphere_equator_field(radius, lam, b_0) / b_0
    moment = sphere_dipole_moment(radius, lam, b_0)

    with right:
        with ui.card("card_sphere_map"):
            ui.card_title("Field map")
            ui.live_chart(visualization.sphere_section(radius, lam, b_0, t), key="sphere_visual")
            st.caption("x-z section: real in-plane field lines. Outside is the uniform field "
                       "plus the London dipole; inside is the exact London field.")
        k = st.columns(3, gap="small")
        k[0].metric("London depth λ(T)", f"{lam * 1e9:.1f} nm", f"λ₀ = {lambda_0 * 1e9:.0f} nm",
                    delta_color="off", delta_arrow="off")
        k[1].metric("Centre field B(0)/B₀", ratio_text(centre_ratio),
                    help="B₀ (R/λ) / sinh(R/λ).")
        k[2].metric("Axial mean ⟨B_z⟩/B₀", ratio_text(mean_ratio), help="B₀ (1 − S).")
        k = st.columns(3, gap="small")
        k[0].metric("Screening factor S", f"{screening:.4f}",
                    help="1 − 3(λ/R) coth(R/λ) + 3λ²/R². Sets the exterior dipole.")
        k[1].metric("Equator field B_eq/B₀", ratio_text(equator_ratio),
                    help="(3/2) B₀ (1 − i₁(R/λ)/sinh(R/λ)). Tends to 1.5 B₀.")
        k[2].metric("Dipole moment |m|", f"{ui.sci(abs(moment))} A m²",
                    help="m = −(2πR³B₀/µ₀) S, negative for B₀ > 0 (diamagnetic).")

    theta = np.linspace(0.0, np.pi, n_points)
    b_r_s, b_theta_s = sphere_field(radius, theta, radius, lam, b_0)
    surface = np.sqrt(b_r_s**2 + b_theta_s**2) / b_0

    c1, c2 = st.columns(2, gap="medium")
    with c1, ui.card("card_sphere_radial"):
        ui.card_title("Equatorial radial profile")
        ui.live_chart(temperature_chart(radius, lambda_0, t_c, b_0, t_k, view, t),
                      key="sphere_radial")
        st.caption("|B_θ| at θ = π/2, where B_r = 0. Faint: T/T_c family at fixed R. "
                   "Bold: the current T. A rise in T raises λ and weakens screening.")
    with c2, ui.card("card_sphere_surface"):
        ui.card_title("Surface field against latitude")
        ui.live_chart(surface_chart(theta, surface, view, b_0, t), key="sphere_surface")
        st.caption("Demagnetization: the surface field peaks at the equator and falls "
                   "towards the poles, so the applied field is not the boundary value.")

    st.caption(f"R/λ = {ratio:.1f}. " + NOTE)
