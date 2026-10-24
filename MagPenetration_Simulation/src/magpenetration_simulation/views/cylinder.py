"""Cylinder page: the axial infinite cylinder in the London model."""

import plotly.graph_objects as go
import streamlit as st

from magpenetration_simulation import ui, visualization
from magpenetration_simulation import ui_config as cfg
from magpenetration_simulation.cylinder import (
    cylinder_center_field,
    cylinder_field_profile,
    cylinder_mean_field,
)
from magpenetration_simulation.temperature import lambda_of_temperature


def family_chart(lam, b_0, size, radius, view, t):
    """The R/lam family, with the current size bold even when it is between family values."""
    fig = ui.cylinder_family_chart(lam, b_0, size, view, t)
    if size not in ui.GEOM_FAMILY:
        yf = 1.0 if view == ui.GEOM_VIEWS[0] else b_0 * 1e3
        r, b = cylinder_field_profile(radius, lam, b_0, ui.TEMP_POINTS)
        fig.add_scatter(x=r / radius, y=b / b_0 * yf, name=f"R/λ = {size:g}", mode="lines",
                        line={"color": t["Shooting"], "width": 2.6}, hovertemplate="%{y:.4f}")
    return fig


def flat_chart(b_0, view, t):
    """Flat field at B0 across the radius when screening has vanished."""
    yf = 1.0 if view == ui.GEOM_VIEWS[0] else b_0 * 1e3
    fig = go.Figure()
    fig.add_scatter(x=[0.0, 1.0], y=[yf, yf], mode="lines", name="B",
                    line={"color": t["analytic"], "width": 2.6}, hovertemplate="%{y:.4f}")
    fig.update_xaxes(title="r / R", range=[0, 1])
    fig.update_yaxes(title="B / B₀" if view == ui.GEOM_VIEWS[0] else "B (mT)", range=[0, cfg.PLOT["y_headroom"] * yf])
    return ui.style(fig, t, ui.CHART_H, uirevision=f"cyl-flat-{view}")


def temperature_chart(radius, lambda_0, t_c, b_0, t_k, view, t):
    """Profiles over the fixed T/Tc family at the current radius, current T in bold."""
    yf = 1.0 if view == ui.GEOM_VIEWS[0] else b_0 * 1e3
    last = len(ui.FAMILY_T) - 1
    fig = go.Figure()
    for i, f in enumerate(ui.FAMILY_T):
        lam = lambda_of_temperature(f * t_c, lambda_0, t_c)
        r, b = cylinder_field_profile(radius, lam, b_0, ui.TEMP_POINTS)
        fig.add_scatter(x=r / radius, y=b / b_0 * yf, mode="lines", hoverinfo="skip",
                        line={"color": ui.alpha(t["Finite difference"], 0.22 + 0.5 * i / last),
                              "width": 1.4})
    lam = lambda_of_temperature(t_k, lambda_0, t_c)
    r, b = cylinder_field_profile(radius, lam, b_0, ui.TEMP_POINTS)
    fig.add_scatter(x=r / radius, y=b / b_0 * yf, name="T", mode="lines",
                    line={"color": t["analytic"], "width": 3}, hovertemplate="%{y:.4f}")
    fig.update_xaxes(title="r / R")
    fig.update_yaxes(title="B / B₀" if view == ui.GEOM_VIEWS[0] else "B (mT)", range=[0, cfg.PLOT["y_headroom"] * yf])
    return ui.style(fig, t, ui.CHART_H, uirevision=f"cyl-temp-{view}")


LIMITS = ("Model limits: axial field, infinite length, local London theory. "
          "The ends of the cylinder are not included.")
MAP_NOTE = ("Longitudinal r–z section: the axial field as a heatmap of |B|/B₀. "
            "The field lines are straight: the model is axial and infinite in length.")


def render():
    """Show the field map, the radial profiles and the R/λ family of an infinite cylinder."""
    t = ui.theme()
    ui.page_title("Cylinder", "Axial field in an infinite cylinder, screened over the London depth λ.")

    left, right = st.columns([0.34, 0.66], gap="medium")
    with left:
        material = ui.material_controls("cylinder", compact=True)
        b_0, lam = material["b_0"], material["lambda"]
        with ui.card("controls"):
            size = ui.size_input("cylinder", cfg.CYLINDER["size"], cfg.CYLINDER["r_nm"],
                                 "Radius R/λ", "Radius R (nm)", lam, t)
            view = st.segmented_control("Units", ui.GEOM_VIEWS, default=ui.GEOM_VIEWS[0],
                                        label_visibility="collapsed",
                                        key="cylinder_units") or ui.GEOM_VIEWS[0]

    radius, ratio = size["size_m"], size["ratio"]

    if material["normal_state"]:
        # R/λ is 0 when set in Physical and R is infinite when set by ratio, so size the map by hand.
        radius = size["size_m"] if size["mode"] == "Physical" else size["ratio"] * material["lambda_0"]
        with right, ui.card("card_cyl_map"):
            ui.card_title("Field map")
            ui.live_chart(visualization.cylinder_section(radius, lam, b_0, t), key="cyl_visual")
            st.caption(MAP_NOTE)
        with right:
            k = st.columns(2)
            k[0].metric("Surface and axis field", ui.field_text(b_0), "flat at B₀",
                        delta_color="off", delta_arrow="off")
            k[1].metric("Screening", "none", "λ → ∞", delta_color="off", delta_arrow="off")

        c1, c2 = st.columns(2, gap="medium")
        with c1, ui.card("card_cyl_profile"):
            ui.card_title("Radial profile")
            ui.live_chart(flat_chart(b_0, view, t), key="cyl_profile")
            st.caption("Normal state: λ → ∞, so the field is flat at B₀ from the wall to the axis.")
        with c2:
            ui.normal_state_infographic(t)
        st.caption(LIMITS)
        return

    center_ratio = cylinder_center_field(radius, lam, b_0) / b_0
    mean_ratio = cylinder_mean_field(radius, lam, b_0) / b_0
    centre_text = f"{center_ratio:.3f}" if center_ratio > 1e-3 else f"{center_ratio:.1e}"

    with right, ui.card("card_cyl_map"):
        ui.card_title("Field map")
        ui.live_chart(visualization.cylinder_section(radius, lam, b_0, t), key="cyl_visual")
        st.caption(MAP_NOTE)
    with right:
        k = (*st.columns(2), *st.columns(2))
        k[0].metric("London depth λ", f"{lam * 1e9:,.0f} nm", f"at T = {material['t']:.2f} K",
                    delta_color="off", delta_arrow="off")
        k[1].metric("Centre field B(0)/B₀", centre_text, help="1 / I₀(R/λ).")
        k[2].metric("Mean field ⟨B⟩/B₀", f"{mean_ratio:.3f}",
                    help="(2λ/R) I₁(R/λ) / I₀(R/λ), volume-averaged.")
        k[3].metric("Screened fraction", f"{1.0 - center_ratio:.4f}",
                    help="1 − B(0)/B₀: the share of the applied field kept out of the axis.")
        if ratio < 1:
            st.warning("Weak screening: R/λ < 1, so the field reaches the axis almost unchanged.")

    c1, c2 = st.columns(2, gap="medium")
    with c1, ui.card("card_cyl_profile"):
        ui.card_title("Radial profile")
        ui.live_chart(temperature_chart(radius, material["lambda_0"], material["t_c"],
                                        b_0, material["t"], view, t), key="cyl_profile")
        st.caption("Faint: T/T_c = " + ", ".join(f"{f:g}" for f in ui.FAMILY_T)
                   + ". Bold: the current T. R is held fixed, so R/λ falls as T rises.")
    with c2, ui.card("card_cyl_family"):
        ui.card_title("R/λ family")
        ui.live_chart(family_chart(lam, b_0, ratio, radius, view, t), key="cyl_family")

    st.caption(LIMITS)
