"""Cylinder page: the axial infinite cylinder in the London model."""

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from magpenetration_simulation import ui
from magpenetration_simulation.cylinder import (
    cylinder_center_field,
    cylinder_field_profile,
    cylinder_mean_field,
)
from magpenetration_simulation.temperature import lambda_of_temperature


def near_wall(radius, lam, n_points):
    """Thick-cylinder wall limit sqrt(R/r) * exp(-s/lam) with s = R - r, as r/R and B/B0.

    Only the part of the grid where the limit stays near the field (value < 1.02)
    is returned; closer to the axis the limit diverges and does not apply.
    """
    r = np.linspace(0.0, radius, n_points)
    with np.errstate(divide="ignore"):
        ref = np.sqrt(radius / r) * np.exp(-(radius - r) / lam)
    keep = ref <= 1.02
    return r[keep] / radius, ref[keep]


def profile_chart(r, b, radius, lam, size, b_0, view, t):
    """Radial field against r/R, with the near-wall thick-cylinder limit when R/lam >= 10."""
    yf = 1.0 if view == ui.GEOM_VIEWS[0] else b_0 * 1e3
    fig = go.Figure()
    fig.add_scatter(x=r / radius, y=b / b_0 * yf, name="London", mode="lines",
                    line={"color": t["analytic"], "width": 2.6}, hovertemplate="%{y:.4f}")
    if size >= 10:
        x, y = near_wall(radius, lam, len(r))
        fig.add_scatter(x=x, y=y * yf, mode="lines", hoverinfo="skip",
                        line={"color": t["Shooting"], "width": 1.6, "dash": "dash"})
    fig.update_xaxes(title="r / R")
    fig.update_yaxes(title="B / B₀" if view == ui.GEOM_VIEWS[0] else "B (mT)", range=[0, 1.05 * yf])
    return ui.style(fig, t, ui.CHART_H, uirevision=f"cyl-profile-{view}")


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
    fig.update_yaxes(title="B / B₀" if view == ui.GEOM_VIEWS[0] else "B (mT)", range=[0, 1.05 * yf])
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
    fig.update_yaxes(title="B / B₀" if view == ui.GEOM_VIEWS[0] else "B (mT)", range=[0, 1.05 * yf])
    return ui.style(fig, t, ui.CHART_H, uirevision=f"cyl-temp-{view}")


LIMITS = ("Model limits: axial field, infinite length, local London theory. "
          "The ends of the cylinder are not included.")


def render(material):
    """Show the radial field of an axial infinite cylinder at the shared lambda and B0."""
    t = ui.theme()
    b_0 = material["b_0"]
    ui.page_title("Cylinder", "Axial field in an infinite cylinder, screened over the London depth λ.")

    with ui.card("controls"):
        c = st.columns([1.3, 1.1, 1], gap="large")
        with c[0]:
            size = float(ui.live_slider("Radius R/λ", 0.5, 20.0, 5.0, 0.5, "c_size", t, decimals=1))
        with c[1]:
            n_points = int(ui.live_slider("Grid points N (odd)", 41, 2001, 401, 20, "c_points", t))
        with c[2]:
            view = st.segmented_control("Units", ui.GEOM_VIEWS, default=ui.GEOM_VIEWS[0],
                                        label_visibility="collapsed", key="c_units") or ui.GEOM_VIEWS[0]

    if material["normal_state"]:
        k = st.columns(2)
        k[0].metric("Surface and axis field", ui.field_text(b_0), "flat at B₀",
                    delta_color="off", delta_arrow="off")
        k[1].metric("Screening", "none", "λ → ∞", delta_color="off", delta_arrow="off")

        with ui.card("card_cyl_profile"):
            ui.card_title("Radial profile")
            ui.live_chart(flat_chart(b_0, view, t), key="cyl_profile")
        st.caption("Normal state: λ → ∞, so the field is flat at B₀ from the wall to the axis.")
        ui.normal_state_infographic(t)
        st.caption(LIMITS)
        return

    lam = material["lambda"]
    radius = size * lam
    r, b = cylinder_field_profile(radius, lam, b_0, n_points)

    center_ratio = cylinder_center_field(radius, lam, b_0) / b_0
    mean_ratio = cylinder_mean_field(radius, lam, b_0) / b_0
    centre_text = f"{center_ratio:.3f}" if center_ratio > 1e-3 else f"{center_ratio:.1e}"

    k = st.columns(4)
    k[0].metric("London depth λ", f"{lam * 1e9:,.0f} nm", f"at T = {material['t']:.2f} K",
                delta_color="off", delta_arrow="off")
    k[1].metric("Centre field B(0)/B₀", centre_text, help="1 / I₀(R/λ).")
    k[2].metric("Mean field ⟨B⟩/B₀", f"{mean_ratio:.3f}",
                help="(2λ/R) I₁(R/λ) / I₀(R/λ), volume-averaged.")
    k[3].metric("Screened fraction", f"{1.0 - center_ratio:.4f}",
                help="1 − B(0)/B₀: the share of the applied field kept out of the axis.")

    if size < 1:
        st.warning("Weak screening: R/λ < 1, so the field reaches the axis almost unchanged.")

    left, right = st.columns(2, gap="medium")
    with left, ui.card("card_cyl_profile"):
        ui.card_title("Radial profile")
        ui.live_chart(profile_chart(r, b, radius, lam, size, b_0, view, t), key="cyl_profile")
        if size >= 10:
            st.caption("Dashed: √(R/r)·exp(−s/λ), the thick-cylinder limit near the wall "
                       "(s = R − r, valid for R/λ ≥ 10).")
    with right, ui.card("card_cyl_family"):
        ui.card_title("R/λ family")
        ui.live_chart(family_chart(lam, b_0, size, radius, view, t), key="cyl_family")

    with ui.card("card_cyl_temp"):
        ui.card_title("Temperature family at the current radius")
        ui.live_chart(temperature_chart(radius, material["lambda_0"], material["t_c"],
                                        b_0, material["t"], view, t), key="cyl_temp")
        st.caption("Faint: T/T_c = " + ", ".join(f"{f:g}" for f in ui.FAMILY_T)
                   + ". Bold: the current T from the sidebar. R is held fixed, so R/λ falls as T rises.")

    st.caption(LIMITS)
