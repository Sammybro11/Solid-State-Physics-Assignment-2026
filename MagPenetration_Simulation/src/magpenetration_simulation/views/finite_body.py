"""Finite body page: the 2D scalar London field on a rectangle."""

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from magpenetration_simulation import ui
from magpenetration_simulation import ui_config as cfg
from magpenetration_simulation.finite_body import (
    rect_center_field,
    rect_field_finite_difference,
    rect_valid_mask,
)

MARGIN = 3.0
LATER = "The 2D field-line view of the same solve is a later slice."


@st.cache_data(show_spinner=False)
def solve(width, height, lam, b_0, n):
    """Cached red-black SOR solve, so a slider rerun reuses the same grid result."""
    return rect_field_finite_difference(width, height, lam, b_0, nx=n, ny=n)


def heatmap_chart(x_nm, y_nm, z, cbar, valid, uirevision, t):
    """B over the rectangle; the dashed outline is the valid zone when valid is given."""
    lo, hi = float(np.min(z)), float(np.max(z))
    if hi - lo <= 1e-12:
        lo = 0.0
    fig = go.Figure(go.Heatmap(
        x=x_nm, y=y_nm, z=z, zmin=lo, zmax=hi,
        colorscale=[[0.0, t["card"]], [1.0, t["Finite difference"]]],
        hovertemplate="x %{x:.0f} nm<br>y %{y:.0f} nm<br>%{z:.4f}<extra></extra>",
        colorbar={"title": {"text": cbar}, "thickness": 12, "outlinewidth": 0,
                  "tickfont": {"size": 11, "color": t["muted"]}}))
    if valid is not None:
        fig.add_trace(go.Contour(
            x=x_nm, y=y_nm, z=valid.astype(float), showscale=False, hoverinfo="skip",
            contours={"coloring": "none", "showlines": True,
                      "start": 0.5, "end": 0.5, "size": 1.0},
            line={"color": t["Relaxation (SOR)"], "width": 1.5, "dash": "dash"}))
    ui.style(fig, t, ui.CHART_H, uirevision=uirevision,
             yaxis={"scaleanchor": "x", "scaleratio": 1})
    fig.update_xaxes(title="x (nm)")
    fig.update_yaxes(title="y (nm)")
    return fig


def caveats():
    """The model limits, shown in both the normal and the superconducting state."""
    with ui.card("card_fb_caveat"):
        ui.card_title("What this solve does not include")
        st.markdown(
            "- Scalar B_z on the sample only. There are no in-plane components.\n"
            "- Dirichlet B = B₀ on all four edges.\n"
            "- No outside-field matching, no demagnetization, and no turning of the screening "
            "current at an edge.\n"
            "- Edge and corner values are grid-dependent: refine N and watch them move.\n"
            f"- The dashed valid zone marks the points at least {MARGIN:g}λ from every edge, "
            "where the 1D slab limit is useful.\n"
            "- Box and disk geometries are later slices."
        )


def normal_view(t, b_0, lambda_0, w_over_lam, h_over_lam, n, units):
    """Flat field at B₀: λ is infinite, so 1/λ² = 0 and the solve is a constant."""
    k = st.columns(2)
    k[0].metric("λ(T)", "∞", "normal state", delta_color="off", delta_arrow="off")
    k[1].metric("Field B/B₀", "1.000", "1/λ² = 0", delta_color="off", delta_arrow="off",
                help="With ∇²B = B/λ² → 0 and B = B₀ on the edges, B is B₀ everywhere.")

    width = w_over_lam * lambda_0
    height = h_over_lam * lambda_0
    x = np.linspace(-width / 2.0, width / 2.0, n)
    y = np.linspace(-height / 2.0, height / 2.0, n)
    yf = b_0 * 1e3 if units == "B (mT)" else 1.0
    cbar = "B (mT)" if units == "B (mT)" else "B / B₀"

    with ui.card("card_fb_field"):
        ui.card_title("Field over the sample")
        ui.live_chart(heatmap_chart(x * 1e9, y * 1e9, np.ones((n, n)) * yf, cbar, None,
                                    f"fb-field-normal-{units}", t), key="fb_field")
    st.caption("Normal state: λ → ∞, so the field is flat at B₀ and the box is sized in λ₀.")
    caveats()
    st.caption(LATER)


def render(material):
    """Show the scalar field over a rectangle of finite width and height."""
    t = ui.theme()
    lam, lambda_0, b_0 = material["lambda"], material["lambda_0"], material["b_0"]

    ui.page_title("Finite body",
                  "2D scalar London field on a rectangle: a five-point relaxation on the sample, "
                  "with all four edges held at B₀.")

    with ui.card("controls"):
        c = st.columns([1, 1, 1, 1], gap="large", vertical_alignment="bottom")
        with c[0]:
            w_cfg = cfg.FINITE_BODY["width"]
            w_over_lam = ui.live_slider("Width w/λ", w_cfg["min"], w_cfg["max"],
                                        w_cfg["default"], w_cfg["step"], "fb_width", t, decimals=1)
        with c[1]:
            h_cfg = cfg.FINITE_BODY["height"]
            h_over_lam = ui.live_slider("Height h/λ", h_cfg["min"], h_cfg["max"],
                                        h_cfg["default"], h_cfg["step"], "fb_height", t, decimals=1)
        with c[2]:
            n_cfg = cfg.FINITE_BODY["n"]
            n = int(ui.live_slider("Grid points N", n_cfg["min"], n_cfg["max"],
                                   n_cfg["default"], n_cfg["step"], "fb_points", t))
        units = c[3].segmented_control("Units", ui.GEOM_VIEWS, default=ui.GEOM_VIEWS[0],
                                       key="fb_units") or ui.GEOM_VIEWS[0]

    if material["normal_state"]:
        normal_view(t, b_0, lambda_0, w_over_lam, h_over_lam, n, units)
        return

    width = w_over_lam * lam
    height = h_over_lam * lam
    solution = solve(width, height, lam, b_0, n)

    centre = rect_center_field(solution) / b_0
    mean = float(solution.b.mean()) / b_0
    valid = rect_valid_mask(solution.x, solution.y, width, height, lam, margin=MARGIN)

    k = st.columns(4)
    k[0].metric("λ(T)", f"{lam * 1e9:,.0f} nm", f"at T = {material['t']:.2f} K",
                delta_color="off", delta_arrow="off")
    k[1].metric("Centre field B(0)/B₀", f"{centre:.3f}", ui.field_text(centre * b_0),
                delta_color="off", delta_arrow="off",
                help="Field at the node nearest the rectangle centre.")
    k[2].metric("Mean field ⟨B⟩/B₀", f"{mean:.3f}", delta_color="off", delta_arrow="off",
                help="Grid average over the whole sample, edge nodes included.")
    k[3].metric("Valid-zone fraction", f"{valid.mean():.3f}",
                f"{MARGIN:g}λ from every edge", delta_color="off", delta_arrow="off",
                help="Share of grid points where the 1D slab limit is useful.")

    k = st.columns(3)
    k[0].metric("Residual", ui.sci(solution.residual), "T/m²", delta_color="off", delta_arrow="off",
                help="Largest finite-difference residual in the interior at the last sweep.")
    k[1].metric("Sweeps", f"{solution.sweeps:,}", "red-black SOR", delta_color="off",
                delta_arrow="off", help="Sweep pairs until the residual met the tolerance.")
    k[2].metric("Status", "converged" if solution.converged else "max sweeps", delta_color="off",
                delta_arrow="off",
                help="Residual at or below 1e-9 of B₀/λ²; max sweeps means it did not get there.")

    yf = b_0 * 1e3 if units == "B (mT)" else 1.0
    with ui.card("card_fb_field"):
        ui.card_title("Field over the sample")
        ui.live_chart(heatmap_chart(
            solution.x * 1e9, solution.y * 1e9, solution.b / b_0 * yf,
            "B (mT)" if units == "B (mT)" else "B / B₀", valid,
            f"fb-field-{units}-{w_over_lam:g}-{h_over_lam:g}", t), key="fb_field")
        st.caption(f"Corners and edges are held at B₀, so the refinement there is grid-dependent. "
                   f"Dashed: the valid zone, at least {MARGIN:g}λ from every edge.")

    caveats()
    st.caption(LATER)
