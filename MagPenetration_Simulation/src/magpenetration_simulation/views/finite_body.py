"""Finite body page: the 2D scalar London field on a rectangle."""

import math

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from magpenetration_simulation import ui, visualization
from magpenetration_simulation import ui_config as cfg
from magpenetration_simulation.finite_body import (
    rect_1d_limit,
    rect_center_field,
    rect_field_finite_difference,
    rect_valid_mask,
)

MARGIN = 3.0


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


def centreline_chart(solution, width, height, lam, b_0, t):
    """B/B₀ up the vertical centreline of the solve, against the 1D slab limit."""
    fig = go.Figure()
    fig.add_scatter(x=solution.y / lam, y=solution.b[:, solution.b.shape[1] // 2] / b_0,
                    mode="lines", name="Rectangle, x = 0",
                    line={"color": t["Finite difference"], "width": 2.4},
                    hovertemplate="%{y:.4f}")
    y, b_slab = rect_1d_limit(width, height, lam, b_0)
    fig.add_scatter(x=y / lam, y=b_slab / b_0, mode="lines", name="1D slab limit",
                    line={"color": t["Relaxation (SOR)"], "width": 1.8, "dash": "dash"},
                    hovertemplate="%{y:.4f}")
    fig.update_xaxes(title="y / λ from the centre")
    fig.update_yaxes(title="B / B₀", range=[0, cfg.PLOT["y_headroom"]])
    return ui.style(fig, t, cfg.CHART_H, uirevision="fb-centre").update_layout(showlegend=True)


def flat_centreline_chart(yf, units, t):
    """Normal state: λ → ∞, so the centreline is flat at B₀."""
    fig = go.Figure()
    fig.add_scatter(x=[-0.5, 0.5], y=[yf, yf], mode="lines", name="B",
                    line={"color": t["Finite difference"], "width": 2.4},
                    hovertemplate="%{y:.4f}")
    fig.update_xaxes(title="y / h from the centre", range=[-0.5, 0.5])
    fig.update_yaxes(title="B (mT)" if units == "B (mT)" else "B / B₀",
                     range=[0, cfg.PLOT["y_headroom"] * yf])
    return ui.style(fig, t, cfg.CHART_H, uirevision=f"fb-centre-normal-{units}")


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


def size_controls(t, ref, ref_label):
    """Width and height from one switch: a ratio to the size length, or nanometres."""
    mode = st.segmented_control("Size by", ["Ratio", "Physical"], default="Physical",
                                key="fb_size_mode") or "Physical"
    if mode == "Physical":
        w_cfg, h_cfg = cfg.FINITE_BODY["width_nm"], cfg.FINITE_BODY["height_nm"]
        width_nm = ui.live_slider("Width w (nm)", w_cfg["min"], w_cfg["max"], w_cfg["default"],
                                  w_cfg["step"], "fb_width_nm", t, decimals=0, unit=" nm")
        height_nm = ui.live_slider("Height h (nm)", h_cfg["min"], h_cfg["max"], h_cfg["default"],
                                   h_cfg["step"], "fb_height_nm", t, decimals=0, unit=" nm")
        width, height = width_nm * 1e-9, height_nm * 1e-9
    else:
        w_cfg, h_cfg = cfg.FINITE_BODY["width"], cfg.FINITE_BODY["height"]
        w_over = ui.live_slider(f"Width w/{ref_label}", w_cfg["min"], w_cfg["max"],
                                w_cfg["default"], w_cfg["step"], "fb_width_ratio", t, decimals=1)
        h_over = ui.live_slider(f"Height h/{ref_label}", h_cfg["min"], h_cfg["max"],
                                h_cfg["default"], h_cfg["step"], "fb_height_ratio", t, decimals=1)
        width, height = w_over * ref, h_over * ref
    st.caption(f"w = {width * 1e9:,.0f} nm · h = {height * 1e9:,.0f} nm · "
               f"w/{ref_label} = {width / ref:.1f} · h/{ref_label} = {height / ref:.1f}")
    return width, height, width / ref, height / ref


def render():
    """Show the scalar field over a rectangle of finite width and height."""
    t = ui.theme()
    ui.page_title("Finite body",
                  "2D scalar London field on a rectangle: a five-point relaxation on the sample, "
                  "with all four edges held at B₀.")

    left, right = st.columns([0.34, 0.66], gap="medium")
    with left:
        material = ui.material_controls("finite_body", compact=True)
        lam, lambda_0, b_0 = material["lambda"], material["lambda_0"], material["b_0"]
        # λ is infinite in the normal state, so a ratio is measured against λ₀ there.
        normal = not math.isfinite(lam)
        ref, ref_label = (lambda_0, "λ₀") if normal else (lam, "λ")

        with ui.card("controls"):
            width, height, w_over_lam, h_over_lam = size_controls(t, ref, ref_label)
            n_cfg = cfg.FINITE_BODY["n"]
            n = int(ui.live_slider("Grid points N", n_cfg["min"], n_cfg["max"],
                                   n_cfg["default"], n_cfg["step"], "fb_points", t))
            units = st.segmented_control("Units", ui.GEOM_VIEWS, default=ui.GEOM_VIEWS[0],
                                         key="fb_units") or ui.GEOM_VIEWS[0]

    yf = b_0 * 1e3 if units == "B (mT)" else 1.0
    cbar = "B (mT)" if units == "B (mT)" else "B / B₀"

    if normal:  # λ → ∞ means 1/λ² = 0, so the field is flat at B₀ and there is nothing to solve
        x = np.linspace(-width / 2.0, width / 2.0, n)
        y = np.linspace(-height / 2.0, height / 2.0, n)
        ratio, centre, mean, valid = np.ones((n, n)), 1.0, 1.0, None
        map_fig = heatmap_chart(x * 1e9, y * 1e9, ratio * yf, cbar, None,
                                f"fb-field-normal-{units}", t)
        centre_fig = flat_centreline_chart(yf, units, t)
    else:
        solution = solve(width, height, lam, b_0, n)
        x, y, ratio = solution.x, solution.y, solution.b / b_0
        centre, mean = rect_center_field(solution) / b_0, float(solution.b.mean()) / b_0
        valid = rect_valid_mask(x, y, width, height, lam, margin=MARGIN)
        map_fig = visualization.finite_body_section(solution, b_0, t)
        centre_fig = centreline_chart(solution, width, height, lam, b_0, t)

    with right:
        with ui.card("card_fb_field"):
            ui.card_title("Field map")
            ui.live_chart(map_fig, key="fb_visual")
            st.caption("x-y cross-section; B is out of plane. The scalar solver has no in-plane "
                       "vector field, so there are no field lines.")
        m = st.columns(2, gap="small")
        m[0].metric("London depth λ(T)", "∞" if normal else f"{lam * 1e9:,.0f} nm",
                    f"at T = {material['t']:.2f} K", delta_color="off", delta_arrow="off")
        m[1].metric("Centre field B(0)/B₀", f"{centre:.3f}", ui.field_text(centre * b_0),
                    delta_color="off", delta_arrow="off",
                    help="Field at the node nearest the rectangle centre.")
        m = st.columns(2, gap="small")
        m[0].metric("Mean field ⟨B⟩/B₀", f"{mean:.3f}", delta_color="off", delta_arrow="off",
                    help="Grid average over the whole sample, edge nodes included.")
        m[1].metric("Valid-zone fraction", "n/a" if valid is None else f"{valid.mean():.3f}",
                    f"{MARGIN:g}λ from every edge", delta_color="off", delta_arrow="off",
                    help="Share of grid points where the 1D slab limit is useful.")
        if not normal:
            st.caption(f"Solver: {solution.sweeps:,} sweeps · residual {ui.sci(solution.residual)} "
                       f"T/m² · {'converged' if solution.converged else 'not converged'}")

    c1, c2 = st.columns(2, gap="medium")
    with c1, ui.card("card_fb_valid"):
        ui.card_title("Valid zone")
        ui.live_chart(heatmap_chart(x * 1e9, y * 1e9, ratio * yf, cbar, valid,
                                    f"fb-valid-{units}-{w_over_lam:g}-{h_over_lam:g}", t),
                      key="fb_valid")
        if normal:
            st.caption("Normal state: λ → ∞, so there is no screened band and no valid zone.")
        else:
            st.caption(f"Dashed: the points at least {MARGIN:g}λ from every edge, where the "
                       "1D slab limit is useful. Edge and corner values are grid-dependent.")
    with c2, ui.card("card_fb_centre"):
        ui.card_title("Centreline comparison")
        ui.live_chart(centre_fig, key="fb_centre")
        if normal:
            st.caption("Normal state: λ → ∞, so the field is flat at B₀ over the sample.")
        else:
            st.caption("B/B₀ up the centreline x = 0 against the 1D slab limit of the same "
                       "height. The two agree only where the width is much longer than λ.")

    caveats()
