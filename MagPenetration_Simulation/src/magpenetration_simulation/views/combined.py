"""Combined page: every implemented geometry at one common half-size."""

import math
from dataclasses import replace

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from magpenetration_simulation import ui
from magpenetration_simulation import ui_config as cfg
from magpenetration_simulation import visualization
from magpenetration_simulation.cylinder import cylinder_center_field, cylinder_mean_field
from magpenetration_simulation.film import (
    film_parallel_center_field,
    film_parallel_mean_field,
    film_regime,
    pearl_ratio,
)
from magpenetration_simulation.finite_body import (
    rect_center_field,
    rect_field_finite_difference,
)
from magpenetration_simulation.sphere import (
    sphere_center_field,
    sphere_equator_field,
    sphere_field,
    sphere_mean_field,
)


def ratios(size):
    """Centre, mean and screened fraction of B/B0 for the slab and the cylinder.

    size is the common half-size s / lambda = d / (2 lambda) = R / lambda. The
    cylinder helpers take R/lam = size at lambda = 1 and b0 = 1, so they return
    the field ratios directly.
    """
    slab_centre = 1.0 / math.cosh(size)
    slab_mean = math.tanh(size) / size
    cyl_centre = cylinder_center_field(size, 1.0, 1.0)
    cyl_mean = cylinder_mean_field(size, 1.0, 1.0)
    return (slab_centre, slab_mean, 1.0 - slab_centre,
            cyl_centre, cyl_mean, 1.0 - cyl_centre)


@st.cache_data(show_spinner=False)
def finite_body_solution(size):
    """Cached fast square solve of side 2*size at lam = b0 = 1; size is the half-side."""
    n = cfg.FINITE_BODY["n"]["min"]
    return rect_field_finite_difference(2.0 * size, 2.0 * size, 1.0, 1.0, nx=n, ny=n)


def finite_body_ratios(size):
    """Centre and mean B/B0 from the cached square solve."""
    solution = finite_body_solution(size)
    return rect_center_field(solution), float(solution.b.mean())


def comparison_rows(size):
    """Centre, mean and note per geometry at the common half-size, in B/B0 units."""
    n = cfg.FINITE_BODY["n"]["min"]
    fb_centre, fb_mean = finite_body_ratios(size)
    return (
        ("Slab", 1.0 / math.cosh(size), math.tanh(size) / size,
         "Flat face at B₀, so the boundary value is B₀."),
        ("Cylinder", cylinder_center_field(size, 1.0, 1.0),
         cylinder_mean_field(size, 1.0, 1.0), "Axial field, radial 1D solution."),
        ("Sphere", sphere_center_field(size, 1.0, 1.0), sphere_mean_field(size, 1.0, 1.0),
         f"Angular field; equator |B| = {sphere_equator_field(size, 1.0, 1.0):.3f} B₀."),
        ("Finite body", fb_centre, fb_mean,
         f"Square side 2s, {n}×{n} scalar solve with Dirichlet edges."),
        ("Film", film_parallel_center_field(2.0 * size, 1.0, 1.0),
         film_parallel_mean_field(2.0 * size, 1.0, 1.0),
         f"Parallel field equals the slab; Pearl Λ/λ = {pearl_ratio(1.0, 2.0 * size):.2f}, "
         f"{film_regime(1.0, 2.0 * size)} regime."),
    )


def comparison_table(rows):
    """Render the per-geometry rows as a five-column table."""
    st.dataframe(
        [{"Geometry": name, "Center B/B₀": f"{centre:.3f}", "Mean B/B₀": f"{mean:.3f}",
          "Screened": f"{1.0 - centre:.1%}", "Note": note}
         for name, centre, mean, note in rows],
        hide_index=True, width="stretch")


def comparison_chart(rows, t):
    """Centre and mean B/B0 side by side, one group per geometry."""
    names = [row[0] for row in rows]
    fig = go.Figure()
    for label, column, color in (("Centre", 1, t["Finite difference"]),
                                 ("Mean", 2, t["Relaxation (SOR)"])):
        fig.add_bar(x=names, y=[row[column] for row in rows], name=label,
                    marker={"color": color}, hovertemplate="%{y:.3f}<extra>" + label + "</extra>")
    fig.update_layout(barmode="group", showlegend=True,
                      legend={"x": 1.0, "y": 1.0, "xanchor": "right", "yanchor": "top",
                              "orientation": "h", "font": {"size": 12}})
    fig.update_yaxes(title="B / B₀", range=[0, 1.25])
    return ui.style(fig, t, ui.CHART_H, uirevision="combined-bars")


def combined_chart(s, slab, s_cyl, cyl, s_sphere, sphere, s_fb, fb, view, b_0, size, temp, t):
    """Every geometry against depth from the nearest face, in the selected units."""
    yf = 1.0 if view == cfg.GEOM_VIEWS[0] else b_0 * 1e3
    # The sphere equator sits above B0 at strong screening, so its peak sets the top.
    top = max(cfg.PLOT["y_headroom"], 1.05 * max(1.0, float(np.max(sphere)))) * yf
    fig = go.Figure()
    for name, xs, ys, color in (
            ("Slab", s, slab, t["Finite difference"]),
            ("Cylinder", s_cyl, cyl, t["Shooting"]),
            ("Sphere (equator)", s_sphere, sphere, t["analytic"]),
            ("Finite body (centerline)", s_fb, fb, t["Relaxation (SOR)"])):
        fig.add_scatter(x=xs, y=ys * yf, name=name, mode="lines",
                        line={"color": color, "width": 2.4},
                        hovertemplate="%{y:.4f}<extra>" + name + "</extra>")
    fig.update_xaxes(title="s / λ from the nearest face")
    fig.update_yaxes(title="B / B₀" if view == cfg.GEOM_VIEWS[0] else "B (mT)", range=[0, top])
    return ui.style(fig, t, ui.CHART_H, uirevision=f"combined-{view}-{size}-{temp}")


def render():
    """Compare every implemented geometry at the common half-size."""
    t = ui.theme()
    ui.page_title(
        "Combined",
        "Slab, cylinder, sphere, finite body and film at one common half-size "
        "s/λ = d/(2λ) = R/λ, measured from the nearest face.")

    left, right = st.columns([0.34, 0.66], gap="medium")
    with left:
        material = ui.material_controls("combined", compact=True)
        lam, b_0 = material["lambda"], material["b_0"]
        normal = material["normal_state"] or not math.isfinite(lam)
        with ui.card("controls"):
            ui.card_title("Common size, grid and units")
            size_spec = ui.size_input("combined", cfg.COMBINED["size"], cfg.COMBINED["size_nm"],
                                      "Common half-size s/λ = d/(2λ) = R/λ",
                                      "Common half-size (nm)", lam, t)
            n_points = int(ui.live_slider("Grid points N (odd)", cfg.GRID["min"], cfg.GRID["max"],
                                          cfg.GRID["default"], cfg.GRID["step"], "combined_n", t))
            view = st.segmented_control("Units", ui.GEOM_VIEWS, default=ui.GEOM_VIEWS[0],
                                        key="combined_units") or ui.GEOM_VIEWS[0]

    s_over, size_m = size_spec["ratio"], size_spec["size_m"]

    if normal:
        # λ is infinite, so the size is read in units of λ₀; all curves are flat at B₀.
        span = s_over if size_spec["mode"] == "Ratio" else size_m / material["lambda_0"]
        s = np.linspace(0.0, span, n_points)
        s_cyl = np.linspace(0.0, span, n_points)
        slab = cyl = np.ones(n_points)
        v = (1.0, 1.0, 0.0, 1.0, 1.0, 0.0)
        section = visualization.combined_section(span, math.inf, b_0, None, t)
        profile_fig = ui.geometry_compare_chart(s, slab, s_cyl, cyl, "Both", view, b_0, t)
        for trace in profile_fig.data:  # both curves are the subject here, so draw both bold
            trace.line.width, trace.line.dash = 2.6, "solid"
        profile_fig.update_xaxes(title="s / (set size)")
    else:
        s, slab, s_cyl, cyl = ui.geometry_curves(s_over, lam, b_0, n_points)
        r = np.linspace(0.0, size_m, n_points)
        _, b_theta = sphere_field(r, math.pi / 2.0, size_m, lam, b_0)
        s_sphere = (size_m - r) / lam
        ratio_sphere = np.abs(b_theta) / b_0
        solution = finite_body_solution(s_over)
        mid = solution.x.size // 2 + 1
        s_fb = solution.x[:mid] - solution.x[0]  # distance from the edge, 0 to size
        ratio_fb = solution.b[solution.y.size // 2, :mid]
        # The solve at λ = 1 carries the same field; rescale x and y to metres for the section.
        section = visualization.combined_section(
            s_over, lam, b_0, replace(solution, x=solution.x * lam, y=solution.y * lam), t)
        v = ratios(s_over)
        profile_fig = combined_chart(s, slab, s_cyl, cyl, s_sphere, ratio_sphere,
                                     s_fb, ratio_fb, view, b_0, s_over, material["t"], t)

    with right:
        with ui.card("card_combined_map"):
            ui.card_title("Field map")
            ui.live_chart(section, key="combined_visual")
            st.caption("Five geometries at the common half-size s/λ. Each panel shows |B|/B₀ "
                       "on a fixed 0–1 scale; the individual pages carry the field lines.")
        k = [*st.columns(3, gap="small"), *st.columns(3, gap="small")]
        labels = ("Slab centre 1/cosh(s)", "Slab mean tanh(s)/s", "Slab screened",
                  "Cylinder centre 1/I₀(s)", "Cylinder mean (2/s)·I₁/I₀", "Cylinder screened")
        for col, label, value in zip(k, labels, v):
            shown = f"{value:.1%}" if "screened" in label else f"{value:.3f} B₀"
            col.metric(label, shown)

    c1, c2 = st.columns(2, gap="medium")
    with c1, ui.card("card_combined"):
        ui.card_title("Screening profile")
        ui.live_chart(profile_fig, key="combined_chart")
        if normal:
            st.caption("Normal state: λ → ∞, so the field is flat at B₀ in both geometries and "
                       "the axis is in units of the set size.")
        else:
            st.caption("The slab and cylinder curves stop at s = size: the slab centre and the "
                       "cylinder axis. Slab, cylinder and finite body start at B₀ at the surface; "
                       "the finite body centre is at s = size. The sphere starts at its equator "
                       "value, which is above B₀ at strong screening. The film parallel profile is "
                       "identical to the slab, so it is not drawn again.")
    with c2, ui.card("card_combined_bars"):
        ui.card_title("Centre and mean")
        if normal:
            st.caption("Normal state: every geometry is flat at B₀, so the centre and the mean "
                       "are both 1.000.")
        else:
            ui.live_chart(comparison_chart(comparison_rows(s_over), t), key="combined_bars")
            st.caption("Centre and mean B/B₀ at the same common half-size. The geometries do "
                       "not share the same boundary value or symmetry.")

    if normal:
        st.caption("Normal state: λ is infinite, so nothing is screened and no geometry is "
                   "solved; the comparison table is not shown.")
        ui.normal_state_infographic(t)
    else:
        with ui.card("card_combined_table"):
            ui.card_title("Same size across geometries")
            comparison_table(comparison_rows(s_over))
        st.caption(f"The finite body is a square of side 2s, solved with the fast "
                   f"N = {cfg.FINITE_BODY['n']['min']} scalar relaxation and Dirichlet edges. "
                   "The film's Pearl length and regime sit in its row.")

    st.caption("What changes with T: a rise in T raises λ and weakens screening, so the centre "
               "and mean fields of every geometry move closer to B₀. Set T and B₀ on the "
               "Material and field card in the left column.")