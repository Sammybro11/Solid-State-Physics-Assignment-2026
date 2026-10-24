"""Combined page: every implemented geometry at one common size."""

import math

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from magpenetration_simulation import ui
from magpenetration_simulation import ui_config as cfg
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
    """Centre, mean and screened fraction of B/B0 for both geometries.

    The cylinder helpers take R/lam = size at lambda = 1 and b0 = 1, so they
    return the field ratios directly.
    """
    slab_centre = 1.0 / math.cosh(size / 2.0)
    slab_mean = 2.0 / size * math.tanh(size / 2.0)
    cyl_centre = cylinder_center_field(size, 1.0, 1.0)
    cyl_mean = cylinder_mean_field(size, 1.0, 1.0)
    return (slab_centre, slab_mean, 1.0 - slab_centre,
            cyl_centre, cyl_mean, 1.0 - cyl_centre)


@st.cache_data(show_spinner=False)
def finite_body_solution(size):
    """Cached fast square scalar solve at lam = b0 = 1, shared by the table and the curve."""
    n = cfg.FINITE_BODY["n"]["min"]
    return rect_field_finite_difference(size, size, 1.0, 1.0, nx=n, ny=n)


def finite_body_ratios(size):
    """Centre and mean B/B0 from the cached square solve."""
    solution = finite_body_solution(size)
    return rect_center_field(solution), float(solution.b.mean())


def comparison_rows(size):
    """Centre, mean and note per geometry at the common size, in B/B0 units."""
    n = cfg.FINITE_BODY["n"]["min"]
    fb_centre, fb_mean = finite_body_ratios(size)
    return (
        ("Slab", 1.0 / math.cosh(size / 2.0), 2.0 / size * math.tanh(size / 2.0),
         "Flat face at B₀, so the boundary value is B₀."),
        ("Cylinder", cylinder_center_field(size, 1.0, 1.0),
         cylinder_mean_field(size, 1.0, 1.0), "Axial field, radial 1D solution."),
        ("Sphere", sphere_center_field(size, 1.0, 1.0), sphere_mean_field(size, 1.0, 1.0),
         f"Angular field; equator |B| = {sphere_equator_field(size, 1.0, 1.0):.3f} B₀."),
        ("Finite body", fb_centre, fb_mean,
         f"Square section, {n}×{n} scalar solve with Dirichlet edges."),
        ("Film", film_parallel_center_field(size, 1.0, 1.0),
         film_parallel_mean_field(size, 1.0, 1.0),
         f"Parallel field equals the slab; Pearl Λ/λ = {pearl_ratio(1.0, size):.2f}, "
         f"{film_regime(1.0, size)} regime."),
    )


def comparison_table(rows):
    """Render the per-geometry rows as a five-column table."""
    st.dataframe(
        [{"Geometry": name, "Center B/B₀": f"{centre:.3f}", "Mean B/B₀": f"{mean:.3f}",
          "Screened": f"{1.0 - centre:.1%}", "Note": note}
         for name, centre, mean, note in rows],
        hide_index=True, width="stretch")


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


def render(material):
    """Compare every implemented geometry at the common size."""
    t = ui.theme()
    lam, b_0 = material["lambda"], material["b_0"]
    normal = material["normal_state"] or not math.isfinite(lam)

    ui.page_title(
        "Combined",
        "Slab, cylinder, sphere, finite body and film at one common size; the profile "
        "below uses s/λ from the nearest face.")

    with ui.card("controls"):
        c = st.columns([1.6, 1.3, 1.1], gap="large")
        with c[0]:
            size_cfg = cfg.COMBINED["size"]
            size = ui.live_slider("Common size d/λ = R/λ", size_cfg["min"], size_cfg["max"],
                                  size_cfg["default"], size_cfg["step"], "combined_size", t,
                                  decimals=1)
        with c[1]:
            n_points = int(ui.live_slider("Grid points N (odd)", cfg.GRID["min"], cfg.GRID["max"],
                                          cfg.GRID["default"], cfg.GRID["step"], "combined_n", t))
        with c[2]:
            view = st.segmented_control("Units", ui.GEOM_VIEWS, default=ui.GEOM_VIEWS[0],
                                        key="combined_units") or ui.GEOM_VIEWS[0]

    if normal:
        s = np.linspace(0.0, size / 2.0, n_points)
        s_cyl = np.linspace(0.0, size, n_points)
        slab = cyl = np.ones(n_points)
    else:
        s, slab, s_cyl, cyl = ui.geometry_curves(size, lam, b_0, n_points)
        d = size * lam
        r = np.linspace(0.0, d, n_points)
        _, b_theta = sphere_field(r, np.pi / 2.0, d, lam, b_0)
        s_sphere = (d - r) / lam
        ratio_sphere = np.abs(b_theta) / b_0
        solution = finite_body_solution(size)
        mid = solution.x.size // 2 + 1
        s_fb = solution.x[:mid] - solution.x[0]
        ratio_fb = solution.b[solution.y.size // 2, :mid]
    with ui.card("card_combined"):
        ui.card_title("Screening profile")
        if normal:
            fig = ui.geometry_compare_chart(s, slab, s_cyl, cyl, "Both", view, b_0, t)
            for trace in fig.data:  # both curves are the subject here, so draw both bold
                trace.line.width, trace.line.dash = 2.6, "solid"
        else:
            fig = combined_chart(s, slab, s_cyl, cyl, s_sphere, ratio_sphere,
                                 s_fb, ratio_fb, view, b_0, size, material["t"], t)
        ui.live_chart(fig, key="combined_chart")
        if normal:
            st.caption("The slab curve stops at its centre s = d/2; the cylinder curve reaches "
                       "its axis s = R. Both curves start at B₀ at the surface.")
        else:
            st.caption("The slab curve stops at its centre s = d/2; the cylinder curve reaches "
                       "its axis s = R. Slab, cylinder, and finite body start at B₀ at the "
                       "surface; the finite body centre is at s = d/2. The sphere starts at "
                       "its equator value, which is above B₀ at strong screening. The film "
                       "parallel profile is identical to the slab; the table lists its Pearl "
                       "length and regime. The sphere and finite body have their own pages "
                       "for angular and 2D views.")

    if normal:
        v = (1.0, 1.0, 0.0, 1.0, 1.0, 0.0)
        st.caption("Normal state: λ is infinite, so the field is flat at B₀ in both "
                   "geometries and nothing is screened.")
        ui.normal_state_infographic(t)
    else:
        v = ratios(size)

    k = st.columns(6, gap="small")
    labels = ("Slab centre 1/cosh(size/2)", "Slab mean (2/size)·tanh(size/2)", "Slab screened",
              "Cylinder centre 1/I₀(size)", "Cylinder mean (2/size)·I₁/I₀", "Cylinder screened")
    for col, label, value in zip(k, labels, v):
        shown = f"{value:.1%}" if "screened" in label else f"{value:.3f} B₀"
        col.metric(label, shown)
    if normal:
        st.caption("Normal-state table is not shown; the London sphere and finite-body "
                   "models are Meissner-state slices.")
    else:
        with ui.card("card_combined_table"):
            ui.card_title("Same size across geometries")
            comparison_table(comparison_rows(size))
        st.caption("The profile chart compares every geometry on s/λ from the nearest face. "
                   "The table uses the same size, but the geometries do not share the same "
                   "boundary value or symmetry; the finite-body centre uses the same fast "
                   f"N = {cfg.FINITE_BODY['n']['min']} scalar solve with Dirichlet edges.")

    st.caption("What changes with T: a rise in T raises λ and weakens screening, so both "
               "centre fields and means move closer to B₀. The sidebar sets T.")
    st.caption("The 2D field-line comparison will land here.")
