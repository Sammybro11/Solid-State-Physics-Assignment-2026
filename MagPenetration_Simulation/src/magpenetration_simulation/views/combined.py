"""Combined page: slab and axial cylinder on one depth axis."""

import math

import numpy as np
import streamlit as st

from magpenetration_simulation import ui
from magpenetration_simulation import ui_config as cfg
from magpenetration_simulation.cylinder import cylinder_center_field, cylinder_mean_field


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


def render(material):
    """Compare the slab and the axial cylinder at the same size."""
    t = ui.theme()
    lam, b_0 = material["lambda"], material["b_0"]
    normal = material["normal_state"] or not math.isfinite(lam)

    ui.page_title(
        "Combined",
        "Slab and axial cylinder on one depth axis, s/λ measured from the nearest face.")

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

    with ui.card("card_combined"):
        ui.card_title("Screening profile")
        fig = ui.geometry_compare_chart(s, slab, s_cyl, cyl, "Both", view, b_0, t)
        for trace in fig.data:  # both curves are the subject here, so draw both bold
            trace.line.width, trace.line.dash = 2.6, "solid"
        ui.live_chart(fig, key="combined_chart")
        st.caption("The slab curve stops at its centre s = d/2; the cylinder curve reaches "
                   "its axis s = R. Both curves start at B₀ at the surface.")

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

    st.caption("What changes with T: a rise in T raises λ and weakens screening, so both "
               "centre fields and means move closer to B₀. The sidebar sets T.")
    st.caption("The 2D field-line comparison will land here.")
