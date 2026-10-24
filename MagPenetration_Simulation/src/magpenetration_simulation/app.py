"""Streamlit shell and page router for the London-equation penetration models.

The header and the top navigation pick one page; only that page computes. The
sidebar holds the material and field settings that every page shares.

Run from MagPenetration_Simulation/ with:
    uv run streamlit run src/magpenetration_simulation/app.py
"""

import math

import streamlit as st

from magpenetration_simulation import ui, ui_config as cfg, views
from magpenetration_simulation.temperature import (
    clamp_temperature,
    lambda_of_temperature,
)

DEFAULT_PAGE = "theory"
VIEWS = {
    "theory": views.theory.render,
    "methods": views.methods.render,
    "slab": views.slab.render,
    "cylinder": views.cylinder.render,
    "combined": views.combined.render,
    "sphere": views.sphere.render,
}


def current_page():
    """The page from the query string, falling back to Theory."""
    page = st.query_params.get("page", DEFAULT_PAGE)
    if isinstance(page, list):
        page = page[0] if page else DEFAULT_PAGE
    return page if page in VIEWS else DEFAULT_PAGE


def goto(page):
    """Switch page without reloading, so the sidebar values survive."""
    if st.query_params.get("page") != page:
        st.query_params["page"] = page
    st.rerun()


def material_state():
    """The shared sidebar settings every page reads, in SI units (m, T, K)."""
    lambda_0 = st.session_state.g_lam0 * 1e-9
    t_c = float(st.session_state.g_tc)
    normal_state = bool(st.session_state.g_normal)
    return {
        "t": float(st.session_state.g_t),
        "t_c": t_c,
        "lambda_0": lambda_0,
        "b_0": st.session_state.g_b0 * 1e-3,
        "lambda": math.inf if normal_state else lambda_of_temperature(
            float(st.session_state.g_t), lambda_0, t_c),
        "normal_state": normal_state,
    }


st.set_page_config(page_title="London penetration model", layout="wide")
if "theme" not in st.session_state:
    st.session_state.theme = "light"

choice = st.segmented_control("Theme", [":material/light_mode:", ":material/dark_mode:"],
                              default=":material/light_mode:" if st.session_state.theme == "light"
                              else ":material/dark_mode:", label_visibility="collapsed", key="theme_pick")
if choice:
    st.session_state.theme = "light" if choice == ":material/light_mode:" else "dark"
t = ui.theme()
ui.inject_css(t)

page = current_page()
if ui.brand("London penetration model"):
    goto("theory")
clicked = ui.nav(page)
if clicked:
    goto(clicked)

with st.sidebar:
    st.markdown('<p class="card-title">Material and field</p>', unsafe_allow_html=True)
    t_c_cfg, t_cfg = cfg.SIDEBAR["t_c"], cfg.SIDEBAR["t"]
    st.number_input("Critical temperature T_c (K)", min_value=t_c_cfg["min"],
                    max_value=t_c_cfg["max"], value=t_c_cfg["default"],
                    step=t_c_cfg["step"], key="g_tc")
    clamp = float(clamp_temperature(st.session_state.g_tc, cfg.TEMPERATURE["clamp_factor"]))
    # A lower T_c leaves the stored T above the new clamp, so pull it down first.
    st.session_state.setdefault("g_t", min(t_cfg["default"], clamp))
    st.session_state.g_t = min(st.session_state.g_t, clamp)
    st.slider("Temperature T (K)", t_cfg["min"], clamp, step=t_cfg["step"], key="g_t")
    lam_cfg, b_cfg = cfg.SIDEBAR["lambda_0"], cfg.SIDEBAR["b_0"]
    st.slider("Penetration depth λ₀ (nm)", lam_cfg["min"], lam_cfg["max"],
              lam_cfg["default"], lam_cfg["step"], key="g_lam0")
    st.slider("Applied field B₀ (mT)", b_cfg["min"], b_cfg["max"],
              b_cfg["default"], b_cfg["step"], key="g_b0")
    st.toggle("Normal state (T > T_c)", value=False, key="g_normal",
              help="Screening vanishes: λ → ∞ and the field is flat at B₀.")

VIEWS[page](material_state())
