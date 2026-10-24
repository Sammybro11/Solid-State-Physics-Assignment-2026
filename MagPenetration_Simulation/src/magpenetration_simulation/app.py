"""Streamlit shell and page router for the London-equation penetration models.

The header and the top navigation pick one page; only that page computes. Each
page owns its own material and field controls, so no sidebar is needed.

Run from MagPenetration_Simulation/ with:
    uv run streamlit run src/magpenetration_simulation/app.py
"""

import streamlit as st

from magpenetration_simulation import ui, views

DEFAULT_PAGE = "theory"
VIEWS = {
    "theory": views.theory.render,
    "methods": views.methods.render,
    "slab": views.slab.render,
    "cylinder": views.cylinder.render,
    "combined": views.combined.render,
    "sphere": views.sphere.render,
    "finite_body": views.finite_body.render,
    "film": views.film.render,
}


def current_page():
    """The page from the query string, falling back to Theory."""
    page = st.query_params.get("page", DEFAULT_PAGE)
    if isinstance(page, list):
        page = page[0] if page else DEFAULT_PAGE
    return page if page in VIEWS else DEFAULT_PAGE


def goto(page):
    """Switch page without reloading, so each page keeps its own values."""
    if st.query_params.get("page") != page:
        st.query_params["page"] = page
    st.rerun()


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

VIEWS[page]()
