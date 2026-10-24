"""Shared presentation helpers for the penetration model pages.

Themes, CSS, chart styling, the live chart and slider components, and the page
chrome (header, navigation, cards). Model maths stays in slab.py, cylinder.py
and temperature.py.
"""

import json

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from plotly.offline import get_plotlyjs
from streamlit.components.v2 import component

from magpenetration_simulation import ui_config as cfg
from magpenetration_simulation.cylinder import cylinder_field_profile
from magpenetration_simulation.slab import slab_field_analytic
from magpenetration_simulation.temperature import (
    clamp_temperature,
    lambda_curve,
    lambda_of_temperature,
    temperature_profile,
)
from magpenetration_simulation.ui_config import (
    CHART_H,
    FAMILY_T,
    GEOM_FAMILY,
    GEOM_VIEWS,
    PAGES,
    TEMP_POINTS,
)

SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
THEMES = {
    "light": {
        "bg": "#f7f8fa", "card": "#ffffff", "border": "#e6e8ec", "text": "#0f172a",
        "muted": "#64748b", "grid": "#eef1f5", "axis": "#cbd5e1", "slab": "rgba(59,130,246,0.06)",
        "analytic": "#0f172a", "hover": "#ffffff", "accent": "#2563eb",
        "Finite difference": "#2563eb", "Relaxation (SOR)": "#d97706", "Shooting": "#059669",
        "pill": "#f8fafc", "pill_on": "#ffffff", "pill_on_border": "#94a3b8",
    },
    "dark": {
        "bg": "#0b0f17", "card": "#121826", "border": "#1f2937", "text": "#e5e7eb",
        "muted": "#94a3b8", "grid": "#1c2433", "axis": "#334155", "slab": "rgba(96,165,250,0.07)",
        "analytic": "#f1f5f9", "hover": "#1e293b", "accent": "#60a5fa",
        "Finite difference": "#60a5fa", "Relaxation (SOR)": "#fbbf24", "Shooting": "#34d399",
        "pill": "#0f1520", "pill_on": "#1e293b", "pill_on_border": "#475569",
    },
}


def theme():
    """The active theme dict, from st.session_state.theme."""
    return THEMES[st.session_state.get("theme", "light")]


# ---------------------------------------------------------------- number formatting

def sci(value, digits=2):
    """Format a number as 7.85×10¹¹."""
    mantissa, exponent = f"{value:.{digits}e}".split("e")
    return f"{mantissa}×10{str(int(exponent)).translate(SUPERSCRIPT)}"


def field_text(tesla):
    """Format a magnetic field with a readable SI prefix, e.g. 0.141 µT or 7.4 pT."""
    for scale, unit in ((1, "T"), (1e-3, "mT"), (1e-6, "µT"), (1e-9, "nT"), (1e-12, "pT")):
        if abs(tesla) >= scale:
            return f"{tesla / scale:.3g} {unit}"
    return f"{tesla:.1e} T"


# ---------------------------------------------------------------- charts

def alpha(hex_color, a):
    """The given hex colour at opacity a, as an rgba() string."""
    h = hex_color.lstrip("#")
    return f"rgba({int(h[0:2], 16)},{int(h[2:4], 16)},{int(h[4:6], 16)},{a})"


def style(fig, t, height, **layout):
    """Shared look for every chart."""
    axis = {"automargin": True, "gridcolor": t["grid"], "zeroline": False, "showline": True, "linecolor": t["axis"],
            "ticks": "outside", "tickcolor": t["axis"], "ticklen": 4, "title_standoff": 8,
            "tickfont": {"size": 12, "color": t["muted"]}, "title_font": {"size": 13, "color": t["muted"]}}
    fig.update_xaxes(**axis)
    fig.update_yaxes(**axis)
    fig.update_layout(template=go.layout.Template(), height=height, autosize=True,
                      margin={"l": 4, "r": 8, "t": 6, "b": 4},
                      paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                      font={"family": "Source Sans Pro, sans-serif", "color": t["text"]},
                      hovermode="x unified",
                      hoverlabel={"bgcolor": t["hover"], "bordercolor": t["border"],
                                  "font": {"color": t["text"], "size": 12}},
                      showlegend=False, **layout)
    return fig


def lambda_chart(t_c, lambda_0_m, t_k, t):
    """Penetration depth against reduced temperature, with the current T marked."""
    factor = cfg.TEMPERATURE["clamp_factor"]
    ratio = np.linspace(0.0, clamp_temperature(t_c, factor) / t_c, TEMP_POINTS)
    lam = lambda_curve(ratio * t_c, lambda_0_m, t_c) / lambda_0_m
    fig = go.Figure()
    fig.add_scatter(x=ratio, y=lam, mode="lines", name="λ(T)",
                    line={"color": t["analytic"], "width": 2.2}, hovertemplate="%{y:.3f}")
    fig.add_scatter(x=[t_k / t_c], y=[lambda_of_temperature(t_k, lambda_0_m, t_c) / lambda_0_m],
                    mode="markers", hoverinfo="skip",
                    marker={"size": 11, "color": t["card"],
                            "line": {"color": t["Finite difference"], "width": 2}})
    fig.add_vline(x=factor, line={"color": t["axis"], "dash": "dot", "width": 1})
    fig.update_xaxes(title="T / T_c")
    fig.update_yaxes(title="λ(T) / λ₀")
    return style(fig, t, CHART_H, uirevision="temp-lambda")


def temperature_profile_chart(t_c, lambda_0_m, d_m, b_0, t_k, units, t):
    """Field profiles at the fixed T/Tc family, faint, plus the current T in bold."""
    x = np.linspace(-d_m / 2.0, d_m / 2.0, TEMP_POINTS)
    yf = b_0 * 1e3 if units == "B (mT)" else 1.0
    fig = go.Figure()
    last = len(FAMILY_T) - 1
    for i, f in enumerate(FAMILY_T):
        y = temperature_profile(x, d_m, lambda_0_m, t_c, f * t_c, b_0) * yf
        fig.add_scatter(x=x * 1e9, y=y, mode="lines", hoverinfo="skip",
                        line={"color": alpha(t["Finite difference"], 0.22 + 0.5 * i / last),
                              "width": 1.4})
    y = temperature_profile(x, d_m, lambda_0_m, t_c, t_k, b_0) * yf
    fig.add_scatter(x=x * 1e9, y=y, mode="lines", name="T",
                    line={"color": t["analytic"], "width": 3}, hovertemplate="%{y:.4f}")
    fig.update_xaxes(title="x (nm)", range=[-d_m * 1e9 / 2.0, d_m * 1e9 / 2.0])
    fig.update_yaxes(title="B (mT)" if units == "B (mT)" else "B / B₀", range=[0, cfg.PLOT["y_headroom"] * yf])
    return style(fig, t, CHART_H, uirevision=f"temp-family-{units}-{d_m:.6e}")


def geometry_curves(size, lam, b0, n_points):
    """Slab and cylinder fields at the same size, d = R = size * lam.

    Returns s/lam for each curve, both measured from the nearest face, with the
    fields as B/B0. The slab curve stops at its centre, s = size/2; the
    cylinder curve reaches its axis, s = size.
    """
    d = size * lam
    s = np.linspace(0.0, size / 2.0, n_points)
    slab = slab_field_analytic(d / 2.0 - s * lam, d, lam, b0) / b0
    r, cyl = cylinder_field_profile(d, lam, b0, n_points)
    return s, slab, size - r[::-1] / lam, cyl[::-1] / b0


def geometry_compare_chart(s, slab, s_cyl, cyl, geometry, view, b0, t):
    """Slab and cylinder B/B0 against depth from the nearest face, selected geometry in bold."""
    yf = 1.0 if view == GEOM_VIEWS[0] else b0 * 1e3
    fig = go.Figure()
    for name, xs, ys, color in (("Slab", s, slab, t["Finite difference"]),
                                ("Cylinder", s_cyl, cyl, t["Shooting"])):
        bold = name == geometry
        fig.add_scatter(x=xs, y=ys * yf, name=name, mode="lines",
                        line={"color": color, "width": 2.6 if bold else 1.6,
                              "dash": "solid" if bold else "dot"},
                        hovertemplate="%{y:.4f}")
    fig.update_xaxes(title="s / λ from the nearest face")
    fig.update_yaxes(title="B / B₀" if view == GEOM_VIEWS[0] else "B (mT)", range=[0, cfg.PLOT["y_headroom"] * yf])
    return style(fig, t, CHART_H, uirevision=f"geom-compare-{view}-{s[-1]:.6e}")


def cylinder_family_chart(lam, b0, size, view, t):
    """Cylinder B/B0 against r/R for the R/lambda family, current size in bold."""
    yf = 1.0 if view == GEOM_VIEWS[0] else b0 * 1e3
    last = len(GEOM_FAMILY) - 1
    fig = go.Figure()
    for i, f in enumerate(GEOM_FAMILY):
        r, b = cylinder_field_profile(f * lam, lam, b0, TEMP_POINTS)
        current = f == size
        fig.add_scatter(
            x=r / (f * lam), y=b / b0 * yf, name=f"R/λ = {f:g}", mode="lines",
            line={"color": t["Shooting"] if current else
                  alpha(t["Shooting"], 0.2 + 0.5 * i / last),
                  "width": 2.6 if current else 1.3},
            hovertemplate="%{y:.4f}")
    fig.update_xaxes(title="r / R")
    fig.update_yaxes(title="B / B₀" if view == GEOM_VIEWS[0] else "B (mT)", range=[0, cfg.PLOT["y_headroom"] * yf])
    return style(fig, t, CHART_H, uirevision=f"geom-family-{view}")


def normal_state_infographic(t):
    """Explain, in three steps, why B is flat above Tc. Colours come from inject_css."""
    steps = (
        ("Superfluid density", "n<sub>s</sub> → 0", "The superconducting carriers vanish."),
        ("London depth", "λ → ∞", "The screening length exceeds the slab."),
        ("Field equation", "∇²B = B/λ² → 0",
         "With B(±d/2) = B₀, the field is flat at B₀. The screening current J is zero."),
    )
    cards = "".join(
        f'<div class="normal-card"><div class="normal-step">Step {i + 1}</div>'
        f'<div class="normal-title">{title}</div><div class="normal-eq">{eq}</div>'
        f'<div class="normal-note">{note}</div></div>'
        for i, (title, eq, note) in enumerate(steps)
    )
    st.markdown(f'<div class="normal-grid">{cards}</div>', unsafe_allow_html=True)


# ---------------------------------------------------------------- custom components
# st.slider only reports on release and st.plotly_chart redraws from scratch on
# every rerun. These two small components report slider values while dragging
# and update the charts in place with Plotly.react, so the page follows the drag.

LIVE_SLIDER_HTML = """<div class="live">
  <div class="row"><span class="label"></span><span class="value"></span></div>
  <input type="range">
</div>"""

LIVE_SLIDER_CSS = """
.live {font-family: "Source Sans Pro", sans-serif; padding: 2px 2px 6px 2px;}
.row {display: flex; justify-content: space-between; align-items: baseline; margin-bottom: 6px;}
.label {font-size: 0.85rem; color: var(--text);}
.value {font-size: 0.85rem; font-weight: 600; color: var(--accent); font-variant-numeric: tabular-nums;}
input[type=range] {-webkit-appearance: none; appearance: none; width: 100%; height: 4px; border-radius: 4px;
  background: linear-gradient(to right, var(--accent) var(--fill), var(--track) var(--fill)); outline: none;}
input[type=range]::-webkit-slider-thumb {-webkit-appearance: none; width: 16px; height: 16px; border-radius: 50%;
  background: var(--accent); border: 3px solid var(--card); box-shadow: 0 0 0 1px var(--accent); cursor: grab;}
input[type=range]::-moz-range-thumb {width: 12px; height: 12px; border-radius: 50%; background: var(--accent);
  border: 3px solid var(--card); cursor: grab;}
"""

LIVE_SLIDER_JS = """
export default function (component) {
  const { data, parentElement, setStateValue } = component;
  const root = parentElement.querySelector(".live");
  const input = root.querySelector("input");
  const valueEl = root.querySelector(".value");
  for (const [k, v] of Object.entries(data.colors)) root.style.setProperty(`--${k}`, v);
  root.querySelector(".label").textContent = data.label;
  if (!input.dataset.ready) {
    Object.assign(input, {min: data.min, max: data.max, step: data.step});
    input.value = data.value;
    input.dataset.ready = "1";
  }
  const show = () => {
    const v = Number(input.value);
    valueEl.textContent = v.toFixed(data.decimals) + (data.unit ? " " + data.unit : "");
    root.style.setProperty("--fill", (100 * (v - data.min)) / (data.max - data.min) + "%");
  };
  let last = 0, pending = null;
  const send = () => { pending = null; last = performance.now(); setStateValue("value", Number(input.value)); };
  input.oninput = () => {
    show();
    const wait = data.throttle_ms - (performance.now() - last);
    if (wait <= 0) send();
    else if (!pending) pending = setTimeout(send, wait);
  };
  input.onchange = () => { clearTimeout(pending); send(); };
  show();
}
"""

LIVE_CHART_HTML = '<div class="live-chart"></div>'

LIVE_CHART_JS = """
export default function (component) {
  const { data, parentElement } = component;
  const root = parentElement.querySelector(".live-chart");
  root.style.height = data.figure.layout.height + "px";
  const draw = () => window.Plotly.react(root, data.figure.data, data.figure.layout, data.config);
  if (window.Plotly) { draw(); return; }
  if (data.plotly_js && !window.__plotlyLoading) {
    window.__plotlyLoading = true;
    const script = document.createElement("script");
    script.textContent = data.plotly_js;
    document.head.appendChild(script);
  }
  const timer = setInterval(() => { if (window.Plotly) { clearInterval(timer); draw(); } }, 30);
  return () => clearInterval(timer);
}
"""

_live_slider = component(
    "live_slider", html=LIVE_SLIDER_HTML, css=LIVE_SLIDER_CSS, js=LIVE_SLIDER_JS)
_live_chart = component(
    "live_chart", html=LIVE_CHART_HTML, js=LIVE_CHART_JS, isolate_styles=False)


def live_slider(label, min_value, max_value, value, step, key, t, decimals=0, unit=""):
    """A slider that reports its value while it is being dragged, not only on release."""
    colors = {"accent": t["Finite difference"], "track": t["border"], "text": t["text"], "card": t["card"]}
    result = _live_slider(
        key=key,
        data={"label": label, "min": min_value, "max": max_value, "step": step, "value": value,
              "decimals": decimals, "unit": unit, "throttle_ms": 60, "colors": colors},
        default={"value": value},
        on_value_change=lambda: None,
    )
    return result.value


def live_chart(fig, key):
    """Draw a Plotly figure that updates in place on every rerun, without flashing.

    Plotly.js is sent to the browser once per session, with the first chart.
    """
    data = {"figure": json.loads(fig.to_json()),
            "config": {"displayModeBar": True, "modeBarButtons": [["resetScale2d"]],
                       "displaylogo": False, "responsive": True}}
    if not st.session_state.get("plotly_js_sent"):
        data["plotly_js"] = get_plotlyjs()
        st.session_state.plotly_js_sent = True
    _live_chart(key=key, data=data, height="content")


# ---------------------------------------------------------------- page chrome

def brand(title):
    """The clickable header. Returns True when the user clicks it."""
    return st.button(title, key="brand_home", type="tertiary", help="Back to Theory")


def nav(active):
    """Top navigation. Returns the page key the user clicked, or None."""
    clicked = None
    with st.container(key="nav"):
        cols = st.columns(len(PAGES), gap="small")
        for col, (key, label) in zip(cols, PAGES):
            if col.button(label, key=f"nav_{key}", width="stretch",
                          type="primary" if key == active else "tertiary"):
                clicked = key
    return clicked


def page_title(title, note=None):
    """The heading of a page, with an optional one-line note under it."""
    st.markdown(f'<div class="page-title">{title}</div>', unsafe_allow_html=True)
    if note:
        st.caption(note)


def card(key):
    """A bordered card, styled by inject_css through its key."""
    return st.container(border=True, key=key)


def card_title(text, where=None):
    """A card heading; pass where to place it inside a column."""
    (where or st).markdown(f'<p class="card-title">{text}</p>', unsafe_allow_html=True)


# ---------------------------------------------------------------- page style

def inject_css(t):
    """Style the page for the chosen theme."""
    dots = "".join(
        f".st-key-pills_profile button:nth-child({i + 1}) p::before, "
        f".st-key-pills_error button:nth-child({i}) p::before {{color: {c};}}"
        for i, c in enumerate([t["analytic"], t["Finite difference"], t["Relaxation (SOR)"], t["Shooting"]])
    )
    st.markdown(f"""<style>
    [data-testid="stHeader"] {{display: none;}}
    .stApp [data-stale="true"] {{opacity: 1 !important; transition: none !important;}}
    .st-key-theme_pick {{position: fixed; top: 0.7rem; right: 1rem; z-index: 999990; width: auto;}}
    .stApp {{background: {t['bg']};}}
    .block-container {{padding: 0.9rem 1.8rem 0.4rem 1.8rem; max-width: 1600px;}}
    [data-testid="stSidebar"] {{background: {t['card']}; border-right: 1px solid {t['border']};}}
    .stApp [data-testid="stWidgetLabel"] p {{color: {t['text']}; font-size: 0.85rem;}}
    .st-key-brand_home button {{background: transparent; border: none; box-shadow: none;
        padding: 0 0 0.1rem 0; width: 100%; justify-content: flex-start; text-align: left;}}
    .st-key-brand_home button p {{font-size: 1.9rem !important; font-weight: 700; letter-spacing: -0.02em;
        color: {t['text']}; line-height: 1.2;}}
    .st-key-brand_home button:hover p {{color: {t['accent']};}}
    .st-key-nav {{border-bottom: 1px solid {t['border']}; margin-bottom: 0.9rem;}}
    .st-key-nav button {{background: transparent; border: none; border-bottom: 2px solid transparent;
        border-radius: 0; box-shadow: none; padding: 0.4rem 0.2rem;}}
    .st-key-nav button p {{font-size: 0.92rem; color: {t['muted']}; font-weight: 500;}}
    .st-key-nav button:hover {{background: {t['pill']};}}
    .st-key-nav button:hover p {{color: {t['text']};}}
    .st-key-nav [data-testid="stBaseButton-primary"] {{border-bottom-color: {t['accent']};}}
    .st-key-nav [data-testid="stBaseButton-primary"] p {{color: {t['text']}; font-weight: 650;}}
    .stApp .page-title {{font-size: 1.35rem; font-weight: 700; color: {t['text']};
        letter-spacing: -0.01em; margin: 0.1rem 0 0.1rem 0;}}
    .stApp .card-title {{font-weight: 650; font-size: 1rem !important; color: {t['text']}; margin: 0;}}
    [class*="st-key-card_"], .st-key-controls {{background: {t['card']}; border-color: {t['border']} !important;}}
    .st-key-controls {{padding-top: 0.6rem; padding-bottom: 0.2rem;}}
    [data-testid="stMetric"] {{background: {t['card']}; border: 1px solid {t['border']};
        border-radius: 0.55rem; padding: 0.4rem 0.8rem;}}
    [data-testid="stMetricLabel"] p {{color: {t['muted']} !important; font-size: 0.78rem !important;}}
    [data-testid="stMetricValue"] {{color: {t['text']}; font-size: 1.15rem !important; font-weight: 600;}}
    [data-testid="stMetricDelta"] {{background: {t['pill']}; color: {t['muted']};}}
    [data-testid="stMetricDelta"] div {{font-size: 0.75rem;}}
    .stApp button[data-variant] {{background: {t['pill']}; border-color: {t['border']};}}
    .stApp button[data-variant] p {{color: {t['muted']}; font-size: 0.85rem;}}
    .stApp button[data-variant][aria-pressed="true"] {{background: {t['pill_on']};
        border-color: {t['pill_on_border']};}}
    .stApp button[data-variant][aria-pressed="true"] p {{color: {t['text']}; font-weight: 600;}}
    .st-key-pills_profile button p::before, .st-key-pills_error button p::before {{
        content: "●"; margin-right: 0.35rem;}}
    .stApp button[data-variant][aria-pressed="false"] p::before {{opacity: 0.35;}}
    {dots}
    .normal-grid {{display: grid; grid-template-columns: repeat(3, 1fr); gap: 0.7rem; margin-top: 0.3rem;}}
    .normal-card {{background: {t['card']}; border: 1px solid {t['border']}; border-radius: 0.55rem;
        padding: 0.6rem 0.8rem;}}
    .normal-step {{font-size: 0.66rem; text-transform: uppercase; letter-spacing: 0.06em;
        color: {t['muted']}; font-weight: 600;}}
    .normal-title {{font-size: 0.9rem; font-weight: 650; color: {t['text']}; margin: 0.15rem 0 0.35rem 0;}}
    .normal-eq {{font-size: 1.05rem; color: {t['Finite difference']}; margin-bottom: 0.35rem;}}
    .normal-note {{font-size: 0.78rem; color: {t['muted']}; line-height: 1.35;}}
    @media (max-width: 900px) {{ .normal-grid {{grid-template-columns: 1fr;}} }}
    </style>""", unsafe_allow_html=True)
