"""Streamlit page for the 1D slab London-equation model.

Run from MagPenetration_Simulation/ with:
    uv run streamlit run src/magpenetration_simulation/app.py
"""

import json

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from plotly.offline import get_plotlyjs
from streamlit.components.v2 import component

from magpenetration_simulation.slab import (
    _grid,
    slab_field_analytic,
    slab_field_finite_difference,
    slab_field_relaxation,
    slab_field_shooting,
)
from magpenetration_simulation.temperature import (
    clamp_temperature,
    guard_metrics,
    lambda_curve,
    lambda_of_temperature,
    temperature_profile,
)

METHODS = {
    "Finite difference": slab_field_finite_difference,
    "Relaxation (SOR)": slab_field_relaxation,
    "Shooting": slab_field_shooting,
}
SHORT = {"Finite difference": "FD", "Relaxation (SOR)": "SOR", "Shooting": "Shooting"}
LABEL = {"Analytic": "Analytic", **SHORT}
TARGET = 1e-3
CHART_H = 480
MU0 = 1.25663706212e-6  # vacuum permeability in H/m
SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

FAMILY_T = (0.0, 0.217, 0.435, 0.652, 0.870, 0.957)
TEMP_POINTS = 401
THEMES = {
    "light": {
        "bg": "#f7f8fa", "card": "#ffffff", "border": "#e6e8ec", "text": "#0f172a",
        "muted": "#64748b", "grid": "#eef1f5", "axis": "#cbd5e1", "slab": "rgba(59,130,246,0.06)",
        "analytic": "#0f172a", "hover": "#ffffff",
        "Finite difference": "#2563eb", "Relaxation (SOR)": "#d97706", "Shooting": "#059669",
        "pill": "#f8fafc", "pill_on": "#ffffff", "pill_on_border": "#94a3b8",
    },
    "dark": {
        "bg": "#0b0f17", "card": "#121826", "border": "#1f2937", "text": "#e5e7eb",
        "muted": "#94a3b8", "grid": "#1c2433", "axis": "#334155", "slab": "rgba(96,165,250,0.07)",
        "analytic": "#f1f5f9", "hover": "#1e293b",
        "Finite difference": "#60a5fa", "Relaxation (SOR)": "#fbbf24", "Shooting": "#34d399",
        "pill": "#0f1520", "pill_on": "#1e293b", "pill_on_border": "#475569",
    },
}


# ---------------------------------------------------------------- model calls

@st.cache_data
def solve(d_over_lam, b0_mT, lam0_nm, n_points):
    """Run all four methods on the shared grid.

    Returns positions in units of lam and fields and errors in units of b0,
    so the charts do not depend on the absolute scales.
    """
    lam, b0 = lam0_nm * 1e-9, b0_mT * 1e-3
    d = d_over_lam * lam
    x = _grid(d, n_points)
    exact = slab_field_analytic(x, d, lam, b0)
    res = {}
    for name, method in METHODS.items():
        try:
            _, b = method(d, lam, b0, n_points)
            res[name] = {"b": b / b0, "signed": (b - exact) / b0,
                         "rel": np.abs(b - exact).max() / b0, "abs": np.abs(b - exact).max()}
        except RuntimeError as err:
            res[name] = {"error": str(err)}
    return {"xt": x / lam, "dt": d_over_lam, "n": n_points, "exact": exact / b0, "m": res,
            "lam": lam, "b0": b0, "d": d}


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

def _alpha(hex_color, a):
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


def x_axis(r, units):
    """Positions and axis title: x/lam when scaled, x in nm when physical."""
    factor = 1.0 if units == "Scaled" else r["lam"] * 1e9
    return r["xt"] * factor, factor, ("x / λ" if units == "Scaled" else "x (nm)")


def slab_backdrop(fig, r, t, xf):
    """Shade the slab and mark one penetration depth in from each face.

    xf converts lengths in units of lam to the chart's x units.
    """
    half = r["dt"] / 2
    fig.add_vrect(x0=-half * xf, x1=half * xf, fillcolor=t["slab"], line_width=0, layer="below")
    if 2 <= r["dt"] <= 40:
        for s in (-1, 1):
            fig.add_vline(x=s * (half - 1) * xf, line={"color": t["axis"], "dash": "dot", "width": 1})
        fig.add_annotation(x=-(half - 1) * xf, y=1, yref="y domain", text="λ", showarrow=False,
                           xshift=-7, yshift=-2, font={"size": 12, "color": t["muted"]})


def profile_chart(r, shown, units, t):
    """Field profile, as B/b0 against x/lam (scaled) or B in mT against x in nm (physical).

    Markers are staggered so methods do not stack; an invisible full-resolution
    line per method drives the hover readout.
    """
    x, xf, x_title = x_axis(r, units)
    yf = 1.0 if units == "Scaled" else r["b0"] * 1e3
    fmt = "%{y:.6f}" if units == "Scaled" else "%{y:.4f} mT"
    half = r["dt"] / 2
    pad = max(0.6, 0.08 * r["dt"])
    fig = go.Figure()
    slab_backdrop(fig, r, t, xf)
    for xs in ([-(half + pad) * xf, -half * xf], [half * xf, (half + pad) * xf]):
        fig.add_scatter(x=xs, y=[yf, yf], mode="lines", hoverinfo="skip",
                        line={"color": t["analytic"], "width": 1.5, "dash": "dot"})
    if "Analytic" in shown:
        fig.add_scatter(x=x, y=r["exact"] * yf, name="Analytic", mode="lines",
                        line={"color": t["analytic"], "width": 2.2}, hovertemplate=fmt)
    methods = [name for name in shown if name in METHODS]
    step = max(4, r["n"] // 20)
    for k, name in enumerate(methods):
        m = r["m"][name]
        if "error" in m:
            continue
        offset = (k - (len(methods) - 1) / 2) * step / max(len(methods), 1)
        s = slice(int(step / 2 + offset), None, step)
        fig.add_scatter(x=x[s], y=m["b"][s] * yf, name=SHORT[name], mode="markers",
                        marker={"size": 8, "color": t["card"], "line": {"color": t[name], "width": 2}},
                        hoverinfo="skip")
        fig.add_scatter(x=x, y=m["b"] * yf, name=SHORT[name], mode="lines",
                        line={"color": t[name], "width": 0}, hovertemplate=fmt)
    fig.update_xaxes(title=x_title, range=[-(half + pad) * xf, (half + pad) * xf])
    fig.update_yaxes(title="B / B₀" if units == "Scaled" else "B (mT)", range=[0, 1.05 * yf])
    return style(fig, t, CHART_H, uirevision=f"profile-{units}")


def error_chart(r, shown, scale, units, t):
    """Pointwise error of each chosen method, signed (linear) or absolute (log),
    in units of b0 (scaled) or in µT (physical)."""
    x, xf, x_title = x_axis(r, units)
    yf = 1.0 if units == "Scaled" else r["b0"] * 1e6
    fig = go.Figure()
    slab_backdrop(fig, r, t, xf)
    for name in shown:
        m = r["m"][name]
        if "error" in m:
            continue
        y = m["signed"] * yf if scale == "Signed" else np.abs(m["signed"]) * yf
        if scale == "Log":
            y = np.where(y > 0, y, np.nan)
        fill = "tozeroy" if scale == "Signed" and len(shown) == 1 else None
        fig.add_scatter(x=x, y=y, name=SHORT[name], mode="lines", fill=fill,
                        fillcolor=_alpha(t[name], 0.12), line={"color": t[name], "width": 2},
                        hovertemplate="%{y:.2e}" + ("" if units == "Scaled" else " µT"))
    if scale == "Log":
        top = np.log10(yf)
        fig.add_hline(y=TARGET * yf, line={"color": "#ef4444", "dash": "dash", "width": 1})
        fig.update_yaxes(type="log", range=[top - 16, top], dtick=2, exponentformat="power",
                         title="|B − B_analytic| / B₀" if units == "Scaled" else "|B − B_analytic| (µT)")
    else:
        fig.add_hline(y=0, line={"color": t["axis"], "width": 1})
        fig.update_yaxes(type="linear", autorange=True, range=None, dtick=None,
                         exponentformat="power",
                         title="(B − B_analytic) / B₀" if units == "Scaled" else "B − B_analytic (µT)")
    fig.update_xaxes(title=x_title)
    return style(fig, t, CHART_H, uirevision=f"error-{units}-{scale}")


def lambda_chart(t_c, lambda_0_m, t_k, theme):
    """Penetration depth against reduced temperature, with the current T marked."""
    ratio = np.linspace(0.0, clamp_temperature(t_c) / t_c, TEMP_POINTS)
    lam = lambda_curve(ratio * t_c, lambda_0_m, t_c) / lambda_0_m
    fig = go.Figure()
    fig.add_scatter(x=ratio, y=lam, mode="lines", name="λ(T)",
                    line={"color": theme["analytic"], "width": 2.2}, hovertemplate="%{y:.3f}")
    fig.add_scatter(x=[t_k / t_c], y=[lambda_of_temperature(t_k, lambda_0_m, t_c) / lambda_0_m],
                    mode="markers", hoverinfo="skip",
                    marker={"size": 11, "color": theme["card"],
                            "line": {"color": theme["Finite difference"], "width": 2}})
    fig.add_vline(x=0.995, line={"color": theme["axis"], "dash": "dot", "width": 1})
    fig.update_xaxes(title="T / T_c")
    fig.update_yaxes(title="λ(T) / λ₀")
    return style(fig, theme, CHART_H, uirevision="temp-lambda")


def temperature_profile_chart(t_c, lambda_0_m, d_m, b_0, t_k, units, theme):
    """Field profiles at the fixed T/Tc family, faint, plus the current T in bold."""
    x = np.linspace(-d_m / 2.0, d_m / 2.0, TEMP_POINTS)
    yf = b_0 * 1e3 if units == "B (mT)" else 1.0
    fig = go.Figure()
    last = len(FAMILY_T) - 1
    for i, f in enumerate(FAMILY_T):
        y = temperature_profile(x, d_m, lambda_0_m, t_c, f * t_c, b_0) * yf
        fig.add_scatter(x=x * 1e9, y=y, mode="lines", hoverinfo="skip",
                        line={"color": _alpha(theme["Finite difference"], 0.22 + 0.5 * i / last),
                              "width": 1.4})
    y = temperature_profile(x, d_m, lambda_0_m, t_c, t_k, b_0) * yf
    fig.add_scatter(x=x * 1e9, y=y, mode="lines", name="T",
                    line={"color": theme["analytic"], "width": 3}, hovertemplate="%{y:.4f}")
    fig.update_xaxes(title="x (nm)", range=[-d_m * 1e9 / 2.0, d_m * 1e9 / 2.0])
    fig.update_yaxes(title="B (mT)" if units == "B (mT)" else "B / B₀", range=[0, 1.05 * yf])
    return style(fig, theme, CHART_H, uirevision=f"temp-family-{units}")


def normal_state_infographic(theme):
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
    .stApp [data-testid="stWidgetLabel"] p {{color: {t['text']}; font-size: 0.85rem;}}
    .stApp .app-title {{text-align: center; font-size: 2rem !important; font-weight: 700; color: {t['text']};
        letter-spacing: -0.02em; line-height: 1.2; margin: 0.1rem 0 0.5rem 0;}}
    .stApp .card-title {{font-weight: 650; font-size: 1rem !important; color: {t['text']}; margin: 0;}}
    .st-key-controls, .st-key-card_profile, .st-key-card_error, .st-key-card_normal, .st-key-card_temp_lambda, .st-key-card_temp_profile {{background: {t['card']};
        border-color: {t['border']} !important;}}
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


# ---------------------------------------------------------------- page

st.set_page_config(page_title="London slab model", layout="wide")
if "theme" not in st.session_state:
    st.session_state.theme = "light"

choice = st.segmented_control("Theme", [":material/light_mode:", ":material/dark_mode:"],
                              default=":material/light_mode:" if st.session_state.theme == "light"
                              else ":material/dark_mode:", label_visibility="collapsed", key="theme_pick")
if choice:
    st.session_state.theme = "light" if choice == ":material/light_mode:" else "dark"
t = THEMES[st.session_state.theme]
inject_css(t)
st.markdown('<div class="app-title">Field penetration into a superconducting slab</div>',
            unsafe_allow_html=True)
slab_tab, temp_tab = st.tabs(["Slab", "Temperature"])

with slab_tab:
    with st.container(border=True, key="controls"):
        c = st.columns(4, gap="large")
        with c[0]:
            d_over_lam = live_slider("Slab thickness d/λ", 0.5, 100.0, 5.0, 0.5, "d_over_lam", t, decimals=1)
        with c[1]:
            b0_mT = live_slider("Applied field B₀", 10, 500, 50, 10, "b0", t, unit="mT")
        with c[2]:
            lam0_nm = live_slider("Penetration depth λ₀", 10, 200, 50, 5, "lam0", t, unit="nm")
        with c[3]:
            n_points = int(live_slider("Grid points N (odd)", 41, 2001, 401, 20, "n_points", t))

    r = solve(d_over_lam, b0_mT, lam0_nm, n_points)

    centre = r["exact"][n_points // 2]
    per_lam = (n_points - 1) / d_over_lam
    dx_nm = r["d"] / (n_points - 1) * 1e9
    surface_current = r["b0"] / (MU0 * r["lam"]) * np.tanh(d_over_lam / 2)

    k = st.columns(7)
    k[0].metric("Centre field B(0)", field_text(centre * r["b0"]),
                f"{centre:.3f} B₀" if centre > 1e-3 else f"{centre:.1e} B₀",
                delta_color="off", delta_arrow="off")
    k[1].metric("Slab thickness d", f"{r['d'] * 1e9:,.0f} nm", f"λ = {lam0_nm} nm",
                delta_color="off", delta_arrow="off")
    k[2].metric("Surface current J", f"{sci(surface_current)} A/m²", "at each face",
                delta_color="off", delta_arrow="off",
                help="Screening current density at each face: J = (1/μ₀) dB/dx = B₀ tanh(d/2λ) / (μ₀λ).")
    k[3].metric("Resolution", f"{per_lam:.0f} pts/λ" if per_lam >= 10 else f"{per_lam:.1f} pts/λ",
                f"Δx = {dx_nm:.3g} nm", delta_color="off", delta_arrow="off")
    for col, name in zip(k[4:], METHODS):
        m = r["m"][name]
        if "error" in m:
            col.metric(f"{SHORT[name]} max error", "failed", help=m["error"])
        else:
            col.metric(f"{SHORT[name]} max error", field_text(m["abs"]), f"{m['rel']:.1e} B₀",
                       delta_color="off", delta_arrow="off",
                       help="Maximum |B − B_analytic| over the grid, absolute and as a fraction of B₀.")

    UNITS = ["Scaled", "Physical"]
    left, right = st.columns(2, gap="medium")
    with left, st.container(border=True, key="card_profile"):
        a, b, c1 = st.columns([0.9, 2, 1.1], vertical_alignment="center")
        a.markdown('<p class="card-title">Field profile</p>', unsafe_allow_html=True)
        shown_p = b.pills("Shown on profile", ["Analytic", *METHODS], format_func=LABEL.__getitem__,
                          selection_mode="multi", default=["Analytic", "Finite difference"],
                          label_visibility="collapsed", key="pills_profile")
        units_p = c1.segmented_control("Profile units", UNITS, default="Scaled",
                                       label_visibility="collapsed", key="units_profile") or "Scaled"
        live_chart(profile_chart(r, shown_p, units_p, t), key="profile_chart")
    with right, st.container(border=True, key="card_error"):
        a, b, c2, c3 = st.columns([0.5, 1.45, 0.85, 1.1], vertical_alignment="center")
        a.markdown('<p class="card-title">Error</p>', unsafe_allow_html=True)
        shown_e = b.pills("Shown on error", list(METHODS), format_func=LABEL.__getitem__, selection_mode="multi",
                          default=["Finite difference"], label_visibility="collapsed", key="pills_error")
        scale = c2.segmented_control("Scale", ["Signed", "Log"], default="Signed",
                                     label_visibility="collapsed", key="scale") or "Signed"
        units_e = c3.segmented_control("Error units", UNITS, default="Scaled",
                                       label_visibility="collapsed", key="units_error") or "Scaled"
        live_chart(error_chart(r, shown_e, scale, units_e, t), key="error_chart")

    bad = [SHORT[n] for n, m in r["m"].items() if m.get("rel", 1) >= TARGET]
    if bad:
        st.warning(f"{', '.join(bad)} above the 10⁻³ target. Increase N or reduce d/λ.")

with temp_tab:
    tc_col, t_col, lam_col, d_col, b_col = st.columns(5, gap="large")
    with tc_col:
        t_c = st.number_input("Critical temperature T_c (K)", min_value=1.0, max_value=20.0,
                              value=9.2, step=0.1, key="temp_tc")
    clamp = clamp_temperature(t_c)
    if "temp_t" in st.session_state and st.session_state.temp_t > clamp:
        st.session_state.temp_t = float(clamp)
    with t_col:
        t_k = st.slider("Temperature T (K)", 0.0, float(clamp), min(4.2, float(clamp)), 0.05, key="temp_t")
    with lam_col:
        temp_lam0_nm = st.slider("Penetration depth λ₀ (nm)", 10, 200, 50, 5, key="temp_lam0")
    with d_col:
        temp_d_nm = st.slider("Slab thickness d (nm)", 50, 600, 200, 10, key="temp_d")
    with b_col:
        temp_b0_mT = st.slider("Applied field B₀ (mT)", 10, 500, 50, 10, key="temp_b0")

    lambda_0_m = temp_lam0_nm * 1e-9
    d_m = temp_d_nm * 1e-9
    b_0 = temp_b0_mT * 1e-3
    g = guard_metrics(lambda_0_m, t_c, d_m)
    lambda_t = lambda_of_temperature(t_k, lambda_0_m, t_c)
    centre_t = 1.0 / np.cosh(d_m / (2.0 * lambda_t))

    j = st.columns(3)
    j[0].metric("Clamp temperature", f"{g['clamp_t']:.2f} K", "0.995 T_c",
                delta_color="off", delta_arrow="off")
    j[1].metric("λ at clamp", f"{g['lambda_clamp'] * 1e9:.1f} nm",
                f"{g['lambda_clamp'] / lambda_0_m:.2f} λ₀", delta_color="off", delta_arrow="off")
    j[2].metric("B(0)/B₀ at clamp", f"{g['center_ratio']:.3f}", delta_color="off", delta_arrow="off")

    st.caption(f"T = {t_k:.2f} K ({t_k / t_c:.3f} T_c), λ = {lambda_t * 1e9:.1f} nm, "
               f"B(0)/B₀ = {centre_t:.3f}.")
    if t_k > 0.99 * t_c:
        st.warning("T is within 1% of Tc. The two-fluid model is outside its reliable range.")

    st.session_state.setdefault("show_normal", False)
    if st.button("Hide normal-state demo" if st.session_state["show_normal"]
                 else "Show T > Tc normal state", key="normal_toggle"):
        st.session_state["show_normal"] = not st.session_state["show_normal"]

    if st.session_state["show_normal"]:
        t_demo = 1.05 * t_c
        xn = np.linspace(-d_m / 2.0, d_m / 2.0, TEMP_POINTS)
        flat = temperature_profile(xn, d_m, lambda_0_m, t_c, t_demo, b_0)
        # Above Tc the family units control still decides the chart units.
        units_n = st.session_state.get("temp_units") or "B/B₀"
        if units_n == "B (mT)":
            y_n = flat * 1e3
            y_title = "B (mT)"
            y_top = 1.05 * b_0 * 1e3
        else:
            y_n = flat / b_0
            y_title = "B / B₀"
            y_top = 1.05
        with st.container(border=True, key="card_normal"):
            st.markdown('<p class="card-title">Normal state above Tc</p>', unsafe_allow_html=True)
            n = st.columns(3)
            n[0].metric("T/Tc", "1.05", delta_color="off", delta_arrow="off")
            n[1].metric("λ", "∞", delta_color="off", delta_arrow="off")
            n[2].metric("B(0)/B0", "1.000", delta_color="off", delta_arrow="off")
            fig = go.Figure()
            fig.add_scatter(x=xn * 1e9, y=y_n, mode="lines", name="B",
                            line={"color": t["analytic"], "width": 3}, hovertemplate="%{y:.4f}")
            fig.update_xaxes(title="x (nm)", range=[-d_m * 1e9 / 2.0, d_m * 1e9 / 2.0])
            fig.update_yaxes(title=y_title, range=[0, y_top])
            live_chart(style(fig, t, 320, uirevision="temp-normal"), key="temp_normal")
            normal_state_infographic(t)
            st.caption("Static normal-state limit. The model ignores the AC skin depth. "
                       "A real sample near Tc can have fluctuations or a mixed state.")

    left_t, right_t = st.columns(2, gap="medium")
    with left_t, st.container(border=True, key="card_temp_lambda"):
        st.markdown('<p class="card-title">Penetration depth</p>', unsafe_allow_html=True)
        live_chart(lambda_chart(t_c, lambda_0_m, t_k, t), key="temp_lambda")
    with right_t, st.container(border=True, key="card_temp_profile"):
        st.markdown('<p class="card-title">Field profiles</p>', unsafe_allow_html=True)
        units_t = st.segmented_control("Temperature units", ["B/B₀", "B (mT)"], default="B/B₀",
                                       label_visibility="collapsed", key="temp_units") or "B/B₀"
        live_chart(temperature_profile_chart(t_c, lambda_0_m, d_m, b_0, t_k, units_t, t),
                   key="temp_family")
        st.caption("Faint curves: T/Tc = 0, 0.217, 0.435, 0.652, 0.870, 0.957. Bold curve: current T.")
