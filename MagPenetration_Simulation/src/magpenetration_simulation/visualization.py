"""Styled Plotly field sections for the geometry pages.

Each section shows a heatmap of |B|/B0 on a fixed 0..1 scale, so the pages can be
compared at a glance. The slab, cylinder, sphere and film sections add the field
lines the scalar 1D models imply; the finite-body solver is scalar with B out of
plane, so its section carries markers instead of fake in-plane lines.

Themes come from the caller as a ui.THEMES dict. No Streamlit here.
"""


import numpy as np
import plotly.graph_objects as go
from scipy.integrate import solve_ivp

from magpenetration_simulation import ui
from magpenetration_simulation import ui_config as cfg
from magpenetration_simulation.cylinder import cylinder_field_analytic
from magpenetration_simulation.finite_body import rect_field_finite_difference
from magpenetration_simulation.slab import slab_field_analytic
from magpenetration_simulation.sphere import (
    MU0,
    sphere_center_field,
    sphere_dipole_moment,
    sphere_field,
)

__all__ = [
    "combined_section",
    "cylinder_section",
    "film_section",
    "finite_body_section",
    "slab_section",
    "sphere_section",
]


# ---------------------------------------------------------------- shared pieces

def _heatmap(x, y, ratio, t, title="|B|/B₀", colorbar=True):
    """A heatmap on the fixed B/B0 scale, coloured from the theme background to its accent."""
    return go.Heatmap(x=x, y=y, z=ratio, zmin=0.0, zmax=1.0,
                      colorscale=[[0.0, t["bg"]], [1.0, t["accent"]]],
                      showscale=colorbar, hovertemplate=title + " = %{z:.3f}<extra></extra>",
                      colorbar={"title": {"text": title}, "thickness": 12, "len": 0.85})


def _vlines(xs, ratio, y0, y1, t):
    """Vertical field lines whose opacity follows the local B/B0."""
    return [go.Scatter(x=[x, x], y=[y0, y1], mode="lines", hoverinfo="skip",
                       opacity=float(np.clip(o, 0.08, 1.0)),
                       line={"color": t["analytic"], "width": 1.2})
            for x, o in zip(xs, ratio)]


def _b0_note(fig, t):
    """Mark the applied-field direction on a section."""
    fig.add_annotation(text="B₀ ↑", x=0.02, y=0.98, xref="paper", yref="paper",
                       showarrow=False, font={"size": 12, "color": t["muted"]})


def _finish(fig, t, x_title, y_title, key):
    """Equal-axis section layout, with hover and zoom keyed to the geometry."""
    fig.update_xaxes(title=x_title)
    fig.update_yaxes(title=y_title, scaleanchor="x", scaleratio=1)
    # style sets hovermode="x unified", so closest is applied after it.
    fig = ui.style(fig, t, cfg.VIS["height"], uirevision=key).update_layout(hovermode="closest")
    return fig.update_xaxes(showgrid=False).update_yaxes(showgrid=False)  # grid lines clutter a heatmap


def _outline(name, ux, uy, t):
    """The geometry edge on a small multiple: sphere surface, rectangle, or the two faces at x = ±1."""
    if name == "Sphere":
        phi = np.linspace(0.0, 2.0 * np.pi, 121)
        x, y = np.cos(phi), np.sin(phi)
    elif name == "Finite body":  # the solved domain is the body, so its extent is the edge
        x = [ux[0], ux[-1], ux[-1], ux[0], ux[0]]
        y = [uy[0], uy[0], uy[-1], uy[-1], uy[0]]
    else:  # slab, film and cylinder: a face each side, None breaks the line between them
        x, y = [-1.0, -1.0, None, 1.0, 1.0], [uy[0], uy[-1], None, uy[0], uy[-1]]
    return go.Scatter(x=x, y=y, mode="lines", hoverinfo="skip",
                      line={"color": t["axis"], "width": 1.4})


def _slab_ratio(xs, d, lam, b0):
    """|B|/b0 across the slab section: the slab field inside, flat at B0 in the vacuum."""
    inside = np.abs(xs) <= d / 2.0
    if not np.isfinite(lam):  # normal state: no screening, so nothing to solve
        return np.ones_like(xs)
    return np.where(inside, slab_field_analytic(np.clip(xs, -d / 2.0, d / 2.0), d, lam, b0) / b0, 1.0)


def _slab_grid(d, lam, b0):
    """xs, zs and |B|/b0 over the slab x-z section, flat at B0 in the vacuum margin."""
    xs = np.linspace(-0.75 * d, 0.75 * d, cfg.VIS["points"])
    zs = np.linspace(-0.8 * d, 0.8 * d, cfg.VIS["points"])
    return xs, zs, np.broadcast_to(_slab_ratio(xs, d, lam, b0)[None, :], (zs.size, xs.size))

def _cylinder_ratio(rs, radius, lam, b0):
    """|B|/b0 across the cylinder section: the wall field inside, flat at B0 outside."""
    inside = np.abs(rs) <= radius
    if not np.isfinite(lam):
        return np.ones_like(rs)
    r = np.minimum(np.abs(rs), radius)
    return np.where(inside, cylinder_field_analytic(r, radius, lam, b0) / b0, 1.0)


def _cylinder_grid(radius, lam, b0):
    """rs, zs and |B|/b0 over the cylinder r-z section, flat at B0 outside the wall."""
    rs = np.linspace(-1.3 * radius, 1.3 * radius, cfg.VIS["points"])
    zs = np.linspace(-1.6 * radius, 1.6 * radius, cfg.VIS["points"])
    return rs, zs, np.broadcast_to(_cylinder_ratio(rs, radius, lam, b0)[None, :], (zs.size, rs.size))


def _sphere_ratio(x, z, radius, lam, b0):
    """B/B0 in the x-z plane: London field inside, uniform field plus dipole outside."""
    x, z = np.broadcast_arrays(np.atleast_1d(np.asarray(x, dtype=float)),
                              np.atleast_1d(np.asarray(z, dtype=float)))
    # The dipole terms are in units of R: mu0 m / (4 pi b0 R**3) is dimensionless.
    xn, zn = x / radius, z / radius
    r = np.hypot(xn, zn)
    inside = r <= 1.0
    m = _m_norm(radius, lam, b0)
    # max(r, 1) keeps the exterior r**-3 and r**-5 terms finite inside, overwritten below.
    r_safe = np.maximum(r, 1.0)
    bx = np.where(inside, 0.0, 3.0 * m * xn * zn / r_safe**5)
    bz = np.where(inside, 0.0, m * (3.0 * zn**2 / r_safe**5 - 1.0 / r_safe**3) + 1.0)
    if np.any(inside):
        th = np.arctan2(xn[inside], zn[inside])
        b_r, b_th = sphere_field(r[inside] * radius, th, radius, lam, b0)
        bx[inside] = (b_r * np.sin(th) + b_th * np.cos(th)) / b0
        bz[inside] = (b_r * np.cos(th) - b_th * np.sin(th)) / b0
    return bx, bz


def _m_norm(radius, lam, b0):
    """mu0 * m / (4 pi b0 R**3), which equals -S/2 for the diamagnetic sphere."""
    return MU0 * sphere_dipole_moment(radius, lam, b0) / (4.0 * np.pi * b0 * radius**3)


def _sphere_grid(radius, lam, b0):
    """xs, zs and |B|/b0 over the square x-z section around the sphere."""
    xs = np.linspace(-1.8 * radius, 1.8 * radius, cfg.VIS["points"])
    bx, bz = _sphere_ratio(xs[None, :], xs[:, None], radius, lam, b0)
    return xs, xs, np.hypot(bx, bz)


def _rect_grid(solution, b0):
    """solution.x, solution.y and B/B0 from the solved rectangle."""
    return solution.x, solution.y, solution.b / b0


# ---------------------------------------------------------------- sections

def slab_section(d, lam, b0, t):
    """Longitudinal x-z section of a slab of thickness d, in a vacuum margin on each side."""
    xs, zs, ratio = _slab_grid(d, lam, b0)
    fig = go.Figure()
    fig.add_trace(_heatmap(xs / d, zs / d, ratio, t))
    seeds = np.linspace(-0.75 * d, 0.75 * d, cfg.VIS["line_seeds"])
    fig.add_traces(_vlines(seeds / d, _slab_ratio(seeds, d, lam, b0), zs[0] / d, zs[-1] / d, t))
    for face in (-0.5, 0.5):  # the two faces only: a slab has no top or bottom edge
        fig.add_vline(x=face, line={"color": t["axis"], "width": 1.5})
    _b0_note(fig, t)
    return _finish(fig, t, "x / d", "z / d", f"slab-{d:.6e}")


def cylinder_section(radius, lam, b0, t):
    """Longitudinal r-z section of an infinite cylinder, in a vacuum margin."""
    rs, zs, ratio = _cylinder_grid(radius, lam, b0)
    fig = go.Figure()
    fig.add_trace(_heatmap(rs / radius, zs / radius, ratio, t))
    seeds = np.linspace(-1.3 * radius, 1.3 * radius, cfg.VIS["line_seeds"])
    fig.add_traces(_vlines(seeds / radius, _cylinder_ratio(seeds, radius, lam, b0),
                           zs[0] / radius, zs[-1] / radius, t))
    for wall in (-1.0, 1.0):  # the two walls only: an infinite cylinder has no end caps
        fig.add_vline(x=wall, line={"color": t["axis"], "width": 1.5})
    _b0_note(fig, t)
    return _finish(fig, t, "r / R", "z / R", f"cylinder-{radius:.6e}")


def film_section(d, lam, b0, t):
    """A thin film section: the slab model, with d as the film thickness and no substrate."""
    return slab_section(d, lam, b0, t)


def sphere_section(radius, lam, b0, t):
    """Real in-plane field lines of the sphere in the x-z plane, in units of R."""
    xs, zs, ratio = _sphere_grid(radius, lam, b0)
    fig = go.Figure()
    fig.add_trace(_heatmap(xs / radius, zs / radius, ratio, t))
    phi = np.linspace(0.0, 2.0 * np.pi, 121)
    fig.add_scatter(x=np.cos(phi), y=np.sin(phi), mode="lines", hoverinfo="skip",
                    line={"color": t["axis"], "width": 1.5})

    def direction(_s, p):
        bx, bz = _sphere_ratio(p[0], p[1], radius, lam, b0)
        n = max(float(np.hypot(bx, bz)[0]), 1e-12)
        return bx[0] / n, bz[0] / n

    for x0 in np.linspace(-1.5, 1.5, cfg.VIS["line_seeds"]):
        path = solve_ivp(direction, (0.0, 12.0), (x0, -cfg.VIS["margin"]), max_step=0.03,
                         rtol=1e-6, atol=1e-9)
        on_plot = np.maximum(np.abs(path.y[0]), np.abs(path.y[1])) <= 1.75
        fig.add_scatter(x=path.y[0][on_plot], y=path.y[1][on_plot], mode="lines", hoverinfo="skip",
                        line={"color": t["analytic"], "width": 1.0})
    _b0_note(fig, t)
    return _finish(fig, t, "x / R", "z / R", f"sphere-{radius:.6e}")


def finite_body_section(solution, b0, t):
    """x-y cross-section of the solved rectangle in nm; B is out of plane, so no in-plane lines."""
    width = (solution.x[-1] - solution.x[0]) * 1e9
    height = (solution.y[-1] - solution.y[0]) * 1e9
    fig = go.Figure()
    fig.add_trace(_heatmap(solution.x * 1e9, solution.y * 1e9, solution.b / b0, t, "B/B₀"))
    fig.add_shape(type="rect", x0=-width / 2, y0=-height / 2, x1=width / 2, y1=height / 2,
                  line={"color": t["axis"], "width": 1.5})
    gx, gy = np.meshgrid(np.linspace(-0.3, 0.3, 3) * width, np.linspace(-0.3, 0.3, 3) * height)
    fig.add_scatter(x=gx.ravel(), y=gy.ravel(), mode="markers", hoverinfo="skip",
                    marker={"symbol": "circle-open", "size": 9, "color": t["text"]})
    fig.add_annotation(text="B is out of plane (⊙)", x=0.5, y=0.5, xref="paper", yref="paper",
                       xanchor="right", yanchor="bottom", showarrow=False,
                       font={"size": 11, "color": t["muted"]})
    return _finish(fig, t, "x (nm)", "y (nm)", f"rect-{width:.6e}-{height:.6e}")


def combined_section(size, lam, b0, solution, t):
    """2 x 3 small multiples of |B|/B0 at the common size s/λ = d/(2λ) = R/λ = size."""
    d = 2.0 * size * lam
    if np.isfinite(lam):
        scale = size * lam
        panels = (("Slab", _slab_grid(d, lam, b0)), ("Cylinder", _cylinder_grid(scale, lam, b0)),
                  ("Sphere", _sphere_grid(scale, lam, b0)), ("Finite body", _rect_grid(solution, b0)),
                  ("Film", _slab_grid(d, lam, b0)))
    else:  # normal state: no screening anywhere, so every panel is flat at B0
        u = np.linspace(-1.5, 1.5, cfg.VIS["points"])
        flat = np.ones((u.size, u.size))
        scale = 1.0
        panels = tuple((name, (u, u, flat)) for name in
                       ("Slab", "Cylinder", "Sphere", "Finite body", "Film"))
    fig = go.Figure()
    fig.set_subplots(rows=2, cols=3, subplot_titles=[name for name, _ in panels] + ["Colour"],
                     horizontal_spacing=0.06, vertical_spacing=0.12)
    for i, (name, (xs, ys, ratio)) in enumerate(panels):
        row, col = i // 3 + 1, i % 3 + 1
        ux, uy = xs / scale, ys / scale
        fig.add_trace(_heatmap(ux, uy, ratio, t, colorbar=False), row=row, col=col)
        fig.add_trace(_outline(name, ux, uy, t), row=row, col=col)
        fig.update_xaxes(range=[ux[0], ux[-1]], showgrid=False, row=row, col=col)
        fig.update_yaxes(range=[uy[0], uy[-1]], showgrid=False, row=row, col=col)
    for suffix in ("", "2", "3", "4", "5", "6"):  # one pixel per unit in every panel
        fig.update_layout({f"xaxis{suffix}": {"scaleanchor": f"y{suffix}", "scaleratio": 1}})
    # The empty sixth cell holds the shared colour bar; opacity keeps its own pixels clear.
    fig.add_trace(go.Heatmap(z=[[0.0, 1.0]], zmin=0.0, zmax=1.0, showscale=True,
                             colorscale=[[0.0, t["bg"]], [1.0, t["accent"]]], opacity=0.01,
                             hoverinfo="skip",
                             colorbar={"title": {"text": "|B|/B₀"}, "x": 0.98, "y": 0.25,
                                       "len": 0.35, "thickness": 12}), row=2, col=3)
    fig.update_xaxes(visible=False, row=2, col=3)
    fig.update_yaxes(visible=False, row=2, col=3)
    fig = ui.style(fig, t, cfg.VIS["height"], uirevision=f"combined-{size:.6e}")
    fig.update_layout(hovermode="closest", height=cfg.VIS["height"] + 60,
                      margin={"l": 12, "r": 12, "t": 36, "b": 12})  # room for the panel titles
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False)
    return fig


# ---------------------------------------------------------------- self-check

def _self_check():
    """Assert the sections reproduce the model values and stay finite."""
    t = ui.THEMES["light"]
    lam, b0 = 50e-9, 0.05
    d = 10 * lam

    _, _, slab = _slab_grid(d, lam, b0)
    mid = slab.shape[1] // 2
    assert slab[0, 0] == 1.0, "slab vacuum margin must be B0"
    assert slab[0, mid] < 1.0, "slab centre must be screened"

    rs, _, cyl = _cylinder_grid(5 * lam, lam, b0)
    wall = int(np.argmin(np.abs(rs - 5 * lam)))
    assert cyl[0, rs.size // 2] < cyl[0, wall], "cylinder axis must be below the wall"

    radius = 5 * lam
    _, _, sphere = _sphere_grid(radius, lam, b0)
    centre = sphere.shape[0] // 2
    assert abs(sphere[centre, centre] - sphere_center_field(radius, lam, b0) / b0) < 1e-9
    assert sphere[0, 0] < 3.0, "the exterior dipole must stay on the B0 scale"
    paths = sphere_section(radius, lam, b0, t).data[2:]  # after the heatmap and the outline
    for path in paths:
        assert np.all(np.isfinite(path.x)) and np.all(np.isfinite(path.y)), "field line is not finite"
    assert any(len(path.x) > 5 for path in paths), "no field line was integrated"

    solution = rect_field_finite_difference(6 * lam, 6 * lam, lam, b0, nx=41, ny=41, tol=1e-6)
    assert np.all(solution.b >= 0.0) and np.all(solution.b <= b0 + 1e-12), "finite body B out of range"
    finite_body_section(solution, b0, t)
    wide = finite_body_section(
        rect_field_finite_difference(6 * lam, 3 * lam, lam, b0, nx=41, ny=41, tol=1e-6), b0, t)
    rect = wide.layout.shapes[0]
    assert abs((rect.x1 - rect.x0) / (rect.y1 - rect.y0) - 2.0) < 1e-9, "finite-body pane must keep its aspect"

    fig = combined_section(5.0, lam, b0, solution, t)
    maps = [trace for trace in fig.data if isinstance(trace, go.Heatmap)]
    assert len(maps) == 6, f"combined must hold five panels plus the colour-bar cell, got {len(maps)}"
    assert all(np.all(np.isfinite(m.z)) for m in maps), "combined heatmap holds non-finite values"
    assert sum(bool(getattr(trace, "showscale", False)) for trace in fig.data) == 1, \
        "combined must carry one shared colour bar"
    assert len(fig.layout.annotations) == 6, "combined must label all six cells"
    wide_fig = combined_section(
        5.0, lam, b0, rect_field_finite_difference(6 * lam, 3 * lam, lam, b0, nx=41, ny=41, tol=1e-6), t)
    x_range, y_range = wide_fig.layout.xaxis4.range, wide_fig.layout.yaxis4.range  # finite-body cell
    assert abs((x_range[1] - x_range[0]) / (y_range[1] - y_range[0]) - 2.0) < 1e-9, \
        "combined finite-body panel must keep its aspect"

    s, s_slab, s_cyl, cyl_curve = ui.geometry_curves(5.0, lam, b0, 161)
    assert abs(s[-1] - 5.0) < 1e-12 and abs(s_cyl[-1] - 5.0) < 1e-12, "curves must end at s = size"
    assert s_slab[-1] < 1.0 and cyl_curve[-1] < 1.0, "both curves must end screened"


if __name__ == "__main__":
    _self_check()
