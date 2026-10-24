"""Tunable UI settings for the penetration model pages.

Edit this file to tune the app. It holds UI ranges, defaults, steps, chart
sizes, families, the temperature clamp, the error target, and plot headroom.
Colours and CSS stay in `ui.py`.
"""

# Each widget group: min, max, default and step. No max means the clamp sets it.
SIDEBAR = {
    "t_c": {"min": 1.0, "max": 20.0, "default": 9.2, "step": 0.1},
    "t": {"min": 0.0, "default": 4.2, "step": 0.05},
    "lambda_0": {"min": 10, "max": 200, "default": 50, "step": 5},
    "b_0": {"min": 10, "max": 500, "default": 50, "step": 10},
}
TEMPERATURE = {"clamp_factor": 0.995}
GRID = {"min": 41, "max": 2001, "default": 401, "step": 20}
METHODS = {
    "d_over_lam": {"min": 0.5, "max": 100.0, "default": 5.0, "step": 0.5},
    "target": 1e-3,
}
SLAB = {
    "d_over_lam": {"min": 0.5, "max": 100.0, "default": 5.0, "step": 0.5},
}
# near_wall_min and near_wall_clip set when the thick-cylinder limit is drawn.
CYLINDER = {
    "size": {"min": 0.5, "max": 20.0, "default": 5.0, "step": 0.5},
    "near_wall_min": 10.0,
    "near_wall_clip": 1.02,
}
COMBINED = {
    "size": {"min": 0.5, "max": 20.0, "default": 5.0, "step": 0.5},
}
# y_headroom and x_pad_* set the chart axis padding; log_decades sets the error axis.
PLOT = {
    "y_headroom": 1.05,
    "x_pad_fraction": 0.08,
    "x_pad_min": 0.6,
    "log_decades": 16,
}
CHART_H = 480
TEMP_POINTS = 401
FAMILY_T = (0.0, 0.217, 0.435, 0.652, 0.870, 0.957)
GEOM_FAMILY = (1.0, 2.0, 5.0, 10.0, 20.0)
GEOM_VIEWS = ["B/B₀", "B (mT)"]
PAGES = (
    ("theory", "Theory"),
    ("methods", "Methods"),
    ("slab", "Slab"),
    ("cylinder", "Cylinder"),
    ("combined", "Combined"),
    ("sphere", "Sphere"),
)
