"""Page views for the penetration model app.

Each module exposes render() with no arguments. A page that needs material
values renders its own controls with ui.material_controls(prefix), so pages do
not share values.
"""

from magpenetration_simulation.views import (
    combined, cylinder, film, finite_body, methods, slab, sphere, theory)

__all__ = [
    "combined", "cylinder", "film", "finite_body", "methods", "slab", "sphere", "theory"]
