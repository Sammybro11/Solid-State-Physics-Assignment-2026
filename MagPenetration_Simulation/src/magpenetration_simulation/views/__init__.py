"""Page views for the penetration model app.

Each module exposes render(material), where material is the dict built by
app.material_state(): t, t_c in K, lambda_0 in m, b_0 in T, lambda in m, and
normal_state as a bool.
"""

from magpenetration_simulation.views import (
    combined, cylinder, film, finite_body, methods, slab, sphere, theory)

__all__ = [
    "combined", "cylinder", "film", "finite_body", "methods", "slab", "sphere", "theory"]
