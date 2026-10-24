"""Theory page: the London model, its length scales, and its limits."""

import math

import streamlit as st

from magpenetration_simulation import ui
from magpenetration_simulation import ui_config as cfg


def render(material):
    """Explain the London equation, the geometry families, and the limits."""
    lam = material["lambda"]
    lam_0 = material["lambda_0"]
    t_c = material["t_c"]
    ratio = material["t"] / t_c
    normal = not math.isfinite(lam)

    ui.page_title("Theory", "The London model behind every page of this app.")

    st.markdown(
        "A superconductor screens a static magnetic field over one length, "
        "the penetration depth λ. The London equation gives the field profile "
        "inside a sample. This page states that equation, the two-fluid "
        "temperature law, and the geometry families that the app solves."
    )

    # ---------------------------------------------------------------- London
    with ui.card("card_london"):
        ui.card_title("1. The London equation")
        st.markdown(
            "The London theory relates the screening current J to the vector "
            "potential **A**. In the London gauge the relation is local."
        )
        st.latex(r"\mathbf{J} = -\frac{\mathbf{A}}{\mu_0\lambda^2}")
        st.markdown(
            "Ampère's law then turns this relation into a field equation for "
            "the flux density **B** inside the sample."
        )
        st.latex(r"\nabla^2\mathbf{B} = \frac{\mathbf{B}}{\lambda^2}")
        st.markdown(
            "The sign on the right is positive. A field decays exponentially "
            "with depth and does not oscillate. The decay length is λ."
        )
        st.markdown(
            "This decay describes the Meissner state. The sample expels flux "
            "from its interior. A screening current flows in a surface layer "
            "of thickness of order λ. That current cancels the applied field "
            "deeper inside."
        )

    # ---------------------------------------------------------- depth and T
    left, right = st.columns(2, gap="medium")
    with left, ui.card("card_depth"):
        ui.card_title("2. The penetration depth")
        st.markdown(
            "λ is the distance over which the screening current reduces the "
            "field by a factor of e. It is a material property and it depends "
            "on temperature."
        )
        st.markdown("Typical zero-temperature values are:")
        st.markdown(
            "- Aluminum: **16 nm**\n"
            "- Niobium: **39 nm**\n"
            "- YBCO: **150 nm** in the CuO₂ planes"
        )
        st.markdown(
            "These values come from Tinkham (1996) and Poole et al. (2007). "
            "An anisotropic material needs one λ per crystal axis. YBCO is "
            "such a material, and its λ applies to the planes only."
        )
        st.markdown("The local London form needs a slow field variation:")
        st.latex(r"q\xi \ll 1, \qquad q = \frac{1}{\lambda}")
        st.markdown(
            "Here ξ is the coherence length. A long coherence length breaks "
            "this condition, and the nonlocal Pippard form applies instead."
        )

    with right, ui.card("card_temp"):
        ui.card_title("3. Temperature")
        st.markdown(
            "The two-fluid model splits the carriers into a normal part and a "
            "superfluid part. The superfluid fraction follows (T/T_c)⁴."
        )
        st.latex(r"\lambda(T) = \frac{\lambda_0}{\sqrt{1 - (T/T_c)^4}}")
        st.markdown(
            "λ starts at λ₀ and grows with temperature. The depth diverges as "
            "T approaches the critical temperature T_c, because the superfluid "
            "density falls to zero."
        )
        st.markdown(
            "At T_c the equation gives a flat field at the applied value B₀. "
            "The material is then normal and the field fills the sample. The "
            "screening current is zero."
        )
        st.markdown(
            f"The temperature slider stops at **{cfg.TEMPERATURE['clamp_factor']:.3f} T_c**. "
            "The clamp keeps λ finite. The slab equation then has a solution."
        )

    with ui.card("card_now"):
        ui.card_title("Current settings")
        k = st.columns(4, gap="medium")
        k[0].metric("λ₀", f"{lam_0 * 1e9:,.0f} nm", delta_color="off", delta_arrow="off")
        k[1].metric("T_c", f"{t_c:.2f} K", delta_color="off", delta_arrow="off")
        k[2].metric("T / T_c", f"{ratio:.3f}", delta_color="off", delta_arrow="off")
        k[3].metric("λ(T)", "infinite" if normal else f"{lam * 1e9:,.0f} nm",
                    delta_color="off", delta_arrow="off",
                    help="The sidebar sets these values. Every page uses them.")

    # ------------------------------------------------------------- geometry
    with ui.card("card_geom"):
        ui.card_title("4. Geometry families")
        st.markdown(
            "Each family below reduces the vector London equation to a form "
            "with a closed solution or a numeric solution."
        )
        a, b = st.columns(2, gap="medium")
        with a:
            st.markdown("**Slab**")
            st.markdown(
                "A wide plate of thickness d sits in a field parallel to its "
                "faces. One Cartesian component remains."
            )
            st.latex(r"B(x) = B_0\,\frac{\cosh(x/\lambda)}{\cosh\left(d/2\lambda\right)}")
            st.markdown(
                "x is the distance from the midplane. The field equals B₀ at "
                "both faces and is lowest at the center."
            )
            st.markdown("**Axial cylinder**")
            st.markdown(
                "A long cylinder of radius R sits in a field along its axis. "
                "The radial equation gives modified Bessel functions."
            )
            st.latex(r"B(r) = B_0\,\frac{I_0(r/\lambda)}{I_0(R/\lambda)}")
            st.markdown(
                "I₀ is the modified Bessel function of order zero. Only this "
                "solution has a finite value at the axis."
            )
        with b:
            st.markdown("**Sphere**")
            st.markdown(
                "A sphere in a uniform field keeps both a radial and an "
                "angular component. The exact London solution uses the "
                "modified spherical Bessel function i₁."
            )
            st.markdown(
                "The Sphere page reports this vector solution only. It draws "
                "no profile yet, and the next geometry slice will add one."
            )
            st.markdown("**Finite bodies**")
            st.markdown(
                "A box, a rectangle, or a disk has no closed London solution. "
                "This geometry needs a numeric solve over the sample and its "
                "surroundings."
            )
            st.markdown(
                "Near an edge the screening current turns, and the one-"
                "dimensional forms fail. A corner peak depends on the grid "
                "size, not on the sample. Only an ellipsoid has an exact "
                "demagnetizing factor."
            )
            st.markdown("**Thin film**")
            st.markdown(
                "Below d = 2λ a film screens in plane over the Pearl length "
                "instead of λ."
            )
            st.latex(r"\Lambda = \frac{2\lambda^2}{d}")
            st.markdown(
                "The slab profile is exact at every thickness when the field "
                "is parallel. A perpendicular field spreads the current over Λ."
            )

    # ------------------------------------------------------------ this app
    with ui.card("card_compute"):
        ui.card_title("5. What this app computes")
        st.markdown(
            "The app solves the slab equation three ways and compares each "
            "answer with the exact cosh profile:"
        )
        st.markdown(
            "- **Finite difference**: a banded linear solve on the shared grid.\n"
            "- **Relaxation (SOR)**: a sweep to a residual tolerance.\n"
            "- **Shooting**: an integration from one face to the other."
        )
        st.markdown(
            "The **Methods** page shows the three error curves against the "
            "exact slab solution. A smaller d/λ or a larger N improves them."
        )
        st.markdown(
            "The **Slab** page plots the analytic profile, the temperature "
            "family, and a thin-film warning for d < 2λ. The **Cylinder** page "
            "plots the Bessel profile over the R/λ family. The **Combined** "
            "page puts the slab and the cylinder on one depth axis, measured "
            "from the nearest face. The **Sphere** page reports the exact "
            "vector solution."
        )
        st.markdown(
            "The navigation bar at the top of the window shows all six pages. "
            "Select a name there to change page."
        )

    # -------------------------------------------------------------- limits
    with ui.card("card_limits"):
        ui.card_title("6. Assumptions and limits")
        st.markdown(
            "- The theory is local London theory. A long coherence length "
            "needs the Pippard kernel.\n"
            "- The profiles describe the Meissner state only. A type-II sample "
            "above its first critical field holds vortices, and this app does "
            "not model them.\n"
            "- The applied field is uniform and static. The model has no time "
            "dependence and no skin effect.\n"
            "- The closed forms need full symmetry. Keep a reported value "
            "several λ away from an edge or an end.\n"
            "- The slab and cylinder boundary values equal B₀. Demagnetization "
            "raises the field at a real boundary, and this app does not "
            "correct for it.\n"
            "- The two-fluid law comes from a fit to measured data. It rounds "
            "off near T_c.\n"
            "- The 2D field-line layer is a separate drawing. It does not "
            "change the numbers on the profile charts."
        )
