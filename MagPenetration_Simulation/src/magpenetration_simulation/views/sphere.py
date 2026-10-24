"""Sphere page: the exact London sphere solution, shown without a solver."""

import streamlit as st

from magpenetration_simulation import ui


def render(material):
    """Placeholder until the Sphere solver lands, with the exact solution on show."""
    ui.page_title(
        "Sphere",
        "A sphere of radius R in a uniform field B0. The field has two components here.",
    )

    st.markdown(
        "**The interactive sphere solver is the next geometry slice.** "
        "This page shows the exact London solution and its limits. No values are computed here."
    )

    with ui.card("card_sphere_solution"):
        ui.card_title("Exact London solution")
        st.markdown(
            "A toroidal potential psi(r) sin(theta) separates the vector problem. "
            "The radial part uses the modified spherical Bessel function of order 1."
        )
        st.latex(
            r"\psi(r) = \frac{3 B_0 R}{2 \sinh(R/\lambda)}\,"
            r"i_1\!\left(\frac{r}{\lambda}\right)"
        )
        st.markdown("The inside field then has these components:")
        st.latex(r"B_r = \frac{2 \psi(r)}{r} \cos\theta")
        st.latex(r"B_\theta = -\frac{\psi(r) + r\,\psi'(r)}{r} \sin\theta")

    with ui.card("card_sphere_moment"):
        ui.card_title("Outside field and scalars")
        st.markdown(
            "The outside field is the applied field plus a dipole. "
            "The screening factor S sets the dipole moment."
        )
        st.latex(r"m = -\frac{2 \pi R^3 B_0}{\mu_0}\, S")
        st.latex(
            r"S = 1 - \frac{3 \lambda}{R} \coth\!\left(\frac{R}{\lambda}\right)"
            r" + \frac{3 \lambda^2}{R^2}"
        )
        st.markdown("Two scalar values follow from the same solution:")
        st.latex(r"B(0) = B_0\, \frac{R/\lambda}{\sinh(R/\lambda)}")
        st.latex(
            r"\langle B_z \rangle = B_0 \left("
            r"\frac{3 \lambda}{R} \coth\!\left(\frac{R}{\lambda}\right)"
            r" - \frac{3 \lambda^2}{R^2}\right)"
        )
        st.markdown(
            "B(0) is the field at the center of the sphere. "
            "The second value is the mean along the polar axis."
        )

    left, right = st.columns(2, gap="medium")

    with left, ui.card("card_sphere_demag"):
        ui.card_title("Demagnetization caveat")
        st.markdown(
            "- The boundary value varies with latitude.\n"
            "- The equator value approaches 3 B0 / 2 in the ideal limit.\n"
            "- The field at the poles approaches zero in the same limit.\n"
            "- Only ellipsoids have exact demagnetizing factors. "
            "The sphere has N = 1/3.\n\n"
            "The applied field is then not the boundary value. "
            "The geometry redistributes the field."
        )

    with right, ui.card("card_sphere_scalar"):
        ui.card_title("Scalar-form caveat")
        st.markdown(
            "Some notes reduce the sphere to a radial scalar form, "
            "with zero angular momentum.\n\n"
            "- The radial scalar form is not the full vector solution.\n"
            "- It gives the center value only.\n"
            "- It misses the angular part, so the equator field departs from it.\n\n"
            "Use it as a rough magnitude guide. Use the vector form for real numbers."
        )

    with ui.card("card_sphere_next"):
        ui.card_title("What will land here")
        st.markdown(
            "- Sphere radial functions\n"
            "- Vector field components\n"
            "- Demagnetization-aware boundary values\n"
            "- A 2D field-line view in the next visualization slice"
        )