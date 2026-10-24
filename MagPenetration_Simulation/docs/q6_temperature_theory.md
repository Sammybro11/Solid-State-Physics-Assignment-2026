# Temperature Dependence of the Magnetic-Field Penetration Depth

**Package** MagPenetration_Simulation · **Issue** #2 · **Author** Researcher 2 · **Supports** Plan.md items 6.4 and 6.7

## 1. Purpose and Scope

This report explains how temperature changes the magnetic-field penetration depth of a superconductor.
The depth $\lambda(T)$ sets the length scale for screening.
A rise in $T$ weakens the screening, and the field then penetrates farther into the metal.

Item 6.4 of Plan.md covers the slider demo.
Item 6.7 covers the results interpretation.
The slab model of item 6.2 consumes $\lambda(T)$ as its input.
The cylinder case and the heat maps stay out of scope.
Section 6 lists the numerical guards that the app needs.

## 2. The Two-Fluid Model

Gorter and Casimir proposed the two-fluid picture in 1934 [2].
The model splits the conduction electrons into two fluids.
The normal component has density $n_n$, and the superfluid component has density $n_s$.
Only the superfluid conducts the dissipationless current and screens the field.
The densities obey one conservation law:

$$n = n_n + n_s$$

Here $n$ is the total conduction-electron density, and $n$ does not vary with temperature.
The model sets the fractions with an empirical law [2]:

$$f_n = \left(\frac{T}{T_c}\right)^4, \qquad f_s = 1 - \left(\frac{T}{T_c}\right)^4, \qquad n_s = f_s\, n$$

The fourth power is a fit to the classic data, and the law is not exact.
Table 1 defines every symbol and unit.

| Symbol | Definition | Unit |
|---|---|---|
| $T$ | temperature | K |
| $T_c$ | critical temperature | K |
| $\lambda(T)$ | penetration depth at temperature $T$ | m |
| $\lambda_0$ | penetration depth at $T = 0$ | m |
| $n$ | total conduction-electron density | $\mathrm{m}^{-3}$ |
| $n_s$ | superfluid electron density | $\mathrm{m}^{-3}$ |
| $n_n$ | normal-fluid electron density | $\mathrm{m}^{-3}$ |
| $f_s$, $f_n$ | superfluid and normal fractions | dimensionless |
| $m$ | electron mass | kg |
| $e$ | elementary charge | C |
| $\mu_0$ | vacuum permeability | $\mathrm{H\,m}^{-1}$ |
| $d$ | slab thickness | m |
| $B_0$ | applied field at $x = \pm d/2$ | T |

Substitution then gives the law of item 6.4:

$$\lambda(T) = \frac{\lambda_0}{\sqrt{1 - (T/T_c)^4}}$$

The ratio $T/T_c$ measures the distance from the transition.
The factor $(T/T_c)^4$ is the normal fraction $f_n$.
The factor $1 - (T/T_c)^4$ is the superfluid fraction $f_s$.
The constant $\lambda_0$ fixes the depth scale.
Textbook tables list values of order 10 nm to 100 nm for the classic elements [4], [5].
At $T_c$ the superfluid fraction hits zero and the law diverges.

## 3. Derivation Sketch

London and London postulated a local relation between the supercurrent and the vector potential [1]:

$$\mathbf{J}_s = -\frac{n_s e^2}{m}\,\mathbf{A}$$

Ampère's law turns this relation into a field equation inside the metal [4], [6]:

$$\nabla^2 \mathbf{B} = \frac{1}{\lambda^2}\,\mathbf{B}, \qquad \lambda^2 = \frac{m}{\mu_0 n_s e^2}$$

The depth $\lambda$ is the length scale in the steady state.
The two-fluid law gives $n_s(T) = n\,[1 - (T/T_c)^4]$.
One substitution then produces the temperature dependence:

$$\lambda^2(T) = \frac{m}{\mu_0 e^2 n\,\left[1 - (T/T_c)^4\right]} = \frac{\lambda_0^2}{1 - (T/T_c)^4}, \qquad \lambda_0^2 = \frac{m}{\mu_0 e^2 n}$$

The square root of both sides returns the law of section 2:

$$\lambda(T) = \frac{\lambda_0}{\sqrt{1 - (T/T_c)^4}}$$

The substitution inherits every approximation of both models:

1. The model takes $n_s(0) = n$, and every conduction electron then joins the superfluid at $T = 0$.
This step is an idealization, and a fit to data fixes $\lambda_0$.
2. The symbol $m$ means the band effective mass $m^*$ of the carriers.
3. The exponent 4 comes from the Gorter-Casimir fit [2] and has no microscopic derivation.
4. The London relation is local.
Nonlocal electrodynamics enters when the coherence length $\xi_0$ is large against $\lambda$ [4].
5. The law covers $0 \le T < T_c$, and it has no value at $T_c$.

## 4. Limits of $\lambda(T)$

Let $t = T/T_c$.

**At $T \to 0$.** Here $t^4 \to 0$, $n_s \to n$, and $\lambda \to \lambda_0$.
The expansion $1/\sqrt{1 - t^4} \approx 1 + t^4/2$ shows how flat the curve is.
At $t = 0.5$ the depth has grown by only 3.3 %.

**Near $T_c$.** Put $t = 1 - \epsilon$ for a small $\epsilon$.
Then $t^4 \approx 1 - 4\epsilon$, and the denominator reduces to $1 - t^4 \approx 4(1 - T/T_c)$:

$$\lambda(T) \approx \frac{\lambda_0}{2\,\sqrt{1 - T/T_c}}$$

The depth diverges with the exponent $-1/2$.
The superfluid fraction vanishes, and the screening fails with it.
At $T = T_c$ the formula has no value, and the exact limit is the uniform field $B(x) = B_0$.

**The flat field at $T \ge T_c$.** At $T_c$ the superfluid density $n_s$ goes to zero.
The London depth $\lambda$ then diverges.
The London equation becomes $\nabla^2 \mathbf{B} = \mathbf{B}/\lambda^2 \to 0$.
With the boundary values $B(\pm d/2) = B_0$, the solution is the uniform field $B(x) = B_0$.
The screening current $\mathbf{J} = \nabla \times \mathbf{B}/\mu_0$ is then zero.
The material is normal.
In the static limit, a nonmagnetic normal metal has a uniform field inside.
At $T_c$ the critical field $B_c(T)$ goes to zero, so any small $B_0$ destroys the superconducting state.
The flat field is the normal-state limit that the toy model uses.
A real normal metal in an AC field has a skin depth.
This model is static, so it has no skin effect.
A real sample near $T_c$ can have fluctuations and can enter a mixed state.
A type-I sample can enter an intermediate state.

**Why the model breaks in the critical region.** Both curves here are mean-field curves, the two-fluid fit and the Ginzburg-Landau curve [7].
Mean-field theory ignores fluctuations of the order parameter.
The fluctuations grow close to $T_c$.
The observed transition then rounds off.
The measured curve bends away from the divergent model curves [5].
The region is narrow for the classic superconductors.
The demo can ignore it.

**The BCS gap.** BCS computes the superfluid density from the energy gap $\Delta(T)$ [3].
Thermal quasiparticles account for the missing density.
Their number contains the factor $\exp(-\Delta/k_B T)$.
The low-$T$ correction to $\lambda$ then contains that exponential, and the power law lacks it.
Weak-coupling BCS places the gap at $\Delta_0 \approx 1.76\,k_B T_c$ [3].
Microscopic theory near $T_c$ gives $n_s \propto 1 - T/T_c$.
The BCS depth then diverges with the same exponent $-1/2$ [3], [5].
Ginzburg and Landau found the same exponent from their mean-field functional [7].
The prefactor differs, because $\lambda_0$ is an empirical scale and the microscopic scale follows from the gap.

## 5. Effect of Temperature on the Slab Profile

Consider a slab of thickness $d$ in a uniform applied field $B_0$ parallel to the slab face.
The coordinate $x$ runs from the mid-plane.
The faces follow the coordinates $x = \pm d/2$.
The static London equation with the listed boundary values gives the solution of item 6.2:

$$B(x) = B_0\,\frac{\cosh(x/\lambda)}{\cosh\left(\frac{d}{2\lambda}\right)}$$

The temperature enters through $\lambda(T)$ alone.
A rise in $T$ grows $\lambda$ and shrinks the ratio $d/(2\lambda)$.
The profile then flattens.
In the limit $\lambda \gg d$ the profile approaches the uniform field $B_0$, and the screening stops.
For a sample thick against $\lambda$, each face decays on its own as $B(x) \approx B_0\,e^{-x/\lambda}$.
Here $x$ is the distance from that face.

A useful demo plots $B(x)/B_0$ against $x$ for several temperatures in one figure.
The dip at the mid-plane fills in as $T$ rises.
The field equals $B_0$ at both faces.
The figure [`figures/temperature_lambda.png`](figures/temperature_lambda.png) plots $\lambda(T)$.
The figure [`figures/temperature_profiles.png`](figures/temperature_profiles.png) plots the profile family.

**Table 2** uses $\lambda_0 = 50$ nm, $T_c = 9.2$ K, and $d = 200$ nm.
The probe depth is $\lambda_0 = 50$ nm below each face, halfway to the mid-plane.
The quantity $1 - B(x)/B_0$ gives the screened fraction at the probe.

| $T$ (K) | $T/T_c$ | $\lambda(T)$ (nm) | $B/B_0$ at the probe | Screened fraction | $B(0)/B_0$ |
|---|---|---|---|---|---|
| 0.0 | 0.000 | 50.0 | 0.410 | 0.590 | 0.266 |
| 2.0 | 0.217 | 50.1 | 0.411 | 0.589 | 0.266 |
| 4.0 | 0.435 | 50.9 | 0.419 | 0.581 | 0.275 |
| 6.0 | 0.652 | 55.2 | 0.458 | 0.542 | 0.319 |
| 8.0 | 0.870 | 76.4 | 0.615 | 0.385 | 0.504 |
| 8.8 | 0.957 | 123.9 | 0.806 | 0.194 | 0.744 |
| 9.2 | 1.000 | diverges | 1.000 | 0.000 | 1.000 |

The table shows three regimes.
Below about $0.5\,T_c$ the depth barely moves.
The profile then differs little from the zero-T profile.
At 8.8 K the depth grows to about $2.5\,\lambda_0$.
The screened fraction then falls from 59 % to 19 %.
Close to $T_c$ the depth grows without bound.
The last row shows the limit with $B(x) = B_0$ everywhere.

At $T = 0$ the exact value 0.410 exceeds the half-space estimate $e^{-1} = 0.368$, because the field leaks in from the opposite face.

The columns cover the core of item 6.7.
A rise in $T$ raises the residual field at every interior point.

## 6. Numerical Design for Items 6.4 and 6.6

The app needs four guards, so a dragged slider cannot break the plots.

**Guard 1: clamp the temperature slider.**
The denominator $1 - (T/T_c)^4$ hits zero at $T = T_c$.
A scalar division then raises `ZeroDivisionError` in Python, and NumPy returns `inf` with a `RuntimeWarning`.
The demo clamps the top of the $T$ slider to $0.995\,T_c$.
At the clamp, $\lambda \approx 7.1\,\lambda_0$.
For the default $\lambda_0 = 50$ nm, $\lambda$ at the clamp is 354.9 nm.
The center field then gives $B(0)/B_0 = 0.9616$ at the default slab.

**Guard 2: branch at $T \ge T_c$.**
The function `lambda_of_temperature` in
[`temperature.py`](../src/magpenetration_simulation/temperature.py) does this branch.
Its signature is `lambda_of_temperature(t, lambda_0, t_c)`.
It returns `float("inf")` for `t >= t_c`, where the law has no value.

The app maps an infinite depth to the limit $B(x) = B_0$, and the panel can then draw the exact flat profile at $T_c$.
The flat branch is the normal-state limit of section 4.
In that limit the superfluid density is zero.
The London equation then reduces to $\nabla^2 \mathbf{B} = 0$, and the screening current is zero.

**Guard 3: stabilize the cosh.**
The argument $d/(2\lambda)$ can grow large.
The `cosh` function overflows past an argument of about 710.
The planned ranges produce a ratio below 50, and plain `cosh` is then safe.
The rewritten form is cheap insurance for future ranges, and every exponent is nonpositive for $x \ge 0$:

$$B(x) = B_0\,e^{-(d/2 - x)/\lambda}\,\frac{1 + e^{-2x/\lambda}}{1 + e^{-d/\lambda}} \qquad (x \ge 0)$$

The mirror form covers $x < 0$.

**Guard 4: one recompute per change.**
A change in $T$, $\lambda_0$, $T_c$, or $d$ recomputes $\lambda$ once, and the panel then rebuilds the profile on a fixed grid:

`x = np.linspace(-d/2, d/2, 401)`

With ipywidgets the code observes the `value` trait of each slider and redraws in one callback.
With Streamlit the slider reruns the script from the top, and the recompute happens by construction.

**Show $\lambda(T)$ beside $B(x)$.**
The design pairs the depth view with the profile view:

1. An annotation on the profile panel reports the current $\lambda$ in nm.
2. A second panel plots $\lambda$ against $T$, with a marker at the current value.
3. A family of profiles at the temperatures of Table 2 shows the trend in one frame.

Table 3 lists the ranges and defaults for the app of item 6.6:

| Control | Range | Default | Step |
|---|---|---|---|
| $T$ | 0 to $0.995\,T_c$ | 4.2 K | 0.05 K |
| $\lambda_0$ | 20 nm to 200 nm | 50 nm | 1 nm |
| $T_c$ | 1 K to 125 K | 9.2 K | 0.1 K |
| $d$ | 20 nm to 2000 nm | 200 nm | 10 nm |
| $B_0$ | 0 mT to 200 mT | 50 mT | 5 mT |

The defaults target niobium with $T_c = 9.2$ K.
When the profile panel plots $B(x)/B_0$, the value of $B_0$ cancels.
The $B_0$ slider then matters only for absolute-field plots.
The normalized view in the app does not change when the user moves the $B_0$ slider.
The $T_c$ range covers niobium and the cuprate scale.

## 7. Assumptions and Caveats

- **Local London limit.**
Pippard's nonlocal theory matters when the coherence length $\xi_0$ is large against $\lambda$ [4].
For niobium the two lengths are close, and small deviations can appear at low $T$.
- **Clean-metal limit.**
The fit assumes a clean metal.
Impurity scattering raises $\lambda_0$ and changes the curve of $\lambda(T)$ [4].
Thin films with many grain boundaries act like dirty samples.
- **Isotropic gap.**
The fit assumes an isotropic $s$-wave gap.
A $d$-wave gap contains nodes.
The low-$T$ correction then follows a power of $T$ rather than an exponential [5].
- **Mean-field validity.**
The two-fluid curve and the Ginzburg-Landau curve are both mean-field curves [7].
Fluctuations round the transition inside a narrow critical region.
For the classic superconductors the region is far narrower than the 0.05 K slider step of Table 3.
- **Geometry.**
The cosh solution covers the slab only.
For films thinner than $\lambda$, the Pearl length $\lambda^2/d$ replaces $\lambda$ as the transverse scale [5].
- **Field strength.**
The report assumes the Meissner state at low field.
Near $T_c$ the critical field approaches zero.
Even a weak field then destroys the state.
- **The empirical status of the fit.**
The two-fluid law is a fit to data.
BCS derives the superfluid fraction from the gap $\Delta(T)$ [3].
BCS and the fit agree on the exponent $-1/2$ near $T_c$.
They disagree at low $T$.
The report must not call the two-fluid formula exact.

## 8. What Was Done and Why

The issue body says "Explain your work."
The bullets below explain the work.

- We collected the two-fluid model from Gorter and Casimir [2] and the field equations from London and London [1].
- We derived $\lambda(T)$ by one substitution into $\lambda^2 = m/(\mu_0 n_s e^2)$.
Section 3 lists the approximations of the step.
- We checked the limits by expansion at both ends.
We connected the divergence at $T_c$ to mean-field theory and to the BCS gap [3], [5], [7].
- We picked the slab of item 6.2 as the demonstration case.
Table 2 comes from a plain-Python computation, and we cross-checked one row against the exponential estimate.
- We designed the four guards, because a user can park the $T$ slider on $T_c$ exactly.
- [`temperature.py`](../src/magpenetration_simulation/temperature.py) now contains the model.
Its `lambda_of_temperature(t, lambda_0, t_c)` function is the two-fluid law of section 2.
- The Temperature tab in [`app.py`](../src/magpenetration_simulation/app.py) uses that module.
The tab shows the sliders of Table 3 and the views of section 6.
- [`scripts/make_figures.py`](../scripts/make_figures.py) writes the two temperature figures:
[`figures/temperature_lambda.png`](figures/temperature_lambda.png) and
[`figures/temperature_profiles.png`](figures/temperature_profiles.png).
- The figure script checks the clamp values, and Guard 1 fixes the two numbers.
The check holds `lambda` at the clamp to 354.9 nm and `B(0)/B_0` at the clamp to 0.9616.
The script exits 1 when a value leaves its bound.

In short: a rise in $T$ grows $\lambda$, and the growth weakens the screening at every interior point.
The profile of $B(x)$ flattens toward $B_0$.
The screened fraction at the probe falls by a factor of three from 0 K to 8.8 K.

## 9. References

[1] F. London and H. London, "The electromagnetic equations of the supraconductor," Proceedings of the Royal Society A, vol. 149, pp. 71-88, 1935.

[2] C. J. Gorter and H. B. G. Casimir, "On supraconductivity I," Physica, vol. 1, pp. 306-320, 1934.

[3] J. Bardeen, L. N. Cooper, and J. R. Schrieffer, "Theory of superconductivity," Physical Review, vol. 108, no. 5, pp. 1175-1204, 1957.

[4] M. Tinkham, Introduction to Superconductivity, 2nd ed., McGraw-Hill, New York, 1996.

[5] C. P. Poole, H. A. Farach, R. J. Creswick, and R. Prozorov, Superconductivity, 2nd ed., Academic Press, Amsterdam, 2007.

[6] C. Kittel, Introduction to Solid State Physics, 8th ed., Wiley, Hoboken, 2005.

[7] V. L. Ginzburg and L. D. Landau, "On the theory of superconductivity," Zh. Eksp. Teor. Fiz., vol. 20, p. 1064, 1950.
