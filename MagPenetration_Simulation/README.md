# Magnetic-field penetration: London-equation models (Q6)

This app shows how a static magnetic field decays inside a superconductor.
It uses the London equation

$$\lambda^2 \frac{d^2B}{dx^2} = B, \qquad B(\pm d/2) = B_0,$$

whose exact slab solution is $B(x) = B_0\cosh(x/\lambda)/\cosh(d/2\lambda)$.
The app shows the slab, the axial cylinder, the sphere, a rectangle, and a thin film.
It also compares the two closed forms.
The code checks three numerical slab methods against that analytic answer.
The techniques and the reasons for each choice are in
[`docs/q6_simulation_techniques.md`](docs/q6_simulation_techniques.md).

## Setup

The project uses [uv](https://docs.astral.sh/uv/) and Python 3.14.
uv downloads the interpreter when it is missing.
Run this command from the project folder:

```bash
uv sync
```

## Interactive app

```bash
uv run streamlit run src/magpenetration_simulation/app.py
```

The page opens at <http://localhost:8501>. Run it from this folder so Streamlit
picks up the theme in `.streamlit/config.toml`.

All UI ranges, defaults, steps, chart sizes, families, the temperature clamp, and
the error target are in `src/magpenetration_simulation/ui_config.py`.

### Navigation

The app has eight pages.
The top navigation bar shows them in this order:

`Theory | Methods | Slab | Cylinder | Sphere | Finite body | Film | Combined`

The large title at the top is a button.
Click the title to return to **Theory**.
Only the active page computes.

### Sidebar

The sidebar provides the material and field controls.
The active page reads these values.

- Critical temperature $T_c$ in K.
- Temperature $T$ in K. The control stops at $0.995\,T_c$.
- Penetration depth at zero temperature $\lambda_0$ in nm.
- Applied field $B_0$ in mT.
- A **Normal state** toggle, which sends $\lambda$ to infinity.
- The light and dark theme switch.

Temperature is no longer a tab.
The sidebar holds $T$, and the geometry pages use $\lambda(T)$.

### Pages

| Page | Content |
|---|---|
| Theory | The London equation, the penetration depth, the two-fluid law, the geometry families, the assumptions, and the limits. |
| Methods | Finite difference, SOR, and shooting against the analytic profile, with error cards. |
| Slab | The analytic slab profile, the temperature family, the $\lambda(T)$ chart, and the Pearl warning. |
| Cylinder | The axial Bessel profile, the $R/\lambda$ family, and the temperature family. |
| Sphere | The exact vector solution, the equatorial and surface profiles, and the center, mean, screening, equator, and dipole cards. |
| Finite body | The 2D scalar field map of a rectangle, the valid-zone mask, and the center, mean, residual, and sweep cards. |
| Film | The Pearl length, the regime, the parallel profile, the temperature family, and the center and mean cards. |
| Combined | The slab and the cylinder on one depth axis, with six metric cards. |

**Methods** defines its own controls: $d/\lambda$, $N$, the method selector, and
the scaled or physical units.
The error chart also has a signed or log scale.
The target is a maximum relative error below $10^{-3}$.

**Slab** gives cards for the center field, the mean field, and the screened
fraction.
It warns when $d/\lambda$ drops below 2, where the Pearl length replaces
$\lambda$.

**Cylinder** gives the same three cards.
It adds a near-wall reference curve for $R/\lambda \ge 10$.

**Sphere** shows the exact vector London solution.
It plots the equatorial radial profile and the boundary field against latitude.
It gives six cards: the center field, the axial mean, the screening factor $S$,
the equator field, and the dipole moment.
This page does not model the normal state.
A 2D field-line view is a later slice.

**Finite body** solves the 2D scalar equation on a rectangle.
It holds $B = B_0$ on all four edges and draws the field as a heatmap.
A dashed mask marks the valid zone, at least $3\lambda$ from every edge.
The cards give the center field, the mean field, the residual, and the sweeps.
The outside field, demagnetization, and edge-current turning stay outside this model.
Box and disk geometries are later slices.

**Film** shows the Pearl length $\Lambda = 2\lambda^2/d$ and the regime.
The regime is `pearl` below $d = 2\lambda$, `crossover` at that thickness, and
`thick` above it.
It plots the exact parallel-field profile and the temperature family.
The cards give the center field and the mean field.
This page does not solve the perpendicular field, its edge currents, or
demagnetization.

**Combined** plots both geometries against the depth from the nearest face,
$s/\lambda$.
A later slice will add the 2D field-line view to that page.

The theme switch applies to each page.

## Smoke check and figures

```bash
uv run python scripts/make_figures.py
```

The script runs every method at the defaults below.
It prints the maximum error of each numerical method.
It prints the temperature guard values and the geometry checks.
It writes nine figures: `docs/figures/slab_profile.png`,
`docs/figures/slab_error.png`, `docs/figures/temperature_lambda.png`,
`docs/figures/temperature_profiles.png`, `docs/figures/geometry_cylinder.png`,
`docs/figures/geometry_comparison.png`, `docs/figures/geometry_sphere.png`,
`docs/figures/geometry_finite_body.png`, and `docs/figures/geometry_film.png`.
It also prints the clamp `lambda` and the clamp `B(0)/B_0`.
It prints the sphere, finite-body, and film checks after the cylinder checks.
The script exits with status 1 if any error is above $10^{-3}$ or a bound
fails. Expected output:

```
Method                 max |dB| (T)   max |dB| / B0   result
Finite difference         1.410e-07       2.819e-06   pass
Relaxation (SOR)          1.411e-07       2.822e-06   pass
Shooting                  7.430e-12       1.486e-10   pass

Temperature defaults: lambda0 = 50 nm, Tc = 9.2 K, d = 200 nm, B0 = 50 mT, T = 4.2 K
lambda at clamp: 354.9 nm (expected 355 +/- 3.55 nm)
B(0)/B0 at clamp: 0.9616 (expected 0.96 +/- 0.005)

Geometry checks: lambda = 50 nm, B0 = 50 mT, N = 401
wall B(R)/B0 - 1, max over the R/lambda family   0   pass
cylinder centre, R = 5 lambda          0.0367109   pass
cylinder mean,   R = 5 lambda          0.357353    pass
cylinder centre, R = 20 lambda        2.29577e-8   pass
cylinder mean,   R = 20 lambda         0.0974671   pass
slab centre,     d = 5 lambda          0.163071    pass
slab mean,       d = 5 lambda          0.394646    pass
...
Sphere/finite-body/film defaults: lambda = 50 nm, B0 = 50 mT, rectangle N = 81
Check                                                        value            expected   result
sphere centre, R = 5 lambda                          0.06738252915       0.06738252915   pass
sphere equator, R = 2.9 lambda                         1.157976253         1.157976253   pass
sphere screening factor, R = 50 lambda                      0.9412              0.9412   pass
finite body centre, 20 x 5 lambda                     0.1631508352        0.1631508352   pass
finite body centre, 5 x 5 lambda                      0.2573616334        0.2573616334   pass
finite body centre, 20 x 1 lambda                      0.886820218         0.886820218   pass
wide rectangle vs slab analytic, d = 5 lambda      0.0004881504792                   0   pass
film Lambda/lambda, d = 2 lambda                                 1                   1   pass
film regime, d = 2 lambda                                crossover           crossover   pass
film Lambda/lambda, d = 1 lambda                                 2                   2   pass
film regime, d = 1 lambda                                    pearl               pearl   pass
film Lambda/lambda, d = 10 lambda                              0.2                 0.2   pass
film regime, d = 10 lambda                                   thick               thick   pass
film parallel centre, d = 5 lambda                    0.1630712319        0.1630712319   pass
film parallel mean, d = 5 lambda                      0.3946457193        0.3946457193   pass
...
All three methods within the 0.001 target.
Temperature guards within bound.
Geometry checks within bound.
```

## Methods

The slab functions live in `src/magpenetration_simulation/slab.py`.
They take SI units (m, T) and use the same grid of `n_points` (odd) from
$-d/2$ to $d/2$.

| Function | Method |
|---|---|
| `slab_field_analytic(x, d, lam, b0)` | exact cosh solution, written so it cannot overflow |
| `slab_field_finite_difference(d, lam, b0, n_points)` | central differences, tridiagonal solve with `scipy.linalg.solve_banded` |
| `slab_field_relaxation(d, lam, b0, n_points, tol=1e-10, max_sweeps=10000)` | red-black successive over-relaxation (SOR) |
| `slab_field_shooting(d, lam, b0, n_points=401, tol=1e-10)` | linear shooting from the slab center with `scipy.integrate.solve_ivp` |

The numerical functions return `(x, b)`.

The cylinder functions live in `src/magpenetration_simulation/cylinder.py`.

| Function | Method |
|---|---|
| `cylinder_field_analytic(r, radius, lam, b0)` | exact Bessel solution $B_0 I_0(r/\lambda)/I_0(R/\lambda)$ |
| `cylinder_field_profile(radius, lam, b0, n_points=401)` | the same solution on a grid from the axis to the wall |
| `cylinder_center_field(radius, lam, b0)` | the axis value $B_0/I_0(R/\lambda)$ |
| `cylinder_mean_field(radius, lam, b0)` | the volume mean $B_0 (2\lambda/R) I_1(R/\lambda)/I_0(R/\lambda)$ |

The cylinder uses `scipy.special.i0e` and `scipy.special.i1e`.
These scaled forms cannot overflow at a large $R/\lambda$.

The sphere functions live in `src/magpenetration_simulation/sphere.py`.

| Function | Method |
|---|---|
| `sphere_potential(r, radius, lam, b0)` | the toroidal potential $\psi(r) = 3B_0R\,i_1(r/\lambda)/(2\sinh(R/\lambda))$ |
| `sphere_field(r, theta, radius, lam, b0)` | the interior components $B_r$ and $B_\theta$ |
| `sphere_center_field(radius, lam, b0)` | the center value $B_0 (R/\lambda)/\sinh(R/\lambda)$ |
| `sphere_mean_field(radius, lam, b0)` | the axial mean $B_0 (1 - S)$ |
| `sphere_screening_factor(radius, lam, b0)` | $S = 1 - 3(\lambda/R)\coth(R/\lambda) + 3\lambda^2/R^2$ |
| `sphere_dipole_moment(radius, lam, b0)` | the exterior dipole $-(2\pi R^3 B_0/\mu_0) S$ |
| `sphere_equator_field(radius, lam, b0)` | the equator value, which approaches $1.5 B_0$ |

The module evaluates $i_1$ and its derivative in a scaled form.
The scaled form cannot overflow at a large $R/\lambda$.

The finite-body functions live in
`src/magpenetration_simulation/finite_body.py`.

| Function | Method |
|---|---|
| `rect_field_finite_difference(width, height, lam, b0, nx=81, ny=81, ...)` | 2D scalar solve with a five-point stencil and red-black SOR, and returns a `RectSolution` |
| `rect_valid_mask(x, y, width, height, lam, margin=3.0)` | the mask of points at least `margin * lam` from every edge |
| `rect_edge_distance(x, y, width, height)` | the distance from each point to the nearest edge |
| `rect_center_field(solution)` | the field at the node nearest the center |
| `rect_1d_limit(width, height, lam, b0, n_points=401)` | the 1D slab profile across the height |

`RectSolution` holds `x`, `y`, `b`, `residual`, `sweeps`, `hx`, `hy`, and
`converged`.
Each edge stays at $B_0$.
The update divides by $2/h_x^2 + 2/h_y^2 + 1/\lambda^2$.
The plus sign is correct, because the minus sign makes the sweep diverge.

The film functions live in `src/magpenetration_simulation/film.py`.

| Function | Method |
|---|---|
| `pearl_length(lam, d)` | the Pearl length $\Lambda = 2\lambda^2/d$ |
| `film_crossover(lam)` | the crossover thickness $2\lambda$ |
| `pearl_ratio(lam, d)` | the ratio $\Lambda/\lambda = 2\lambda/d$ |
| `inplane_screening_length(lam, d)` | the in-plane length $\max(\Lambda, \lambda)$ |
| `film_regime(lam, d)` | `pearl`, `crossover`, or `thick` |
| `film_parallel_profile(x, d, lam, b0)` | the exact 1D slab profile at any thickness |
| `film_parallel_center_field(d, lam, b0)` | the center value $B_0/\cosh(d/2\lambda)$ |
| `film_parallel_mean_field(d, lam, b0)` | the mean $B_0 (2\lambda/d)\tanh(d/2\lambda)$ |

## Defaults

| Quantity | Value |
|---|---|
| Penetration depth $\lambda_0$ | 50 nm |
| Critical temperature $T_c$ | 9.2 K |
| Temperature $T$ | 4.2 K |
| Applied field $B_0$ | 50 mT |
| Slab thickness $d$ | $5\lambda$ |
| Cylinder radius $R$ | $5\lambda$ |
| Sphere radius $R$ | $5\lambda$ |
| Finite-body size | $5\lambda \times 5\lambda$ |
| Finite-body grid $N$ | 81 |
| Film thickness $d$ | $1\lambda$ |
| Grid points $N$ | 401 |
| Temperature clamp | $0.995\,T_c$ |
| Acceptance | max relative error below $10^{-3}$ for every numerical method |

The geometry functions take $\lambda$ directly.
The sidebar temperature controls and the geometry pages supply $\lambda(T)$.
The theory is in
[`docs/q6_temperature_theory.md`](docs/q6_temperature_theory.md).
