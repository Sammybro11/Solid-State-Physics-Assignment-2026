# Magnetic-field penetration: London-equation slab model (Q6)

Simulates how a static magnetic field decays inside a superconducting slab,
using the London equation

$$\lambda^2 \frac{d^2B}{dx^2} = B, \qquad B(\pm d/2) = B_0,$$

whose exact solution is $B(x) = B_0\cosh(x/\lambda)/\cosh(d/2\lambda)$.
The code checks three numerical methods against that analytic answer.
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

The page has two tabs.

**Slab** holds the slab controls, the field profile, and the error charts.

- Sliders for slab thickness $d/\lambda$, applied field $B_0$, penetration depth
  $\lambda_0$ and grid points $N$. The graphs update while you drag a slider.
- **Field profile** and **pointwise error** graphs, each with its own selector for
  which methods to show. A **Scaled / Physical** switch sets the axes to
  $x/\lambda$ and $B/B_0$, or to nm and mT. In scaled units the profile depends
  only on $d/\lambda$.
- Cards for the center field, slab thickness, surface screening current, and the
  grid resolution. Further cards give the maximum absolute error and the maximum
  relative error of each method.

**Temperature** holds the two-fluid demo.

- Sliders for $T$, $T_c$, $\lambda_0$, $d$, and $B_0$.
- The $\lambda(T)$ curve, with a marker at the current $T$.
- The profile family at several temperatures, in one frame.
- Guard metrics, and a warning when $T$ comes near $T_c$.
- A normal-state demo, which shows the flat field at $T \ge T_c$.

Both tabs share the light and dark themes (switch at the top right).

## Smoke check and figures

```bash
uv run python scripts/make_figures.py
```

The script runs every method at the defaults below.
It prints the maximum error of each numerical method and the two temperature guard values.
It writes four figures: `docs/figures/slab_profile.png`, `docs/figures/slab_error.png`,
`docs/figures/temperature_lambda.png`, and `docs/figures/temperature_profiles.png`.
It also prints the clamp `lambda` and the clamp `B(0)/B_0`.
The script exits with status 1 if any error is above $10^{-3}$ or a temperature guard
leaves its bound. Expected output:

```
Method                 max |dB| (T)   max |dB| / B0   result
Finite difference         1.410e-07       2.819e-06   pass
Relaxation (SOR)          1.411e-07       2.822e-06   pass
Shooting                  7.430e-12       1.486e-10   pass

Temperature defaults: lambda0 = 50 nm, Tc = 9.2 K, d = 200 nm, B0 = 50 mT, T = 4.2 K
lambda at clamp: 354.9 nm (expected 355 +/- 3.55 nm)
B(0)/B0 at clamp: 0.9616 (expected 0.96 +/- 0.005)
...
All three methods within the 0.001 target.
Temperature guards within bound.
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

## Defaults

| Quantity | Value |
|---|---|
| Penetration depth $\lambda_0$ | 50 nm |
| Slab thickness $d$ | $5\lambda$ |
| Applied field $B_0$ | 50 mT |
| Grid points $N$ | 401 |
| Acceptance | max relative error below $10^{-3}$ for every numerical method |

The slab functions take $\lambda$ directly.
The `temperature` module and the Temperature tab supply $\lambda(T)$.
The theory is in
[`docs/q6_temperature_theory.md`](docs/q6_temperature_theory.md).
