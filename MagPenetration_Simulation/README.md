# Magnetic-field penetration: London-equation slab model (Q6)

Simulates how a static magnetic field decays inside a superconducting slab,
using the London equation

$$\lambda^2 \frac{d^2B}{dx^2} = B, \qquad B(\pm d/2) = B_0,$$

whose exact solution is $B(x) = B_0\cosh(x/\lambda)/\cosh(d/2\lambda)$.
Three independent numerical methods are checked against that analytic answer.
The techniques and the reasons for each choice are in
[`docs/q6_simulation_techniques.md`](docs/q6_simulation_techniques.md).

## Setup

The project uses [uv](https://docs.astral.sh/uv/) and Python 3.14 (uv downloads
the interpreter if needed). From this folder:

```bash
uv sync
```

## Interactive app

```bash
uv run streamlit run src/magpenetration_simulation/app.py
```

The page opens at <http://localhost:8501>. Run it from this folder so Streamlit
picks up the theme in `.streamlit/config.toml`.

- Sliders for slab thickness $d/\lambda$, applied field $B_0$, penetration depth
  $\lambda_0$ and grid points $N$. The graphs update while a slider is dragged.
- **Field profile** and **pointwise error** graphs, each with its own selector for
  which methods to show, and a **Scaled / Physical** switch ($x/\lambda$ and $B/B_0$,
  or nm and mT). In scaled units the profile depends only on $d/\lambda$.
- Cards for the centre field, slab thickness, surface screening current, grid
  resolution, and each method's maximum absolute and relative error.
- Light and dark themes (switch at the top right).

## Smoke check and figures

```bash
uv run python scripts/make_figures.py
```

Runs every method at the defaults below, prints the maximum error of each
numerical method, writes `docs/figures/slab_profile.png` and
`docs/figures/slab_error.png`, and exits with status 1 if any error is above
$10^{-3}$. Expected output:

```
Method                 max |dB| (T)   max |dB| / B0   result
Finite difference         1.410e-07       2.819e-06   pass
Relaxation (SOR)          1.411e-07       2.822e-06   pass
Shooting                  7.430e-12       1.486e-10   pass
```

## Methods

All functions live in `src/magpenetration_simulation/slab.py`, take SI units
(m, T) and use the same grid of `n_points` (odd) from $-d/2$ to $d/2$.

| Function | Method |
|---|---|
| `slab_field_analytic(x, d, lam, b0)` | exact cosh solution, written so it cannot overflow |
| `slab_field_finite_difference(d, lam, b0, n_points)` | central differences, tridiagonal solve with `scipy.linalg.solve_banded` |
| `slab_field_relaxation(d, lam, b0, n_points, tol=1e-10, max_sweeps=10000)` | red-black successive over-relaxation (SOR) |
| `slab_field_shooting(d, lam, b0, n_points=401, tol=1e-10)` | linear shooting from the slab centre with `scipy.integrate.solve_ivp` |

The numerical functions return `(x, b)`.

## Defaults

| Quantity | Value |
|---|---|
| Penetration depth $\lambda_0$ | 50 nm |
| Slab thickness $d$ | $5\lambda$ |
| Applied field $B_0$ | 50 mT |
| Grid points $N$ | 401 |
| Acceptance | max relative error below $10^{-3}$ for every numerical method |

Temperature dependence of $\lambda$ is handled separately (PR #10); here
$\lambda_0$ is used directly as $\lambda$.
