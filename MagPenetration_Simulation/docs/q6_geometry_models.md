# Geometry models for magnetic-field penetration (Q6, issue 3)

This report answers issue 3 with the usable geometry models.
Each entry gives the solution, the assumptions, the valid range, and the failure modes.
The last two sections guide the plots of issue 5 and the interface of issue 6.

## Shared framework

The local theory of London and London (1935) gives the electrodynamics of the supercurrent.
In the London gauge, the supercurrent follows from the vector potential:

$$\mathbf{J} = -\frac{\mathbf{A}}{\mu_0 \lambda^2}$$

The Maxwell curl equation turns this into a field equation inside the sample:

$$\nabla^2 \mathbf{B} = \frac{\mathbf{B}}{\lambda^2}$$

Here $\lambda$ is the London penetration depth with units of length.
Typical clean values are 16 nm for aluminum, 39 nm for niobium, and 150 nm for YBCO.
The YBCO value applies to the CuO2 planes.
These values come from the tables in Tinkham (1996) and Poole et al. (2007).

When the field keeps one Cartesian component, the vector equation reduces to one scalar equation for that component.
The scalar reduction gives the slab and the axial cylinder solutions below.
The sphere, the transverse fields, and all finite bodies keep the full vector structure.

### Symbols and units

- $\lambda$ is the London penetration depth in meters.
- $B_0$ is the applied field magnitude in tesla.
- $d$ is the slab or film thickness in meters.
- $R$ is the cylinder or sphere radius in meters.
- $x$ is the distance from the midplane in meters.
- $s$ is the depth from the nearest face in meters.
- $I_0$ and $I_1$ denote modified Bessel functions of orders 0 and 1.
- $i_1$ is the modified spherical Bessel function of order 1, dimensionless.
- $m$ is a dipole moment in A $\mathrm{m^2}$.
- $\Lambda$ is the Pearl length in meters.

### Limits of the local London form

The local form takes the current density at one point from the potential at that point.
The field must vary over scales well past the coherence length $\xi$.
The condition reads $q\xi \ll 1$, with $q = 1/\lambda$.
Pippard (1953) found the nonlocal correction in the opposite regime.

High-$\kappa$ type-II materials lie deep in the local regime (Abrikosov, 1957).
Aluminum with $\lambda = 16$ nm and $\xi = 1.6$ micrometers lies far outside it.
Niobium lies near the border between the two types.
Anisotropic materials replace the scalar $\lambda$ by a tensor.
Thin films add the Pearl crossover of the film section.
For one vortex, the London equation gives the decay profile $K_0(r/\lambda)$.
Above the first critical field, vortices populate the sample (Abrikosov, 1957).
The uniform-field profiles apply to the Meissner state only.

## Model 1: one-dimensional slab

The slab occupies $-d/2 \le x \le d/2$, with the midplane at $x = 0$.
The applied field $B_0$ runs parallel to the broad faces.
With one component only, the field equation reads:

$$\frac{d^2 B}{dx^2} = \frac{B}{\lambda^2}$$

The even solution that takes the value $B_0$ at both faces is:

$$B(x) = B_0\,\frac{\cosh(x/\lambda)}{\cosh(d/(2\lambda))}$$

The field at both faces remains $B_0$ at every $d/\lambda$.
For $x$ near the face at $d/2$, write the depth as $s = d/2 - x$.

### Thick slab, $d \gg \lambda$

A surface layer of thickness of order $\lambda$ holds the field:

$$B(x) \approx B_0\,e^{-s/\lambda}$$

For $d = 20\lambda$ the midplane value is $9 \times 10^{-5} B_0$.

### Thin slab, $d \ll \lambda$

The sample screens weakly, and the profile is nearly flat:

$$B(x) \approx B_0 \left(1 + \frac{x^2/2 - d^2/8}{\lambda^2}\right)$$

The field at both faces remains $B_0$, and the midplane dips below it.
The relative dip is $d^2/(8\lambda^2)$.
Integration over the thickness gives the volume mean:

$$\langle B \rangle = B_0\,\frac{2\lambda}{d}\,\tanh\!\left(\frac{d}{2\lambda}\right)$$

For $d \ll \lambda$ the mean is $B_0 (1 - d^2/(12\lambda^2))$.
The screening stays under one percent while $d < 0.35\lambda$.

### Caveats

- The field must remain parallel to the faces.
- A tilted or perpendicular field needs the transverse treatment.
- The faces must be much wider than $\lambda$.
- The model takes the local London limit.

## Model 2: infinite cylinder with an axial field

The cylinder has radius $R$ and infinite length.
The applied field $B_0$ runs along the axis, and the field keeps one component $B_z(r)$.
In cylindrical coordinates the London equation becomes the modified Bessel equation of order zero:

$$\frac{1}{r}\frac{d}{dr}\!\left(r\frac{dB}{dr}\right) = \frac{B}{\lambda^2}$$

The independent solutions are $I_0(r/\lambda)$ and $K_0(r/\lambda)$.
$K_0$ diverges at the axis.
Only $I_0$ remains finite there, and the boundary condition at $r = 0$ picks it.
With no transport current, the London theory sets the boundary value at the applied value:

$$B(R) = B_0$$

The solution reads:

$$B(r) = B_0\,\frac{I_0(r/\lambda)}{I_0(R/\lambda)}$$

### Thick cylinder, $R \gg \lambda$

Near the wall of the cylinder the profile approaches:

$$B(r) \approx B_0\sqrt{\frac{R}{r}}\;e^{-s/\lambda}$$

The curvature factor $\sqrt{R/r}$ departs from 1 by about $s/(2R)$.
For $R = 20\lambda$ the factor lifts the $s = \lambda$ value by 2.6 percent.
Further in, the slab exponential applies.

### Thin cylinder, $R \ll \lambda$

The sample barely screens the field.
The value at the axis is $B_0/I_0(R/\lambda) \approx B_0 (1 - R^2/(4\lambda^2))$.
The radial mean of the solution gives:

$$\langle B \rangle = B_0\,\frac{2\lambda}{R}\,\frac{I_1(R/\lambda)}{I_0(R/\lambda)} \approx B_0 \left(1 - \frac{R^2}{8\lambda^2}\right)$$

The flux reduction scales like $R^2$, weaker than the thin-slab factor $d^2/(12\lambda^2)$ by the ratio $2/3$.

### Scope and caveats

- This model covers an axial field and an infinite length only.
- In a transverse field the angular part changes, and the radial functions come from the order-one equation.
- In the ideal limit the field at the equator rises to $2 B_0$.
- Real cylinders have ends.
  End corrections enter when the length $L$ is not many times $R$ (Prozorov, 2021).

## Model 3: sphere

A sphere of radius $R$ sits in a uniform field $B_0$ along the polar axis.
Here the vector equation does not reduce, and the field keeps both components $B_r$ and $B_\theta$.
The London theory separates the problem with a toroidal potential $\psi(r)\sin\theta$ in the $\phi$ direction.
The radial function obeys:

$$\psi'' + \frac{2}{r}\psi' - \frac{2}{r^2}\psi = \frac{\psi}{\lambda^2}$$

The regular solution uses the modified spherical Bessel function of order 1:

$$\psi(r) = \frac{3 B_0 R}{2\,\sinh(R/\lambda)}\,i_1\!\left(\frac{r}{\lambda}\right), \qquad i_1(x) = \frac{x \cosh x - \sinh x}{x^2}$$

Continuity of both components at the boundary of the sphere fixes that prefactor.
The interior field reads:

$$B_r = \frac{2\psi(r)}{r}\,\cos\theta, \qquad B_\theta = -\frac{\psi(r) + r\psi'}{r}\,\sin\theta$$

Outside the sphere, the applied field joins a dipole with moment:

$$m = -\frac{2\pi R^3 B_0}{\mu_0}\,S, \qquad S = 1 - \frac{3\lambda}{R}\,\coth\!\left(\frac{R}{\lambda}\right) + \frac{3\lambda^2}{R^2}$$

In the ideal limit $S$ approaches 1.
Then $m$ approaches the perfect diamagnetic value (London, 1950).
For $R \ll \lambda$, $S$ falls to about $R^2/(15\lambda^2)$, and the sphere barely screens.

The same solution gives two useful scalar values:

$$B(0) = B_0\,\frac{R/\lambda}{\sinh(R/\lambda)}, \qquad \langle B_z \rangle = B_0 \left(\frac{3\lambda}{R}\,\coth\!\frac{R}{\lambda} - \frac{3\lambda^2}{R^2}\right)$$

$B(0)$ is the field at the center.
$\langle B_z \rangle$ is the mean along the polar axis.
For $R = 5\lambda$ the center value is $0.067 B_0$.

### Demagnetization caveat

The boundary value varies strongly with latitude.
At the equator it comes to $\frac{3}{2} B_0 (1 - i_1/\sinh)$ and approaches $\frac{3}{2} B_0$ in the ideal limit.
At the poles it approaches zero in the same limit.
The applied field is then not the boundary value, and the geometry redistributes the field.
Only ellipsoids have one exact demagnetizing factor per axis (Osborn, 1945).
The sphere has $N = 1/3$.

### Scalar-form caveat

Some notes reduce the sphere to a radial magnitude with zero angular momentum:

$$B_{\mathrm{scal}}(r) = B_0\,\frac{R}{r}\,\frac{\sinh(r/\lambda)}{\sinh(R/\lambda)}$$

That scalar form solves the spherically symmetric reduction only.
It gives the correct value at the center.
It misses the angular part, and the field at the equator departs from it.
Use it as a rough magnitude guide only, and use the vector form for real numbers.

## Model 4: finite rectangle, box, and disk

Finite samples have no closed-form London solution in general.
Inside the sample, each Cartesian component still obeys the scalar equation.
The hard part is the boundary.
The matching conditions couple all components, and they need the outside field.
The outside field solves vacuum magnetostatics with the boundary currents.

A working numerical recipe follows:

- Put the grid over the sample and its surroundings.
- Write the London equation per Cartesian component with the five-point stencil.
  The interior update reads $B = (B_E + B_W + B_N + B_S)/(4 - h^2/\lambda^2)$ in a 2D grid of step $h$.
- Sweep until the residual falls under the chosen tolerance.
  Successive over-relaxation accelerates the sweep.
- Take the boundary data of the interior from the continuity conditions.

Near an edge, the Meissner current turns, and the 1D and cylindrical profiles miss that region.
Keep the closed forms several $\lambda$ away from any edge.
Corners are sharper.
The corner field rises sharply, and the value depends on the grid.
A reported corner peak is a property of the grid, not the sample.
Round the corners in the model, or state the cutoff.

Only ellipsoids have exact demagnetizing factors (Osborn, 1945).
Rectangles, boxes, and finite cylinders use tabulated or computed factors, and the factor depends on the aspect ratio.
For finite cylinders in an axial field, Prozorov (2021) gives the susceptibility with the end corrections.

The closed forms of Models 1 and 2 fail near edges and ends for two reasons.
They assume straight field lines and full translation or rotation symmetry.
A finite body distorts the outside field, and the boundary value departs from $B_0$.

## Model 5: thin films and the Pearl regime

Pearl (1964) set up the film problem.
The characteristic length is:

$$\Lambda = \frac{2\lambda^2}{d}$$

Here $d$ is the film thickness.
The crossover occurs at $d = 2\lambda$, where $\Lambda$ equals $\lambda$.
For $d < 2\lambda$, $\Lambda$ exceeds $\lambda$, and the in-plane fields and currents spread over $\Lambda$ instead of $\lambda$.
The London penetration depth alone is then not enough for in-plane profiles.
Pearl (1964) also gives the current around one fluxoid in the film.
The current decays only like $1/r$ out to the Pearl length.

For a parallel field, the 1D slab profile remains exact at every thickness.
For a perpendicular field, the demagnetizing factor approaches 1.
Then the mean field in the film barely drops, and the currents spread over $\Lambda$.

## Comparison table

| Model | Geometry | Governing equation | Closed-form solution | Main assumptions | Valid regime | Main caveat |
|---|---|---|---|---|---|---|
| 1D slab | infinite plate | $B'' = B/\lambda^2$ | $\cosh(x/\lambda)/\cosh(d/2\lambda)$ | field parallel | every $d/\lambda$ | wide faces |
| Axial infinite cylinder | radius $R$, endless | $B'' + B'/r = B/\lambda^2$ | $B_0 I_0(r/\lambda)/I_0(R/\lambda)$ | field axial | every $R/\lambda$ | end effects |
| Sphere | radius $R$ | vector London form | $i_1$ radial factor | local London | every $R/\lambda$ | boundary value varies |
| Box, rectangle, disk | finite body | vector London form | none | local London | numeric only | needs the outside field |
| Thin film | film thickness $d$ | Pearl model | none | thin limit | $d < 2\lambda$ | $\Lambda = 2\lambda^2/d$ replaces $\lambda$ |

## Guidance for the plots of issue 5

- Draw the slab profile with the full cosh form at every $d/\lambda$.
  The exponential $B_0 e^{-s/\lambda}$ is a dashed thick-limit reference only.
- Draw the axial cylinder profile at several $R/\lambda$ values on its own panel.
- Map $B(x)$ over the slab thickness with `imshow` or `contourf`.
  With a parallel field, the field outside the plate remains $B_0$, and the map omits an exterior computation.
- Label the depth axis in nm and in units of $\lambda(T)$.
- Label the field axis or the colorbar as $B/B_0$, and print the $B_0$ value in tesla.
- Normalize lengths with $\lambda(T) = \lambda_0/\sqrt{1 - (T/T_c)^4}$ from the two-fluid model (Poole et al., 2007).
- Provide a log color option, since interior values fall several decades for thick plates.
- Keep boundary values above $B_0$ visible.
  Demagnetization raises the field just inside the boundary above $B_0$.
  The sphere equator gives $1.16 B_0$ for $R = 2.9\lambda$.
- Do not compare geometries on one depth axis.
  Plot every curve against the depth from the nearest face, in units of $\lambda$.
  State the boundary value of every curve in the caption.
- Mark the valid zone on rectangle maps.
  The 1D midplane profile applies several $\lambda$ away from an edge.

## Guidance for the interface of issue 6

- Provide sliders for $T$, $\lambda_0$, $B_0$, the geometry, and $d$ or $R$.
- Update $\lambda$ from the two-fluid form above, and show the value in nm.
- Warn when $d < 2\lambda$, and show $\Lambda = 2\lambda^2/d$.
- Warn when the user picks a rectangle or box geometry.
  State the 1D limit of that map in the caption.
- Keep units on every axis, and keep the normalization $B/B_0$ on every plot.
- For geometry comparisons, plot $B/B_0$ against $s/\lambda$, and give the boundary value per curve.

## What was done and why

Issue 3 calls for a list of usable geometry models with caveats.
The list below states each family and its limits.

- The 1D slab gives the closed form of issue 6.2.
  It needs a parallel field and faces wider than $\lambda$.
- The axial cylinder gives the Bessel form of issue 6.3.
  It remains exact at every radius.
  Transverse fields and finite lengths need extra work.
- The sphere has an exact London solution, and this report states it in full.
  The demagnetization caveat belongs with every sphere curve.
- Finite rectangles, boxes, and disks need the numeric path.
  Issue 6.5 can map the midplane profile, and the map must mark the edge zone where the profile fails.
- Thin films cross into the Pearl regime below $d = 2\lambda$.
  The interface must switch the length scale to $\Lambda$ there.

Every closed form follows from the local London equation.
The numeric path inherits the same limits.
Results outside the local regime need the Pippard kernel, an anisotropic $\lambda$, or a thin-film treatment.

## References

F. London and H. London published "The electromagnetic equations of the supra-conductor" in 1935.
It appeared in Proc. R. Soc. A 149, 71.

F. London wrote the monograph "Superfluids" in 1950.
It appeared at Wiley, New York, with the subtitle "Macroscopic Theory of Superconductivity."

A. A. Abrikosov published "On the magnetic properties of superconductors of the second group" in 1957.
It appeared in Sov. Phys. JETP 5, 1174.

J. A. Osborn published "Demagnetizing factors of the general ellipsoid" in 1945.
It appeared in Phys. Rev. 67, 351.

A. B. Pippard proposed the nonlocal modification in 1953.
The exact title reads "An experimental and theoretical study of the relation between magnetic field and current in a superconductor".
It appeared in Proc. R. Soc. A 216, 547.

J. Pearl published "Current distribution in superconducting films carrying quantized fluxoids" in 1964.
It appeared in Appl. Phys. Lett. 5, 65.

M. Tinkham published "Introduction to Superconductivity" in 1996.
It appeared in the 2nd edition at McGraw-Hill, New York.

C. P. Poole, H. A. Farach, R. J. Creswick, and R. Prozorov published "Superconductivity" in 2007.
It appeared in the 2nd edition at Academic Press, Amsterdam.

R. Prozorov published "Meissner-London susceptibility of superconducting right circular cylinders in an axial magnetic field" in 2021.
It appeared in Phys. Rev. Applied 16, 024014.
