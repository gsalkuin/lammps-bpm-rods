# Magnetic Helical Rods

**Sec. 3.3** — Based on [Sano (2022)](https://doi.org/10.1088/1751-8121/ac4de2)

A helical rod magnetized along its screw axis contracts under dipole–dipole interactions.
An axial-gradient magnetic field stretches the helix, and quasi-static loading–unloading cycles produce mechanical hysteresis.

## Parameters

| Parameter | Value |
|-----------|-------|
| Centerline length $L$ | 103 mm |
| Rod diameter $d$ | 2 mm |
| Young's modulus $E$ | 1.16 MPa |
| Poisson's ratio $\nu$ | 0.5 |
| Density $\rho$ | 1.23 g/cm$^3$ |
| Number of turns | 1.64 |
| Normalized pitch $\tilde{\psi}$ | 0.96 |
| Number of bonds $N$ | 103 (bond length $\ell_0 = 1$ mm) |
| Dipole strength $\gamma$ | varies |
| Field strength $\lambda_m$ | varies |

## Setup

The helix is initialized in its stress-free configuration and the topmost atom at $z \approx 0$ is clamped.
The centerline is parametrized by arclength $s$:

$$\mathbf{r}(s) = \langle R\cos(Ks),\; R\sin(Ks),\; s\cos\psi \rangle,$$

where $K = R^{-1}\sin\psi$ and $\tilde{\psi} = \psi/(\pi/2) = 0.96$.

All dipoles are initialized along the $+z$ screw axis.
Dipole–dipole interactions are controlled by $\gamma = \mu_0 M^2 / E$.
Since the problem is governed by dimensionless parameters, the value of $E$ in the input script can be decreased to speed up the simulation.
No gravity is applied because the experiments were performed in a glycerol bath.

The external gradient field is defined by the magnetic scalar potential

$$\psi = -\frac{b}{\mu_0}\left[\frac{x^2 + y^2}{4} - \frac{z^2}{2}\right],$$

which gives $\mathbf{B} = -\mu_0\nabla\psi = b\langle x/2,\; y/2,\; -z\rangle$. The prefactor $b$ is related to the dimensionless parameters by

$$b = \sqrt{\frac{\mu_0}{\gamma E}}\,\frac{EI}{A}\left(\frac{\lambda_m}{L}\right)^3,$$

where $\lambda_m = L\,[MbA/(EI)]^{1/3}$ measures the external field strength relative to elastic stiffness.

The field strength is ramped up from $\lambda_m = 0$ to $3.5$ in increments of $0.1$, then ramped back down.

### Key LAMMPS settings

- **Units**: `units cgs`
- **Atom style**: `atom_style hybrid bpm/sphere dipole`
- **Contact**: `pair_style gran/hertz/history`; `special_bonds lj 0 1 1 coul 0 1 1` disables contact and dipole–dipole forces for bonded (1–2) neighbors
- **Dipole–dipole interactions**: `pair_style lj/cut/dipole/cut` with the LJ energy set to zero, combined with contact through `pair_style hybrid/overlay`
- **Magnetic loading**: `fix efield/lepton`
- **Time integration**: `fix nve/bpm/sphere update dipole`
- **Damping**: `fix langevin` for quasi-static equilibration

## Files

### Data file generator

The generator writes `helix-N103.lam`.
It stores zero dipoles in the data file; the LAMMPS input script sets the dipole magnitude from the selected $\gamma$.

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">create_helical_rod.py</span>
  <a class="download-btn" href="create_helical_rod.py" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```python
--8<-- "docs/examples/helical-rod/create_helical_rod.py"
```

</details>

### LAMMPS input script

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">deform-helix.in</span>
  <a class="download-btn" href="deform-helix.in" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```lammps
--8<-- "docs/examples/helical-rod/deform-helix.in"
```

</details>

## References

- Sano, T.G. (2022). Reduced theory for hard magnetic rods with dipole–dipole interactions. *J. Phys. A: Math. Theor.*, **55**(10), 104002. [doi:10.1088/1751-8121/ac4de2](https://doi.org/10.1088/1751-8121/ac4de2)
- Sano, T.G., Pezzulla, M., Reis, P.M. (2022). A Kirchhoff-like theory for hard magnetic rods under geometrically nonlinear deformation in three dimensions. *J. Mech. Phys. Solids*, **160**, 104739. [doi:10.1016/j.jmps.2021.104739](https://doi.org/10.1016/j.jmps.2021.104739)

---

[← Back to Examples](../../index.md#examples){.back-btn}
