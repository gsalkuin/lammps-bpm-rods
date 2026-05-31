# Magnetic Beams

**Based on [Yan et al. (2022)](https://doi.org/10.1016/j.ijsolstr.2021.111319)**

These scripts simulate the three configurations in Fig. 6(b) of our paper (cf. Fig. 9 in Yan et al.).
For Fig. 6(a) (cf. Fig. 8(a)), add gravity along $+x$ to Configuration B1 using `fix gravity` (see the [twisting rod example](../twisting-rod/index.md)) and center the field at $(L/2, 0, 0)$ by replacing $x$ with $x-L/2$ in the potential.

## Parameters

| Parameter | Value |
|-----------|-------|
| Beam geometry ($L \times W \times H$) | 25.8 mm $\times$ 1.21 mm $\times$ 0.49 mm |
| Young's modulus $E$ | 1.16 MPa |
| Poisson's ratio $\nu$ | 0.49 |
| Density $\rho$ | 2.01 g/cm$^3$ |
| Magnetization magnitude $M$ | 94.1 kA/m |
| Number of bonds $N$ | 258 (bond length $\ell_0 = 0.1$ mm) |

## Setup

All three cases use the same straight rectangular beam.
The leftmost atom is clamped, motion is restricted to 2D, and the beam is equilibrated quasi-statically with Langevin dynamics after each magnetic loading increment.
The quantity of interest is the normalized tip deflection $\delta_y/L$, measured at the 
rightmost atom as in the [cantilever beam demo](../cantilever-beam/index.md).

The field strength per unit $\lambda_m^C$ is $B_1 = EI/(MAL^2)$.
For gradient-field configs (B1, B2), the gradient is scaled by $\lambda_m^\nabla \approx 0.36\,\lambda_m^C$ so that all three scripts report loading on the same $\lambda_m^C$ scale.

### Configuration A (`config-A.in`)

The beam is magnetized in the $+x$ direction and a uniform magnetic field is applied in the $+y$ direction:

$$\mathbf{B} = B \hat{\mathbf{y}}.$$

The simulation ramps the dimensionless parameter $\lambda_m^C = MBAL^2/(EI)$.
Here `fix efield` is used by analogy with magnetic dipoles; 
it applies the dipole torque $\boldsymbol{\tau} = \mathbf{m} \times \mathbf{B}$.

### Configuration B1 (`config-B1.in`)

The beam is magnetized in the $+y$ direction and a constant-gradient magnetic field is applied.
The magnetic field is defined using the scalar potential

$$\psi = -\frac{b}{\mu_0}\left[\frac{y^2}{2} - \frac{x^2 + z^2}{4}\right].$$

This is applied with `fix efield/lepton`.
Unlike `fix efield`, which only applies the dipole torque, `fix efield/lepton` also captures the force from the field gradient through the scalar potential.
This yields $\mathbf{B} = -\mu_0 \nabla \psi = b \langle -x/2, y, -z/2 \rangle$ and $\nabla \mathbf{B} = b\,\mathrm{diag}(-1/2, 1, -1/2)$, where $\lambda_m^\nabla \approx 0.36 \lambda_m^C$.

!!! note
    `fix efield/lepton` is used through the electric-field analogy and expects an electric potential expression.
    The equivalent potential passed to LAMMPS is therefore written to obtain the desired gradient prefactor `dB`; it is not the actual magnetic scalar potential, so no explicit $\mu_0$ appears in the script.
    Since these simulations are 2D with $z=0$, the input scripts only include the $x$ and $y$ terms.

### Configuration B2 (`config-B2.in`)

The beam is magnetized in the $+x$ direction with the same field as B1, with the addition of a small uniform $+y$ perturbation to kick the beam out of the initially straight equilibrium.

### Key LAMMPS settings

- **Units**: `units lj` using our [milli unit system](../../misc/unit-system.md)
- **Atom style**: `atom_style hybrid bpm/sphere dipole`
- **Pair style**: `pair_style zero`; no contact or dipole-dipole interactions
- **Magnetic loading**: `fix efield` for the uniform and `fix efield/lepton` for the non-uniform fields
- **Time integration**: `fix nve/bpm/sphere update dipole`
- **Damping**: `fix langevin` for quasi-static equilibration

## Files

### Data file generator

The generator writes `magnetic-beam-N258.lam`, which is read by all three input scripts.

!!! note
    The data file generator writes the dipole as zero.
    Each configuration script sets the dipole direction: $+x$ for A and B2, $+y$ for B1.

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">create_magnetic_beam.py</span>
  <a class="download-btn" href="create_magnetic_beam.py" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```python
--8<-- "docs/examples/magnetic-beam/create_magnetic_beam.py"
```

</details>

### LAMMPS input scripts

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">config-A.in</span>
  <a class="download-btn" href="config-A.in" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```lammps
--8<-- "docs/examples/magnetic-beam/config-A.in"
```

</details>

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">config-B1.in</span>
  <a class="download-btn" href="config-B1.in" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```lammps
--8<-- "docs/examples/magnetic-beam/config-B1.in"
```

</details>

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">config-B2.in</span>
  <a class="download-btn" href="config-B2.in" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```lammps
--8<-- "docs/examples/magnetic-beam/config-B2.in"
```

</details>

## References

- Yan, D., Abbasi, A., Reis, P.M. (2022). A comprehensive framework for hard-magnetic beams: reduced-order theory, 3D simulations, and experiments. *Int. J. Solids Struct.*, **257**, 111319. [doi:10.1016/j.ijsolstr.2021.111319](https://doi.org/10.1016/j.ijsolstr.2021.111319)

---

[← Back to Examples](../../index.md#examples){.back-btn}
