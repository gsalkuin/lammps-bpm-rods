# Magnetically Actuated Cilia Array

**Based on [Gu et al. (2020)](https://doi.org/10.1038/s41467-020-16458-4)**

An $8 \times 8$ array of magnetically actuated cilia beats in an antiplectic metachronal wave pattern, pumping fluid along a channel.
The cilia magnetization varies along the array, with cilia at the same $x$ position sharing the same dipole orientation, and the array is driven by a rotating external magnetic field.
The phase offset between neighboring cilia produces a traveling wave whose wavelength $\lambda$ controls the pumping efficiency.

## Parameters

| Parameter | Value |
|-----------|-------|
| Cilium length $L$ | 4 mm |
| Cilium diameter $d$ | 0.8 mm |
| Young's modulus $E$ | 185 kPa |
| Shear modulus $G$ | 61.6 kPa |
| Solid density $\rho_s$ | 2.39 g/cm$^3$ |
| Number of bonds $N$ | 10 (bond length $\ell_0 = 0.4$ mm) |
| Magnetization $M$ | 20 kA/m |
| Applied field $B$ | 80 mT |
| Field angular velocity $\omega_B$ | −30 deg/s |
| Fluid kinematic viscosity $\nu_f$ | 0.913 mm$^2$/ms |
| Fluid density $\rho_f$ | 1.26 mg/mm$^3$ |
| Domain ($x \times y \times z$) | 90 mm $\times$ 36 mm $\times$ 45 mm |
| Metachronal wavelength $\lambda$ | varies |

## Setup

The cilia are immersed in a viscous fluid (99% glycerol) and driven by a magnetic field rotating in the $x$–$z$ plane:

$$\mathbf{B}(t) = B\bigl[\cos\phi(t)\,\hat{\mathbf{x}} + \sin\phi(t)\,\hat{\mathbf{z}}\bigr], \qquad \phi(t) = \frac{\pi}{2} + \omega t.$$

The magnetization pattern varies with position and is defined in the generated data file.
One thousand tracer particles are randomly distributed in the channel to visualize the flow and measure the net fluid transport per cycle.

### CoupLB

The fluid is periodic in $x$, bounded by no-slip walls in $y$, and no-slip at the bottom ($z = 0$) with free-slip at the top.
The LBM grid has spacing $\Delta x_\text{LB} = 1.0$ mm ($90 \times 36 \times 45$ nodes).
The BPM bond length is shorter than the grid spacing ($\mathrm{Gr} = \Delta x_\text{BP}/\Delta x_\text{LB} = 0.4$).
No subcycling is used; the LAMMPS and LBM timesteps are both 0.008 ms.
The data file generator includes a callable function `lbm_analysis(dt)` to check LBM and BPM stability.

### Key LAMMPS settings

- **Units**: `units lj` using our [milli unit system](../../misc/unit-system.md)
- **Atom style**: `atom_style hybrid molecular bpm/sphere dipole`
- **Boundary**: `p f f` — periodic in $x$, fixed walls in $y$ and $z$
- **Pair style**: `pair_style zero` — no contact
- **Magnetic actuation**: `fix efield v_bx 0 v_bz` applies the rotating field
- **FSI**: `fix couplb` with `wall_y 1 1` and `wall_z 1 3`
- **Time integration**: `fix nve/bpm/sphere update dipole` for cilia, `fix nve` for tracers
- **Confinement**: `fix wall/region` to prevent lost atoms
- **Damping**: no explicit damping — all dissipation comes from the fluid coupling

!!! note "Atom types and molecule IDs"
    The atom types are: 1 (plate), 2 (cilia), and 3 (tracers).
    The molecule IDs are: 1 (floor = plate + cilia base), 2 (free cilia atoms).
    Each cilium has $N+1$ atoms, where the base atom is excluded from the dynamics.
    The plate atoms are only for visualization.

## Files

### Data file generator

The generator writes a data file for the $8 \times 8$ cilia array using the normalized wavelength set by the `lamda` variable.
Call `lbm_analysis(dt)` to check LBM and BPM stability for a given timestep.

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">create_magnetic_cilia_array.py</span>
  <a class="download-btn" href="create_magnetic_cilia_array.py" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```python
--8<-- "docs/examples/fsi-magnetic-cilia-array/create_magnetic_cilia_array.py"
```

</details>

### LAMMPS input script

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">magnetic-cilia-array.in</span>
  <a class="download-btn" href="magnetic-cilia-array.in" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```lammps
--8<-- "docs/examples/fsi-magnetic-cilia-array/magnetic-cilia-array.in"
```

</details>

## References

- Gu, H., Boehler, Q., Cui, H., Secchi, E., Savorana, G., De Marco, C., Gervasoni, S., Peyron, Q., Huang, T.-Y., Pane, S., Hirt, A.M., Ahmed, D., Nelson, B.J. (2020). Magnetic cilia carpets with programmable metachronal waves. *Nat. Commun.*, **11**, 2637. [doi:10.1038/s41467-020-16458-4](https://doi.org/10.1038/s41467-020-16458-4)

---

[← Back to Examples](../../index.md#examples){.back-btn}
