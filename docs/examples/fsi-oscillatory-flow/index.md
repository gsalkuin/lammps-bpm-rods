# Oscillatory Channel Flow

**Sec. 4.2** — Based on [Pinelli et al. (2017)](https://doi.org/10.1007/s11012-016-0513-0)

Five elastic filaments are immersed in an oscillatory channel flow between no-slip walls.

## Lattice Boltzmann and Immersed Boundary Method

The fluid–structure interaction (FSI) is handled by the [CoupLB](https://github.com/tengzhang48/CoupLB) package, a custom, lightweight LBM implementation in LAMMPS.
Momentum is exchanged between the Lagrangian BPM particles and the Eulerian fluid lattice through an immersed boundary method (IBM).
A penalty force drives each particle toward the local fluid velocity, and the equal and opposite reaction is spread back onto the fluid grid.

Because the stiff BPM springs require a smaller timestep than the LBM, multiple LAMMPS steps are performed per LB step.
This subcycling is controlled by the `md_per_lb` keyword in `fix couplb`.
The fluid velocity is interpolated once at the start of each LBM cycle, and the penalty force is recomputed at every substep using each particle's current velocity.

## Parameters

| Parameter | Value |
|-----------|-------|
| Filament length $L$ | 10 mm |
| Filament diameter $d$ | 1 mm |
| Young's modulus $E$ | 1.23 MPa |
| Poisson's ratio $\nu_s$ | 0.3 |
| Filament density $\rho_s$ | 0.97 g/cm$^3$ |
| Fluid kinematic viscosity $\nu_f$ | 0.1 mm$^2$/ms |
| Fluid density $\rho_f$ | 1.26 g/cm$^3$ |
| Peak velocity $u_\text{max}$ | 0.6 m/s |
| Oscillation frequency $f$ | 1 Hz |
| Reynolds number $\mathrm{Re}$ | 60 |
| Domain ($x \times y \times z$) | $10L \times 6L \times 5L$ |

## Setup

Five straight filaments are clamped to the lower wall ($y = 0$) of the channel and arranged in a row along $x$ at intervals of $L/2$, with the middle filament at the center.
The fluid is driven by a sinusoidal body-force acceleration

$$g_x(t) = g_0 \sin(2\pi f\, t),$$

where the amplitude is $g_0 = u_\text{max} \cdot 2\pi f$, determined from the analytical solution for oscillatory channel flow in the limit where the Stokes layer thickness is much smaller than the half-channel height.
No-slip walls bound the channel in $y$; the remaining directions are periodic.

### Discretization

The LBM grid has spacing $\Delta x_\text{LB} = 0.625$ mm ($160 \times 96 \times 80$ nodes).
The BPM bond length matches the LBM spacing ($\mathrm{Gr} = \Delta x_\text{BP}/\Delta x_\text{LB} = 1$), giving 16 bonds per filament.
The stiff BPM springs require subcycling at 21 LAMMPS steps per LB step (`md_per_lb 21`).

The data file generator includes a callable function `lbm_analysis()` that documents how the LBM parameters were selected.

### Key LAMMPS settings

- **Units**: `units lj` using our [milli unit system](../../misc/unit-system.md)
- **Atom style**: `atom_style hybrid molecular bpm/sphere`
- **Boundary**: `p f p` — periodic in $x$ and $z$, fixed in $y$ (walls)
- **Pair style**: `pair_style zero` — no contact
- **FSI coupling**: `fix couplb` with `wall_y 1 1` for no-slip walls and `gravity v_gx_osc 0 0` for the oscillating body force
- **Time integration**: `fix nve/bpm/sphere`
- **Damping**: no explicit damping — all dissipation comes from the fluid coupling

!!! note "Atom types"
    Three atom types are used: type 1 (clamped base), type 2 (interior), and type 3 (tip).
    The base atoms are excluded from the dynamics.

## Files

### Data file generator

The generator writes `filaments.lam` with 5 filaments ($\mathrm{Gr} = 1$ by default).
Run `lbm_analysis()` to print the full LBM parameter analysis.

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">create_filaments.py</span>
  <a class="download-btn" href="create_filaments.py" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```python
--8<-- "docs/examples/fsi-oscillatory-flow/create_filaments.py"
```

</details>

### LAMMPS input script

The input script is configured for 45 MPI ranks (`processors 5 3 3`).

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">flow.in</span>
  <a class="download-btn" href="flow.in" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```lammps
--8<-- "docs/examples/fsi-oscillatory-flow/flow.in"
```

</details>

## References

- Pinelli, A., Omidyeganeh, M., Brücker, C., Revell, A., Sarkar, A., Alinovi, E. (2017). The PELskin project: part IV — control of bluff body wakes using hairy filaments. *Meccanica*, **52**, 1503–1514. [doi:10.1007/s11012-016-0513-0](https://doi.org/10.1007/s11012-016-0513-0)
- Agrawal, V., Kulachenko, A., Scapin, N., Tammisola, O., Brandt, L. (2024). An efficient isogeometric/finite-difference immersed boundary method for the fluid–structure interactions of slender flexible structures. *Comput. Methods Appl. Mech. Engrg.*, **418**, 116495. [doi:10.1016/j.cma.2023.116495](https://doi.org/10.1016/j.cma.2023.116495)

---

[← Back to Examples](../../index.md#examples){.back-btn}
