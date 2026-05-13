# Extreme Twisting & Plectoneme Formation

**Sec. 3.1** — Based on [Lazarus et al. (2013)](https://doi.org/10.1016/j.jmps.2013.04.002)

A heavy elastic rod is axially compressed until it buckles under gravity, then twisted up to 9 full rotations.
The simulation captures writhing instabilities and plectoneme formation with self-contact via Hertz–Mindlin interactions.
Two cases are considered: an initially straight rod and a naturally curved rod.

## Parameters

| Parameter | Value |
|-----------|-------|
| Rod length $L$ | 300 mm |
| Rod diameter $d$ | 3.1 mm |
| Density $\rho$ | 1.2 g/cm$^3$ |
| Young's modulus $E$ | 1.3 MPa |
| Shear modulus $G$ | 433 kPa |
| Number of bonds $N$ | 300 (bond length $\ell_0 \approx 1$ mm) |
| Curvature (curved rod) | $\kappa = 44.84$ m$^{-1}$ (2.14 coils) |

## Setup

Both cases share the same material properties and boundary conditions.
The left end is clamped (position and orientation fixed); the right end is displaced and rotated using `fix move`.

### Straight rod (`twist-straight.in`)

1. **Compress**: the right end is displaced by 80 mm at a constant velocity.
2. **Relax**: the rod equilibrates under Langevin dynamics with gravity.
3. **Twist**: the right end is rotated up to 9 full turns at a constant rate using `fix move rotate`.

### Curved rod (`twist-curved.in`)

The rod is initialized in its stress-free coiled configuration and bonds are created in this state.
An additional stage precedes the compression:

0. **Straighten**: a pure bending moment $M = \kappa EI$ is applied at the free end via `fix addtorque/atom` to uncoil the rod. Motion is constrained to 2D using `fix setforce NULL NULL 0`. Contact and gravity are disabled during this stage.
1. **Compress**: similar to the straight rod.
2. **Relax**: similar to the straight rod.
3. **Twist**: similar to the straight rod, with an additional frictionless wall at the rotating end (see below).

**Atom types and pair style.**
The left and right end atoms are assigned different atom types (`set group left type 2`, `set group right type 3`).
A `pair_style hybrid zero ... granular` is used so that during straightening (Stage 0), no granular contact is active — the coiled rod overlaps itself and contact would interfere with the uncoiling.
The `pair_coeff 2 3 granular ...` assignment exists only to satisfy LAMMPS's requirement that every sub-style in a hybrid pair style is referenced by at least one pair_coeff; types 2 and 3 are at opposite ends of the rod and never contact each other.
After straightening, type 1–1 granular contact is enabled for the compression and twist stages.

**Wall at the rotating end.**
During twisting, the plectoneme loop can migrate toward the rotating end and slip off, causing the rod to untwist.
A frictionless wall (`fix wall/gran ... xplane 0 220`) is placed near the rotating end to prevent this.
The straight rod does not need this wall because its plectoneme forms farther from the end.

### Key LAMMPS settings

- **Units**: `units lj` with *milli* units (see [Unit System](../../misc/unit-system.md))
- **Atom style**: `atom_style hybrid bpm/sphere dipole` — a "dipole" is assigned to each atom to visualize the twist. 
                   Note that the dipole vector does not capture the full 3D rotation.
- **Bond style**: `bond_style bpm/rotational` with `smooth no`, `break no`
- **Contact**: `pair_style granular` with Hertz–Mindlin; atom diameter set to $\approx 3 \ell_0$ for smooth contact surface; `special_bonds lj 0 0 0` disables pair forces for bonded neighbors up to 3 bonds away. This is necessary because the atom diameter causes bonded neighbors to overlap.
- **Time integration**: `fix nve/bpm/sphere update dipole`
- **Damping**: `fix langevin` for quasi-static equilibration

!!! note
    The atom diameter is independent of the rod diameter $d$. The latter is only used to compute the BPM parameters.

## Files

### Data file generator

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">create_straight_rod.py</span>
  <a class="download-btn" href="create_straight_rod.py" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>
```python
--8<-- "docs/examples/twisting-rod/create_straight_rod.py"
```
</details>

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">create_curved_rod.py</span>
  <a class="download-btn" href="create_curved_rod.py" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>
```python
--8<-- "docs/examples/twisting-rod/create_curved_rod.py"
```
</details>

### LAMMPS input scripts

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">twist-straight.in</span>
  <a class="download-btn" href="twist-straight.in" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>
```lammps
--8<-- "docs/examples/twisting-rod/twist-straight.in"
```
</details>

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">twist-curved.in</span>
  <a class="download-btn" href="twist-curved.in" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>
```lammps
--8<-- "docs/examples/twisting-rod/twist-curved.in"
```
</details>

## References

- Lazarus, A., Miller, J., Reis, P.M. (2013). Continuation of equilibria and stability of slender elastic rods using an asymptotic numerical method. *J. Mech. Phys. Solids*, **61**(8), 1712–1736. [doi:10.1016/j.jmps.2013.04.002](https://doi.org/10.1016/j.jmps.2013.04.002)
- Lazarus, A., Miller, J.T., Metlitz, M.M., Reis, P.M. (2013). Contorting a heavy and naturally curved elastic rod. *Soft Matter*, **9**(34), 8274–8281. [doi:10.1039/c3sm50873k](https://doi.org/10.1039/c3sm50873k)
- Grange, S., Bertrand, D. (2024). Co-rotational 3D beam element using quaternion algebra. *Int. J. Solids Struct.*, **293**, 112746. [doi:10.1016/j.ijsolstr.2024.112746](https://doi.org/10.1016/j.ijsolstr.2024.112746)

---

[← Back to Examples](../../index.md#examples){.back-btn}