# Demo: Cantilever Beams

!!! note
    This demo assumes familiarity with running LAMMPS simulations (input scripts, data files, fixes, and computes). 
    It is aimed at users experienced with molecular dynamics who want to learn the BPM workflow for simulating elastic rods.

This demo walks through the dynamic cantilever beam simulation based on the example from
[Le et al. 2011](https://doi.org/10.1007/s00466-011-0585-6).
By the end, you'll learn:

- How to generate a data file with Python
- How to define the BPM parameters
- How to set up a simple simulation

Additionally, the scripts for the static loading examples are provided at the end.

## Cantilever beam under dynamic loading

A rectangular beam with one end clamped has its free end subject to a sinusoidal transverse force
$P(t) = P_0 \sin(\omega t)$, where $P_0 = 10$ MN and $\omega = 50$ rad/s.
The following table lists the parameters.

| Parameter | Value |
|-----------|-------|
| Beam geometry ($L \times w \times h$) | $10$ m $\times$ $0.5$ m $\times$ $0.25$ m |
| Young's modulus $E$ | 210 GPa |
| Poisson's ratio $\nu$ | 0.3 |
| Density $\rho$ | 7850 kg/m$^3$ |
| Tip force amplitude $P_0$ | 10 MN |
| Angular frequency $\omega$ | 50 rad/s |

### Part 1: Data File Generation

Although it is possible to create the data file within LAMMPS, especially for simple 1D systems, 
we prefer to do this using a Python script for flexibility.

For a straight beam, atoms are placed at evenly spaced positions along the x-axis with sequential IDs. 
Bonds connect consecutive atom pairs (1–2, 2–3, ..., N–N+1), giving $N$ bonds and $N+1$ atoms.
A Python script is provided below.

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">create_beam.py</span>
  <a class="download-btn" href="create_beam.py" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```python
--8<-- "docs/examples/cantilever-beam/create_beam.py"
```

</details>

Important details:

- **`atom_style bpm/sphere`** uses per-atom mass set via density and not per-type mass.
  The effective density is back-calculated from the target mass and the chosen diameter.
- **Mass**: each bond with rest length $\ell_0 = L/N$ represents a rectangular segment of volume $\ell_0 \times w \times h$. 
            The masses are split between the two atoms forming the bond.  
- **Diameter**: the atom diameter is a geometric parameter for LAMMPS, not the physical rod diameter.
                Here we set it to $h = 0.25$ m.

### Part 2: BPM Parameters

!!! note
    Although the BPM assumes cylindrical bond segments, the four spring constants can be defined independently. For planar bending without twisting, $K_r$, $K_s$, and $K_b$ can be consistently defined for any cross-section. The twist stiffness $K_t$ does not affect the result in 2D.

For a rectangular beam with width $w$ and depth $h$:

- Area: $A = wh$
- Second moment of area: $I = wh^3/12$
- Polar second moment: $J = wh(w^2 + h^2)/12$

The four spring constants are as follows:

| Stiffness | Formula | Role |
|-----------|---------|------|
| $K_r$ | $EA / \ell_0$ | Axial stretching |
| $K_s$ | $12EI / \ell_0^3$ | Transverse shear |
| $K_t$ | $GJ / \ell_0$ | Twisting |
| $K_b$ | $EI / \ell_0$ | Bending |

The $K_s$ and $K_b$ formulas follow from Euler–Bernoulli beam theory [(Chen et al. 2022)](https://doi.org/10.1007/s10035-021-01187-2).
All damping coefficients are set to zero since this is a fully dynamic, undamped simulation to test the natural frequency response.

#### Atom diameter and moment of inertia

The atom diameter is set to $h$ in the data file. LAMMPS computes the moment of inertia of each `bpm/sphere` atom as $I_\text{sphere} = \tfrac{1}{10} m d^2$. With $d = h$, this is close enough to the Timoshenko rotary inertia $\tfrac{1}{12} m h^2$.

An alternative is to adjust the density and the diameter of the atom to produce the desired mass and moment of inertia.

!!! note
    The `fix nve/bpm/sphere` also has a `disc` keyword, which only changes the default moment of inertia formula. 
    Since the diameter and density in the data file assume spherical atoms, do *not* add the `disc` keyword.

### Part 3: Simulation Setup

The `left` and `right` groups are defined using regions that only capture the leftmost and rightmost atoms, respectively (assuming $N < 100$).
The clamped boundary condition on the left end is imposed by excluding `left` from time integration.
This means the leftmost atom's position and orientation are never updated.
The sinusoidal tip force is applied to the rightmost atom via `fix addforce`.
The total simulation time is $1.5$ s.

#### Timestep selection

The timestep $\Delta t$ for explicit time integration must be smaller than the critical timestep set by the highest-frequency mode.
For translational modes, $\omega = \sqrt{k/m}$; for rotational modes, $\omega = \sqrt{k/I_\text{sphere}}$.

Assuming $N = 100$ ($\ell_0 = 0.1$ m, $m \approx 98$ kg, $I \approx 0.61$ kg m$^2$):

| Mode | Stiffness | $\omega$ [rad/s] | $\Delta t_c = 2/\omega$ [s] |
|------|-----------|-------------------|-------------------|
| Radial | $K_r = EA/\ell_0 \approx 2.6 \times 10^{11}$ | $5.2 \times 10^{4}$ | $3.9 \times 10^{-5}$ |
| Shear | $K_s = 12EI/\ell_0^3 \approx 1.6 \times 10^{12}$ | $1.3 \times 10^{5}$ | $1.5 \times 10^{-5}$ |
| Twist | $K_t = GJ/\ell_0 \approx 2.6 \times 10^{9}$ | $6.5 \times 10^{4}$ | $3.1 \times 10^{-5}$ |
| Bend | $K_b = EI/\ell_0 \approx 1.4 \times 10^{9}$ | $4.7 \times 10^{4}$ | $4.2 \times 10^{-5}$ |

The shear spring is the stiffest, giving $\Delta t_c \approx 1.5 \times 10^{-5}$ s. 
In molecular dynamics, a good rule of thumb is to set $\Delta t < 0.1 \Delta t_c$.
Here, we use $\Delta t = 10^{-8}$ s, which is well below the threshold.

#### Notes:

- **`special_bonds`, `newton` and `comm_modify` settings** documented in [*How to BPM*](https://docs.lammps.org/Howto_bpm.html).
- **Pair style**: `pair_style zero` — no contact needed for this problem
- **`fix addforce`** to apply the loading
- **Time integration**: `fix nve/bpm/sphere` (no damping, no dipole update)
- **Output**: tip deflection $\delta_x / L$ and $\delta_y / L$ printed to file at regular intervals. The tip displacement is tracked using `compute reduce ave` on the `right` group.

The LAMMPS simulation script is provided below.
<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">bend-dynamic-point-load.in</span>
  <a class="download-btn" href="bend-dynamic-point-load.in" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```lammps
--8<-- "docs/examples/cantilever-beam/bend-dynamic-point-load.in"
```

</details>

## Cantilever beam under static loading

These examples use a different beam geometry: a cylindrical rod with $L = 100$, $d = 1$, $E = 1$, $\nu = 0.5$ in dimensionless units (`units lj`). The applied loads are chosen such that the maximum tip deflection is 10% of the beam length, and the simulated shapes are compared against the exact Euler–Bernoulli deflection curves (Hibbeler 2018).

#### Key differences from the dynamic example:

- **Damping**: `fix langevin` is used to reach quasi-static equilibrium. The load is first ramped from zero, then held constant while the beam relaxes.
- **Loading**: point load (`fix addforce` on `right`) or uniform load (`fix addforce` on all free atoms, with total load distributed evenly).
- **Energy verification**: the scripts compare the LAMMPS bond energy against the theoretical strain energy ($U = P^2 L^3 / (6EI)$ for point load; $U = w^2 L^5 / (40EI)$ for uniform load).

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">create_beam_for_static.py</span>
  <a class="download-btn" href="create_beam_for_static.py" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```python
--8<-- "docs/examples/cantilever-beam/create_beam_for_static.py"
```

</details>

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">bend-static-point-load.in</span>
  <a class="download-btn" href="bend-static-point-load.in" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```lammps
--8<-- "docs/examples/cantilever-beam/bend-static-point-load.in"
```

</details>

<details class="code-file" markdown>
<summary class="code-file-header">
  <span class="filename">bend-static-uniform-load.in</span>
  <a class="download-btn" href="bend-static-uniform-load.in" download onclick="event.stopPropagation()">⬇ Download</a>
</summary>

```lammps
--8<-- "docs/examples/cantilever-beam/bend-static-uniform-load.in"
```

</details>


## References

- Le, T.-N., Battini, J.-M., Hjiaj, M. (2011). Efficient formulation for dynamics of corotational 2D beams. *Computational Mechanics*, **48**(2), 153–161. [doi:10.1007/s00466-011-0585-6](https://doi.org/10.1007/s00466-011-0585-6)
- Chen, X., Peng, D., Morrissey, J.P., Ooi, J.Y. (2022). A comparative assessment and unification of bond models in DEM simulations. *Granular Matter*, **24**(1), 29. [doi:10.1007/s10035-021-01187-2](https://doi.org/10.1007/s10035-021-01187-2)
- Hibbeler, R.C. (2018). *Mechanics of Materials*. Pearson.

[← Back to Examples](../../index.md#examples){.back-btn}
