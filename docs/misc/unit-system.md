The unit systems in LAMMPS can be found [here](https://docs.lammps.org/units.html).
In addition to SI units, we also use the two unit systems described below.

The `cgs` units assume ESU, but it can also be used for magnetostatics.
The CGS EMU units are provided below. Note that there are many different conventions for
Gaussian units. We use [this](https://www.nist.gov/system/files/documents/pml/electromagnetics/magnetics/magnetic_units.pdf).

### CGS EMU

| Quantity | Unit | SI |
|----------|------|----|
| Length | cm | $10^{-2}$ m |
| Time | s | 1 s |
| Mass | g | $10^{-3}$ kg |
| Force | dyne | $10^{-5}$ N |
| Energy | erg | $10^{-7}$ J |
| Density | g/cm$^3$ | $10^3$ kg/m$^3$ |
| Stress | Ba | $10^{-1}$ Pa |
| Magnetic dipole | emu = erg/G | $10^{-3}$ A m$^2$ |
| Magnetic field | G | $10^{-4}$ T |
| Kinematic viscosity | cm$^2$/s | $10^{-4}$ m$^2$/s |
| Temperature | K | 1 K |

A custom *milli* unit system (mm, mg, ms, $\mu_0/4\pi = 1$) is also used in some examples. 
All other units except temperature can be derived from the four fundamental units.

!!! Note
    The **temperature** in a DEM simulation is only used as a proxy for the kinetic energy 
    of the particles, relevant only when a Langevin thermostat is applied for quasi-static equilibration. 
    Its scale is set by $\{k_B\}$, the numerical value of the Boltzmann constant in the chosen LAMMPS unit system 
    ($\{k_B\} = 1$ in `units lj`; $\{k_B\} = 1.38 \times 10^{-16}$ in `units cgs`).

### milli

| Quantity | Unit | SI |
|----------|------|----|
| Length | mm | $10^{-3}$ m |
| Time | ms | $10^{-3}$ s |
| Mass | mg | $10^{-6}$ kg |
| Force | mN | $10^{-3}$ N |
| Energy | $\mu$J | $10^{-6}$ J |
| Density | mg/mm$^3$ | $10^3$ kg/m$^3$ |
| Stress | kPa | $10^3$ Pa |
| Magnetic dipole | -- | $10^{-4}$ A m$^2$ |
| Magnetic field | -- | $10^{-2}$ T |
| Kinematic viscosity | mm$^2$/ms | $10^{-3}$ m$^2$/s |
| Temperature | -- | $10^{-6}/\{k_B\}$ K |
