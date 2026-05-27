import numpy as np

lamda = 1   # metachronal wavelength (in units of Ncx * cilia_spacing)

# Physical parameters (milli units: mm, ms, mg, mu0 = 4pi)
L_plate = 36.           # plate side length
W_plate = 36.
t_plate = 0.8           # plate thickness

L = 4.                  # cilium length
D = 0.8                 # cilium diameter
rho_s = 2.39            # solid density

E = 185.                # Young's modulus (kPa)
G = 61.6                # shear modulus (kPa)
nu_s = E / (2 * G) - 1

# Magnetic properties
M_mag = 0.2             # magnetization (1e5 A/m)
B_ext = 8.0             # applied field strength

# Fluid (99% glycerol)
nuf = 0.913             # kinematic viscosity (mm^2/ms)
rho_f = 1.26            # fluid density (mg/mm^3)

# Actuation
omega_B = -30 * np.pi / 180 / 1000  # -30 deg/s in rad/ms
u_phys = L * np.abs(omega_B)        # tip speed estimate

# Domain
xbox = 90.
ybox = 36.
zbox = 45.

# Cilia array
Ncx = 8                 # cilia in x
Ncy = 8                 # cilia in y
Na = 11                 # atoms per cilium
cilia_spacing = 4.0     # mm between cilia centers

# Plate discretization
Nx_plate = 45
Ny_plate = 45

# LBM grid
dx_LB = 1.0

# BPM discretization
A_cs = np.pi * D**2 / 4
I_cs = np.pi * D**4 / 64
J_cs = np.pi * D**4 / 32
l0 = L / (Na - 1)

# =========================================================================
# Metachronal wave
# =========================================================================

mu_mag = M_mag * A_cs * l0

def theta_dipole(x_pos, lamda, theta0=0, sign=1):
    """Dipole phase angle as function of x position."""
    return sign * (2 * np.pi * x_pos / (lamda * Ncx * cilia_spacing) + theta0)

# =========================================================================
# Create geometry
# =========================================================================

Natoms = 0
V_xyz, V_type, V_molid, V_diam, V_rho, V_mu = [], [], [], [], [], []
bond_atoms = []

# --- Plate atoms (type 1, molecule 1) ---
dx_plate = L_plate / Nx_plate
dy_plate = W_plate / Ny_plate

plate_xyz = []
for i in range(Nx_plate):
    for j in range(Ny_plate):
        x = dx_plate * (i + 0.5)
        y = dy_plate * (j + 0.5)
        z = D / 2
        plate_xyz.append([x, y, z])

plate_xyz = np.array(plate_xyz)
n_plate = len(plate_xyz)

V_xyz.append(plate_xyz)
V_type.append(np.ones((n_plate, 1)))
V_molid.append(np.ones((n_plate, 1)))
V_diam.append(t_plate * np.ones((n_plate, 1)))
V_rho.append(rho_s * np.ones((n_plate, 1)))
V_mu.append(np.zeros((n_plate, 3)))
Natoms += n_plate

# --- Cilia (type 2, molecule 2) ---
def create_cilium(x0, y0, mol_id=2):
    global Natoms

    xyz = np.zeros((Na, 3))
    xyz[:, 0] = x0
    xyz[:, 1] = y0
    xyz[:, 2] = np.linspace(0, L, Na) + D / 2

    types = 2 * np.ones((Na, 1))
    molids = mol_id * np.ones((Na, 1))
    molids[0] = 1       # base atom belongs to plate molecule

    th = theta_dipole(x0, lamda=lamda, theta0=np.pi/2, sign=1)   # antiplectic
    mu = np.zeros((Na, 3))
    mu[:, 0] = mu_mag * np.cos(th)
    mu[:, 2] = mu_mag * np.sin(th)

    diams = D * np.ones((Na, 1))
    rhos = rho_s * np.ones((Na, 1))

    bonds = np.zeros((Na - 1, 2))
    for k in range(Na - 1):
        bonds[k, 0] = Natoms + k + 1
        bonds[k, 1] = Natoms + k + 2

    V_xyz.append(xyz)
    V_type.append(types)
    V_molid.append(molids)
    V_diam.append(diams)
    V_rho.append(rhos)
    V_mu.append(mu)
    bond_atoms.append(bonds)
    Natoms += Na

for i in range(1, Ncx + 1):
    for j in range(1, Ncy + 1):
        create_cilium(cilia_spacing * i, cilia_spacing * j)

# Concatenate
V_xyz = np.concatenate(V_xyz, axis=0)
V_type = np.concatenate(V_type, axis=0)
V_molid = np.concatenate(V_molid, axis=0)
V_diam = np.concatenate(V_diam, axis=0)
V_rho = np.concatenate(V_rho, axis=0)
V_mu = np.concatenate(V_mu, axis=0)
bond_atoms = np.concatenate(bond_atoms, axis=0)
Nbonds = bond_atoms.shape[0]

# Center structure in box
x_center = xbox / 2
y_center = ybox / 2
x_shift = x_center - V_xyz[:n_plate, 0].mean()
y_shift = y_center - V_xyz[:n_plate, 1].mean()
V_xyz[:, 0] += x_shift
V_xyz[:, 1] += y_shift

# =========================================================================
# Write data file
# =========================================================================

filename = f'magnetic-cilia-array-8x8-antiplectic-lambda-{lamda}.lam'

# atom_style hybrid molecular bpm/sphere dipole
# Columns: id type x y z molid diameter density q mux muy muz
with open(filename, 'w') as f:
    f.write('LAMMPS DATA FILE\n\n')

    f.write(f'{Natoms} atoms\n')
    f.write(f'{Nbonds} bonds\n\n')

    f.write('3 atom types\n')
    f.write('1 bond types\n\n')

    f.write(f'{0:.6f} {xbox:.6f} xlo xhi\n')
    f.write(f'{0:.6f} {ybox:.6f} ylo yhi\n')
    f.write(f'{0:.6f} {zbox:.6f} zlo zhi\n\n')

    f.write('Atoms\n\n')
    for i in range(Natoms):
        f.write('%i %i %f %f %f %i %f %f 0 %f %f %f\n' %
                (i+1, int(V_type[i, 0]),
                 V_xyz[i, 0], V_xyz[i, 1], V_xyz[i, 2],
                 int(V_molid[i, 0]), V_diam[i, 0], V_rho[i, 0],
                 V_mu[i, 0], V_mu[i, 1], V_mu[i, 2]))
    f.write('\n')

    f.write('Bonds\n\n')
    for i in range(Nbonds):
        f.write(f'{i+1} 1 {int(bond_atoms[i, 0])} {int(bond_atoms[i, 1])}\n')

print(f"Written: {filename}")
print(f"  {Natoms} atoms ({n_plate} plate + {Ncx*Ncy*Na} cilia)")
print(f"  {Nbonds} bonds")


def lbm_analysis(dt):
    """Check LBM and BPM stability for the given timestep (no subcycling)."""

    cs2 = 1.0 / 3.0
    nx_lb = int(xbox / dx_LB)
    ny_lb = int(ybox / dx_LB)
    nz_lb = int(zbox / dx_LB)

    U_star = u_phys * dt / dx_LB
    nu_LB = nuf * dt / dx_LB**2
    tau = 0.5 + nu_LB / cs2
    Ma = U_star / (1.0 / np.sqrt(3))
    Re = u_phys * L / nuf

    # BPM stability
    Kr = E * A_cs / l0
    Ks = 12 * E * I_cs / l0**3
    Kt = G * J_cs / l0
    Kb = E * I_cs / l0
    m_seg = rho_s * A_cs * l0
    I_seg = m_seg * l0**2
    omega_max = max(np.sqrt(Kr/m_seg), np.sqrt(Ks/m_seg),
                    np.sqrt(Kt/I_seg), np.sqrt(Kb/I_seg))
    dt_crit = np.pi / omega_max

    print(f"\nRe ≈ {Re:.4f}")
    print(f"LBM grid: {nx_lb} x {ny_lb} x {nz_lb}  (dx = {dx_LB})")
    print(f"BPM grid: Na = {Na}, l0 = {l0:.4f}")

    print(f"\n  dt      = {dt}")
    print(f"  U*      = {U_star:.6e}")
    print(f"  nu*     = {nu_LB:.6f}")
    print(f"  tau*    = {tau:.6f}", end="")
    if tau <= 0.5:      print("  *** UNSTABLE ***")
    elif tau < 0.55:    print("  *** WARNING: marginal (close to 0.5) ***")
    elif tau > 2.0:     print("  *** WARNING: over-damped ***")
    else:               print("  OK")

    print(f"  Ma      = {Ma:.6e}", end="")
    if Ma > 0.3:        print("  *** WARNING: high ***")
    else:               print("  OK")

    print(f"  dt_crit = {dt_crit:.6f}")
    print(f"  dt/dt_crit = {dt/dt_crit:.4f}", end="")
    if dt > dt_crit:    print("  *** UNSTABLE ***")
    else:               print("  OK")


if __name__ == '__main__':
    lbm_analysis(dt=0.008)
