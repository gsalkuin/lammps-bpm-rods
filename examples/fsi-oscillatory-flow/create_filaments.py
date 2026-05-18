import numpy as np

Gr = 1.0  # grid ratio dx_BP/dx_LB

# Physical parameters (milli units: mm, ms, mg)
L = 10.             # filament length
D = 1.              # filament diameter
rho = 0.97          # filament density

E = 1230.           # Young's modulus (kPa)
nu_s = 0.3          # Poisson's ratio
G = E / (2 * (1 + nu_s))

nuf = 0.1           # fluid kinematic viscosity (mm^2/ms)
rho_f = 1.26        # fluid density (mg/mm^3)
u_phys = 0.6        # characteristic velocity (mm/ms)

freq = 1e-3         # oscillation frequency (1/ms = 1 Hz)

# Domain
x_length = 10*L
y_length = 6*L
z_length = 5*L
H = y_length        # channel height

# Cilia layout
Ncilia = 5
sep = L / 2         # spacing between filaments

# LBM grid spacing
dx_LB = 0.625

# BPM discretization
dx_BP = dx_LB * Gr
Nb = int(np.round(L / dx_BP))      # bonds per filament
Na = Nb + 1                         # atoms per filament
l0 = L / Nb

# Cross-section
A = np.pi * D**2 / 4
I = np.pi * D**4 / 64
J = np.pi * D**4 / 32

# =========================================================================
# Write data file
# =========================================================================

x_center = x_length / 2
z_center = z_length / 2
x_positions = x_center + np.linspace(-(Ncilia-1)*sep/2, (Ncilia-1)*sep/2, Ncilia)

Natoms = Ncilia * Na
Nbonds = Ncilia * Nb

# atom_style hybrid molecular bpm/sphere
# Columns: id type x y z molid diameter density
with open('filaments.lam', 'w') as f:
    f.write('LAMMPS DATA FILE\n\n')

    f.write(f'{Natoms} atoms\n')
    f.write(f'{Nbonds} bonds\n\n')

    f.write('3 atom types\n')
    f.write('1 bond types\n\n')

    f.write(f'{0:.6f} {x_length:.6f} xlo xhi\n')
    f.write(f'{0:.6f} {y_length:.6f} ylo yhi\n')
    f.write(f'{0:.6f} {z_length:.6f} zlo zhi\n\n')

    f.write('Atoms\n\n')
    atom_id = 0
    for c in range(Ncilia):
        for k in range(Na):
            atom_id += 1
            x = x_positions[c]
            y = D / 2 + k * l0
            z = z_center

            if k == 0:
                atype = 1       # bottom (clamped)
            elif k == Na - 1:
                atype = 3       # tip
            else:
                atype = 2       # interior

            f.write('%i %i %f %f %f %i %f %f\n' %
                    (atom_id, atype, x, y, z, c+1, D, rho))
    f.write('\n')

    f.write('Bonds\n\n')
    bond_id = 0
    for c in range(Ncilia):
        for k in range(Nb):
            bond_id += 1
            a1 = c * Na + k + 1
            a2 = c * Na + k + 2
            f.write('%i 1 %i %i\n' % (bond_id, a1, a2))
    f.write('\n')


# =========================================================================
# LBM parameter selection and stability analysis
#
# Call lbm_analysis() to print the full report.
# This documents how the parameters in flow.in were chosen.
#
# Option A (Krueger textbook): fix dx_LB, choose U*
#   U*   = u_phys * dt / dx_LB       =>  dt = U* * dx_LB / u_phys
#   nu*  = nu * dt / dx_LB^2
#   tau* = 0.5 + nu* / cs2           (cs2 = 1/3)
#   Ma   = U* / cs                   (cs = 1/sqrt(3))
#
# Subcycling (CoupLB md_per_lb):
#   Multiple LAMMPS steps per LB step when dt_LB > dt_critical for BPM
# =========================================================================

def lbm_analysis(U_star=0.08):
    """Print LBM parameter selection, stability analysis, and LAMMPS values."""

    cs = 1.0 / np.sqrt(3)
    cs2 = 1.0 / 3.0
    nx = int(x_length / dx_LB)
    ny = int(y_length / dx_LB)
    nz = int(z_length / dx_LB)

    # Derived LBM parameters
    dt_LB = U_star * dx_LB / u_phys
    nu_LB = nuf * dt_LB / dx_LB**2
    tau = 0.5 + nu_LB / cs2
    Ma = U_star / cs

    l_star = L / dx_LB
    Re_phys = u_phys * L / nuf
    Re_LB = U_star * l_star / nu_LB

    vel_scale = dx_LB / dt_LB
    force_scale = rho_f * dx_LB**4 / dt_LB**2

    print(f"Re = {Re_phys:.2f}")
    print(f"\nLBM grid:  {nx} x {ny} x {nz}  (dx_LB = {dx_LB})")
    print(f"BPM grid:  Nb = {Nb}, Na = {Na}, l0 = {l0:.6f}, Gr = {l0/dx_LB:.4f}")

    print(f"\n{'=' * 70}")
    print("LATTICE BOLTZMANN PARAMETERS")
    print(f"{'=' * 70}")
    print(f"\n  Chosen:  U* = {U_star}")
    print(f"\n  Derived:")
    print(f"    dt_LB  = {dt_LB:.6f}")
    print(f"    nu*    = {nu_LB:.6f}")
    print(f"    tau*   = {tau:.6f}")
    print(f"    Ma     = {Ma:.6f}")
    print(f"    l*     = {l_star:.1f}")
    print(f"\n  Re check:  Re_phys = {Re_phys:.4f}   Re_LB = {Re_LB:.4f}   "
          f"{'MATCH' if abs(Re_phys - Re_LB)/max(Re_phys, 1e-15) < 1e-6 else 'MISMATCH'}")
    print(f"\n  Conversion scales:")
    print(f"    vel_scale   = {vel_scale:.6f}")
    print(f"    force_scale = {force_scale:.6e}")

    # Body force amplitude (oscillatory Stokes flow)
    #
    # For g(t) = g0 sin(wt) between walls separated by H, the centerline
    # velocity is u_max ~ g0/w in the thick-channel limit (H >> delta_s).
    # Stokes layer: delta_s = sqrt(2*nu/w).  Here delta_s ~ 5.6 mm and
    # H/2 = 30 mm, so the approximation is valid (error < 1%).
    omega = 2 * np.pi * freq
    delta_s = np.sqrt(2 * nuf / omega)
    G0 = u_phys * omega

    print(f"\n{'=' * 70}")
    print("BODY FORCE AMPLITUDE (oscillatory Stokes flow)")
    print(f"{'=' * 70}")
    print(f"\n  omega     = 2*pi*f = {omega:.6e}")
    print(f"  delta_s   = sqrt(2*nu/omega) = {delta_s:.4f} mm")
    print(f"  H/(2*delta_s) = {H/(2*delta_s):.2f}  (>> 1: thick-channel limit)")
    print(f"\n  G0 = u_max * omega = {G0:.6e}  (mm/ms^2)")

    # Stability analysis
    print(f"\n{'=' * 70}")
    print("STABILITY ANALYSIS")
    print(f"{'=' * 70}")

    all_ok = True

    print(f"\n  1. RELAXATION TIME:  tau* = {tau:.6f}", end="  ")
    if tau <= 0.5:
        print("UNSTABLE"); all_ok = False
    elif tau < 0.55:
        print("MARGINAL")
    elif tau > 2.0:
        print("HIGH (over-damped)")
    else:
        print("GOOD")

    print(f"  2. MACH NUMBER:     Ma   = {Ma:.6f}", end="  ")
    if Ma > 0.3:
        print("HIGH"); all_ok = False
    elif Ma > 0.1:
        print("ACCEPTABLE")
    else:
        print("EXCELLENT")

    print(f"  3. GRID RESOLUTION: l*   = {l_star:.1f}", end="    ")
    if l_star < 10:
        print("COARSE")
    elif l_star < 20:
        print("MODERATE")
    else:
        print("GOOD")

    # Sensitivity: U*
    print(f"\n  4. SENSITIVITY (varying U* at dx_LB = {dx_LB}):")
    print(f"     {'U*':<8} {'dt_LB':<12} {'nu*':<10} {'tau*':<10} {'Ma':<10} {'Status'}")
    print(f"     {'-'*66}")
    for U_t in [0.02, 0.05, 0.08, 0.10, 0.15, 0.20]:
        dt_t = U_t * dx_LB / u_phys
        nu_t = nuf * dt_t / dx_LB**2
        tau_t = 0.5 + 3 * nu_t
        Ma_t = U_t / cs
        st = ("UNSTABLE" if tau_t <= 0.5 else "marginal tau" if tau_t < 0.55
              else "high Ma" if Ma_t > 0.3 else "over-damped" if tau_t > 2.0 else "OK")
        mk = " <--" if abs(U_t - U_star) < 1e-10 else ""
        print(f"     {U_t:<8.3f} {dt_t:<12.6f} {nu_t:<10.6f} {tau_t:<10.6f} {Ma_t:<10.4f} {st}{mk}")

    # Sensitivity: dx_LB
    print(f"\n  5. SENSITIVITY (varying dx_LB at U* = {U_star}):")
    print(f"     {'dx_LB':<10} {'dt_LB':<12} {'nu*':<10} {'tau*':<10} {'l*':<8} {'Nx':<6} {'Status'}")
    print(f"     {'-'*72}")
    for dx_t in [0.3125, 0.625, 1.0, 1.25, 2.0]:
        dt_t = U_star * dx_t / u_phys
        nu_t = nuf * dt_t / dx_t**2
        tau_t = 0.5 + 3 * nu_t
        l_t = L / dx_t
        nx_t = int(x_length / dx_t)
        st = ("UNSTABLE" if tau_t <= 0.5 else "marginal tau" if tau_t < 0.55
              else "over-damped" if tau_t > 2.0 else "too coarse" if l_t < 10 else "OK")
        mk = " <--" if abs(dx_t - dx_LB) < 1e-10 else ""
        print(f"     {dx_t:<10.4f} {dt_t:<12.6f} {nu_t:<10.6f} {tau_t:<10.6f} {l_t:<8.1f} {nx_t:<6} {st}{mk}")

    # Timestep synchronization
    Kr = E * A / l0
    Ks = 12 * E * I / l0**3
    Kt = G * J / l0
    Kb = E * I / l0

    m_seg = rho * A * l0
    I_seg = m_seg * l0**2

    omega_r = np.sqrt(Kr / m_seg)
    omega_s = np.sqrt(Ks / m_seg)
    omega_t = np.sqrt(Kt / I_seg)
    omega_b = np.sqrt(Kb / I_seg)
    omega_max = max(omega_r, omega_s, omega_t, omega_b)
    f_max = omega_max / (2 * np.pi)

    dt_crit = np.pi / omega_max
    dt_safe = 0.1 / f_max

    dt_lammps_target = min(0.5 * dt_safe, 0.5 * dt_crit)
    lammps_per_lb = max(1, int(np.ceil(dt_LB / dt_lammps_target)))
    dt_lammps = dt_LB / lammps_per_lb

    print(f"\n{'=' * 70}")
    print("TIMESTEP SYNCHRONIZATION (CoupLB md_per_lb subcycling)")
    print(f"{'=' * 70}")
    print(f"\n  Spring frequencies:")
    print(f"    omega_r = {omega_r:.4f}  omega_s = {omega_s:.4f}  "
          f"omega_t = {omega_t:.4f}  omega_b = {omega_b:.4f}")
    print(f"    omega_max = {omega_max:.4f}  f_max = {f_max:.4f}")
    print(f"\n  Stability limits:")
    print(f"    dt_critical = {dt_crit:.6f}    dt_safe = {dt_safe:.6f}")
    print(f"\n  Timesteps:")
    print(f"    dt_LB     = {dt_LB:.6f}")
    print(f"    dt_LAMMPS = {dt_lammps:.6f}")
    print(f"    md_per_lb = {lammps_per_lb}")

    r_crit = dt_lammps / dt_crit
    r_safe = dt_lammps / dt_safe
    print(f"\n  Safety:  dt/dt_crit = {r_crit:.4f}   dt/dt_safe = {r_safe:.4f}")

    # Simulation time
    T_osc = 1.0 / freq                 # period = 1000 ms = 1 s
    n_cycles = 10
    sim_time = n_cycles * T_osc        # 10 s
    n_steps = int(np.ceil(sim_time / dt_lammps))
    n_lb = n_steps // lammps_per_lb

    print(f"\n{'=' * 70}")
    print("SIMULATION TIME")
    print(f"{'=' * 70}")
    print(f"  Target:     {sim_time:.0f} ms ({sim_time/1000:.0f} s, {n_cycles} cycles at {1/T_osc*1000:.0f} Hz)")
    print(f"  LAMMPS:     {n_steps:,} steps @ dt = {dt_lammps:.6f}")
    print(f"  LB:         {n_lb:,} steps @ dt = {dt_LB:.6f}")
    print(f"  md_per_lb = {lammps_per_lb}")

    # Summary
    print(f"\n{'=' * 70}")
    print("VALUES FOR flow.in")
    print(f"{'=' * 70}")
    print(f"  timestep    {dt_lammps:.6f}")
    print(f"  nuf         {nuf:.6f}")
    print(f"  rhof        {rho_f:.4f}")
    print(f"  G0          {G0:.6e}")
    print(f"  freq        {freq:.6e}")
    print(f"  grid        {nx} {ny} {nz}")
    print(f"  md_per_lb   {lammps_per_lb}")
    print(f"  run         {n_steps}")

    if all_ok:
        print(f"\nALL CHECKS PASSED")
    else:
        print(f"\nISSUES DETECTED -- REVIEW ABOVE")


if __name__ == '__main__':
    lbm_analysis()
