import numpy as np

N = 103  # number of bonds

# Units: CGS EMU

# Physical parameters
L = 10.3
d = 0.2
rho = 1.23

# Helix geometry
R = 1.0              # coil radius
psi_tilde = 0.96     # normalized pitch angle
psi = psi_tilde * np.pi / 2
K = np.sin(psi) / R
nturns = K * L / (2 * np.pi)
cospsi = np.cos(psi)

# Discretization
ds = L / N
A = np.pi * d**2 / 4
mass = rho * ds * A

diam = d
m = mass * np.ones(N+1)
dens = m / (np.pi * diam**3 / 6)

# Helix positions: z goes from -H0 to 0, clamped atom at z = 0 (top)
s = np.linspace(0, L, N+1)
x = R * np.cos(-K * (s - L))
y = R * np.sin(-K * (s - L))
z = cospsi * (s - L)

with open('helix-N{}.lam'.format(int(N)), 'w') as file:
    file.write('LAMMPS DATA FILE \n\n')

    file.write(str(N+1) + ' atoms \n')
    file.write(str(N) + ' bonds \n\n')

    file.write('1 atom types \n')
    file.write('1 bond types \n\n')

    file.write('%f %f xlo xhi \n' % (-2*R, 2*R))
    file.write('%f %f ylo yhi \n' % (-2*R, 2*R))
    file.write('%f %f zlo zhi \n\n' % (z[0]-d, z[-1]+d))

    # atom_style hybrid bpm/sphere dipole
    file.write('Atoms \n\n')
    for k in range(N+1):
        # atom-id type x y z molid diameter density q mux muy muz
        file.write('%i 1 %f %f %f 1 %f %f 0 0 0 0\n' %
                   (k+1, x[k], y[k], z[k], diam, dens[k]))
    file.write('\n')

    file.write('Bonds \n\n')
    for k in range(N):
        file.write('%i %i %i %i \n' % (k+1, 1, k+1, k+2))
    file.write('\n')
