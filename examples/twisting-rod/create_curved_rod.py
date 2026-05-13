import numpy as np

N = 300  # number of bonds

# Units: (mm, mg, ms)

# Physical parameters
L = 300.
d = 3.1
rho = 1.2

# Curvature
R = 1/44.84 * 1000
theta = L/R

# Discretization
t = np.linspace(0, theta, N+1)
x = R * np.cos(t-np.pi/2)
y = R * (1 + np.sin(t-np.pi/2))

dx = L/N
A = np.pi * d**2 / 4
mass = rho * dx * A

diam = 3 * dx # atom diameter for contact != rod diameter
m = mass * np.ones(N+1)
m[0] = m[-1] = 0.5 * mass
dens = m / (np.pi * diam**3 / 6)

with open('curved-rod-N{}.lam'.format(int(N)), 'w') as file:
    file.write('LAMMPS DATA FILE \n\n')

    file.write(str(N+1) + ' atoms \n')
    file.write(str(N) + ' bonds \n\n')

    file.write('3 atom types \n')
    file.write('1 bond types \n\n')

    # Create box
    file.write('%f %f xlo xhi \n'%(-1.1*L, 1.1*L))
    file.write('%f %f ylo yhi \n'%(-1.1*L, 1.1*L))
    file.write('%f %f zlo zhi \n\n'%(-0.6*L, 0.6*L))

    # atom_style hybrid bpm/sphere dipole
    file.write('Atoms \n\n')
    for k in range(N+1):
        # atom-id atom-type x y z molid diameter density q mux muy muz
        file.write('%i 1 %f %f 0 1 %f %f 0 0 1 0\n'%(k+1, x[k], y[k], diam, dens[k]))
    file.write('\n')

    file.write('Bonds \n\n')
    for k in range(N):
        file.write('%i %i %i %i \n'%(k+1, 1, k+1, k+2))
    file.write('\n')
