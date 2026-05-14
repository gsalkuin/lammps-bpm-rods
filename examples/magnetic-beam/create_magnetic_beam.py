import numpy as np

N = 258  # number of bonds

# Physical properties (milli units: mm, ms, mg)
L = 25.8
W = 1.21
H = 0.49
density = 2.01  # mg/mm^3

# Discretization
d0 = L / N
diam = d0  # arbitrary
mass = density * d0 * W * H
dens = mass / (np.pi / 6 * diam**3)

# atom_style hybrid bpm/sphere dipole
with open('magnetic-beam-N{}.lam'.format(int(N)), 'w') as file:
    file.write('LAMMPS DATA FILE \n\n')

    file.write(str(N+1) + ' atoms \n')
    file.write(str(N) + ' bonds \n\n')

    file.write('1 atom types \n')
    file.write('1 bond types \n\n')

    file.write('%f %f xlo xhi \n' % (-1, 27))
    file.write('%f %f ylo yhi \n' % (-27, 27))
    file.write('%f %f zlo zhi \n\n' % (-0.1, 0.1))

    file.write('Atoms \n\n')
    x = 0.
    for k in range(N+1):
        # atom-id type x y z molid diameter density q mux muy muz
        file.write('%i 1 %f 0 0 1 %f %f 0 0 0 0 \n' % (k+1, x, diam, dens))
        x += d0
    file.write('\n')

    file.write('Bonds \n\n')
    for k in range(N):
        file.write('%i %i %i %i \n' % (k+1, 1, k+1, k+2))
    file.write('\n')