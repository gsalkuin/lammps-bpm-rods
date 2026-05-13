import numpy as np

N = 8  # number of bonds

# Dimensionless units
L = 100.
dx = L / N

with open('beam-{}.lam'.format(int(N)), 'w') as file:
    file.write('LAMMPS DATA FILE \n\n')

    file.write(str(N+1) + ' atoms \n')
    file.write(str(N) + ' bonds \n\n')

    file.write('1 atom types \n')
    file.write('1 bond types \n\n')

    file.write('%f %f xlo xhi \n' % (-1, 101))
    file.write('%f %f ylo yhi \n' % (-101, 101))
    file.write('%f %f zlo zhi \n\n' % (-0.1, 0.1))

    # atom_style bpm/sphere
    file.write('Atoms \n\n')
    x = 0.
    for k in range(N+1):
        # atom-id atom-type molid diameter density x y z
        file.write('%i 1 1 1.0 1.0 %f 0 0 \n' % (k+1, x))
        x += dx
    file.write('\n')

    file.write('Bonds \n\n')
    for k in range(N):
        file.write('%i %i %i %i \n' % (k+1, 1, k+1, k+2))
    file.write('\n')
