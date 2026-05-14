import numpy as np

N = 64  # number of bonds

# Units: SI

# Physical parameters
L = 10.
w = 0.5
h = 0.25
rho = 7850.

# Discretization
x = np.linspace(0, L, N+1)
dx = L/N

diam = h                    # set atom diameter to h; not required
mass = rho * dx * w * h     # mass of each bond...
m = mass * np.ones(N+1)    
m[0] = m[-1] = 0.5 * mass   # ...evenly split between the atoms
dens = m / (np.pi * diam**3 / 6)

with open('beam-{}.lam'.format(int(N)), 'w') as file:
    file.write('LAMMPS DATA FILE \n\n')

    file.write(str(N+1) + ' atoms \n')
    file.write(str(N) + ' bonds \n\n')

    file.write('1 atom types \n')
    file.write('1 bond types \n\n')

    # Create box
    file.write('%f %f xlo xhi \n'%(-1, 11))
    file.write('%f %f ylo yhi \n'%(-11, 11))
    file.write('%f %f zlo zhi \n\n'%(-0.1, 0.1))

    # atom_style bpm/sphere
    file.write('Atoms \n\n')
    for k in range(N+1):
        # atom-id atom-type molid diameter density x y z
        file.write('%i 1 1 %f %f %f 0 0 \n'%(k+1, diam, dens[k], x[k]))
    file.write('\n')

    file.write('Bonds \n\n')
    for k in range(N):
        file.write('%i %i %i %i \n'%(k+1, 1, k+1, k+2))
    file.write('\n')