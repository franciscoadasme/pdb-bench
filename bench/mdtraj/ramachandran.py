import sys
from shared import timeit

import mdtraj as md


def ramachandran(traj):
    md.compute_phi(traj)
    md.compute_psi(traj)


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

traj = md.load(pdbfile)

print(timeit(ramachandran, traj, repeats=repeats))
