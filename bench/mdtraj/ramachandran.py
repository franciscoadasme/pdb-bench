import sys
from shared import timeit

import numpy as np
import mdtraj as md
from mdtraj.geometry.dihedral import indices_phi, indices_psi


def ramachandran(traj):
    md.compute_phi(traj)
    md.compute_psi(traj)


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

traj = md.load(pdbfile)

print(timeit(ramachandran, traj, repeats=repeats))
