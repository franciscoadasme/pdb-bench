import sys
from shared import timeit

import mdtraj as md


def align(traj):
    traj.superpose(traj[0])


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

traj = md.load(pdbfile)

print(timeit(align, traj, repeats=repeats))
