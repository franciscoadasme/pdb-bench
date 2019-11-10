import sys
from shared import timeit

import mdtraj as md


def count(traj):
    return sum(1 for r in traj.topology.residues if r.name == "ALA")


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

traj = md.load(pdbfile)

print(timeit(count, traj, repeats=repeats))
