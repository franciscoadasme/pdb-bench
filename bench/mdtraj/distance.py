import sys
from shared import timeit

import mdtraj as md


def distance(traj):
    da, _ = md.compute_contacts(traj, [(49, 59)], scheme="closest")
    return da.min()


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

traj = md.load(pdbfile)

print(timeit(distance, traj, repeats=repeats))
