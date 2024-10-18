import sys
from shared import timeit

from chemfiles import Selection, Trajectory


def distance(f):
    phi_angles = []
    for i in range(2, len(f.topology.residues) + 1):
        idxs = Selection(
            f"(resid {i - 1} and name C) or (resid {i} and name N CA C)"
        ).evaluate(f)

        if len(idxs) % 4 != 0:
            continue
        for j in range(0, len(idxs), 4):
            phi_angles.append(f.dihedral(*idxs[j : j + 4]))

    psi_angles = []
    for i in range(1, len(f.topology.residues)):
        idxs = Selection(
            f"(resid {i} and name N CA C) or (resid {i + 1} and name N)"
        ).evaluate(f)

        if len(idxs) % 4 != 0:
            continue
        for j in range(0, len(idxs), 4):
            psi_angles.append(f.dihedral(*idxs[j : j + 4]))

    return phi_angles, psi_angles


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

f = Trajectory(pdbfile).read()
print()

print(timeit(distance, f, repeats=repeats))
