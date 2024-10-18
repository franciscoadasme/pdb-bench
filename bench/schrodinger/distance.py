import sys
from shared import timeit

from schrodinger.structure import Structure
from schrodinger.structutils.measure import measure_distance as measure_distance


def distance(struc):
    r1 = struc.findResidue("A:50")
    r2 = struc.findResidue("A:60")

    min_dist = float("inf")
    for a in r1.atom:
        for b in r2.atom:
            d = measure_distance(a, b)
            if d < min_dist:
                min_dist = d


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

struc = Structure.read(pdbfile)

print(timeit(distance, struc, repeats=repeats))
