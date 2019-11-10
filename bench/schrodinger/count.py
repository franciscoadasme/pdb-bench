import sys
from shared import timeit

from schrodinger.structure import Structure


def count(st):
    return sum(1 for r in st.residue if r.pdbres == "ALA")


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

st = Structure.read(pdbfile)

print(timeit(count, st, repeats=repeats))
