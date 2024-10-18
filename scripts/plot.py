# /usr/bin/env python3

import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def parse_data(path: str) -> Tuple[List[str], List[str], Dict[str, Dict[str, float]]]:
    data: Dict[str, Dict[str, float]] = {}
    programs: Set[str] = set()
    tests: Set[str] = set()
    with open(path) as file:
        for line in file:
            if not line.strip():
                continue
            line, time = line.split()
            program, test, protein = line.split("/")
            if test == "parse":
                test = f"parse/{protein}"
            if test not in data:
                data[test] = {}
            try:
                value = float(time)
            except ValueError:
                value = float("nan")
            data[test][program] = value
            programs.add(program)
            tests.add(test)
    programs = sorted(programs)
    tests = sorted(data.keys(), key=lambda x: (0 if "parse" in x else 1, x))
    return tests, programs, data


path = Path(sys.argv[1])
tests, programs, data = parse_data(path)

colormap = plt.get_cmap("tab10")
colors = {program: colormap(i) for i, program in enumerate(programs)}
width = 0.5

fig, ax = plt.subplots(figsize=(8, 4))
for i, test in enumerate(tests):
    for program in programs:
        time = np.log10(data[test][program]) if program in data[test] else np.nan
        ax.plot([i - width / 2, i + width / 2], [time, time], c=colors[program], lw=3)

for i, [left, right] in enumerate(zip(tests[:-1], tests[1:])):
    for program in programs:
        rhs = np.log10(data[left][program]) if program in data[left] else np.nan
        lhs = np.log10(data[right][program]) if program in data[right] else np.nan
        ax.plot(
            [i + width / 2, i + 1 - width / 2],
            [rhs, lhs],
            c=colors[program],
            ls=":",
            lw=1,
        )

ax.set_xlim(-2, len(data) + 1)
ax.set_xticks(range(len(tests)))
ax.set_xticklabels(tests, rotation=60, ha="right", family="monospace")
ax.set_ylim(-7, 2)
ax.yaxis.set_major_formatter("$10^{{{x:.0f}}}$")
ax.grid(True, axis="y", c="0.9", lw=1, ls="--")
ax.set_ylabel("Runtime / s", size="large")

ax.legend(
    handles=[
        ax.plot([], c=colors[program], marker="s", label=program, lw=0)[0]
        for program in programs
    ],
    loc="center left",
    ncol=1,
    frameon=False,
    handletextpad=0,
    bbox_to_anchor=(1, 0.5),
)

fig.tight_layout()
fig.savefig(path.with_suffix(".png"))
