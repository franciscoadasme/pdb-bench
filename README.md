# PDB benchmarks

A series of tests to measure performance in handling PDB files. This was
originally developed to compare
[chem.cr](https://github.com/franciscoadasme/chem.cr) capabilities with
widely-used libraries and programs.

The benchmark is designed as follows:

* The tests are implemented using the functionality documented by each library
  in tutorials, examples, etc. Optimized versions may be faster but require
  advanced (possibly undocumented) usage.
* Tests are run ten times (except for 1HTQ, 3 times) and the elapsed time for
  each run is averaged.
* Parsing PDB files
  * [1CRN](http://www.rcsb.org/pdb/explore/explore.do?structureId=1crn) -
    hydrophobic protein (327 atoms).
  * [3JYV](http://www.rcsb.org/pdb/explore/explore.do?structureId=3jyv) - 80S
    rRNA (57,327 atoms).
  * [1HTQ](http://www.rcsb.org/pdb/explore/explore.do?structureId=1htq) -
    multicopy glutamine synthetase (10 models of 97,872 atoms).
* Counting the number of alanine residues in adenylate kinase (1AKE, 3816
  atoms).
* Calculating the distance between residues 50 and 60 of chain A in adenylate
  kinase (1AKE, 3816 atoms).
* Calculating the Ramachandran phi/psi angles in adenylate kinase (1AKE, 3816
  atoms).

Heavily inspired by [PDB benchmarks](https://github.com/jgreener64/pdb-benchmarks).

## Installation

Either clone or download this repository. Libraries and programs listed below
must be installed and available in the executable path. Please refer to the
official webpages for installation instructions.

## Usage

First, please do the following:

* Add `pdb-bench` to the `PYTHONPATH` environment variable
* Set `SCHRODINGER` environment variable to the Schrodinger installation folder
* Create the symbolic link `pdb-bench/lib` pointing to the chem.cr folder
* Ensure that the `vmd` executable is accessible

Then, simply run the `run.py` script to start the benchmark.

```shell
$ cd /path/to/pdb-bench
$ python run.py
biopython/count/1ake                    0.000191
chem/count/1ake                         0.000009
chemfiles/count/1ake                    0.000300
mdanalysis/count/1ake                   0.000036
mdtraj/count/1ake                       0.000075
schrodinger/count/1ake                  0.026060
vmd/count/1ake                          0.000173

biopython/distance/1ake                 0.000162
chem/distance/1ake                      0.000000
chemfiles/distance/1ake                 0.000998
mdanalysis/distance/1ake                0.000378
mdtraj/distance/1ake                    0.000936
schrodinger/distance/1ake               0.041747
vmd/distance/1ake                       0.000381

...
```

You can select which test to run by using the `-t/--tests` option:

```shell
$ python run.py -t distance,ramachandran
biopython/distance/1ake                 0.000162
chem/distance/1ake                      0.000000
chemfiles/distance/1ake                 0.000998
mdanalysis/distance/1ake                0.000378
mdtraj/distance/1ake                    0.000936
schrodinger/distance/1ake               0.041747
vmd/distance/1ake                       0.000381

...
```

Similarly, you can select which libraries to test by using the `-l/--libraries`
option:

```shell
$ python run.py -l chem
chem/count/1ake                         0.000008
chem/distance/1ake                      0.000000
chem/parse/1htq                         1.679115
chem/parse/1ake                         0.008246
chem/parse/3jyv                         0.094727
chem/parse/1crn                         0.001152
chem/ramachandran/1ake                  0.000596
```


## Comparison

**IMPORTANT**: direct comparison of parsing times should be taken with a grain
of salt because each library does something slightly different, e.g., error
checking. Some of this functionality is listed below. Nonetheless, these results
gives an overall picture in terms of the expected performance.

|                      | Biopython | chem.cr | Chemfiles | MDAnalysis | MDTraj | schrodinger |   VMD |
| -------------------- | --------: | ------: | --------: | ---------: | -----: | ----------: | ----: |
| Parse 1CRN [ms]      |     6.521 |   1.028 |     1.668 |      5.059 | 11.923 |      45.497 | 2.285 |
| Parse 3JYV [s]       |     0.837 |   0.086 |     0.199 |      0.404 |  1.490 |       0.766 | 0.162 |
| Parse 1HTQ [s]       |    16.146 |   1.673 |     2.540 |      1.387 | 18.969 |      11.997 | 0.236 |
| Count [ms]           |     0.210 |   0.009 |     0.322 |      0.041 |  0.079 |      25.997 | 0.165 |
| Distance [ms]        |     0.172 |   0.000 |     1.016 |      0.382 |  0.990 |      43.101 | 0.379 |
| Ramachandran [ms]    |   110.450 |   0.607 |         - |    690.201 |  4.947 |      68.758 | 1.814 |
|                      |           |         |           |            |        |             |       |
| License              | Biopython |     MIT |       BSD |      GPLv2 |   LGPL | Proprietary |   VMD |
| Parse Header         |       yes |     yes |       yes |         no |     no |          no |    no |
| Parse CONECT         |        no |     yes |       yes |         no |    yes |         yes |   yes |
| Guess bonds          |        no |      no |       yes |         no |    yes |         yes |   yes |
| Hybrid36             |        no |     yes |        no |        yes |     no |          no |    no |
| Hierarchical parsing |       yes |     yes |        no |         no |     no |          no |    no |
| Supports disorder    |       yes |     yes |        no |         no |    yes |         yes |    no |

Latest update: 2019-11-10

## Environment

#### PC specs

* CPU: Intel(R) Core(TM) i7-5930K CPU @ 3.50GHz
* RAM: 16 GB 2400 MHz DDR4
* SSD: 480GB Western Digital Green SATA3
* OS: Ubuntu 18.04.2

#### Software

| Library     | Version | Language   | Language version |
| ----------- | ------: | ---------- | ---------------: |
| Biopython   |    1.74 | Python     |            3.6.8 |
| chem.cr     |   0.1.0 | Crystal    |           0.31.0 |
| Chemfiles   |   0.9.1 | C++/Python |            3.6.8 |
| MDAnalysis  |  0.20.1 | Python     |            3.6.8 |
| MDTraj      |   1.9.3 | Python     |            3.6.8 |
| Schrodinger |  2018-4 | C++/Python |            3.6.1 |
| VMD         |   1.9.4 | C/Tcl      |                  |

**NOTE**: Python libraries here listed internally use C-powered code for some
performance-critical code such as distance calculation. In the other hand,
C++/Python libraries are implemented in pure C++ and offer access from Python
through bindings.