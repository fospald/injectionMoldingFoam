<p align="center">
  <a href="LICENSE" alt="GPLv3 license"><img src="https://img.shields.io/badge/license-GPLv3-brightgreen.svg" /></a>
  <a href="#" alt="no warranty"><img src="https://img.shields.io/badge/warranty-no-red.svg" /></a>
</p>

# injectionMoldingFoam

A solver for OpenFOAM for the simulation of fiber reinforced injection molding. It is a modification of the interFoam solver distributed with OpenFOAM. The solver is not validated and generally not suitable for industrial use.


## Requirements

The code was developed using [OpenFOAM](https://www.openfoam.com/)-2.3.x.
It might or might not work with newer versions of OpenFOAM (please test and let me know).


## Installation

1. download source
```
git clone https://github.com/fospald/injectionMoldingFoam.git
```
2. run Allwmake, which uses the wmake tool of OpenFOAM
```
sh Allwmake
```


## Run Demo

The demo consists of two joint boxes with two injection points.
The injection parameters may be specified 

```
cd demo
sh Allrun
```
This should run the solver, convert the results (20 snapshots) to VTK and open them in ParaView.
You should add a clip filter on the scalar `alpha.phase1` with value `0.5` in order to visualize the injected material.
Then use the time controls to step forward in time.


## Similar Projects

- [openInjMoldSim](https://github.com/krebeljk/openInjMoldSim)


## Acknowledgements

[Felix Ospald](https://www.tu-chemnitz.de/mathematik/part_dgl/people/ospald) gratefully acknowledges financial support by the [German Research Foundation](http://www.dfg.de/en/) (DFG), [Federal Cluster of Excellence EXC 1075](https://www.tu-chemnitz.de/MERGE/) "MERGE Technologies for Multifunctional Lightweight Structures". Many thanks to [Matti Schneider](https://www.itm.kit.edu/cm/287_3957.php) for his helpful introduction to FFT-based homogenization and ideas regarding the ACG distribution.

