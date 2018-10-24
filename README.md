<p align="center">
  <a href="LICENSE" alt="GPLv3 license"><img src="https://img.shields.io/badge/license-GPLv3-brightgreen.svg" /></a>
  <a href="#" alt="no warranty"><img src="https://img.shields.io/badge/warranty-no-red.svg" /></a>
</p>

# injectionMoldingFoam

A solver for OpenFOAM for the simulation of injection molding. It is a modification of the interFoam solver distributed with OpenFOAM. The solver is not validated and generally not suitable for industrial use.


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

The demo consists of a simple 1x1x2 box. The points and distance field names are specified in `constant/fastMarchingDict`.
For this demo a distance filed `d` is generated, which gives the shortest distance to the point `(0.0 0.2 0.5)`. The `maxDist` parameter can be used perform the distance calculation only up to distance `maxDist` (in order to save computation time). Larger distances are then set to `maxDist`.
```
cd demo
sh Allrun
```
This should run the solver, convert the results to VTK and open them in ParaView.


## Similar Projects

- [openInjMoldSim](https://github.com/krebeljk/openInjMoldSim)


## Acknowledgements

[Felix Ospald](https://www.tu-chemnitz.de/mathematik/part_dgl/people/ospald) gratefully acknowledges financial support by the [German Research Foundation](http://www.dfg.de/en/) (DFG), [Federal Cluster of Excellence EXC 1075](https://www.tu-chemnitz.de/MERGE/) "MERGE Technologies for Multifunctional Lightweight Structures". Many thanks to [Matti Schneider](https://www.itm.kit.edu/cm/287_3957.php) for his helpful introduction to FFT-based homogenization and ideas regarding the ACG distribution.

