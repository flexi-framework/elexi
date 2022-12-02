# FLEXI-PARTICLE

[![license](https://img.shields.io/github/license/flexi-framework/flexi.svg?maxAge=2592000)]()

FLEXI is a high-order numerical framework for solving PDEs, with a focus on Computational Fluid Dynamics. FLEXI is based on the Discontinuous Galerkin Spectral Element Method (DGSEM), which allows for high-order of accuracy and fully unstructured hexahedral meshes.
The solver is parallelized very efficiently and scales up to hundreds of thousand cores.

For the main development branch, see [flexi-framework/flexi](https://github.com/flexi-framework/flexi). This repository contains an extension of FLEXI towards particle-laden flows. The particle tracking approach operates in physical space and is based on methods from ray-tracing to handle intersections with possibly curved boundaries. Particles are implemented with a Lagrangian point-particle ansatz and can either be one- or two-way coupled to fluid phase. Particle-wall interactions are handled via a hard-sphere approach, where multiple models are available. 

For details on the available models and implementation, see our manuscript on [arxiv](https://arxiv.org/abs/2211.05458).

FLEXI has been developed by the [Numerics Research Group (NRG)][nrg] led by Prof. Claus-Dieter Munz at the Institute of Aerodynamics and Gasdynamics at the University of Stuttgart, Germany.

This is a scientific project. If you use FLEXI for publications or presentations in science, please support the project by citing our publications given in [references](REFERENCE.md).

## Installation / Documentation

For installation instructions see [install](INSTALL.md).

See the full documentation of the main development branch including usage instructions and tutorial for FLEXI [here][flexi]. Examples of particle-laden flow setups are given in the [regressioncheck folder](regressioncheck/checks/particle). A full documentation of all available parameters is given by calling `flexi --help`.
 
In case you have questions regarding FLEXI, want to report bugs or contribute to the project, feel free to open an issue or pull request.

## License
FLEXI is Copyright (C) 2016, Prof. Claus-Dieter Munz and is released under the terms of the GNU General Public License v3.0. For the full license terms see the included license file [license](LICENSE.md).

## List of Contributors
Numerous people have worked on and with FLEXI over the last years. We would like to thank all these [contributors](CONTRIBUTORS.md) for their efforts they spent on building FLEXI.

## Used libraries

FLEXI uses several external libraries as well as auxiliary functions from open source projects, including:
* [HDF5](https://www.hdfgroup.org/)
* [MPI](http://www.mcs.anl.gov/research/projects/mpi/)
* [LAPACK](http://www.netlib.org/lapack/)
* [PAPI](http://icl.cs.utk.edu/papi/)
* [OpenMP](http://www.openmp.org/)
* [FFTW](http://www.fftw.org/)

[nrg]:  https://www.iag.uni-stuttgart.de/arbeitsgruppen/numerische-methoden/
[flexi]: https://www.flexi-project.org/

## Regression check

For information about the regression checks, see [reggie](REGGIE.md).

