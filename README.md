[![logo](https://numericsresearchgroup.org/images/icons/elexi.svg "ƎLexi")][flexi]

[![license](https://img.shields.io/github/license/flexi-framework/flexi.svg?maxAge=2592000 "GPL-3.0 License")](LICENSE.md)
[![doi](https://img.shields.io/badge/DOI-10.1016/j.cpc.2023.108762-blue "DOI")](https://doi.org/10.1016/j.cpc.2023.108762)
[![youtube](https://img.shields.io/badge/YouTube-red?logo=youtube "YouTube")](https://www.youtube.com/@nrgiag8633)
[![userguide](https://img.shields.io/badge/Userguide-silver "Userguide")][userguide]
[![gallery](https://img.shields.io/badge/Gallery-teal "Gallery")][gallery]

# About

[ƎLexi][elexi] is a high-order numerical Eulerian-Lagrangian framework for solving PDEs, with a special focus on Computational Fluid Dynamics. It is an extension of [FLEXI][flexi] which is based on the Discontinuous Galerkin Spectral Element Method (DGSEM). DGSEM allows for high-order of accuracy and fully unstructured hexahedral meshes. The solver is parallelized very efficiently for large-scale applications and scales to 500,000+ cores. Moreover, [ƎLexi][elexi] comes with a capable pre- and post-processing suite that enables complex simulation setups up to the finished visualization.

For the main development branch, see [flexi-framework/flexi](https://github.com/flexi-framework/flexi). This repository contains an extension of FLEXI towards particle-laden flows. The particle tracking approach operates in physical space and is based on methods from ray-tracing to handle intersections with possibly curved boundaries. Particles are implemented with a Lagrangian point-particle ansatz and can either be one- or two-way coupled to the fluid phase. Particle-wall interactions are handled via a hard-sphere approach, where multiple models are available. 

For details on the available models and implementation, see our paper in [Computer Physics Communications](https://doi.org/10.1016/j.cpc.2023.108762) or the manuscript on [arxiv](https://arxiv.org/abs/2211.05458).

[ƎLexi][elexi]/[FLEXI][flexi] has been developed by the [Numerics Research Group (NRG)][nrg] founded by Prof. Claus-Dieter Munz and currently lead by Prof. Andrea Beck at the Institute of Aerodynamics and Gasdynamics at the University of Stuttgart, Germany.

You can find detailed installation instructions, the extensive documentation and several tutorial cases for FLEXI [here][flexi].

[ƎLexi][elexi]/[FLEXI][flexi] is Copyright (C) 2016, Prof. Claus-Dieter Munz and is released under the **GNU General Public License v3.0**. For the full license terms see the included [license file](LICENSE.md).

Numerous people have worked on and with [ƎLexi][elexi]/[FLEXI][flexi] over the last years. We would like to thank all these [contributors](CONTRIBUTORS.md) for their efforts they spent on building [ƎLexi][elexi]/[FLEXI][flexi].
 
In case you have questions regarding ƎLexi or want to contribute yourself by either reporting bugs, requesting features or adding somthing different to the project, feel free to open an issue or pull request.

# Cite
ƎLexi is a scientific project. If you use ƎLexi for publications or presentations in science, please support the project by citing it. As general reference, please cite
```
Kopper, P., Schwarz, A., Copplestone, S. M., Ortwein, P., Staudacher, S., and Beck, A.
A framework for high-fidelity particle tracking on massively parallel systems
Computer Physics Communications, 289, (2023), 108762
```
or use the following Bibtex entry

    @article{elexi,
     title     = {A framework for high-fidelity particle tracking on massively parallel systems},
     author    = {Kopper, Patrick and Schwarz, Anna and Copplestone, Stephen M. and Ortwein, Philip and Staudacher, Stephan and Beck, Andrea},
     journal   = {Computer Physics Communications},
     year      = {2023},
     month     = aug,
     pages     = {108762},
     doi       = {10.1016/j.cpc.2023.108762},
     publisher = {Elsevier},
    }

To refer to specific applications and features, you can also cite the appropriate paper from this [list][publications].

# Quick Start Guide
For a more detailed installation instructions, please see the documention [here][userguide].

[ƎLexi][elexi] is tested for various Linux distributions including Ubuntu, OpenSUSE, CentOS, or Arch. ƎLexi also runs on macOS. For the installation, you require the following dependencies:

| Package          | Required | Installed by ƎLexi |
|:-----------------|:--------:|:------------------:|
| Git              |      x   |                    |
| CMake            |      x   |                    |
| C/C++ Compiler   |      x   |                    |
| Fortran Compiler |      x   |                    |
| LAPACK/OpenBLAS  |      x   |      x             |
| HDF5             |      x   |      x             |
| MPI              |     (x)  |                    |

The MPI library is only required for running parallel simulations on multiple ranks. The HDF5 and LAPACK libraries can are optionally built and locally installed during the [FLEXI][flexi] build process. The names of the packages and the package manager might differ depending on the specific distribution used.

### Getting the code
Open a terminal, download [ƎLexi][elexi] via git

    git clone git@github.com:flexi-framework/elexi.git

### Compiling the code
Enter the [ƎLexi][elexi] directory, create a build directory and use CMake to configure and compile the code

    cd elexi
    cmake -B build
    cmake --build build

The executable `elexi` is now contained in the ƎLexi directory in `build/bin/`. Custom configurations of the compiler options, dependencies, and code features can be set using

    ccmake -B build

### Running the code
Navigate to the directory of the tutorial **cavity** and run [ƎLexi][elexi]

    cd tutorials/cavity
    flexi parameter_flexi.ini

# Used libraries
[ƎLexi][elexi] uses several external libraries as well as auxiliary functions from open source projects, including:
* [CMake](https://cmake.org)
* [FFTW](https://www.fftw.org)
* [HDF5](https://www.hdfgroup.org)
* [LAPACK](https://www.netlib.org/lapack)
* [MPI](https://www.mcs.anl.gov/research/projects/mpi)
* [OpenMP](https://www.openmp.org)
* [OpenBLAS](https://www.openblas.net)
* [PAPI](https://icl.cs.utk.edu/papi)
* [Reggie2.0](https://github.com/reggie-framework/reggie2.0)

[nrg]:          https://numericsresearchgroup.org/index.html
[flexi]:        https://numericsresearchgroup.org/flexi_index.html
[elexi]:        https://numericsresearchgroup.org/codes.html#codes_particle
[publications]: https://numericsresearchgroup.org/publications.html#services
[userguide]:    https://numericsresearchgroup.org/userguide/userguide.pdf
[gallery]:      https://numericsresearchgroup.org/gallery.html#portfolio
[youtube]:      https://www.youtube.com/@nrgiag8633 
