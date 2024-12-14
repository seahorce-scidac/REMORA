---
title: 'REMORA: Regional Modeling of Oceans Refined Adaptively (built on AMReX)'

tags:
  - C++
  - ocean modeling
  - mesoscale

authors:
  - name: Hannah Klion
    orcid: 0000-0003-2095-4293
    corresponding: true
    affiliation: 1
  - name: Jean Sexton
    orcid: 0000-0003-2551-1678
    corresponding: false
    affiliation: 1
  - name: Ann Almgren
    orcid: 0000-0003-2103-312X
    corresponding: false
    affiliation: 1
  - name: Robert Hetland
    orcid: 0000-0001-9531-2119
    corresponding: false
    affiliation: 2
  - name: Kyle Hinson
    orcid: 0000-0002-2737-2379
    corresponding: false
    affiliation: 2
  - name: Iulian Grindeanu
    orcid: 0000-0002-0264-8236
    corresponding: false
    affiliation: 3
  - name: Vijay Mahadevan
    orcid: 0000-0002-3337-2607
    corresponding: false
    affiliation: 3

affiliations:
 - name: Lawrence Berkeley National Laboratory
   index: 1
 - name: Pacific Northwest Energy Laboratory
   index: 2
 - name: Argonne National Laboratory
   index: 3

date: January 2025

bibliography: paper.bib
---

# Summary

The REMORA code is a new model that simulates the mesoscale and microscale
dynamics of the ocean using the latest high-performance computing architectures.
It employs hierarchical parallelism using an MPI+X model, where X may be OpenMP on
multicore CPU-only systems, or CUDA, HIP, or SYCL on GPU-accelerated systems.
It is able to be built and run in both single and double precision.
REMORA is built on AMReX [@AMReX:JOSS; @AMReX:IJHPCA],
a block-structured adaptive mesh refinement (AMR) software framework that
provides the underlying performance-portable software infrastructure for block-structured mesh operations.
REMORA, like ROMS, is a regional model, meaning that it simulates the ocean dynamics on
a less-than-global scale, and as such requires boundary conditions derived from data
or from larger-scale models.
The REMORA development is funded by the US Department of Energy's Office of Science
through the Science Discovery through Advanced Computing (SciDAC) partnership program.

# REMORA Features

### Hydrodynamic Evolution

REMORA solves the incompressible time-dependent Navier-Stokes equation with the Boussinesq and hydrostatic approximations.
Temperature, salinity, and a passive scalar are also advected and diffused.
The density is calculated from a linear equation of state. The strength of vertical diffusion and viscosity is parametrized either by a spatially-varying analytical function or a Generic Length Scale (GLS) model [@umlauf:03].

### Time and Space Discretization and Terrain

Like ROMS, REMORA uses a split-explicit time-stepping scheme, where several fast barotropic (2D) steps take place within each baroclinic (3D) update.
In the barotropic steps, the code solves depth-averaged versions of the 3D evolution equations.
These vertically-averaged solutions are used to calculate the sea surface height and vertical-mean velocity.
Full 3D equations are then evolved for for velocities and scalars. 
Specifically, REMORA uses the same time integration as Rutgers ROMS.
That is, each barotropic step consists of a leapfrog predictor followed by a three-time Adams-Moulton corrector.
The 3D momenta are updated with a third-order Adams-Bashforth scheme, and scalars are advanced with a leapfrog step with a trapezoidal correction.

The spatial discretization in REMORA uses the classic Arakawa C-grid with
scalar quantities at cell centers and normal velocities at cell faces.
Bathymetry and sea-surface height are defined at the centers of the cells of the 2D grid.
Horizontally, the evolution equations are discretized over a boundary-following, orthogonal curvilinear grid, specified by metric terms.
This formulation allows for grids that, for example, conform to coastlines.
Land areas can be included in the domain and are represented by masks on cell centers and edges.
Fluxes, velocities, and tracer values are set to zero where the land mask is true.
The advection terms may be calculated using second- through fourth-order accurate
spatial discretizations, including both centered difference and upwind
schemes.

Vertically, the domain is discretized using a stretched, terrain-following vertical coordinate.
There are the same number of vertical levels everywhere; a spatially-varying water column depth (bathymetry and sea-surface height) is captured by cells of different thickness.
Cell thicknesses are determined by a non-linear transformation function that has parameters to control the distribution of levels.

### Physical Forcings and Boundary Conditions

Physical forcings include Coriolis, wind stress forcing, and bottom drag.
Lateral boundary conditions can be specified as periodic, inflow/outflow, radiation (following @orlanski:76), 
or time-varying values read in from external files in NetCDF format.
The solution at the boundary can either be clamped to the value specified from file, or deviations from the specified value can be radiated out.
For the barotropic variables in REMORA, this radiation occurs at the speed of external gravity waves, using the schemes of @flather:76 and @chapman:85 for momenta and sea-sea surface height, respectively.
The Orlanski radiation boundary condition radiates deviations in the 3D momenta and tracers at the local normal phase velocity.
REMORA uses the mixed radiation-nudging boundary condition of @marchesiello:01, where the Orlanski radiation condition is used for cells
where there is outflow, and nudging to a known exterior value is used where there is inflow.
The initial data can be specified by the user analytically or read from NetCDF files.

# Statement of need

Most widely used ocean modeling codes today do not have the
ability to use GPU acceleration, which limits their ability to
efficiently utilize current and next-generation high performance computing
architectures.  REMORA provides an ocean modeling capability that runs on the latest high-performance
computing architectures, from laptops to supercomputers,
whether CPU-only or GPU-accelerated.  In addition, REMORA is based on AMReX,
a modern, well-supported adaptive mesh refinement (AMR) library,
which provides a performance portable interface that shields REMORA
from most of the detailed changes needed to adapt to new systems.
The active and large developer community contributing to AMReX helps ensure
that REMORA will continue to run efficiently as architectures and operating systems
evolve.

# Acknowledgements

Funding for this work was provided by the U.S. Department of Energy
Office of Science.
We acknowledge the help of the AMReX team
in developing and supporting new AMReX features needed by REMORA.
The work at LBNL was supported by the U.S. Department of Energy
under contract No. DE-AC02-05CH11231.
The work at PNNL was supported by the U.S. Department of Energy
under contract No.
The work at ANL was supported by the U.S. Department of Energy
under contract No.

# References
