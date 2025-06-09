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
  - name: Robert Hetland
    orcid: 0000-0001-9531-2119
    corresponding: false
    affiliation: 2
  - name: Jean Sexton
    orcid: 0000-0003-2551-1678
    corresponding: false
    affiliation: 1
  - name: Ann Almgren
    orcid: 0000-0003-2103-312X
    corresponding: false
    affiliation: 1
  - name: Iulian Grindeanu
    orcid: 0000-0002-0264-8236
    corresponding: false
    affiliation: 3
  - name: Kyle Hinson
    orcid: 0000-0002-2737-2379
    corresponding: false
    affiliation: 2
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

date: May 2025

bibliography: paper.bib
---

# Summary

The Regional Model of the Ocean Refined Adaptively (REMORA) is a new implementation
of an existing community standard ocean model,
the Regional Ocean Modeling System
[ROMS, @haidvogel.ea:08; @shchepetkin.mcwilliams:05] that simulates estuarine and oceanic
dynamics using the latest high-performance computing architectures.
REMORA employs hierarchical parallelism using an MPI+X model, where X may be OpenMP on
multicore CPU-only systems, or CUDA, HIP, or SYCL on GPU-accelerated systems.
It is able to be built and run in both single and double precision.
REMORA is built on AMReX [@AMReX:JOSS; @AMReX:IJHPCA],
a block-structured adaptive mesh refinement (AMR) software framework that
provides the underlying performance-portable software infrastructure for block-structured mesh operations.
REMORA, like ROMS, is a regional model, meaning that it is generally used for limited domains,
and as such requires boundary conditions derived analytically, or from larger-scale models.

# Statement of need

Most widely used ocean modeling codes today do not have the
ability to use GPU acceleration, which limits their ability to
efficiently utilize current and next-generation high performance computing
architectures. Oceananigans (Julia-based; @oceananigans) and Veros (Python-based; @veros) have both been developed as flexible,
GPU-enabled models that can be used to simulate regional and global applications.
REMORA is a C++-based alternative. Its ocean modeling capability is directly based on ROMS (a proven Fortran code
that runs efficiently on CPUs), and is able to run on all of the latest high-performance
computing architectures, from laptops to supercomputers, CPU-only or GPU-accelerated.
REMORA is based on AMReX,
a modern, well-supported adaptive mesh refinement (AMR) library,
which provides a performance portable interface that shields REMORA
from most of the detailed changes needed to adapt to new systems. It will also provide a straightforward
avenue for incorporating AMR into REMORA.
The active and large developer community contributing to AMReX helps ensure
that REMORA will continue to run efficiently as architectures and operating systems
evolve.

# REMORA Features

REMORA re-implements the core functionality of ROMS with the AMReX C++ framework.
Unless otherwise indicated, the equations and solvers used are the same in REMORA and ROMS.
REMORA solves the incompressible time-dependent Navier-Stokes equation with the Boussinesq
and hydrostatic approximations [@haidvogel.ea:08; @shchepetkin.mcwilliams:05].
The equations are solved on a curvilinear Arakawa C-grid with a stretched, terrain-following vertical s-coordinate.
We use a split-explicit time-stepping scheme, where several fast barotropic (2D) steps take place within each baroclinic (3D) update [@shchepetkin.mcwilliams:05].
Baroclinic steps are advanced with a third-order Adams-Bashforth scheme, and barotropic steps use a leapfrog predictor followed by a three-time Adams-Moulton corrector.
Scalars are advanced with a leapfrog step and a trapezoidal correction.
Momentum is advected using a third-order upwind (U3) momentum advection scheme.
Tracer advection either uses U3 or a center-difference, fourth-order (C4) scheme.

We have implemented the nonlinear equation of state based on @jackett.macdougall:97.
The user is also provided the option of specifying vertical diffusivity and viscosity analytically, or using the Generic Length Scale (GLS) turbulence closure model [@umlauf:03; @warner:05].
A bulk flux parametrization [@fairall:96; @fairall:03] can optionally be used to calculate surface momentum stress from winds, as well as surface heat flux and effects of evaporation-precipitation.
Bottom drag can be calculated using linear, quadratic, or log-layer prescriptions.
We also provide options for specifying land masking and Coriolis forcing.

REMORA implements many of the ROMS boundary conditions, including periodic and zero-gradient.
Baroclinic variables have the additional option of a radiation boundary condition (e.g. @orlanski:76).
Boundary data provided by file can be used to clamp the solution on the boundaries, or be incorporated using a Chapman/Flather [@flather:76; @chapman:85] condition (in the case of barotropic variables).
Boundary data can also be used for nudging, based on @marchesiello:01.
The solution may also be nudged towards climatology data in an interior sponge region.

Additionally, REMORA provides support for serial I/O with PnetCDF and parallel I/O with AMReX plotfiles (unique to REMORA).

### Next development steps

Subsequent releases of REMORA will include parallel I/O with PnetCDF and point sources and sinks to simulate, e.g. rivers. We will also implement full adaptive mesh refinement. Currently, AMR is tested in REMORA for simple problems such as scalar advection over flat bathymetry. This functionality for more complex problems is a work in progress.

# Acknowledgements

REMORA development is a component of the Study for Exascale Advances in a High-resolution Ocean using ROMS Coupled to E3SM (SEAHORÇE) project funded through the U.S. Department of Energy, Office of Science and Office of Advanced Scientific Computing Research Scientific Discovery through Advanced Computing (SciDAC) program.

We acknowledge the help of the AMReX team
in developing and supporting new AMReX features needed by REMORA.
The work at LBNL was supported by the U.S. Department of Energy
under contract No. DE-AC02-05CH11231.
The work at PNNL was supported by the U.S. Department of Energy
under contract No. DE-AC05-76RL01830.
The work at ANL was supported by the U.S. Department of Energy
under contract No. DE-AC02-06CH11357.

# References
