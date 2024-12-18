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

The Regional Model of the Ocean Refined Adaptively (REMORA) is a new implimentation 
of an existing community standard ocean model, 
the Regional Ocean Modeling System 
(ROMS, [@shchepetkin.mcwilliams:05], [@haidvogel.ea:08]) that simulates esutuarine and oceanic 
dynamics using the latest high-performance computing architectures.
REMORA employs hierarchical parallelism using an MPI+X model, where X may be OpenMP on
multicore CPU-only systems, or CUDA, HIP, or SYCL on GPU-accelerated systems.
It is able to be built and run in both single and double precision.
REMORA is built on AMReX [@AMReX:JOSS; @AMReX:IJHPCA],
a block-structured adaptive mesh refinement (AMR) software framework that
provides the underlying performance-portable software infrastructure for block-structured mesh operations.
REMORA, like ROMS, is a regional model, meaning that it is generally used for limited domains, 
and as such requires boundary conditions derived analytically, or from larger-scale models.

# REMORA Features

Like ROMS, REMORA:
 - solves the incompressible time-dependent Navier-Stokes equation with the Boussinesq and hydrostatic approximations (see [@shchepetkin.mcwilliams:05], [@haidvogel.ea:08])
 - uses a curvilinear Arakawa C-grid
 - uses a streched, terrain-following vertical s-coordinate
 - uses a split-explicit time-stepping scheme, where several fast barotropic (2D) steps take place within each baroclinic (3D) update (see [@shchepetkin.mcwilliams:05])
 - baroclinic steps use a third-order Adams-Bashforth scheme
 - barotropic steps use a leapfrog predictor followed by a three-time Adams-Moulton corrector
 - scalars are advanced with a leapfrog step with a trapezoidal correction
 - uses a nonlinear equation of state based on [@jackett.macdougall:97]
 - uses a third-order (U3) upwind momentum advection scheme
 - uses U3 or center-difference, fourth-order (C4) tracer advection
 - uses analytical vertical diffusivity or viscosity, or uses the Generic Length Scale (GLS) turbulence closure model ([@umlauf:03], [@warner.ea:05]).
 - specified land masking, Coriolis force, and ealistic wind stress
 - uses quadratic or log-layer bottom drag
 - uses periodic, radiation (e.g., @orlanski:76), or clamped time-varying baroclinic lateral boundary conditions
 - uses radiation, chapman/flather (@flather:76, @chapman:85), or clamped barotropic lateral boundary conditions
 - uses optional boundary nudging based on @marchesiello:01
 - uses parallel I/O with netcdf (using the pnetcdf library), or plotfiles (based on yt)

### Next development steps
 - Surface heat flux parameterizations
 - Rivers
 - Evaporation-Precipitation fluxes
 - Nudging to climatology
 - Adaptive mesh refinement


# Statement of need

Most widely used ocean modeling codes today do not have the
ability to use GPU acceleration, which limits their ability to
efficiently utilize current and next-generation high performance computing
architectures.  REMORA provides an ocean modeling capability, based on a proven FORTRAN code
that runs efficiently on CPUs, that is able to run on all of the latest high-performance
computing architectures, from laptops to supercomputers, CPU-only or GPU-accelerated.  
In addition, REMORA is based on AMReX,
a modern, well-supported adaptive mesh refinement (AMR) library,
which provides a performance portable interface that shields REMORA
from most of the detailed changes needed to adapt to new systems.
The active and large developer community contributing to AMReX helps ensure
that REMORA will continue to run efficiently as architectures and operating systems
evolve.

# Acknowledgements

REMORA development is a component of the Study for Exascale Advances in a High-resolution Ocean using ROMS Coupled to E3SM (SEAHORÇE) project funded through the U.S. Department of Energy, Office of Science and Office of Advanced Scientific Computing Research Scientific Discovery through Advanced Computing (SciDAC) program.

We acknowledge the help of the AMReX team
in developing and supporting new AMReX features needed by REMORA.
The work at LBNL was supported by the U.S. Department of Energy
under contract No. DE-AC02-05CH11231.
The work at PNNL was supported by the U.S. Department of Energy
under contract No. DE-AC05-76RL01830
The work at ANL was supported by the U.S. Department of Energy
under contract No. DE-AC02-06CH11357

# References
