 .. role:: cpp(code)
    :language: c++

.. _ROMSXvsWRF:

ROMSX vs WRF
===============

The following comparison is based on the WRF Version 4 Technical Report, titled
"A Description of the Advancd Research WRF Model Version 4"

Similarities
--------------------

**Equations**: both ROMSX and WRF solve the fully-compressible, Eulerian nonhydrostatic equations, and conserve
dry air mass and scalar mass.  ROMSX does not have a hydrostatic option.

**Prognostic Variables**: velocity components (u,v,w); perturbation moist potential temperature.  Optionally,
turbulent kinetic energy and any number of scalars such as water vapor mixing ratio, rain/snow mixing ratio,
cloud water / ice mixing ratio.

**Horizontal grid**: both ROMSX and WRF use Arakawa C-grid staggering.

**Time Integration**: Time-split integration using 3rd-order Runge-Kutta scheme with smaller time step for
acoustic and gravity wave modes.  Variable time step capability.

**Spatial Discretization**: 2nd- to 6th-order advection options in horizontal and vertical

**Turbulent Mixing**: Sub-grid scale turbulence formulation.  Vertically implicit acoustic step off-centering.

**Initial conditions**: both ROMSX and WRF have the ability to initialize problems from
3D "real" data (output of real.exe), "ideal" data (output of ideal.exe) and from 1D input soundings.

**Lateral boundary conditions**: Periodic, open, symmetric and specified (in wrfbdy* files).

**Bottom boundary conditions**: Frictional or free-slip

**Earth's Rotation**: Coriolis terms in ROMSX controlled by run-time input flag

**Mapping to Sphere**: ROMSX supports the use of map scale factors for isotropic projections (read in from
wrfinput files).

**Nesting**: One-way or two-way.  Multiple levels and integer ratios.



Key Differences
--------------------

**Vertical coordinates**: Unlike WRF, ROMSX uses a terrain-following height-based vertical coordinate,
with vertical grid stretching permitted.

**Time Integration**: ROMSX supports using a 3rd-order Runge-Kutta scheme with no substepping as alternative to RK3 with acoustic substepping.

**Initial conditions**: ROMSX has an additional mode of "custom" initialization in which
the user writes the initialization routine.

ROMSX does **not* have the capability for global simulation

