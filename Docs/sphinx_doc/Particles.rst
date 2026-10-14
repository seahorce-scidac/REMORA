
 .. role:: cpp(code)
    :language: c++

 .. _Particles:

Particles
=========

REMORA has the option to include Lagrangian particles in addition to the mesh-based solution.  Currently the
particle functionality is very simple -- the particles are initialized one per mesh cell
in a particular plane, and are advected by the velocity field.

However, the AMReX particle data structure is very general and particles may take on a number of
different roles in future.

To enable the use of particles, set

::

   USE_PARTICLES = TRUE

in the GNUmakefile if using gmake, or add

::

   -DREMORA_ENABLE_PARTICLES:BOOL=ON \

to the cmake command if using cmake.  (See, e.g., ``Build/cmake_with_particles.sh``)

One must also set

::

   remora.use_tracer_particles = true

in the inputs file or on the command line at runtime.

Currently, by default, the particles are initialized at cell centers, one per cell when the cell index is
(3,j,k), with zero initial velocity.  They are advanced in time every time step.

Particles are advected in the terrain-following coordinates of the grid.  Each particle carries a
fractional vertical index :math:`s \in [0, N]` relative to the :math:`w`-levels of its column, with
:math:`s = 0` the bed and :math:`s = N` the free surface.  Horizontally, particles move with the
interpolated :math:`u` and :math:`v` at fixed :math:`s`, so they follow the levels over sloping
bathymetry; vertically, they move across the levels at :math:`\Omega / H_z`, where :math:`\Omega` is
the velocity across the :math:`w`-levels that the model computes from the horizontal mass fluxes.
The time integration is the midpoint method, and :math:`s` is reflected off the bed and the free
surface, so particles always stay in the water column.  The particle ``zvel`` attribute is the rate of
change of the particle's height over the step.

Caveat: the particle information is currently output when using the AMReX-native plotfile format, but not
when using netcdf.  Writing particles into the netcdf files is a WIP.

To see an example of using the particle functionality, build the executable using gmake in Exec and cd into ParticleAdvectionFlat to run.

To visualize the number of particles per cell as a mesh-based variable, add ``tracer_particle_count`` to the line in the inputs file

::

   remora.plot_vars_1 =

