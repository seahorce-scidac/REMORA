
 .. role:: cpp(code)
    :language: c++

 .. _Particles:

Particles
=========

REMORA has the option to include Lagrangian particles in addition to the mesh-based solution.  Currently the
particle functionality is very simple -- the particles are initialized a fixed number per mesh cell
inside a box, and can be advected by the velocity field.

However, the AMReX particle data structure is very general and particles may take on a number of
different roles in future.

To enable the use of particles, set

::

   USE_PARTICLES = TRUE

in the GNUmakefile if using gmake, or add

::

   -DREMORA_ENABLE_PARTICLES:BOOL=ON \

to the cmake command if using cmake.  (See, e.g., ``Build/cmake_with_particles.sh``)

One must also turn on one or both of the particle species at runtime, in the inputs file or on the
command line:

::

   remora.use_tracer_particles = true
   remora.use_hydro_particles  = true

The species are named ``tracer_particles`` and ``hydro_particles``. Each species is configured with
parameters prefixed by ``remora.<species>.``, for example ``remora.tracer_particles.particle_box_lo``:

+--------------------------------+-----------------------------------+-------------------+-------------------+
| Parameter                      | Definition                        | Acceptable        | Default           |
|                                |                                   | Values            |                   |
+================================+===================================+===================+===================+
| **initial_distribution_type**  | how particles are initialized     | ``box``           | ``box``           |
+--------------------------------+-----------------------------------+-------------------+-------------------+
| **particle_box_lo**            | lower corner (x, y, z) of the     | 3 Reals           | domain lower      |
|                                | region in which particles are     |                   | corner            |
|                                | initialized                       |                   |                   |
+--------------------------------+-----------------------------------+-------------------+-------------------+
| **particle_box_hi**            | upper corner (x, y, z) of the     | 3 Reals           | domain upper      |
|                                | region in which particles are     |                   | corner            |
|                                | initialized                       |                   |                   |
+--------------------------------+-----------------------------------+-------------------+-------------------+
| **place_randomly_in_cells**    | place particles at random         | true / false      | true              |
|                                | positions within each cell; if    |                   |                   |
|                                | false, at a fixed position in     |                   |                   |
|                                | the cell (for regression tests)   |                   |                   |
+--------------------------------+-----------------------------------+-------------------+-------------------+
| **initial_particles_per_cell** | number of particles placed in     | Integer >= 0      | 1                 |
|                                | each level-0 cell whose center    |                   |                   |
|                                | lies inside the box               |                   |                   |
+--------------------------------+-----------------------------------+-------------------+-------------------+
| **advect_with_flow**           | whether the particles are moved   | true / false      | true for          |
|                                | by the flow velocity              |                   | tracer particles, |
|                                |                                   |                   | false for hydro   |
|                                |                                   |                   | particles         |
+--------------------------------+-----------------------------------+-------------------+-------------------+

Particles are initialized on level 0 with zero initial velocity and are advanced every time step.
The default box spans ``prob_lo`` to ``prob_hi``, but a cell is tested by the physical height of its
center, so with bathymetry deeper than the domain's ``prob_lo`` in z, or a free surface above
``prob_hi``, the cells outside that range get no particles unless the box is widened.

Caveat: the particle information is currently output when using the AMReX-native plotfile format, but not
when using netcdf.  Writing particles into the netcdf files is a WIP.

To see an example of using the particle functionality, build the executable using gmake in Exec and cd into ParticleAdvectionFlat to run.

To visualize the number of particles per cell as a mesh-based variable, add ``<species>_count`` (for
example ``tracer_particles_count``) to ``remora.plot_vars_3d`` in the inputs file. The particle mass
density is available the same way as ``<species>_mass_density``:

::

   remora.plot_vars_3d = salt temp x_velocity y_velocity z_velocity tracer_particles_count
