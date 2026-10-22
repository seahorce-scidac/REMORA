
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

How the particles move is described in :ref:`ParticleAdvection` below.

Caveat: the particle information is currently output when using the AMReX-native plotfile format, but not
when using netcdf.  Writing particles into the netcdf files is a WIP.

To see an example of using the particle functionality, build the executable using gmake in Exec and cd into ParticleAdvectionFlat to run.

To visualize the number of particles per cell as a mesh-based variable, add ``<species>_count`` (for
example ``tracer_particles_count``) to ``remora.plot_vars_3d`` in the inputs file. The particle mass
density is available the same way as ``<species>_mass_density``:

::

   remora.plot_vars_3d = salt temp x_velocity y_velocity z_velocity tracer_particles_count

.. _ParticleAdvection:

Particle advection
------------------

Particles are advected in the terrain-following coordinates of the grid (see :ref:`VerticalSCoord`).

Vertical coordinate
~~~~~~~~~~~~~~~~~~~

Let :math:`z_k(x,y)`, :math:`k = 0, \dots, N`, be the heights of the :math:`w`-levels of the column at
horizontal position :math:`(x,y)`, where :math:`N` is the number of vertical cells,
:math:`z_0 = -h` is the bed and :math:`z_N = \zeta` is the free surface. Each level height is the
bilinear interpolation of the heights at the cell corners:

.. math::
   z_k(x,y) = \sum_{a,b \in \{0,1\}} \lambda_a(\xi)\, \lambda_b(\eta)\, z^{nd}_{i+a,\,j+b,\,k},
   \qquad \lambda_0(t) = 1-t, \quad \lambda_1(t) = t,

where :math:`i = \lfloor (x - x_{lo})/\Delta x \rfloor`, :math:`\xi = (x - x_{lo})/\Delta x - i`, and
likewise :math:`j` and :math:`\eta` in :math:`y`.

A particle's vertical position is described by a fractional :math:`w`-level index
:math:`s = k + f \in [0, N]`, :math:`0 \le f \le 1`, which places it the fraction :math:`f` of the way
from level :math:`k` to level :math:`k+1`:

.. math::
   z = Z(x,y,s) = z_k(x,y) + f\,\bigl(z_{k+1}(x,y) - z_k(x,y)\bigr).

So :math:`s = 0` is the bed and :math:`s = N` is the free surface.

Equations of motion
~~~~~~~~~~~~~~~~~~~

In these coordinates each particle obeys

.. math::
   \frac{dx}{dt} = u, \qquad \frac{dy}{dt} = v, \qquad \frac{ds}{dt} = \frac{\Omega}{H_z},

where :math:`u` and :math:`v` are the horizontal velocities, :math:`H_z = z_{k+1}(x,y) - z_k(x,y)` is the
thickness of the particle's cell, and :math:`\Omega` is the velocity across the :math:`w`-levels (ROMS's
``omega``), obtained by integrating the divergence of the horizontal mass fluxes up the column.
:math:`\Omega` is built from the same time-averaged mass fluxes that the tracers are advected with, so
the particles move with the same discrete, volume-conserving flow as the tracers.

The physical vertical velocity follows from applying the chain rule to :math:`z = Z(x,y,s,t)`:

.. math::
   w = \frac{dz}{dt} = \Omega
     + u \left.\frac{\partial z}{\partial x}\right|_s
     + v \left.\frac{\partial z}{\partial y}\right|_s
     + \left.\frac{\partial z}{\partial t}\right|_s .

Advancing :math:`s`, rather than integrating :math:`z` with :math:`w`, has two consequences:

- The along-level terms :math:`u\,\partial_x z|_s + v\,\partial_y z|_s` never have to be computed.
  Moving a particle horizontally at fixed :math:`s` and re-evaluating :math:`Z` at its new position
  raises or lowers it by exactly the change in level height.
- :math:`\Omega = 0` at the bed (:math:`s = 0`) and at the free surface (:math:`s = N`), so the
  continuous motion cannot leave :math:`[0, N]`. Integrating :math:`z` with :math:`w` instead keeps
  particles in the water column only if the interpolated :math:`w` exactly cancels the flow up the
  sloping bed, :math:`u \cdot \nabla h`, which it does not after discretization. Particles then leak
  through the bed over sloping bathymetry.

Time integration
~~~~~~~~~~~~~~~~

With :math:`\mathbf{q} = (x, y, s)` and :math:`\mathbf{F}(\mathbf{q}) = (u,\ v,\ \Omega / H_z)`, each
particle is advanced with the midpoint method:

.. math::
   \mathbf{q}^{*} = R\left(\mathbf{q}^n + \tfrac{\Delta t}{2}\,\mathbf{F}(\mathbf{q}^n)\right),
   \qquad
   \mathbf{q}^{n+1} = R\left(\mathbf{q}^n + \Delta t\,\mathbf{F}(\mathbf{q}^{*})\right).

Both stages use the velocities at the end of the step, so the midpoint stage is in space only.
:math:`R` acts only on :math:`s`. It reflects :math:`s` off the bed and the free surface,

.. math::
   s \mapsto -s \quad (s < 0), \qquad s \mapsto 2N - s \quad (s > N),

and then clamps it to :math:`[0, N]`. Because :math:`\Omega` vanishes at both boundaries, only the time
discretization can carry a particle past them. After each stage the particle's height is set to
:math:`z = Z(x, y, s)`.

The velocities are interpolated to the particle trilinearly in index space, in :math:`x`, :math:`y` and
:math:`s`, using each field's own staggering. Values defined at cell centers in the vertical
(:math:`u` and :math:`v`) are held constant over the half cells next to the bed and the free surface,
and :math:`\Omega`, defined on the :math:`w`-levels, is interpolated linearly between them.

After the step, the particle's ``xvel`` and ``yvel`` attributes hold :math:`u` and :math:`v` at the
midpoint, and ``zvel`` holds the rate of change of its height over the step,
:math:`\bigl(z^{n+1} - z^n\bigr) / \Delta t`.

Limitation
~~~~~~~~~~

A particle stores its height :math:`z`, not :math:`s`, so :math:`s` is recovered from :math:`z` at the
start of every step, using the current grid. While the free surface moves between steps, a particle
therefore keeps its height rather than its position relative to the levels. This drops the
level-stretching term

.. math::
   \left.\frac{\partial z}{\partial t}\right|_s = \frac{\partial \zeta}{\partial t}\,\frac{z + h}{\zeta + h}

from its motion. Storing :math:`s` as a particle attribute would remove this, at the cost of a
different checkpoint layout.
