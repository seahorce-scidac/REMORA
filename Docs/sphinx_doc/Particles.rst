
 .. role:: cpp(code)
    :language: c++

 .. _Particles:

Particles
=========

REMORA can carry Lagrangian particles alongside the mesh-based solution. Particles are placed a fixed
number per cell inside a box and can be advected by the flow. The underlying AMReX particle data
structure is general, so particles may take on other roles in future.

Building and enabling particles
-------------------------------

Build with particles by setting ``USE_PARTICLES = TRUE`` in the GNUmakefile, or by adding
``-DREMORA_ENABLE_PARTICLES:BOOL=ON`` to the cmake command (see, e.g.,
``Build/cmake_with_particles.sh``).

Then turn on one or both particle species at runtime:

::

   remora.use_tracer_particles = true
   remora.use_hydro_particles  = true

Each species, ``tracer_particles`` or ``hydro_particles``, takes parameters prefixed by
``remora.<species>.``, for example ``remora.tracer_particles.particle_box_lo``:

+--------------------------------+-----------------------------------+-------------------+-------------------+
| Parameter                      | Definition                        | Acceptable        | Default           |
|                                |                                   | Values            |                   |
+================================+===================================+===================+===================+
| **initial_distribution_type**  | how particles are initialized     | ``box``           | ``box``           |
+--------------------------------+-----------------------------------+-------------------+-------------------+
| **particle_box_lo**            | lower corner (x, y, z) of the     | 3 Reals           | ``prob_lo``       |
|                                | box in which particles are placed |                   |                   |
+--------------------------------+-----------------------------------+-------------------+-------------------+
| **particle_box_hi**            | upper corner (x, y, z) of the     | 3 Reals           | ``prob_hi``       |
|                                | box in which particles are placed |                   |                   |
+--------------------------------+-----------------------------------+-------------------+-------------------+
| **initial_particles_per_cell** | number of particles placed in     | Integer >= 0      | 1                 |
|                                | each level-0 cell whose center    |                   |                   |
|                                | lies inside the box               |                   |                   |
+--------------------------------+-----------------------------------+-------------------+-------------------+
| **place_randomly_in_cells**    | place particles at random within  | true / false      | true              |
|                                | each cell; if false, at a fixed   |                   |                   |
|                                | point in the cell (for regression |                   |                   |
|                                | tests)                            |                   |                   |
+--------------------------------+-----------------------------------+-------------------+-------------------+
| **advect_with_flow**           | whether the flow moves the        | true / false      | true for tracer   |
|                                | particles                         |                   | particles, false  |
|                                |                                   |                   | for hydro         |
|                                |                                   |                   | particles         |
+--------------------------------+-----------------------------------+-------------------+-------------------+

Particles start on level 0 with zero velocity and are advanced every time step.

A cell receives particles if the physical height of its center lies inside the box. The water column,
from the bed :math:`-h` to the free surface :math:`\zeta`, can extend beyond the domain's z range:
the bathymetry may be deeper than ``prob_lo`` in z, or the free surface may rise above ``prob_hi``.
Cells outside the default box then get no particles. To seed the whole column, give the box a z range
larger than the water column, for example:

::

   remora.tracer_particles.particle_box_lo = 0.     0.     -10000.
   remora.tracer_particles.particle_box_hi = 41000. 80000.   1000.

Particles are still placed only in water cells; the box only selects which cells qualify. The box's
x and y are in the geometry's uniform coordinates (see :ref:`ParticleAdvection`), which on a stretched
or curvilinear grid are not physical metres.

Output
------

Particles are written to AMReX-native plotfiles, but not yet to NetCDF output. To plot the number of
particles per cell as a mesh variable, add ``<species>_count`` (for example
``tracer_particles_count``) to ``remora.plot_vars_3d``; the particle mass density is available the
same way as ``<species>_mass_density``:

::

   remora.plot_vars_3d = salt temp x_velocity y_velocity z_velocity tracer_particles_count

Plotfiles hold the particles' physical positions: x and y are the physical coordinates of the cell
corners (the same ones the plotfile's nodal displacement gives a viewer, see :ref:`sec:Plotfiles`),
interpolated bilinearly to each particle, and z is the physical height. A viewer that draws the mapped
grid therefore draws each particle inside the cell that holds it. Checkpoints instead keep the
geometry's uniform coordinates that the particles move in (see :ref:`ParticleAdvection`), which a
restart needs. The ``xvel`` and ``yvel`` attributes are physical velocities in m/s.

For an example, build the executable in ``Exec`` and run the ``ParticleAdvectionFlat`` case.

.. _ParticleAdvection:

Particle advection
------------------

Particles are advected in the terrain-following coordinates of the grid (see :ref:`VerticalSCoord`).

Horizontal coordinates
~~~~~~~~~~~~~~~~~~~~~~

A particle's horizontal position :math:`(x, y)` is in the uniform coordinates of the AMReX geometry,
with cell size :math:`\Delta x = (x_{hi} - x_{lo})/N_x` and :math:`\Delta y = (y_{hi} - y_{lo})/N_y`
(from ``prob_lo``, ``prob_hi`` and ``n_cell``). These are index coordinates scaled by a constant:
:math:`x` lies in cell :math:`i = \lfloor (x - x_{lo})/\Delta x \rfloor`. The physical width of that cell
is :math:`1/p_m`, where the grid metrics :math:`p_m` and :math:`p_n` vary in space on stretched or
curvilinear grids: those read from a NetCDF grid file, or set analytically with
``remora.grid_scale_type = analytic``. On a uniform grid (``grid_scale_type = constant``, the
default), :math:`p_m = 1/\Delta x` and :math:`p_n = 1/\Delta y`.

Vertical coordinate
~~~~~~~~~~~~~~~~~~~

Let :math:`z_k(x,y)`, :math:`k = 0, \dots, N`, be the heights of the :math:`w`-levels of the column at
:math:`(x,y)`, where :math:`N` is the number of vertical cells, :math:`z_0 = -h` is the bed and
:math:`z_N = \zeta` is the free surface. Each level height is interpolated bilinearly from the heights
at the cell corners:

.. math::
   z_k(x,y) = \sum_{a,b \in \{0,1\}} \lambda_a(\xi)\, \lambda_b(\eta)\, z^{nd}_{i+a,\,j+b,\,k},
   \qquad \lambda_0(t) = 1-t, \quad \lambda_1(t) = t,

where :math:`i = \lfloor (x - x_{lo})/\Delta x \rfloor`, :math:`\xi = (x - x_{lo})/\Delta x - i`, and
likewise :math:`j` and :math:`\eta` in :math:`y`.

A particle's vertical position is a fractional :math:`w`-level index :math:`s = k + f \in [0, N]`,
:math:`0 \le f \le 1`, which places it the fraction :math:`f` of the way from level :math:`k` to level
:math:`k+1`:

.. math::
   z = Z(x,y,s) = z_k(x,y) + f\,\bigl(z_{k+1}(x,y) - z_k(x,y)\bigr).

So :math:`s = 0` is the bed and :math:`s = N` is the free surface.

Equations of motion
~~~~~~~~~~~~~~~~~~~

Each particle obeys

.. math::
   \frac{dx}{dt} = u\, p_m \Delta x, \qquad \frac{dy}{dt} = v\, p_n \Delta y, \qquad
   \frac{ds}{dt} = \frac{\Omega}{H_z},

where

- :math:`u` and :math:`v` are the physical horizontal velocities (m/s) along the grid directions;
- :math:`p_m` and :math:`p_n` are interpolated bilinearly to the particle;
- :math:`H_z = z_{k+1}(x,y) - z_k(x,y)` is the thickness of the particle's cell;
- :math:`\Omega` is the velocity across the :math:`w`-levels (ROMS's ``omega``), from integrating the
  divergence of the horizontal mass fluxes up the column.

The factors :math:`p_m \Delta x` and :math:`p_n \Delta y` convert physical metres to the geometry's
uniform coordinates, so particles cover the right physical distance on stretched or curvilinear grids;
on a uniform grid both are 1. Because :math:`u` and :math:`v` are along the grid directions, rotated
grids need no further correction. :math:`\Omega` comes from the same time-averaged mass fluxes that
advect the tracers, so particles move with the same discrete, volume-conserving flow.

Applying the chain rule to :math:`z = Z(x,y,s,t)` gives the physical vertical velocity:

.. math::
   w = \frac{dz}{dt} = \Omega
     + u \left.\frac{\partial z}{\partial x}\right|_s
     + v \left.\frac{\partial z}{\partial y}\right|_s
     + \left.\frac{\partial z}{\partial t}\right|_s .

Advancing :math:`s` instead of integrating :math:`z` with :math:`w` has two advantages:

- The along-level terms :math:`u\,\partial_x z|_s + v\,\partial_y z|_s` never need computing. Moving a
  particle horizontally at fixed :math:`s` and re-evaluating :math:`Z` raises or lowers it by exactly
  the change in level height.
- :math:`\Omega = 0` at the bed (:math:`s = 0`) and the free surface (:math:`s = N`), so the
  continuous motion stays in :math:`[0, N]`. Integrating :math:`z` with :math:`w` keeps particles in the
  water column only if the interpolated :math:`w` exactly cancels the flow up the sloping bed,
  :math:`u \cdot \nabla h`. After discretization it does not, and particles leak through the bed over
  sloping bathymetry.

Time integration
~~~~~~~~~~~~~~~~

With :math:`\mathbf{q} = (x, y, s)` and
:math:`\mathbf{F}(\mathbf{q}) = (u\, p_m \Delta x,\ v\, p_n \Delta y,\ \Omega / H_z)`, each particle is
advanced with the midpoint method:

.. math::
   \mathbf{q}^{*} = R\left(\mathbf{q}^n + \tfrac{\Delta t}{2}\,\mathbf{F}(\mathbf{q}^n)\right),
   \qquad
   \mathbf{q}^{n+1} = R\left(\mathbf{q}^n + \Delta t\,\mathbf{F}(\mathbf{q}^{*})\right).

Both stages use the end-of-step velocities, so the midpoint stage is in space only. :math:`R` acts on
:math:`s` alone: it reflects :math:`s` off the bed and the free surface,

.. math::
   s \mapsto -s \quad (s < 0), \qquad s \mapsto 2N - s \quad (s > N),

and then clamps it to :math:`[0, N]`. Since :math:`\Omega` vanishes at both boundaries, only the time
discretization can carry a particle past them. After each stage, the particle's height is set to
:math:`z = Z(x, y, s)`.

Velocities are interpolated to the particle trilinearly in index space (:math:`x`, :math:`y`,
:math:`s`), using each field's own staggering. :math:`u` and :math:`v`, defined at vertical cell
centers, are held constant over the half cells next to the bed and the free surface; :math:`\Omega`,
defined on the :math:`w`-levels, is interpolated linearly between them.

After the step, the ``xvel`` and ``yvel`` attributes hold :math:`u` and :math:`v` (m/s, along the grid
directions) at the midpoint, and ``zvel`` holds the rate of change of the particle's height,
:math:`\bigl(z^{n+1} - z^n\bigr) / \Delta t`.

Limitation
~~~~~~~~~~

A particle stores its height :math:`z`, not :math:`s`, so :math:`s` is recovered from :math:`z` on the
current grid at the start of every step. When the free surface moves between steps, a particle
therefore keeps its height rather than its position relative to the levels, dropping the
level-stretching term

.. math::
   \left.\frac{\partial z}{\partial t}\right|_s = \frac{\partial \zeta}{\partial t}\,\frac{z + h}{\zeta + h}.

Storing :math:`s` as a particle attribute would remove this, at the cost of a different checkpoint
layout.
