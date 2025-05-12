.. role:: cpp(code)
  :language: c++

.. _sec:Inputs:

******
Inputs
******
.. toctree::
   :maxdepth: 1

The REMORA executable reads run-time information from an “inputs” file which you put on the command line.
This section describes the inputs which can be specified either in the inputs file or on the command line.
If a value is specified on the command line, that value will override a value specified in the inputs file.

.. _geometry-parameters:

Problem Geometry
================

The problem geometry can be specified either by a NetCDF grid file or in the inputs.
Instructions for setting grid, initial, and Dirichelet boundary conditions from NetCDF file can be found
:ref:`here<icbc-parameters>`. If one of these quantities is specified from a NetCDF file, they all must be. Even if the grid is specified in the NetCDF file, the geometry and boundary parameters below must still be set. ``remora.prob_lo`` and ``remora.prob_hi`` do not need to agree with the file in this case.

The z-component of ``remora.prob_lo`` should be more negative than the deepest bathymetry, and the z-compoonent of ``remora.prob_hi`` should be 0.

List of Parameters
------------------

+--------------------------+-----------------+-------------------+-------------+
| Parameter                | Definition      | Acceptable        | Default     |
|                          |                 | Values            |             |
+==========================+=================+===================+=============+
| **remora.prob_lo**       | physical        | [Real Real -Real] | must be set |
|                          | location of low |                   |             |
|                          | corner of the   |                   |             |
|                          | domain          |                   |             |
+--------------------------+-----------------+-------------------+-------------+
| **remora.prob_hi**       | physical        | [Real Real 0]     | must be set |
|                          | location of     |                   |             |
|                          | high corner of  |                   |             |
|                          | the domain      |                   |             |
+--------------------------+-----------------+-------------------+-------------+
| **remora.is_periodic**   | is the domain   | 0 if false, 1     | 0 0 0       |
|                          | periodic in     | if true.          |             |
|                          | this direction  | Z-component must  |             |
|                          |                 | be zero           |             |
+--------------------------+-----------------+-------------------+-------------+

Examples of Usage
-----------------

-  **remora.prob_lo** = 0 0 -200
   defines the low corner of the domain at (0 m,0 m,-200 m) in physical space.

-  **remora.prob_hi** = 1.e8 2.e8 0
   defines the high corner of the domain at (1.e8 m, 2.e8 m, 0 m) in
   physical space.

-  **remora.is_periodic** = 0 1 0
   says the domain is periodic in the y-direction only.

Domain Boundary Conditions
==========================

Instructions for how to specify domain boundary conditions, with usage examples can be found in :ref:`Domain Boundary Conditions <sec:domainBCs>`.

.. _icbc-parameters:

Imposing Boundary and Initial Conditions from NetCDF File
=========================================================

Grid, initial, and time-dependent boundary data can be specified using NetCDF files, as in ROMS. REMORA expects files in the same format as ROMS, in NetCDF classic format (32- or 64-bit). If the file format is incorrect, REMORA will exit with the error:

.. code:: shell

   NetCDF: Attempt to use feature that was not turned on when netCDF was built.

Other versions of NetCDF files can be converted to 64-bit NetCDF classic by running the command:

.. code:: shell

    ncks -5 old_file.nc converted_file.nc

The utility ``ncks`` is part of the `NCO <https://nco.sourceforge.net>`_ suite.

Currently, if one initial, grid, or boundary files are specified, they all must be. Boundary condition options with NetCDF boundary data are equivalent to ROMS clamped, Chapman-Flather, and Orlanski + Nudging boundary conditions. Options and examples can be found in the section on :ref:`Domain Boundary Conditions <sec:domainBCs>`.

List of Parameters
------------------

+-------------------------------+-----------------------------------+-------------+---------------------------+
| Parameter                     | Definition                        | Acceptable  | Default                   |
|                               |                                   | Values      |                           |
+===============================+===================================+=============+===========================+
| **remora.ic_bc_type**         | read initial, grid, and           |             |                           |
|                               | boundary data from NetCDF         | true/false  | false                     |
|                               | files                             |             |                           |
+-------------------------------+-----------------------------------+-------------+---------------------------+
| **remora.nc_init_file_0**     | initial data NetCDF file name     | string      | must be set               |
|                               |                                   |             | if ``remora.ic_bc_type``  |
|                               |                                   |             | is true                   |
+-------------------------------+-----------------------------------+-------------+---------------------------+
| **remora.nc_grid_file_0**     | grid data NetCDF file name        | string      | must be set               |
|                               |                                   |             | if ``remora.ic_bc_type``  |
|                               |                                   |             | is true                   |
+-------------------------------+-----------------------------------+-------------+---------------------------+
| **remora.nc_bdry_file_0**     | boundary data NetCDF file         | string      | must be set               |
|                               | name                              |             | if ``remora.ic_bc_type``  |
|                               |                                   |             | is true                   |
+-------------------------------+-----------------------------------+-------------+---------------------------+
| **remora.nc_frc_file**        | forcing data NetCDF file name     | string      | must be set               |
|                               |                                   |             | if ``remora.wind_type``   |
|                               |                                   |             | or ``remora.smflux_type`` |
|                               |                                   |             | equal ``netcdf``          |
+-------------------------------+-----------------------------------+-------------+---------------------------+
| **remora.bdy_time_varname**   | name of time variable in boundary | string      | ``ocean_time``            |
|                               | file                              |             |                           |
+-------------------------------+-----------------------------------+-------------+---------------------------+
| **remora.frc_time_varname**   | name of time variable in forcing  | string      | ``wind_time`` for wind,   |
|                               | file                              |             | ``sms_time`` for surface  |
|                               |                                   |             | momentum stress           |
+-------------------------------+-----------------------------------+-------------+---------------------------+

Resolution and Tiling
=====================

The REMORA gridding and load balancing strategy is based on that in AMReX.
See the `Gridding`_ section of the AMReX documentation for details.

.. _`Gridding`: https://amrex-codes.github.io/amrex/docs_html/ManagingGridHierarchy_Chapter.html

.. _list-of-parameters-2:

List of Parameters
------------------

+---------------------------+-----------------+---------------------------+-----------------------+
| Parameter                 | Definition      | Acceptable                | Default               |
|                           |                 | Values                    |                       |
+===========================+=================+===========================+=======================+
| **remora.n_cell**         | number of cells | Triplet of integers > 0   | must be set           |
|                           | in each         |                           |                       |
|                           | direction at    | {x,y,z}                   |                       |
|                           | the coarsest    |                           |                       |
|                           | level           |                           |                       |
+---------------------------+-----------------+---------------------------+-----------------------+
| **remora.omp_tile_size**  | target tile     | Triplet of integers > 0   | CPU: 8 8 1024         |
|                           | size            |                           |                       |
|                           |                 | {x,y,z}                   | GPU: 1024 1024 1024   |
+---------------------------+-----------------+---------------------------+-----------------------+

Notes
-----

-  **remora.omp_tile_size** is an alias for **fabarray.mfiter_tile_size**, which controls
   the distribution of work by OpenMP.

-  The domain may not be tiled in the z direction. That is, the last component of **remora.omp_tile_size**
   must be greater than or equal to the number of vertical levels.

.. _examples-of-usage-2:

Examples of Usage
-----------------

-  **remora.n_cell** = 32 64 64

   would define the domain to have 32 cells in the x-direction, 64 cells
   in the y-direction, and 64 cells in the z-direction *at the coarsest level*.

Mesh Refinement and (Re)gridding
================================

Overview
--------

The user defines how to tag individual cells at a given level for refinement.
This list of tagged cells is sent to a grid generation routine, which uses the
Berger–Rigoutsos algorithm to create rectangular grids that contain the
tagged cells.

See :ref:`MeshRefinement` for more details on how to specify regions for
refinement.

Note that because these arguments are primarily used by AMReX and pertain to adaptive mesh
refinement, they use the prefix ``amr``.

.. _list-of-parameters-4:

List of Parameters
------------------

+----------------------------+-----------------+----------------+----------------+
| Parameter                  | Definition      | Acceptable     | Default        |
|                            |                 | Values         |                |
+============================+=================+================+================+
| **amr.max_level**          | number of       | Integer >= 0   | must be set    |
|                            | levels of       |                |                |
|                            | refinement      |                |                |
|                            | above the       |                |                |
|                            | coarsest level  |                |                |
+----------------------------+-----------------+----------------+----------------+
| **amr.ref_ratio**          | ratio of coarse | 2 / 3 / 4      | 2 for all      |
|                            | to fine grid    | (one per level)| levels         |
|                            | spacing between |                |                |
|                            | subsequent      |                |                |
|                            | levels          |                |                |
+----------------------------+-----------------+----------------+----------------+
| **amr.ref_ratio_vect**     | ratio of coarse | 3 integers     | 2 for all      |
|                            | to fine grid    | (one per dir)  | directions     |
|                            | spacing between | 2 / 3 / 4      |                |
|                            | subsequent      |                |                |
|                            | levels          |                |                |
+----------------------------+-----------------+----------------+----------------+
| **amr.regrid_int**         | how often to    | Integer > 0    | must be set    |
|                            | regrid          |                |                |
+----------------------------+-----------------+----------------+----------------+
| **amr.regrid_on_restart**  | should we       | 0 or 1         | 0              |
|                            | regrid          |                |                |
|                            | immediately     |                |                |
|                            | after           |                |                |
|                            | restarting      |                |                |
+----------------------------+-----------------+----------------+----------------+
| **amr.regrid_file**        | name of file    | text           | no file        |
|                            | from which to   |                |                |
|                            | read the grids  |                |                |
+----------------------------+-----------------+----------------+----------------+
| **amr.grid_eff**           | grid            | Real > 0, < 1  | 0.7            |
|                            | efficiency at   |                |                |
|                            | coarse level    |                |                |
|                            | at which grids  |                |                |
|                            | are created     |                |                |
+----------------------------+-----------------+----------------+----------------+
| **amr.n_error_buf**        | radius of       | Integer >= 0   | 1              |
|                            | additional      |                |                |
|                            | tagging around  |                |                |
|                            | already tagged  |                |                |
|                            | cells           |                |                |
+----------------------------+-----------------+----------------+----------------+
| **amr.blocking_factor**    | grid size must  | Integer > 0    | 2              |
|                            | be a multiple   |                |                |
|                            | of this         |                |                |
+----------------------------+-----------------+----------------+----------------+
| **amr.refine_grid_layout** | refine grids    | 0 if false, 1  | 1              |
|                            | more if # of    | if true        |                |
|                            | processors      |                |                |
|                            | :math:`>` # of  |                |                |
|                            | grids           |                |                |
+----------------------------+-----------------+----------------+----------------+
| **amr.do_substep**         | whether to sub- | 0 if false, 1  | 0              |
|                            | step finer      | if true        |                |
|                            | levels in time  |                |                |
|                            |                 | NOTE: true     |                |
|                            |                 | will trigger   |                |
|                            |                 | Assert failure |                |
+----------------------------+-----------------+----------------+----------------+

.. _notes-2:

Notes
-----

-  if **amr.max_level** = 0 then you do not need to set
   **amr.ref_ratio** or **amr.regrid_int**.

-  **amr.n_error_buf**, **remora.max_grid_size** and
   **amr.blocking_factor** can be read in as a single value which is
   assigned to every level, or as multiple values, one for each level

-  **amr.max_grid_size** at every level must be even

-  **amr.blocking_factor** at every level must be a power of 2

-  the domain size **remora.n_cell** must be a multiple of
   **amr.blocking_factor** at level 0

-  **amr.max_grid_size** must be a multiple of **amr.blocking_factor**
   at every level

-  the substepping turned on by **amr.do_substep** is NOT implemented yet so will trigger an Assert.

.. _examples-of-usage-3:

Examples of Usage
-----------------

-  | **amr.max_level** = 2
   | would allow a maximum of 2 refined levels in addition to the coarse
     level. Note that these additional levels will only be created only
     if the tagging criteria are such that cells are flagged as needing
     refinement. The number of refined levels in a calculation must be
     :math:`\leq` **amr.max_level**, but can change in time and need not
     always be equal to **amr.max_level**.

-  | **amr.ref_ratio** = 2 3
   | would set factor 2 refinement between levels 0 and 1, and factor 3
     refinement between levels 1 and 2. Note that you must have at least
     **amr.max_level** values of **amr.ref_ratio** (Additional values
     may appear in that line and they will be ignored).

-  | **amr.ref_ratio_vect** = 2 4 3
   | would set factor {2 in x-dir, 4 in y-dir, 3 in z-dir} refinement between
     all adjacent levels.    Note that you must specify 3 values, one for
     each coordinate direction.

-  | **amr.regrid_int** = 2 2
   | tells the code to regrid every 2 steps. Thus in this example, new
     level-1 grids will be created every 2 level-0 time steps, and new
     level-2 grids will be created every 2 level-1 time steps.

-  | **amr.regrid_file** = *fixed_grids*
   | In this case the list of grids at each fine level are contained in
     the file *fixed_grids*, which will be read during the gridding
     procedure. These grids must not violate the **amr.max_grid_size**
     criterion. The rest of the gridding procedure described below will
     not occur if **amr.regrid_file** is set.

-  | **amr.grid_eff** = 0.9
   | During the grid creation process, at least 90% of the cells in each
     grid at the level at which the grid creation occurs must be tagged
     cells. Note that this is applied at the coarsened level at which
     the grids are actually made, and before **amr.max_grid_size** is
     imposed.

-  | **amr.max_grid_size** = 64
   | The final grids will be no longer than 64 cells on a side at every
     level.

-  | **amr.max_grid_size** = 64 32 16
   | The final grids will be no longer than 64 cells on a side at level
     0, 32 cells on a side at level 1, and 16 cells on a side at level
     2.

-  | **amr.blocking_factor** = 32
   | The dimensions of all the final grids will be multiples of 32 at
     all levels.

-  | **amr.blocking_factor** = 32 16 8
   | The dimensions of all the final grids will be multiples of 32 at
     level 0, multiples of 16 at level 1, and multiples of 8 at level 2.

Simulation Time
===============

.. _list-of-parameters-5:

List of Parameters
------------------

+------------------------+---------------------------+--------------+---------+
| Parameter              | Definition                | Acceptable   | Default |
|                        |                           | Values       |         |
+========================+===========================+==============+=========+
| **remora.max_step**    | maximum number of level 0 | Integer >= 0 | -1      |
|                        | time steps                |              |         |
+------------------------+---------------------------+--------------+---------+
| **remora.stop_time**   | final simulation          | Real >= 0    | -1.0    |
|                        | time                      |              |         |
+------------------------+---------------------------+--------------+---------+
| **remora.start_time**  | initial simulation        | Real >= 0    | 0.0     |
|                        | time                      |              |         |
+------------------------+---------------------------+--------------+---------+

.. _notes-3:

Notes
-----

To control the number of time steps, you can limit by the maximum number
of level-0 time steps (**remora.max_step**), or the final simulation time
(**remora.stop_time**), or both. The code will stop at whichever criterion
comes first. Note that if the code reaches **remora.stop_time** then the final
time step will be shortened so as to end exactly at **remora.stop_time**, not
pass it.

.. _examples-of-usage-4:

Examples of Usage
-----------------

-  **remoraw.max_step** = 1000

-  **remora.stop_time** = 1.0

will end the calculation when either the simulation time reaches 1.0 or
the number of level-0 steps taken equals 1000, whichever comes first.

Time Step
=========

.. _list-of-parameters-6:

List of Parameters for Single-Rate
----------------------------------

+-------------------------------+-----------------+----------------+----------------------------+
| Parameter                     | Definition      | Acceptable     | Default                    |
|                               |                 | Values         |                            |
+===============================+=================+================+============================+
| **remora.cfl**                | CFL number for  | Real > 0 and   | 0.8                        |
|                               | hydro           | <= 1           |                            |
+-------------------------------+-----------------+----------------+----------------------------+
| **remora.fixed_dt**           | set level 0 dt  | Real > 0       | unused if not              |
|                               | as this value   |                | set                        |
|                               | regardless of   |                |                            |
|                               | cfl or other    |                |                            |
|                               | settings        |                |                            |
+-------------------------------+-----------------+----------------+----------------------------+
| **remora.fixed_fast_dt**      | set fast dt     | real > 0       | inferred from **fixed_dt** |
|                               | as this value   |                | and **fixed_ndfast_ratio** |
|                               |                 |                | if not set                 |
+-------------------------------+-----------------+----------------+----------------------------+
| **remora.fixed_ndfast_ratio** | set fast dt     | int            | inferred from **fixed_dt** |
|                               | as slow dt /    |                | and **fixed_fast_dt**      |
|                               | this ratio      |                | if not set                 |
+-------------------------------+-----------------+----------------+----------------------------+
| **remora.change_max**         | factor by which | Real >= 1      | 1.1                        |
|                               | dt can grow     |                |                            |
|                               | in subsequent   |                |                            |
|                               | steps           |                |                            |
+-------------------------------+-----------------+----------------+----------------------------+

.. _examples-of-usage-5:

Examples of Usage
-----------------

-  | **remora.cfl** = 0.9
   | defines the timestep as dt = cfl \* dx / (u+c).  Only relevant if **fixed_dt** not set

-  | **remora.change_max** = 1.1
   | allows the time step to increase by no more than 10% in this case.
     Note that the time step can shrink by any factor; this only
     controls the extent to which it can grow.

-  | **remora.fixed_dt** = 1.e-4
   | sets the level-0 time step to be 1.e-4 for the entire simulation,
     ignoring the other timestep controls.

Restart Capability
==================

See :ref:`sec:Checkpoint` for how to control the checkpoint/restart capability.

PlotFiles
===============================

See :ref:`sec:Plotfiles` for how to control the types and frequency of plotfile
generation.


Screen Output
=============

.. _list-of-parameters-10:

List of Parameters
------------------

+----------------------------+------------------+----------------------------------+----------------+
| Parameter                  | Definition       | Acceptable                       | Default        |
|                            |                  | Values                           |                |
+============================+==================+==================================+================+
| **amr.v**                  | verbosity of     | 0 or 1                           | remora.v       |
|                            | Amr.cpp          |                                  |                |
+----------------------------+------------------+----------------------------------+----------------+
| **remora.v**               | verbosity of     | - 0: none                        | 0              |
|                            | REMORA functions | - 1: integrated/max quantities   |                |
|                            |                  | - 2: print boxes                 |                |
+----------------------------+------------------+----------------------------------+----------------+
| **remora.sum_interval**    | how often (in    |                                  |                |
|                            | level-0 time     |                                  |                |
|                            | steps)           |                                  |                |
|                            | to compute       | Integer                          | -1             |
|                            | integral         |                                  |                |
|                            | quantities       |                                  |                |
+----------------------------+------------------+----------------------------------+----------------+

.. _examples-of-usage-9:

Examples of Usage
-----------------

-  | **remora.sum_interval** = 2
   | if **remora.sum_interval** :math:`> 0` then the code computes and
     prints certain integral quantities, such as total mass, momentum
     and energy in the domain every **remora.sum_interval** level-0 steps.
     In this example the code will print these quantities every two
     coarse time steps. The print statements have the form
   | ``TIME= 1.91717746 MASS= 1.792410279e+34``
   | for example. If this line is commented out or if **remora.v** :math:`<= 0`
     then it will not compute and print these quantities.

Included terms
==============

.. _list-of-parameters-14:

List of Parameters
------------------

+----------------------------------+-----------------------------+-------------------+-------------+
| Parameter                        | Definition                  | Acceptable        | Default     |
|                                  |                             | Values            |             |
+==================================+=============================+===================+=============+
| **remora.use_coriolis**          | Include Coriolis terms.     | true / false      | false       |
+----------------------------------+-----------------------------+-------------------+-------------+
| **remora.flat_bathymetry**       | Use flat bathymetry.        | true / false      | false       |
+----------------------------------+-----------------------------+-------------------+-------------+
| **remora.use_prestep**           | Do prestep terms. Only for  |  true / false     | true        |
|                                  | debugging purposes.         |                   |             |
+----------------------------------+-----------------------------+-------------------+-------------+
| **remora.use_uv3dmix**           | Include harmonic viscosity. | true / false      | true        |
|                                  | Only for debugging purposes.|                   |             |
+----------------------------------+-----------------------------+-------------------+-------------+
| **remora.use_barotropic**        | Include 2d barotropic step. | true / false      | true        |
|                                  | Only for debugging purposes.|                   |             |
+----------------------------------+-----------------------------+-------------------+-------------+

Physics Parameters
==================

.. _list-of-parameters-15:

List of Parameters
------------------

+-----------------------------------+----------------------------------------+-------------------+----------------+
| Parameter                         | Definition                             | Acceptable        | Default        |
|                                   |                                        | Values            |                |
+===================================+========================================+===================+================+
| **remora.ggrav**                  | Gravitational field strength           | Real number       | 9.81           |
|                                   | [kg m/s^2]                             |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.eos_type**               | Which equation of state to use.        | Linear or         | Linear         |
|                                   | Nonlinear corresponds to UNESCO        | Nonlinear         |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.R0**                     | Background density [kg/m^3]            | Real number       | 1028           |
|                                   | used in Linear Equation of             |                   |                |
|                                   | State. May be used in setup            |                   |                |
|                                   | of some problems.                      |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.S0**                     | Background salinity                    | Real number       | 35             |
|                                   | (nondimensional) used in               |                   |                |
|                                   | Linear Equation of State               |                   |                |
|                                   | State. May be used in setup            |                   |                |
|                                   | of some problems.                      |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.T0**                     | Background temperature                 | Real number       | 5              |
|                                   | (Celsius) used in                      |                   |                |
|                                   | Linear Equation of State               |                   |                |
|                                   | State. May be used in setup            |                   |                |
|                                   | of some problems.                      |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.Tcoef**                  | Linear EOS parameter                   | Real number       | 1.7e-4         |
|                                   | (1/Celsius)                            |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.Scoef**                  | Linear EOS parameter                   | Real number       | 0.0            |
|                                   | (nondimensional)                       |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.rho0**                   | Mean density (kg/m^3) used             | Real number       | 1025           |
|                                   | when Boussinesq approx is              |                   |                |
|                                   | inferred                               |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.coriolis_type**          | Type of Coriolis forcing.              | ``beta_plane`` /  | ``beta_plane`` |
|                                   | ``beta_plane`` uses a linear           | ``analytic`` /    |                |
|                                   | approximation. ``analytic`` is         | ``netcdf``        |                |
|                                   | calculated from a function in          |                   |                |
|                                   | ``prob.cpp``, and ``netcdf`` is        |                   |                |
|                                   | read from the netcdf grid file         |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.coriolis_f0**            | f-plane constant for                   | Real number       | 0.0            |
|                                   | Coriolis param                         |                   |                |
|                                   | :math:`f = f_0 + \beta y`              |                   |                |
|                                   | when using beta plane                  |                   |                |
|                                   | Coriolis type                          |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.coriolis_beta**          | beta-plane constant for                | Real number       | 0.0            |
|                                   | Coriolis param                         |                   |                |
|                                   | :math:`f = f_0 + \beta y`              |                   |                |
|                                   | when using beta plane                  |                   |                |
|                                   | Coriolis type                          |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.horizontal_mixing_type** | Horizontal mixing type. ``analytic``   | ``analytic`` /    | ``analytic``   |
|                                   | function is specified in               | ``constant``      |                |
|                                   | ``prob.cpp``.                          |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.visc2**                  | Constant horizontal viscosity,         | Real number       | 0.0            |
|                                   | everywhere. Needed when                |                   |                |
|                                   | ``horizontal_mixing_type`` is          |                   |                |
|                                   | ``constant``.                          |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.tnu2_salt**              | Constant horizontal diffusivity,       | Real number       | 0.0            |
|                                   | everywhere for salt. Needed when       |                   |                |
|                                   | ``horizontal_mixing_type`` is          |                   |                |
|                                   | ``constant``.                          |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.tnu2_temp**              | Constant horizontal diffusivity,       | Real number       | 0.0            |
|                                   | everywhere for temperature. Needed     |                   |                |
|                                   | when ``horizontal_mixing_type`` is     |                   |                |
|                                   | ``constant``.                          |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.tnu2_scalar**            | Constant horizontal diffusivity,       | Real number       | 0.0            |
|                                   | everywhere for passive scalar. Needed  |                   |                |
|                                   | when ``horizontal_mixing_type`` is     |                   |                |
|                                   | ``constant``.                          |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.vertical_mixing_type**   | Vertical mixing type. ``analytic``     | ``analytic`` /    | ``analytic``   |
|                                   | function is specified in               | ``GLS``           |                |
|                                   | ``prob.cpp``.                          |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.gls_stability_type**     | Stability function to use for GLS      | ``Canuto_A`` /    | ``Canuto_A``   |
|                                   |                                        | ``Canuto_B`` /    |                |
|                                   |                                        | ``Galperin``      |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.Akv_bak**                | Minimum/initial value of Akv           | Real number       | 5.0e-6         |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.Akt_bak**                | Minimum/initial value of Akt           | Real number       | 1.0e-6         |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.bulk_fluxes**            | Whether to use bulk fluxes             | true / false      | false          |
|                                   | parametrization                        |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+

.. _list-of-parameters-drag:

List of drag-related parameters
--------------------------------

+-----------------------------------+------------------------------------------+-------------------+----------------+
| Parameter                         | Definition                               | Acceptable        | Default        |
|                                   |                                          | Values            |                |
+===================================+==========================================+===================+================+
| **remora.bottom_stress_type**     | Bottom drag type                         | linear            | linear         |
|                                   |                                          | / quadratic       |                |
|                                   |                                          | / logarithmic     |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+
| **remora.rdrag**                  | Linear drag coefficient (used if         | Positive real     | 3.0e-4         |
|                                   | remora.bottom_stress_type = linear)      | number            |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+
| **remora.rdrag2**                 | Quadratic drag coefficient (used if      | Positive real     | 3.0e-3         |
|                                   | remora.bottom_stress_type = quadratic)   | number            |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+
| **remora.Zob**                    | Bottom roughness [m] (used if            | Positive real     | 0.02           |
|                                   | remora.bottom_stress_type = logarithmic  | number            |                |
|                                   | and for GLS                              |                   |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+
| **remora.Zos**                    | Surface roughness [m] (used for GLS)     | Positive real     | 0.02           |
|                                   |                                          | number            |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+
| **remora.Cdb_min**                | Minimum threshold for transfer           | Positive real     | 1.0e-6         |
|                                   | coefficient of momentum                  | number            |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+
| **remora.Cdb_max**                | Maximum threshold for transfer           | Positive real     | 0.5            |
|                                   | coefficient of momentum                  | number            |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+


.. _list-of-parameters-gls:

List of GLS-specific parameters
-------------------------------

+----------------------------------+--------------------------------------+-------------------+----------------+
| Parameter                        | Definition                           | Acceptable        | Default        |
|                                  |                                      | Values            |                |
+==================================+======================================+===================+================+
| **remora.gls_P**                 |                                      | Real number       | 3.0            |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.gls_M**                 |                                      | Real number       | 1.5            |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.gls_N**                 |                                      | Real number       | -1.0           |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.gls_Kmin**              |                                      | Real number       | 7.6e-6         |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.gls_Pmin**              |                                      | Real number       | 1.0e-12        |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.gls_cmu0**              |                                      | Real number       | 0.5477         |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.gls_c1**                |                                      | Real number       | 1.44           |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.gls_c2**                |                                      | Real number       | 1.92           |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.gls_c3m**               |                                      | Real number       | -0.4           |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.gls_c3p**               |                                      | Real number       | 1.0            |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.gls_sigk**              |                                      | Real number       | 1.0            |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.gls_sigp**              |                                      | Real number       | 1.3            |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.Akk_bak**               | Initial/minimum value of Akk         | Real number       | 5.0e-6         |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.Akp_bak**               | Initial/minimum value of Akp         | Real number       | 5.0e-6         |
|                                  |                                      |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+

.. _list-of-parameters surface-forcing:

List of surface forcing parameters
----------------------------------

+----------------------------------+-----------------------------------------+-------------------+----------------+
| Parameter                        | Definition                              | Acceptable        | Default        |
|                                  |                                         | Values            |                |
+==================================+=========================================+===================+================+
| **remora.smflux_type**           | Input format for surface momentum flux, | ``analytic`` or   | ``analytic``   |
|                                  | if using. ``analytic`` specified in     | ``netcdf``        |                |
|                                  | ``prob.cpp``; ``netcdf`` from file.     |                   |                |
+----------------------------------+-----------------------------------------+-------------------+----------------+
| **remora.wind_type**             | Input format for surface wind speed,    | ``analytic`` or   | ``analytic``   |
|                                  | if using. ``analytic`` specified in     | ``netcdf``        |                |
|                                  | ``prob.cpp``; ``netcdf`` from file.     |                   |                |
+----------------------------------+-----------------------------------------+-------------------+----------------+

.. _list-of-parameters-bulk-fluxes:

List of Bulk Fluxes parameters
------------------------------

+----------------------------------+----------------------------------------+-------------------+----------------+
| Parameter                        | Definition                             | Acceptable        | Default        |
|                                  |                                        | Values            |                |
+==================================+========================================+===================+================+
| **remora.air_temperature**       | Air temperature [C]                    | Real number       | 23.567         |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.air_humidity**          | Relative humidity of air               | Real number from  | 0.776          |
|                                  |                                        | 0 to 1            |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.air_pressure**          | Air pressure [hPa]                     | Real number       | 1013.48        |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.blk_ZQ**                | Height [m] of atmospheric humidity     | Real number       | 10.0           |
|                                  | memasurements for bulk fluxes          |                   |                |
|                                  | parametrization                        |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.blk_ZT**                | Height [m] of atmospheric temperature  | Real number       | 10.0           |
|                                  | memasurements for bulk fluxes          |                   |                |
|                                  | parametrization                        |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.blk_ZW**                | Height [m] of atmospheric wind         | Real number       | 10.0           |
|                                  | memasurements for bulk fluxes          |                   |                |
|                                  | parametrization                        |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.cloud**                 | Cloud cover fraction (0=clear sky,     | Real number from  | 0.0            |
|                                  | 1=overcast)                            | 0 to 1            |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.rain**                  | Precipitation rate [kg/m^2/s]          | Real number       | 0.0            |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.eminusp**               | Whether to do E-P prescription for     | true / false      | false          |
|                                  | evaporation/precipiation               |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.eminusp_correct_ssh**   | Whether to adjust sea surface          | true / false      | false          |
|                                  | height for amount of evaporation       |                   |                |
|                                  | and precipitation                      |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+

Numerical Algorithms
====================

.. _list-of-parameters-16:

List of Parameters
------------------

+-----------------------------------------------+-----------------------------+-------------------+-------------+
| Parameter                                     | Definition                  | Acceptable        | Default     |
|                                               |                             | Values            |             |
+===============================================+=============================+===================+=============+
| **remora.tracer_horizontal_advection_scheme** | Scheme for horizontal       | upstream3,        | upstream3   |
|                                               | advection of tracers        | centered4         |             |
+-----------------------------------------------+-----------------------------+-------------------+-------------+
| **remora.uv_horizontal_advection_scheme**     | Scheme for horizontal       | upstream3,        | upstream3   |
|                                               | advection of momenta        | centered2         |             |
+-----------------------------------------------+-----------------------------+-------------------+-------------+

Vertical Stretch parameters
===========================

.. _list-of-parameters-17:

List of Parameters
------------------

+---------------------------------------+-----------------------------+----------------------------------+-------------+
| Parameter                             | Definition                  | Acceptable                       | Default     |
|                                       |                             | Values                           |             |
+=======================================+=============================+==================================+=============+
| **remora.theta_s**                    | Stretching parameter for    | :math:`0 \leq \theta_S \leq 10`  | 3.0         |
|                                       | surface refinement of       |                                  |             |
|                                       | vertical S-grid             |                                  |             |
+---------------------------------------+-----------------------------+----------------------------------+-------------+
| **remora.theta_b**                    | Stretching parameter for    | :math:`0 \leq \theta_B \leq 4`   | 0.0         |
|                                       | bottom refinement of        |                                  |             |
|                                       | vertical S-grid             |                                  |             |
+---------------------------------------+-----------------------------+----------------------------------+-------------+
| **remora.tcline**                     | Surface/bottom layer width  | Positive number                  | 150         |
|                                       | (m) in vertical S-grid      |                                  |             |
+---------------------------------------+-----------------------------+----------------------------------+-------------+

These parameters are used to calculate the vertical S-grid stretch/transform functions detailed in
:ref:`Vertical S-Coordinate<VerticalSCoord>`.

Climatology parameters
======================

.. _climatology-parameters:

List of Parameters
------------------

+---------------------------------------+-----------------------------+--------------+---------------------------+
| Parameter                             | Definition                  | Acceptable   | Default                   |
|                                       |                             | Values       |                           |
+=======================================+=============================+==============+===========================+
| **remora.do_m3_clim_nudg**            | Whether to nudge 3D         | true / false | false                     |
|                                       | momentum variables          |              |                           |
|                                       | (u and v) to climatology    |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.do_m2_clim_nudg**            | Whether to nudge 2D         | true / false | false                     |
|                                       | momentum variables (ubar    |              |                           |
|                                       | and vbar) to climatology    |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.do_temp_clim_nudg**          | Whether to nudge            | true / false | false                     |
|                                       | temperature                 |              |                           |
|                                       | to climatology              |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.do_salt_clim_nudg**          | Whether to nudge            | true / false | false                     |
|                                       | salinity                    |              |                           |
|                                       | to climatology              |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.nc_clim_his_file**           | NetCDF file name for        | string       | must be set if one of the |
|                                       | climatology data            |              | ``remora.do_*_clim_nudg`` |
|                                       |                             |              | flags is true             |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.nc_clim_coeff_file**         | NetCDF file name for        | string       | must be set if one of the |
|                                       | climatology nudging         |              | ``remora.do_*_clim_nudg`` |
|                                       | coefficients                |              | flags is true             |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.clim_ubar_time_varname**     | name of time variable for   | string       | ``ocean_time``            |
|                                       | ubar climatology            |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.clim_vbar_time_varname**     | name of time variable for   | string       | ``ocean_time``            |
|                                       | vbar climatology            |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.clim_u_time_varname**        | name of time variable for   | string       | ``ocean_time``            |
|                                       | u climatology               |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.clim_v_time_varname**        | name of time variable for   | string       | ``ocean_time``            |
|                                       | v climatology               |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.clim_temp_time_varname**     | name of time variable for   | string       | ``ocean_time``            |
|                                       | temperature climatology     |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.clim_salt_time_varname**     | name of time variable for   | string       | ``ocean_time``            |
|                                       | salinity climatology        |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+

..
  include:: InputsPhysics.rst


Runtime Error Checking
======================

Through `AMReX functionality <https://amrex-codes.github.io/amrex/docs_html/Debugging.html>`_,
REMORA supports options to raise errors when NaNs, division by zero, and overflow errors
are detected. These checks are activated at runtime using the input parameters below.

.. note::

   When running on Macs using the Apple-Clang compilers with optimization
   (``DEBUG = FALSE`` in the ``GNUmakefile``), these checks may lead to false positives
   due to optimizations performed by the compiler and the flags should be turned off.
   It is still possible to run with these error checks with Apple-Clang debug builds.

List of Parameters
------------------

+-----------------------------+---------------------------+-------------------+------------+
| Parameter                   | Definition                | Acceptable Values | Default    |
+=============================+===========================+===================+============+
| **amrex.fpe_trap_invalid**  | Raise errors for NaNs     |  0 / 1            | 0          |
+-----------------------------+---------------------------+-------------------+------------+
| **amrex.fpe_trap_zero**     | Raise errors for divide   |  0 / 1            | 0          |
|                             | by zero                   |                   |            |
+-----------------------------+---------------------------+-------------------+------------+
| **amrex.fpe_trap_overflow** | Raise errors for overflow |  0 / 1            | 0          |
+-----------------------------+---------------------------+-------------------+------------+
