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
|                          |                 |                   |             |
|                          | corner of the   |                   |             |
|                          | domain          |                   |             |
+--------------------------+-----------------+-------------------+-------------+
| **remora.prob_hi**       | physical        | [Real Real 0]     | must be set |
|                          | location of     |                   |             |
|                          |                 |                   |             |
|                          | high corner of  |                   |             |
|                          | the domain      |                   |             |
+--------------------------+-----------------+-------------------+-------------+
| **remora.is_periodic**   | is the domain   | 0 if false, 1     | 0 0 0       |
|                          | periodic in     | if true.          |             |
|                          |                 |                   |             |
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

Grid, initial, and time-dependent boundary data can be specified using NetCDF files, as in ROMS. REMORA expects files in the same format as ROMS, in NetCDF classic format (32- or 64-bit). NetCDF4-classic does not work. If the file format is incorrect, REMORA will exit with the error:

.. code:: shell

   NetCDF: Attempt to use feature that was not turned on when netCDF was built.

Other versions of NetCDF files can be converted to 64-bit NetCDF classic by running the command:

.. code:: shell

    ncks -5 old_file.nc converted_file.nc

The utility ``ncks`` is part of the NCO suite.

Currently, if initial or grid files are specified, they both must be. Boundary condition options with NetCDF boundary data are equivalent to ROMS clamped, Chapman-Flather, and Orlanski + Nudging boundary conditions. Options and examples can be found in the section on :ref:`Domain Boundary Conditions <sec:domainBCs>`.

List of Parameters
------------------

+-----------------------------------+-----------------------------------+-----------------+---------------------------------+
| Parameter                         | Definition                        | Acceptable      | Default                         |
|                                   |                                   |                 |                                 |
|                                   |                                   | Values          |                                 |
+===================================+===================================+=================+=================================+
| **remora.ic_type**                | read initial and grid             | ``analytic``    |                                 |
|                                   |                                   |                 |                                 |
|                                   | data from files                   | or ``netcdf``   | ``analytic``                    |
+-----------------------------------+-----------------------------------+-----------------+---------------------------------+
| **remora.nc_init_file_0**         | initial data NetCDF               | string          | must be set                     |
|                                   |                                   |                 |                                 |
|                                   |                                   |                 | if ``remora.ic_type``           |
|                                   |                                   |                 |                                 |
|                                   | file name                         |                 | is true                         |
+-----------------------------------+-----------------------------------+-----------------+---------------------------------+
| **remora.nc_grid_file_0**         | grid data NetCDF                  | string          | must be set                     |
|                                   |                                   |                 |                                 |
|                                   |                                   |                 | if ``remora.ic_type``           |
|                                   |                                   |                 |                                 |
|                                   | file name                         |                 | is true                         |
+-----------------------------------+-----------------------------------+-----------------+---------------------------------+
| **remora.nc_bdry_file**           | boundary data NetCDF              | string or list  | must be set if                  |
|                                   |                                   |                 |                                 |
|                                   |                                   |                 | any boundary                    |
|                                   | file names(s)                     | of strings      |                                 |
|                                   |                                   |                 | condition requires it:          |
|                                   |                                   |                 |                                 |
|                                   |                                   |                 | ``orlanski_rad_nudg``,          |
|                                   |                                   |                 |                                 |
|                                   |                                   |                 | ``clamped``, ``chapman``,       |
|                                   |                                   |                 |                                 |
|                                   |                                   |                 | or ``flather``.                 |
+-----------------------------------+-----------------------------------+-----------------+---------------------------------+
| **remora.nc_frc_file**            | forcing data NetCDF               | string          | must be set                     |
|                                   |                                   |                 |                                 |
|                                   |                                   |                 | if ``remora.wind_type``         |
|                                   | file name                         |                 |                                 |
|                                   |                                   |                 | or ``remora.smflux_type``       |
|                                   |                                   |                 |                                 |
|                                   |                                   |                 | equal ``netcdf``                |
+-----------------------------------+-----------------------------------+-----------------+---------------------------------+
| **remora.bdy_time_varname**       | default name of time              | string          | ``ocean_time``                  |
|                                   |                                   |                 |                                 |
|                                   | variable in                       |                 |                                 |
|                                   |                                   |                 |                                 |
|                                   | boundary file                     |                 |                                 |
+-----------------------------------+-----------------------------------+-----------------+---------------------------------+
| **remora.bdy_{var}_time_varname** | name of time variable             | string          | None                            |
|                                   |                                   |                 |                                 |
|                                   | for variable {var} (one           |                 |                                 |
|                                   |                                   |                 |                                 |
|                                   | of ``temp``, ``salt``, ``u``,     |                 |                                 |
|                                   |                                   |                 |                                 |
|                                   | ``v``, ``ubar``, ``vbar``,        |                 |                                 |
|                                   |                                   |                 |                                 |
|                                   | ``zeta``)                         |                 |                                 |
+-----------------------------------+-----------------------------------+-----------------+---------------------------------+
| **remora.frc_time_varname**       | name of time variable             | string          | ``wind_time`` for wind,         |
|                                   |                                   |                 |                                 |
|                                   | in forcing file                   |                 | ``sms_time`` for surface        |
|                                   |                                   |                 |                                 |
|                                   |                                   |                 | momentum stress                 |
+-----------------------------------+-----------------------------------+-----------------+---------------------------------+

Notes
-----

-  ``nc_bdry_file`` must either be a string or a space-separated list of strings of boundary data files. They must be in time series order.

-  The time variables in the boundary files may be different for each boundary variable. Any that are not individually specified with
   ``remora.bdy_{var}_time_varname`` will default to the variable name given by ``bdy_time_varname``.

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
|                           |                 | {x,y,z}                   |                       |
|                           | direction at    |                           |                       |
|                           | the coarsest    |                           |                       |
|                           |                 |                           |                       |
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

+----------------------------------+------------------+--------------------+-------------------+
| Parameter                        | Definition       | Acceptable         | Default           |
|                                  |                  | Values             |                   |
+==================================+==================+====================+===================+
| **amr.max_level**                | number of        | Integer >= 0       | must be set       |
|                                  | levels of        |                    |                   |
|                                  |                  |                    |                   |
|                                  | refinement       |                    |                   |
|                                  | above the        |                    |                   |
|                                  |                  |                    |                   |
|                                  | coarsest level   |                    |                   |
+----------------------------------+------------------+--------------------+-------------------+
| **amr.ref_ratio_vect**           | ratio of coarse  | 3 integers         | 2 for all         |
|                                  |                  |                    |                   |
|                                  | to fine grid     |                    | directions        |
|                                  |                  | (one per dir)      |                   |
|                                  |                  |                    |                   |
|                                  | spacing between  | per refinement     |                   |
|                                  |                  |                    |                   |
|                                  | subsequent       | level. Refinement  |                   |
|                                  | levels           |                    |                   |
|                                  |                  | in z must be 1     |                   |
|                                  |                  |                    |                   |
+----------------------------------+------------------+--------------------+-------------------+
| **amr.regrid_int**               | how often to     | Integer            | -1 (don't regrid) |
|                                  | regrid           |                    |                   |
+----------------------------------+------------------+--------------------+-------------------+
| **amr.regrid_on_restart**        | should we        | 0 or 1             | 0                 |
|                                  | regrid           |                    |                   |
|                                  |                  |                    |                   |
|                                  | immediately      |                    |                   |
|                                  |                  |                    |                   |
|                                  | after            |                    |                   |
|                                  | restarting       |                    |                   |
+----------------------------------+------------------+--------------------+-------------------+
| **amr.regrid_file**              | name of file     | text               | no file           |
|                                  |                  |                    |                   |
|                                  | from which to    |                    |                   |
|                                  |                  |                    |                   |
|                                  | read the grids   |                    |                   |
+----------------------------------+------------------+--------------------+-------------------+
| **amr.grid_eff**                 | grid             | Real > 0, < 1      | 0.7               |
|                                  | efficiency at    |                    |                   |
|                                  |                  |                    |                   |
|                                  | coarse level     |                    |                   |
|                                  |                  |                    |                   |
|                                  | at which grids   |                    |                   |
|                                  |                  |                    |                   |
|                                  | are created      |                    |                   |
+----------------------------------+------------------+--------------------+-------------------+
| **amr.n_error_buf**              | radius of        | Integer >= 0       | 1                 |
|                                  | additional       |                    |                   |
|                                  |                  |                    |                   |
|                                  | tagging around   |                    |                   |
|                                  |                  |                    |                   |
|                                  | already tagged   |                    |                   |
|                                  | cells            |                    |                   |
+----------------------------------+------------------+--------------------+-------------------+
| **amr.n_error_buf_{x,y,z}**      | radius of        | Integer >= 0;      | 1                 |
|                                  | additional       |                    |                   |
|                                  |                  | Can specify up     |                   |
|                                  | tagging around   |                    |                   |
|                                  |                  |                    |                   |
|                                  | already tagged   | to one per         |                   |
|                                  |                  | ref. level         |                   |
|                                  | cells in x, y,   |                    |                   |
|                                  | or z             |                    |                   |
+----------------------------------+------------------+--------------------+-------------------+
| **amr.blocking_factor**          | grid size must   | Integer > 0        | 2                 |
|                                  |                  |                    |                   |
|                                  | be a multiple    |                    |                   |
|                                  |                  |                    |                   |
|                                  | of this          |                    |                   |
+----------------------------------+------------------+--------------------+-------------------+
| **amr.blocking_factor_{x,y}**    | grid size in     | Integer > 0;       | 1                 |
|                                  |                  |                    |                   |
|                                  | {x,y} must be    | Can specify up to  |                   |
|                                  |                  |                    |                   |
|                                  | multiple of this | one per ref. level |                   |
+----------------------------------+------------------+--------------------+-------------------+
| **amr.blocking_factor_z**        | grid size in     | Integer > 0;       | 1                 |
|                                  |                  |                    |                   |
|                                  | z must be        | In AMR problems,   |                   |
|                                  |                  |                    |                   |
|                                  | multiple of this | recommend it equal |                   |
|                                  |                  |                    |                   |
|                                  |                  | ``n_cell`` in z    |                   |
+----------------------------------+------------------+--------------------+-------------------+
| **amr.refine_grid_layout**       | refine grids     | 0 if false, 1      | 1                 |
|                                  | more if          | if true            |                   |
|                                  |                  |                    |                   |
|                                  | # of processors  |                    |                   |
|                                  |                  |                    |                   |
|                                  | :math:`>` # of   |                    |                   |
|                                  | grids            |                    |                   |
+----------------------------------+------------------+--------------------+-------------------+
| **amr.do_substep**               | whether to       | 0 if false, 1      | 0                 |
|                                  | sub-step finer   | if true            |                   |
|                                  |                  |                    |                   |
|                                  | levels in time   |                    |                   |
|                                  |                  | NOTE: true         |                   |
|                                  |                  | will               |                   |
|                                  |                  |                    |                   |
|                                  |                  | trigger Assert     |                   |
|                                  |                  |                    |                   |
|                                  |                  | failure            |                   |
+----------------------------------+------------------+--------------------+-------------------+

.. _notes-2:

Notes
-----

-  if **amr.max_level** = 0 then you do not need to set
   **amr.ref_ratio_vect** or **amr.regrid_int**.

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

-  | **amr.ref_ratio** = 2 6
   | would set factor 2 refinement between levels 0 and 1, and factor 3
     refinement between levels 1 and 2 (6 between levels 0 and 2). Note
     that you must have at least **amr.max_level** values of
     **amr.ref_ratio** (Additional values
     may appear in that line and they will be ignored).

-  | **amr.ref_ratio_vect** = 2 4 1
   | would set factor {2 in x-dir, 4 in y-dir, 1 in z-dir} refinement between
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

-  **remora.max_step** = 1000

-  **remora.stop_time** = 1.0

will end the calculation when either the simulation time reaches 1.0 or
the number of level-0 steps taken equals 1000, whichever comes first.

Time Step
=========

.. _list-of-parameters-6:

List of Parameters for Single-Rate
----------------------------------

+-------------------------------+-------------------------+----------------+----------------------------+
| Parameter                     | Definition              | Acceptable     | Default                    |
|                               |                         |                |                            |
|                               |                         | Values         |                            |
+===============================+=========================+================+============================+
| **remora.cfl**                | CFL number for          | Real > 0       | 0.8                        |
|                               | hydro                   |                |                            |
|                               |                         | and <= 1       |                            |
+-------------------------------+-------------------------+----------------+----------------------------+
| **remora.fixed_dt**           | set level 0 dt          | Real > 0       | unused if                  |
|                               |                         |                |                            |
|                               | as this value           |                | not set                    |
|                               |                         |                |                            |
|                               | value regardless of     |                |                            |
|                               | cfl                     |                |                            |
|                               |                         |                |                            |
|                               | or other                |                |                            |
|                               | settings                |                |                            |
+-------------------------------+-------------------------+----------------+----------------------------+
| **remora.fixed_fast_dt**      | set fast dt as this     | real > 0       | inferred from **fixed_dt** |
|                               | value                   |                |                            |
|                               |                         |                | and **fixed_ndfast_ratio** |
|                               |                         |                |                            |
|                               |                         |                | if not set                 |
+-------------------------------+-------------------------+----------------+----------------------------+
| **remora.fixed_ndfast_ratio** | set fast dt as          | int            | inferred from **fixed_dt** |
|                               |                         |                |                            |
|                               |                         |                | and **fixed_fast_dt**      |
|                               |                         |                |                            |
|                               | slow dt /               |                |                            |
|                               | this ratio              |                | if not set                 |
+-------------------------------+-------------------------+----------------+----------------------------+
| **remora.change_max**         | factor by which         | Real >= 1      | 1.1                        |
|                               | dt can grow             |                |                            |
|                               |                         |                |                            |
|                               | in subsequent           |                |                            |
|                               | steps                   |                |                            |
+-------------------------------+-------------------------+----------------+----------------------------+

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
|                            | REMORA           | - 1: integrated/max quantities   |                |
|                            |                  | - 2: print boxes                 |                |
|                            | functions        |                                  |                |
+----------------------------+------------------+----------------------------------+----------------+
| **remora.sum_interval**    | how often (in    |                                  |                |
|                            | level-0          |                                  |                |
|                            |                  |                                  |                |
|                            | time steps) to   |                                  |                |
|                            | to compute       | Integer                          | -1             |
|                            |                  |                                  |                |
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
|                                  |                             |                   |             |
|                                  | debugging purposes.         |                   |             |
+----------------------------------+-----------------------------+-------------------+-------------+
| **remora.use_uv3dmix**           | Include harmonic viscosity. | true / false      | true        |
|                                  |                             |                   |             |
|                                  | Only for debugging purposes.|                   |             |
+----------------------------------+-----------------------------+-------------------+-------------+
| **remora.use_barotropic**        | Include 2d barotropic step. | true / false      | true        |
|                                  |                             |                   |             |
|                                  | Only for debugging purposes.|                   |             |
+----------------------------------+-----------------------------+-------------------+-------------+

Physics Parameters
==================

.. _list-of-parameters-15:

List of Parameters
------------------

+-----------------------------------+----------------------------------------+-------------------+----------------+
| Parameter                         | Definition                             | Acceptable        | Default        |
|                                   |                                        |                   |                |
|                                   |                                        | Values            |                |
+===================================+========================================+===================+================+
| **remora.ggrav**                  | Gravitational field strength           | Real number       | 9.81           |
|                                   |                                        |                   |                |
|                                   | [kg m/s^2]                             |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.eos_type**               | Which equation of state to use.        | Linear or         | Linear         |
|                                   |                                        |                   |                |
|                                   | Nonlinear is UNESCO                    | Nonlinear         |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.R0**                     | Background density [kg/m^3]            | Real number       | 1028           |
|                                   |                                        |                   |                |
|                                   | used in Linear Equation of             |                   |                |
|                                   |                                        |                   |                |
|                                   | State. May be used in setup            |                   |                |
|                                   |                                        |                   |                |
|                                   | of some problems.                      |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.S0**                     | Background salinity                    | Real number       | 35             |
|                                   |                                        |                   |                |
|                                   | (nondimensional) used in               |                   |                |
|                                   |                                        |                   |                |
|                                   | Linear Equation of State               |                   |                |
|                                   |                                        |                   |                |
|                                   | State. May be used in setup            |                   |                |
|                                   |                                        |                   |                |
|                                   | of some problems.                      |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.T0**                     | Background temperature                 | Real number       | 5              |
|                                   |                                        |                   |                |
|                                   | (Celsius) used in                      |                   |                |
|                                   |                                        |                   |                |
|                                   | Linear Equation of State               |                   |                |
|                                   |                                        |                   |                |
|                                   | State. May be used in setup            |                   |                |
|                                   |                                        |                   |                |
|                                   | of some problems.                      |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.Tcoef**                  | Linear EOS parameter                   | Real number       | 1.7e-4         |
|                                   |                                        |                   |                |
|                                   | (1/Celsius)                            |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.Scoef**                  | Linear EOS parameter                   | Real number       | 0.0            |
|                                   |                                        |                   |                |
|                                   | (nondimensional)                       |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.rho0**                   | Mean density (kg/m^3) used             | Real number       | 1025           |
|                                   |                                        |                   |                |
|                                   | when Boussinesq approx is              |                   |                |
|                                   |                                        |                   |                |
|                                   | inferred                               |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.coriolis_type**          | Type of Coriolis forcing.              | ``beta_plane`` /  | ``beta_plane`` |
|                                   |                                        |                   |                |
|                                   | ``beta_plane`` uses a linear           | ``analytic`` /    |                |
|                                   |                                        |                   |                |
|                                   | approximation. ``analytic`` is         | ``netcdf``        |                |
|                                   |                                        |                   |                |
|                                   | calculated from a function in          |                   |                |
|                                   |                                        |                   |                |
|                                   | ``prob.cpp``, and ``netcdf`` is        |                   |                |
|                                   |                                        |                   |                |
|                                   | read from the netcdf grid file         |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.coriolis_f0**            | f-plane constant for                   | Real number       | 0.0            |
|                                   |                                        |                   |                |
|                                   | Coriolis param                         |                   |                |
|                                   |                                        |                   |                |
|                                   | :math:`f = f_0 + \beta y`              |                   |                |
|                                   |                                        |                   |                |
|                                   | when using beta plane                  |                   |                |
|                                   |                                        |                   |                |
|                                   | Coriolis type                          |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.coriolis_beta**          | beta-plane constant for                | Real number       | 0.0            |
|                                   |                                        |                   |                |
|                                   | Coriolis param                         |                   |                |
|                                   |                                        |                   |                |
|                                   | :math:`f = f_0 + \beta y`              |                   |                |
|                                   |                                        |                   |                |
|                                   | when using beta plane                  |                   |                |
|                                   |                                        |                   |                |
|                                   | Coriolis type                          |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.horizontal_mixing_type** | Horizontal mixing type.                | ``analytic`` /    | ``analytic``   |
|                                   |                                        |                   |                |
|                                   | ``analytic`` means function is         | ``constant``      |                |
|                                   |                                        |                   |                |
|                                   | specified in ``prob.cpp``.             | / ``scaled_to_grid`` |             |
|                                   |                                        |                   |                |
|                                   | ``scaled_to_grid`` scales harmonic     |                   |                |
|                                   | viscosity/diffusivity by the grid      |                   |                |
|                                   | cell area. Equivalent to ``DIFF_GRID`` |                   |                |
|                                   | and ``VISC_GRID`` in ROMS.             |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.scaled_to_grid_amr_scaling** | AMR scaling behavior for                | ``none`` /        | ``none``       |
|                                   | ``scaled_to_grid`` and ``constant``     | ``linear``        |                |
|                                   | coefficients on refined levels.         |                   |                |
|                                   | ``linear`` decreases coefficients in    |                   |                |
|                                   | proportion to the horizontal refinement |                   |                |
|                                   | ratio.                                  |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.visc2**                  | Constant horizontal viscosity,         | Real number       | 0.0            |
|                                   |                                        |                   |                |
|                                   | everywhere. Needed when                |                   |                |
|                                   |                                        |                   |                |
|                                   | ``horizontal_mixing_type`` is          |                   |                |
|                                   |                                        |                   |                |
|                                   | ``constant`` or ``scaled_to_grid``     |                   |                |
|                                   | (in this case, it is the maximum       |                   |                |
|                                   | viscosity over the domain).            |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.tnu2_salt**              | Constant horizontal diffusivity,       | Real number       | 0.0            |
|                                   |                                        |                   |                |
|                                   | everywhere for salt.                   |                   |                |
|                                   |                                        |                   |                |
|                                   | Needed when                            |                   |                |
|                                   |                                        |                   |                |
|                                   | ``horizontal_mixing_type`` is          |                   |                |
|                                   |                                        |                   |                |
|                                   | ``constant`` or ``scaled_to_grid``     |                   |                |
|                                   | (in this case, it is the maximum       |                   |                |
|                                   | salt diffusivity over the domain).     |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.tnu2_temp**              | Constant horizontal diffusivity,       | Real number       | 0.0            |
|                                   |                                        |                   |                |
|                                   | everywhere for temperature.            |                   |                |
|                                   |                                        |                   |                |
|                                   | Needed when                            |                   |                |
|                                   |                                        |                   |                |
|                                   | ``horizontal_mixing_type``             |                   |                |
|                                   |                                        |                   |                |
|                                   | is ``constant`` or ``scaled_to_grid``  |                   |                |
|                                   | (in this case, it is the maximum       |                   |                |
|                                   | temperature diffusivity over the       |                   |                |
|                                   | domain).                               |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.tnu2_scalar**            | Constant horizontal diffusivity,       | Real number       | 0.0            |
|                                   |                                        |                   |                |
|                                   | everywhere for passive scalar.         |                   |                |
|                                   |                                        |                   |                |
|                                   | Needed when                            |                   |                |
|                                   |                                        |                   |                |
|                                   | ``horizontal_mixing_type``             |                   |                |
|                                   |                                        |                   |                |
|                                   | is ``constant`` or ``scaled_to_grid``  |                   |                |
|                                   | (in this case, it is the maximum       |                   |                |
|                                   | scalar diffusivity over the domain).   |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.vertical_mixing_type**   | Vertical mixing type. ``analytic``     | ``analytic`` /    | ``analytic``   |
|                                   |                                        |                   |                |
|                                   | function is specified in               | ``GLS``           |                |
|                                   |                                        |                   |                |
|                                   | ``prob.cpp``.                          |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.gls_stability_type**     | Stability function to use for GLS      | ``Canuto_A`` /    | ``Canuto_A``   |
|                                   |                                        |                   |                |
|                                   |                                        | ``Canuto_B`` /    |                |
|                                   |                                        |                   |                |
|                                   |                                        | ``Galperin``      |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.Akv_bak**                | Minimum/initial value of Akv           | Real number       | 5.0e-6         |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.Akt_bak**                | Minimum/initial value of Akt           | Real number       | 1.0e-6         |
+-----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.bulk_fluxes**            | Whether to use bulk fluxes             | true / false      | false          |
|                                   |                                        |                   |                |
|                                   | parametrization                        |                   |                |
+-----------------------------------+----------------------------------------+-------------------+----------------+

Scaled-to-grid horizontal mixing
--------------------------------

If ``remora.horizontal_mixing_type = "scaled_to_grid"``, REMORA follows the ROMS-style approach of
scaling horizontal harmonic mixing coefficients by the grid cell area.

.. math::

   G(i,j) = \sqrt{\frac{1}{pm(i,j)\,pn(i,j)}}

.. math::

   \nu(i,j) = \nu_0 \frac{G(i,j)}{\max(G)}, \qquad
   \kappa_n(i,j) = \kappa_{0,n} \frac{G(i,j)}{\max(G)}

where :math:`\\nu_0` is ``remora.visc2`` and :math:`\\kappa_{0,n}` are the tracer diffusivities
(``remora.tnu2_temp``, ``remora.tnu2_salt``, ``remora.tnu2_scalar``). This ensures the *maximum*
coefficient over the normalization region equals the user-specified value, while varying spatially
with grid size. Note that if the largest cell area occurs over land, then the maximum over *wet*
cells (and thus what you see after applying ``mask_rho`` in post-processing) may be smaller than
the user-specified value.

Implementation details
^^^^^^^^^^^^^^^^^^^^^^

- The normalization ``max(G)`` is computed as a global maximum over the level-0 grid (rho points at
  the surface, i.e. ``k=0``) and does not use land/sea masks. Equivalently, it uses the maximum grid-cell
  area :math:`A(i,j) = 1/(pm\,pn)` via :math:`\\max(G)=\\sqrt{\\max(A)}`.

- AMR refinement scaling: if ``remora.scaled_to_grid_amr_scaling = "linear"``, then on AMR level
  :math:`\\ell` the coefficients are additionally scaled by the cumulative horizontal refinement ratio,
  :math:`1/\\prod_{m<\\ell}\\sqrt{r_x(m)\,r_y(m)}`. For example, with a refinement ratio of ``5 5 1``,
  level 1 coefficients are reduced by a factor of 5 relative to level 0.

- Ghost cells for the coefficient fields are filled using the same periodic/foextrap boundary fill
  used elsewhere in REMORA. This is done so stencil-based operations (e.g., the psi-point averaging for
  ``visc2_p``) have valid neighbor values near domain boundaries and coarse/fine interfaces. Output
  (plotfiles / NetCDF) writes only the valid region.

When this option is active, the spatially varying coefficients are automatically included in plot output
as 2D (vertically homogeneous) fields:

- ``visc2`` (horizontal viscosity at rho points)
- ``diff2_temp``, ``diff2_salt``, ``diff2_tracer`` (horizontal diffusivities at rho points)

.. _list-of-parameters-drag:

List of drag-related parameters
--------------------------------

+-----------------------------------+------------------------------------------+-------------------+----------------+
| Parameter                         | Definition                               | Acceptable        | Default        |
|                                   |                                          | Values            |                |
+===================================+==========================================+===================+================+
| **remora.bottom_stress_type**     | Bottom drag type                         | ``linear``        | ``linear``     |
|                                   |                                          |                   |                |
|                                   |                                          | / ``quadratic``   |                |
|                                   |                                          |                   |                |
|                                   |                                          | / ``logarithmic`` |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+
| **remora.rdrag**                  | Linear drag coefficient (used if         | Positive real     | 3.0e-4         |
|                                   |                                          |                   |                |
|                                   | ``remora.bottom_stress_type`` =          | number            |                |
|                                   |                                          |                   |                |
|                                   | ``linear``)                              |                   |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+
| **remora.rdrag2**                 | Quadratic drag coefficient (used if      | Positive real     | 3.0e-3         |
|                                   |                                          |                   |                |
|                                   | ``remora.bottom_stress_type`` =          | number            |                |
|                                   |                                          |                   |                |
|                                   | ``quadratic``)                           |                   |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+
| **remora.Zob**                    | Bottom roughness [m] (used if            | Positive real     | 0.02           |
|                                   |                                          |                   |                |
|                                   | ``remora.bottom_stress_type`` =          | number            |                |
|                                   |                                          |                   |                |
|                                   | ``logarithmic`` and for GLS)             |                   |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+
| **remora.Zos**                    | Surface roughness [m]                    | Positive real     | 0.02           |
|                                   |                                          |                   |                |
|                                   | (used for GLS)                           | number            |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+
| **remora.Cdb_min**                | Minimum threshold for transfer           | Positive real     | 1.0e-6         |
|                                   |                                          |                   |                |
|                                   | coefficient of momentum                  | number            |                |
+-----------------------------------+------------------------------------------+-------------------+----------------+
| **remora.Cdb_max**                | Maximum threshold for transfer           | Positive real     | 0.5            |
|                                   |                                          |                   |                |
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
| **remora.smflux_type**           | Input format for surface momentum       | ``analytic`` or   | ``analytic``   |
|                                  |                                         |                   |                |
|                                  | flux if using. ``analytic`` specified   | ``netcdf``        |                |
|                                  |                                         |                   |                |
|                                  | in ``prob.cpp``; ``netcdf`` from file.  |                   |                |
+----------------------------------+-----------------------------------------+-------------------+----------------+
| **remora.wind_type**             | Input format for surface wind speed,    | ``analytic`` or   | ``analytic``   |
|                                  |                                         |                   |                |
|                                  | if using. ``analytic`` specified in     | ``netcdf``        |                |
|                                  |                                         |                   |                |
|                                  | ``prob.cpp``; ``netcdf`` from file.     |                   |                |
+----------------------------------+-----------------------------------------+-------------------+----------------+

.. _list-of-parameters-bulk-fluxes:

List of Bulk Fluxes parameters
------------------------------

+----------------------------------+----------------------------------------+-------------------+----------------+
| Parameter                        | Definition                             | Acceptable        | Default        |
|                                  |                                        |                   |                |
|                                  |                                        | Values            |                |
+==================================+========================================+===================+================+
| **remora.air_temperature**       | Air temperature [C] (used as           | Real number       | 23.567         |
|                                  |                                        |                   |                |
|                                  | uniform value if                       |                   |                |
|                                  |                                        |                   |                |
|                                  | **Tair_from_netcdf** is false)         |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.air_humidity**          | Relative humidity of air               | Real number       | 0.776          |
|                                  |                                        |                   |                |
|                                  | (used as uniform value if              | from 0 to 1       |                |
|                                  |                                        |                   |                |
|                                  | **qair_from_netcdf** is false)         |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.air_pressure**          | Air pressure [hPa] (used as            | Real number       | 1013.48        |
|                                  |                                        |                   |                |
|                                  | uniform value if                       |                   |                |
|                                  |                                        |                   |                |
|                                  | **Pair_from_netcdf** is false)         |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.surface_radiation_flux**| Shortwave radiation flux [W/m^2]       | Real number       | 0.0            |
|                                  |                                        |                   |                |
|                                  | (used as uniform value if              |                   |                |
|                                  |                                        |                   |                |
|                                  | **srflx_from_netcdf** is false)        |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.Tair_from_netcdf**      | Load air temperature from NetCDF       | true / false      | false          |
|                                  |                                        |                   |                |
|                                  | file (spatially varying)               |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.qair_from_netcdf**      | Load air humidity from NetCDF          | true / false      | false          |
|                                  |                                        |                   |                |
|                                  | file (spatially varying)               |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.qair_is_percent**       | Convert qair from percentage           | true / false      | false          |
|                                  |                                        |                   |                |
|                                  | (0-100) to fraction (0-1).             |                   |                |
|                                  |                                        |                   |                |
|                                  | Only used if                           |                   |                |
|                                  |                                        |                   |                |
|                                  | **qair_from_netcdf** is true           |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.Pair_from_netcdf**      | Load air pressure from NetCDF          | true / false      | false          |
|                                  |                                        |                   |                |
|                                  | file (spatially varying)               |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.srflx_from_netcdf**     | Load shortwave radiation flux          | true / false      | false          |
|                                  |                                        |                   |                |
|                                  | from NetCDF file                       |                   |                |
|                                  |                                        |                   |                |
|                                  | (spatially varying)                    |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.longwave_down**         | Use file-provided downward longwave    | true / false      | false          |
|                                  | radiation to compute net longwave      |                   |                |
|                                  | (``Lnet = Ldown - sigma*epsilon*T^4``) |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.longwave_down_from_netcdf** | Load longwave field from NetCDF    | true / false      | false          |
|                                  | file (spatially varying)               |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.longwave_netcdf_is_net**| Interpret the NetCDF longwave field    | true / false      | false          |
|                                  | as net longwave (use as-is). If false, |                   |                |
|                                  | interpret as downward longwave and     |                   |                |
|                                  | compute net in the bulk flux routine   |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.longwave_netcdf_varname** | Name of the NetCDF longwave variable | String            | ``lwrad``      |
|                                  | in ``remora.nc_frc_file``              |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.blk_ZQ**                | Height [m] of atmospheric              | Real number       | 10.0           |
|                                  |                                        |                   |                |
|                                  | humidity memasurements for             |                   |                |
|                                  |                                        |                   |                |
|                                  | bulk fluxes parametrization            |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.blk_ZT**                | Height [m] of atmospheric              | Real number       | 10.0           |
|                                  |                                        |                   |                |
|                                  | temperature memasurements              |                   |                |
|                                  |                                        |                   |                |
|                                  | bulk fluxes parametrization            |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.blk_ZW**                | Height [m] of atmospheric wind         | Real number       | 10.0           |
|                                  |                                        |                   |                |
|                                  | memasurements for bulk fluxes          |                   |                |
|                                  |                                        |                   |                |
|                                  | parametrization                        |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.cloud**                 | Cloud cover fraction (0=clear sky,     | Real number       | 0.0            |
|                                  |                                        |                   |                |
|                                  | 1=overcast) (used as uniform           | from 0 to 1       |                |
|                                  |                                        |                   |                |
|                                  | value if **cloud_from_netcdf**         |                   |                |
|                                  |                                        |                   |                |
|                                  | is false)                              |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.rain**                  | Precipitation rate [kg/m^2/s]          | Real number       | 0.0            |
|                                  |                                        |                   |                |
|                                  | (used as uniform value if              |                   |                |
|                                  |                                        |                   |                |
|                                  | **rain_from_netcdf** is false)         |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.cloud_from_netcdf**     | Load cloud cover from NetCDF           | true / false      | false          |
|                                  |                                        |                   |                |
|                                  | file (spatially varying)               |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.rain_from_netcdf**      | Load precipitation rate from           | true / false      | false          |
|                                  |                                        |                   |                |
|                                  | NetCDF file (spatially varying)        |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.EminusP_from_netcdf**   | Load evaporation minus                 | true / false      | false          |
|                                  |                                        |                   |                |
|                                  | precipitation from NetCDF file         |                   |                |
|                                  |                                        |                   |                |
|                                  | (spatially varying)                    |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.eminusp**               | Whether to do E-P prescription for     | true / false      | false          |
|                                  |                                        |                   |                |
|                                  | evaporation/precipiation               |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+
| **remora.eminusp_correct_ssh**   | Whether to adjust sea surface          | true / false      | false          |
|                                  |                                        |                   |                |
|                                  | height for amount of evaporation       |                   |                |
|                                  |                                        |                   |                |
|                                  | and precipitation                      |                   |                |
+----------------------------------+----------------------------------------+-------------------+----------------+

.. note::

   When loading atmospheric forcing variables from NetCDF files (by setting
   **Tair_from_netcdf**, **qair_from_netcdf**, **Pair_from_netcdf**,
   **srflx_from_netcdf**, **longwave_down_from_netcdf**,
   **rain_from_netcdf**, **cloud_from_netcdf**, or **EminusP_from_netcdf**
   to true), these variables are read from the file
   specified by **remora.nc_frc_file** (see :ref:`list-of-parameters surface-forcing`).
   The NetCDF file must contain variables named ``Tair``, ``qair``, ``Pair``,
   ``swrad``, ``rain``, ``cloud``, and ``EminusP`` respectively, with the same
   spatial dimensions as the model grid.
   Time interpolation is performed automatically based on the simulation time.

   For longwave forcing, the variable name is controlled by
   **remora.longwave_netcdf_varname** (default ``lwrad``).
   If **remora.longwave_netcdf_is_net** is true, that variable is treated as
   net longwave radiation and used directly. If false, it is treated as
   downward longwave radiation and net longwave is computed in the bulk-flux
   routine using sea-surface temperature and emissivity.

   The **qair_is_percent** flag should be set to true if the relative humidity
   in the NetCDF file is stored as a percentage (0-100) rather than as a
   fraction (0-1). This conversion is applied after loading and includes proper
   ghost cell synchronization.

Numerical Algorithms
====================

.. _list-of-parameters-16:

List of Parameters
------------------

+-----------------------------------------------+-----------------------------+-------------------+-------------+
| Parameter                                     | Definition                  | Acceptable        | Default     |
|                                               |                             |                   |             |
|                                               |                             | Values            |             |
+===============================================+=============================+===================+=============+
| **remora.tracer_horizontal_advection_scheme** | Scheme for                  | upstream3,        | upstream3   |
|                                               |                             |                   |             |
|                                               | horizontal advection        | centered4         |             |
|                                               |                             |                   |             |
|                                               | of tracers                  |                   |             |
+-----------------------------------------------+-----------------------------+-------------------+-------------+
| **remora.uv_horizontal_advection_scheme**     | Scheme for                  | upstream3,        | upstream3   |
|                                               |                             |                   |             |
|                                               | horizontal advection        | centered2         |             |
|                                               |                             |                   |             |
|                                               | of momenta                  |                   |             |
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
|                                       |                             |                                  |             |
|                                       | surface refinement of       |                                  |             |
|                                       |                             |                                  |             |
|                                       | vertical S-grid             |                                  |             |
+---------------------------------------+-----------------------------+----------------------------------+-------------+
| **remora.theta_b**                    | Stretching parameter for    | :math:`0 \leq \theta_B \leq 4`   | 0.0         |
|                                       |                             |                                  |             |
|                                       | bottom refinement of        |                                  |             |
|                                       |                             |                                  |             |
|                                       | vertical S-grid             |                                  |             |
+---------------------------------------+-----------------------------+----------------------------------+-------------+
| **remora.tcline**                     | Surface/bottom layer width  | Positive number                  | 150         |
|                                       |                             |                                  |             |
|                                       | [m] in vertical S-grid      |                                  |             |
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
|                                       |                             |              |                           |
|                                       |                             | Values       |                           |
+=======================================+=============================+==============+===========================+
| **remora.do_m3_clim_nudg**            | Whether to nudge 3D         | true / false | false                     |
|                                       |                             |              |                           |
|                                       | momentum variables          |              |                           |
|                                       |                             |              |                           |
|                                       | to climatology              |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.do_m2_clim_nudg**            | Whether to nudge 2D         | true / false | false                     |
|                                       |                             |              |                           |
|                                       | momentum variables          |              |                           |
|                                       |                             |              |                           |
|                                       | to climatology              |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.do_temp_clim_nudg**          | Whether to nudge            | true / false | false                     |
|                                       |                             |              |                           |
|                                       | temperature                 |              |                           |
|                                       |                             |              |                           |
|                                       | to climatology              |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.do_salt_clim_nudg**          | Whether to nudge            | true / false | false                     |
|                                       |                             |              |                           |
|                                       | salinity                    |              |                           |
|                                       |                             |              |                           |
|                                       | to climatology              |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.nc_clim_his_file**           | NetCDF file name for        | string       | must be set if one of     |
|                                       |                             |              |                           |
|                                       | climatology data            |              | ``do_*_clim_nudg``        |
|                                       |                             |              |                           |
|                                       |                             |              | flags is true             |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.nc_clim_coeff_file**         | NetCDF file name for        | string       | must be set if one of     |
|                                       |                             |              |                           |
|                                       | climatology nudging         |              | ``do_*_clim_nudg``        |
|                                       |                             |              |                           |
|                                       | coefficients                |              | flags is true             |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.clim_ubar_time_varname**     | name of time variable       | string       | ``ocean_time``            |
|                                       |                             |              |                           |
|                                       | for ubar climatology        |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.clim_vbar_time_varname**     | name of time variable       | string       | ``ocean_time``            |
|                                       |                             |              |                           |
|                                       | for vbar climatology        |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.clim_u_time_varname**        | name of time variable       | string       | ``ocean_time``            |
|                                       |                             |              |                           |
|                                       | for u climatology           |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.clim_v_time_varname**        | name of time variable       | string       | ``ocean_time``            |
|                                       |                             |              |                           |
|                                       | for v climatology           |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.clim_temp_time_varname**     | name of time variable       | string       | ``ocean_time``            |
|                                       |                             |              |                           |
|                                       | for temperature             |              |                           |
|                                       |                             |              |                           |
|                                       | climatology                 |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.clim_salt_time_varname**     | name of time variable       | string       | ``ocean_time``            |
|                                       |                             |              |                           |
|                                       | for salinity                |              |                           |
|                                       |                             |              |                           |
|                                       | climatology                 |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+

Rivers (point sources)
======================

These parameters are used to configure NetCDF-specified river-like point sources and sinks.

+-----------------------------+----------------------------------+--------------+-----------------------------------+
| Parameter                   | Definition                       | Acceptable   | Default                           |
|                             |                                  |              |                                   |
|                             |                                  | Values       |                                   |
+=============================+==================================+==============+===================================+
| **remora.do_rivers**        | Whether to do rivers.            | true / false | false                             |
|                             |                                  |              |                                   |
|                             | Equiavlent to ``LuvSrc`` in      |              |                                   |
|                             |                                  |              |                                   |
|                             | ROMS. Sources always             |              |                                   |
|                             |                                  |              |                                   |
|                             | apply momentum                   |              |                                   |
+-----------------------------+----------------------------------+--------------+-----------------------------------+
| **remora.nc_river_file**    | NetCDF file for river sources    | string       | must be set if                    |
|                             |                                  |              |                                   |
|                             |                                  |              | ``do_rivers``                     |
+-----------------------------+----------------------------------+--------------+-----------------------------------+
| **remora.riv_time_varname** | Name of time variable            | string       | ``river_time``                    |
+-----------------------------+----------------------------------+--------------+-----------------------------------+
| **remora.do_rivers_temp**   | Whether rivers are               | true / false | true; only used                   |
|                             |                                  |              |                                   |
|                             | temperature sources              |              | if ``do_rivers``                  |
+-----------------------------+----------------------------------+--------------+-----------------------------------+
| **remora.do_rivers_salt**   | Whether rivers are               | true / false | true; only used                   |
|                             |                                  |              |                                   |
|                             | salinity sources                 |              | if ``do_rivers``                  |
+-----------------------------+----------------------------------+--------------+-----------------------------------+
| **remora.do_rivers_scalar** | Whether rivers are passive       | true / false | false; only used                  |
|                             |                                  |              |                                   |
|                             | scalar sources                   |              | if ``do_rivers``                  |
+-----------------------------+----------------------------------+--------------+-----------------------------------+

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
