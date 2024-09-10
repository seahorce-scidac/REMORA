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
:ref:`here<icbc-parameters>`. If one of these quantities is specified from a NetCDF file, they all must be. Even if the grid is specified in the NetCDF file, the geometry and boundary parameters below must still be set. ``geometry.prob_lo`` and ``geometry.prob_hi`` do not need to agree with the file in this case.

The z-component of ``geometry.prob_lo`` should be more negative than the deepest bathymetry, and the z-compoonent of ``geometry.prob_hi`` should be 0.

List of Parameters
------------------

+--------------------------+-----------------+-------------------+-------------+
| Parameter                | Definition      | Acceptable        | Default     |
|                          |                 | Values            |             |
+==========================+=================+===================+=============+
| **geometry.prob_lo**     | physical        | [Real Real -Real] | must be set |
|                          | location of low |                   |             |
|                          | corner of the   |                   |             |
|                          | domain          |                   |             |
+--------------------------+-----------------+-------------------+-------------+
| **geometry.prob_hi**     | physical        | [Real Real 0]     | must be set |
|                          | location of     |                   |             |
|                          | high corner of  |                   |             |
|                          | the domain      |                   |             |
+--------------------------+-----------------+-------------------+-------------+
| **geometry.is_periodic** | is the domain   | 0 if false, 1     | 0 0 0       |
|                          | periodic in     | if true.          |             |
|                          | this direction  | Z-component must  |             |
|                          |                 | be zero           |             |
+--------------------------+-----------------+-------------------+-------------+

Examples of Usage
-----------------

-  **geometry.prob_lo** = 0 0 -200
   defines the low corner of the domain at (0 m,0 m,-200 m) in physical space.

-  **geometry.prob_hi** = 1.e8 2.e8 0
   defines the high corner of the domain at (1.e8 m, 2.e8 m, 0 m) in
   physical space.

-  **geometry.is_periodic** = 0 1 0
   says the domain is periodic in the y-direction only.

Domain Boundary Conditions
==========================

.. _list-of-parameters-1:

List of Parameters
------------------

+---------------+---------------------------------+-------------------+-----------------------------+
| Parameter     | Definition                      | Acceptable Values | Default                     |
+===============+=================================+===================+=============================+
| **xlo.type**  | boundary type of xlo face       | see below         | must be set if not periodic |
+---------------+---------------------------------+-------------------+-----------------------------+
| **xhi.type**  | boundary type of xhi face       | see below         | must be set if not periodic |
+---------------+---------------------------------+-------------------+-----------------------------+
| **ylo.type**  | boundary type of ylo face       | see below         | must be set if not periodic |
+---------------+---------------------------------+-------------------+-----------------------------+
| **yhi.type**  | boundary type of yhi face       | see below         | must be set if not periodic |
+---------------+---------------------------------+-------------------+-----------------------------+
| **zlo.type**  | boundary type of zlo face       | slipwall          | must be set                 |
+---------------+---------------------------------+-------------------+-----------------------------+
| **zhi.type**  | boundary type of zhi face       | slipwall          | must be set                 |
+---------------+---------------------------------+-------------------+-----------------------------+

Currently available type of boundary conditions are
``inflow``, ``outflow``, ``slipwall``, ``noslipwall``, or ``symmetry``
(Spelling of the type matters; capitalization does not.) Z-boundaries are always treated as the seafloor
and surface, and boundary type selection here will not affect program behavior. Domain boundary types
are specified on a per-side basis, rather than a per-variable basis. Inflow, outflow, etc

Usage examples can be found :ref:`here<sec:domainBCs>`.

.. _icbc-parameters:

Imposing Boundary and Initial Conditions from NetCDF File
=========================================================

Grid, initial, and time-dependent boundary data can be specified using NetCDF files, as in ROMS. REMORA expects files in the same format as ROMS.
Currently, if one of these are specified in a file, they all must be. Boundary conditions are currently applied as Dirichelet conditions, as in the ROMS ``clamped`` boundary type. More sophisticated boundary conditions like nudging are a work in progress.

The ``outflow`` option must be selected for ``xlo.type``, ``xhi.type``, ``ylo.type``, and ``yhi.type`` to read the boundary data.

List of Parameters
------------------

+---------------------------+-------------------------------+-------------+--------------------------+
| Parameter                 | Definition                    | Acceptable  | Default                  |
|                           |                               | Values      |                          |
+===========================+===============================+=============+==========================+
| **remora.ic_bc_type**     | read initial, grid, and       |             |                          |
|                           | boundary data from NetCDF     | true/false  | false                    |
|                           | files                         |             | is true                  |
+---------------------------+-------------------------------+-------------+--------------------------+
| **remora.nc_init_file_0** | initial data NetCDF file name | string      | must be set              |
|                           |                               |             | if ``remora.ic_bc_type`` |
|                           |                               |             | is true                  |
+---------------------------+-------------------------------+-------------+--------------------------+
| **remora.nc_grid_file_0** | grid data NetCDF file name    | string      | must be set              |
|                           |                               |             | if ``remora.ic_bc_type`` |
|                           |                               |             | is true                  |
+---------------------------+-------------------------------+-------------+--------------------------+
| **remora.nc_bdry_file_0** | boundary data NetCDF file     | string      | must be set              |
|                           | name                          |             | if ``remora.ic_bc_type`` |
|                           |                               |             | is true                  |
+---------------------------+-------------------------------+-------------+--------------------------+

Resolution
==========

.. _list-of-parameters-2:

List of Parameters
------------------

+---------------------------+-----------------+-----------------+-------------+
| Parameter                 | Definition      | Acceptable      | Default     |
|                           |                 | Values          |             |
+===========================+=================+=================+=============+
| **amr.n_cell**            | number of cells | Integer > 0     | must be set |
|                           | in each         |                 |             |
|                           | direction at    |                 |             |
|                           | the coarsest    |                 |             |
|                           | level           |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.max_level**         | number of       | Integer >= 0    | must be set |
|                           | levels of       |                 |             |
|                           | refinement      |                 |             |
|                           | above the       |                 |             |
|                           | coarsest level  |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.ref_ratio**         | ratio of coarse | 2 / 3 / 4       | 2 for all   |
|                           | to fine grid    | (one per level) | levels      |
|                           | spacing between |                 |             |
|                           | subsequent      |                 |             |
|                           | levels          |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.ref_ratio_vect**    | ratio of coarse | 3 integers      | 2 for all   |
|                           | to fine grid    | (one per dir)   | directions  |
|                           | spacing between | 2 / 3 / 4       |             |
|                           | subsequent      |                 |             |
|                           | levels          |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.regrid_int**        | how often to    | Integer > 0     | must be set |
|                           | regrid          |                 |             |
+---------------------------+-----------------+-----------------+-------------+
| **amr.regrid_on_restart** | should we       | 0 or 1          | 0           |
|                           | regrid          |                 |             |
|                           | immediately     |                 |             |
|                           | after           |                 |             |
|                           | restarting      |                 |             |
+---------------------------+-----------------+-----------------+-------------+

Note: if **amr.max_level** = 0 then you do not need to set
**amr.ref_ratio** or **amr.regrid_int**.

.. _examples-of-usage-2:

Examples of Usage
-----------------

-  **amr.n_cell** = 32 64 64

   would define the domain to have 32 cells in the x-direction, 64 cells
   in the y-direction, and 64 cells in the z-direction *at the coarsest level*.

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

Regridding
==========

Overview
--------

The user defines how to tag individual cells at a given level for refinement.
This list of tagged cells is sent to a grid generation routine, which uses the
Berger–Rigoutsos algorithm to create rectangular grids that contain the
tagged cells.

See :ref:`MeshRefinement` for more details on how to specify regions for
refinement.

.. _list-of-parameters-4:

List of Parameters
------------------

+----------------------------+----------------+----------------+----------------+
| Parameter                  | Definition     | Acceptable     | Default        |
|                            |                | Values         |                |
+============================+================+================+================+
| **amr.regrid_file**        | name of file   | text           | no file        |
|                            | from which to  |                |                |
|                            | read the grids |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.grid_eff**           | grid           | Real > 0, < 1  | 0.7            |
|                            | efficiency at  |                |                |
|                            | coarse level   |                |                |
|                            | at which grids |                |                |
|                            | are created    |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.n_error_buf**        | radius of      | Integer >= 0   | 1              |
|                            | additional     |                |                |
|                            | tagging around |                |                |
|                            | already tagged |                |                |
|                            | cells          |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.max_grid_size**      | maximum size   | Integer > 0    | 32             |
|                            | of a grid in   |                |                |
|                            | any direction  |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.max_grid_size**      | maximum size   | Integer        | 32             |
+----------------------------+----------------+----------------+----------------+
| **amr.blocking_factor**    | grid size must | Integer > 0    | 2              |
|                            | be a multiple  |                |                |
|                            | of this        |                |                |
+----------------------------+----------------+----------------+----------------+
| **amr.refine_grid_layout** | refine grids   | 0 if false, 1  | 1              |
|                            | more if # of   | if true        |                |
|                            | processors     |                |                |
|                            | :math:`>` # of |                |                |
|                            | grids          |                |                |
+----------------------------+----------------+----------------+----------------+

.. _notes-2:

Notes
-----

-  **amr.n_error_buf**, **amr.max_grid_size** and
   **amr.blocking_factor** can be read in as a single value which is
   assigned to every level, or as multiple values, one for each level

-  **amr.max_grid_size** at every level must be even

-  **amr.blocking_factor** at every level must be a power of 2

-  the domain size **amr.n_cell** must be a multiple of
   **amr.blocking_factor** at level 0

-  **amr.max_grid_size** must be a multiple of **amr.blocking_factor**
   at every level

.. _examples-of-usage-3:

Examples of Usage
-----------------

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

.. _subsec:grid-generation:

Gridding and Load Balancing
---------------------------

The REMORA gridding and load balancing strategy is based on that in AMReX.
See the `Gridding`_ section of the AMReX documentation for details.

.. _`Gridding`: https://amrex-codes.github.io/amrex/docs_html/ManagingGridHierarchy_Chapter.html

Simulation Time
===============

.. _list-of-parameters-5:

List of Parameters
------------------

+------------------------+---------------------------+--------------+---------+
| Parameter              | Definition                | Acceptable   | Default |
|                        |                           | Values       |         |
+========================+===========================+==============+=========+
| **max_step**           | maximum number of level 0 | Integer >= 0 | -1      |
|                        | time steps                |              |         |
+------------------------+---------------------------+--------------+---------+
| **stop_time**          | final simulation          | Real >= 0    | -1.0    |
|                        | time                      |              |         |
+------------------------+---------------------------+--------------+---------+
| **remora.start_time**  | initial simulation        | Real >= 0    | 0.0     |
|                        | time                      |              |         |
+------------------------+---------------------------+--------------+---------+

.. _notes-3:

Notes
-----

To control the number of time steps, you can limit by the maximum number
of level-0 time steps (**max_step**), or the final simulation time
(**stop_time**), or both. The code will stop at whichever criterion
comes first. Note that if the code reaches **stop_time** then the final
time step will be shortened so as to end exactly at **stop_time**, not
pass it.

.. _examples-of-usage-4:

Examples of Usage
-----------------

-  **max_step** = 1000

-  **stop_time** = 1.0

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
     ignoring the other timestep controls. Note that if
     **remora.init_shrink** :math:`\neq 1` then the first time step will in
     fact be **remora.init_shrink** \* **remora.fixed_dt**.

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
| **amr.v**                  | verbosity of     | 0 or 1                           | 0              |
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
| **remora.flat_bathymetry**       | Use flat bathymetry.        | true / false      | true        |
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

+----------------------------------+--------------------------------------+-------------------+----------------+
| Parameter                        | Definition                           | Acceptable        | Default        |
|                                  |                                      | Values            |                |
+==================================+======================================+===================+================+
| **remora.ggrav**                 | Gravitational field strength         | Real number       | 9.81           |
|                                  | [kg m/s^2]                           |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.R0**                    | Background density [kg/m^3]          | Real number       | 1028           |
|                                  | used in Linear Equation of           |                   |                |
|                                  | State. May be used in setup          |                   |                |
|                                  | of some problems.                    |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.S0**                    | Background salinity                  | Real number       | 35             |
|                                  | (nondimensional) used in             |                   |                |
|                                  | Linear Equation of State             |                   |                |
|                                  | State. May be used in setup          |                   |                |
|                                  | of some problems.                    |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.T0**                    | Background temperature               | Real number       | 5              |
|                                  | (Celsius) used in                    |                   |                |
|                                  | Linear Equation of State             |                   |                |
|                                  | State. May be used in setup          |                   |                |
|                                  | of some problems.                    |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.Tcoef**                 | Linear EOS parameter                 | Real number       | 1.7e-4         |
|                                  | (1/Celsius)                          |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.Scoef**                 | Linear EOS parameter                 | Real number       | 0.0            |
|                                  | (nondimensional)                     |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.rho0**                  | Mean density (kg/m^3) used           | Real number       | 1025           |
|                                  | when Boussinesq approx is            |                   |                |
|                                  | inferred                             |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.coriolis_type**         | Type of Coriolis forcing.            | ``beta_plane`` /  | ``beta_plane`` |
|                                  | ``beta_plane`` uses a linear         | ``custom`` /      |                |
|                                  | approximation. ``custom`` is         | ``real``          |                |
|                                  | calculated from a function in        |                   |                |
|                                  | ``prob.cpp``, and ``real`` is        |                   |                |
|                                  | read from the netcdf grid file       |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.coriolis_f0**           | f-plane constant for                 | Real number       | 0.0            |
|                                  | Coriolis param                       |                   |                |
|                                  | :math:`f = f_0 + \beta y`            |                   |                |
|                                  | when using beta plane                |                   |                |
|                                  | Coriolis type                        |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.coriolis_beta**         | beta-plane constant for              | Real number       | 0.0            |
|                                  | Coriolis param                       |                   |                |
|                                  | :math:`f = f_0 + \beta y`            |                   |                |
|                                  | when using beta plane                |                   |                |
|                                  | Coriolis type                        |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.vertical_mixing_type**  | Vertical mixing type. ``analytical`` | ``analytical`` /  | ``analytical`` |
|                                  | function is specified in             | ``GLS``           |                |
|                                  | ``prob.cpp``.                        |                   |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.gls_stability_type**    | Stability function to use for GLS    | ``Canuto_A`` /    | ``Canuto_A``   |
|                                  |                                      | ``Canuto_B`` /    |                |
|                                  |                                      | ``Galperin``      |                |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.Akv_bak**               | Minimum/initial value of Akv         | Real number       | 5.0e-6         |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.Akt_bak**               | Minimum/initial value of Akt         | Real number       | 1.0e-6         |
+----------------------------------+--------------------------------------+-------------------+----------------+
| **remora.rdrag**                 | Bottom drag                          | Real number       | 3.0e-4         |
+----------------------------------+--------------------------------------+-------------------+----------------+

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
