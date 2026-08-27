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
Biology options are documented in :ref:`sec:Fennel`. That section covers the
``remora.fennel.*`` parameters, biology initial-condition selection
(``remora.biology_ic_type``), and the options used to validate the port against
ROMS (``remora.use_biology_cpp_answer``, ``remora.biology_debug``).

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

By default, bathymetry is specified at level 0 and interpolated to the finer levels, like with any other variable. Bathymetry may also be specified at level ``remora.hires_grid_level > 0`` in file ``remora.nc_grid_file_hires``. Bathymetry on levels ``< remora.hires_grid_level`` is set by averaging down the given bathymetry. Bathymetry on higher levels is set by interpolating. High resolution bathymetry data must be given with a number of grow cells equal to the cumulative refinement ratio between level 0 and ``remora.hires_grid_level``. That is, the refined grid must fully cover the level 0 grid plus one level 0 grow cell. For example, in a problem with ``hires_grid_level = 2``, a refinement ratio of 2 between levels 0 and 1, and 3 between levels 1 and 2, ``nc_grid_file_hires`` must have 6 grow cells on each side of the domain.

Both levels must be greater than 0 (``-1`` means the data is specified at level 0 as usual); 0 is rejected, since the full-domain arrays it would need exist only for levels above 0. A high-resolution file smaller than the rule above is rejected by name rather than read partially, and note that on restart the bathymetry is restored from the checkpoint but ``pm``/``pn`` are re-read, so the high-resolution grid file must still be present.

Initial data may be specified on a high-resolution level the same way, at level ``remora.hires_init_level > 0`` in file ``remora.nc_init_file_hires``, and averaged down to the levels below. The grow cell rule is the same as for bathymetry, applied to this level: the file must have as many grow cells as the cumulative refinement ratio between level 0 and ``remora.hires_init_level``. The fields it supplies are the ones a level 0 initial file supplies, on the same terms: ``temp``, ``salt``, ``u``, ``v``, and ``zeta``, plus each dye the run carries (optional, as at level 0) and the biology tracers when ``remora.biology_ic_type`` selects NetCDF. High-resolution initialization requires ``remora.ic_type = netcdf``; it is not implemented for analytic initial conditions.

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
| **remora.nc_grid_file_hires**     | high-resolution grid data NetCDF  | string          | must be set if                  |
|                                   |                                   |                 |                                 |
|                                   |                                   |                 | ``remora.hires_grid_level``     |
|                                   |                                   |                 |                                 |
|                                   | file name                         |                 | is valid (greater than -1)      |
+-----------------------------------+-----------------------------------+-----------------+---------------------------------+
| **remora.hires_grid_level**       | level where high-resolution       | integer         | -1, meaning grid data will      |
|                                   |                                   |                 |                                 |
|                                   | grid data is specified, either    |                 | be specified at level 0         |
|                                   |                                   |                 |                                 |
|                                   | in NetCDF file or analytically    |                 |                                 |
+-----------------------------------+-----------------------------------+-----------------+---------------------------------+
| **remora.nc_init_file_hires**     | high-resolution initial data      | string          | must be set if                  |
|                                   |                                   |                 |                                 |
|                                   | NetCDF file name                  |                 | ``remora.hires_init_level``     |
|                                   |                                   |                 |                                 |
|                                   |                                   |                 | is valid (greater than -1)      |
+-----------------------------------+-----------------------------------+-----------------+---------------------------------+
| **remora.hires_init_level**       | level where high-resolution       | integer         | -1, meaning initial data        |
|                                   |                                   |                 |                                 |
|                                   | initial data is specified         |                 | will be specified at level 0    |
|                                   |                                   |                 |                                 |
|                                   | in a NetCDF file                  |                 |                                 |
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
| **remora.nc_frc_file**            | forcing data NetCDF               | string or list  | must be set                     |
|                                   |                                   |                 |                                 |
|                                   |                                   |                 | if ``remora.wind_type``         |
|                                   | file name(s)                      | of strings      |                                 |
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
|                                   | of ``u``, ``v``, ``ubar``,        |                 |                                 |
|                                   |                                   |                 |                                 |
|                                   | ``vbar``, ``zeta``, or any        |                 |                                 |
|                                   |                                   |                 |                                 |
|                                   | tracer name)                      |                 |                                 |
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

-  ``nc_frc_file`` may be either a string or a space-separated list of strings of forcing data files. They must be in time series order.

-  The time variables in the boundary files may be different for each boundary variable. Any that are not individually specified with
   ``remora.bdy_{var}_time_varname`` will default to the variable name given by ``bdy_time_varname``.

-  Every cell-centered tracer has its own boundary variable, named for the tracer itself: ``temp``, ``salt``, and then either the
   active biology tracer names or ``tracer``, ``tracer_1``, ... So a run carrying nitrate reads ``NO3_west`` and friends from the
   boundary file and accepts ``remora.bdy_NO3_time_varname``. See :ref:`sec:bc-per-tracer`.

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
| **amr.n_error_buf_{x,y}**        | radius of        | Integer >= 0;      | 1                 |
|                                  | additional       |                    |                   |
|                                  |                  | Can specify up     |                   |
|                                  | tagging around   |                    |                   |
|                                  |                  |                    |                   |
|                                  | already tagged   | to one per         |                   |
|                                  |                  | ref. level         |                   |
|                                  | cells in x or y  |                    |                   |
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

-  The blocking factor in the z-direction will be forced to a large value automatically to
   guarantee the domain will not be decomposed in the z-direction

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
| **remora.time_ref**    | reference date of the     | ``yyyymmdd`` | 0.0     |
|                        |                           |              |         |
|                        | model clock, and with it  | ``.dd``, or  |         |
|                        |                           |              |         |
|                        | the calendar. See below.  | 0, -1, -2    |         |
+------------------------+---------------------------+--------------+---------+

.. _calendar:

Reference date and calendar
---------------------------

``remora.time_ref`` is the reference date the model clock is measured from, and
it also selects the calendar, exactly as ROMS ``TIME_REF`` does:

+----------------+-----------------------+-----------------------+-------------+
| ``time_ref``   | calendar              | epoch                 | year length |
+================+=======================+=======================+=============+
| ``yyyymmdd.dd``| proleptic Gregorian   | the date given        | 365.2425 d  |
+----------------+-----------------------+-----------------------+-------------+
| ``0``          | proleptic Gregorian   | 0001-01-01 00:00:00   | 365.2425 d  |
+----------------+-----------------------+-----------------------+-------------+
| ``-1``         | 360_day: twelve       | 0000-12-30 00:00:00   | 360 d       |
|                | 30-day months, no     |                       |             |
|                | leap years            |                       |             |
+----------------+-----------------------+-----------------------+-------------+
| ``-2``         | Gregorian, as a       | 1968-05-23 00:00:00   | 365.25 d    |
|                | truncated Julian day  |                       |             |
+----------------+-----------------------+-----------------------+-------------+

The fractional part of a positive value is a time of day, so
``remora.time_ref = 20020115.5`` is 15 January 2002 at 12:00. A value below
``-2`` names no calendar and is an error.

Model time is elapsed time since that epoch, so ``remora.start_time`` offsets
the run within the calendar the way ROMS ``DSTART`` does. With the default
``time_ref = 0``, a run starting on 1 January 2020 sets
``remora.start_time = 63713433600``.

Only features that need a date consult this: at present the time-dependent
atmospheric CO2 options of the Fennel biology model, see :ref:`sec:Fennel`. The
conversion is ``remora_caldate`` in ``Source/Utils/REMORA_DateClock.H``, a port
of ROMS ``caldate``.

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
| **remora.use_curvilinear_grid**  | Add curvilinear grid terms  | true / false      | false       |
|                                  |                             |                   |             |
|                                  | for advection               |                   |             |
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

+------------------------------------------+----------------------------------------+------------------------+----------------+
| Parameter                                | Definition                             | Acceptable             | Default        |
|                                          |                                        |                        |                |
|                                          |                                        | Values                 |                |
+==========================================+========================================+========================+================+
| **remora.ggrav**                         | Gravitational field strength           | Real number            | 9.81           |
|                                          |                                        |                        |                |
|                                          | [kg m/s^2]                             |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.eos_type**                      | Which equation of state to use.        | Linear or              | Linear         |
|                                          |                                        |                        |                |
|                                          | Nonlinear is UNESCO                    | Nonlinear              |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.R0**                            | Background density [kg/m^3]            | Real number            | 1028           |
|                                          |                                        |                        |                |
|                                          | used in Linear Equation of             |                        |                |
|                                          |                                        |                        |                |
|                                          | State. May be used in setup            |                        |                |
|                                          |                                        |                        |                |
|                                          | of some problems.                      |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.S0**                            | Background salinity                    | Real number            | 35             |
|                                          |                                        |                        |                |
|                                          | (nondimensional) used in               |                        |                |
|                                          |                                        |                        |                |
|                                          | Linear Equation of State               |                        |                |
|                                          |                                        |                        |                |
|                                          | State. May be used in setup            |                        |                |
|                                          |                                        |                        |                |
|                                          | of some problems.                      |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.T0**                            | Background temperature                 | Real number            | 5              |
|                                          |                                        |                        |                |
|                                          | (Celsius) used in                      |                        |                |
|                                          |                                        |                        |                |
|                                          | Linear Equation of State               |                        |                |
|                                          |                                        |                        |                |
|                                          | State. May be used in setup            |                        |                |
|                                          |                                        |                        |                |
|                                          | of some problems.                      |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.Tcoef**                         | Linear EOS parameter                   | Real number            | 1.7e-4         |
|                                          |                                        |                        |                |
|                                          | (1/Celsius)                            |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.Scoef**                         | Linear EOS parameter                   | Real number            | 0.0            |
|                                          |                                        |                        |                |
|                                          | (nondimensional)                       |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.rho0**                          | Mean density (kg/m^3) used             | Real number            | 1025           |
|                                          |                                        |                        |                |
|                                          | when Boussinesq approx is              |                        |                |
|                                          |                                        |                        |                |
|                                          | inferred                               |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.coriolis_type**                 | Type of Coriolis forcing.              | ``beta_plane`` /       | ``beta_plane`` |
|                                          |                                        |                        |                |
|                                          | ``beta_plane`` uses a linear           | ``analytic`` /         |                |
|                                          |                                        |                        |                |
|                                          | approximation. ``analytic`` is         | ``netcdf``             |                |
|                                          |                                        |                        |                |
|                                          | calculated from a function in          |                        |                |
|                                          |                                        |                        |                |
|                                          | ``prob.cpp``, and ``netcdf`` is        |                        |                |
|                                          |                                        |                        |                |
|                                          | read from the netcdf grid file         |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.coriolis_f0**                   | f-plane constant for                   | Real number            | 0.0            |
|                                          |                                        |                        |                |
|                                          | Coriolis param                         |                        |                |
|                                          |                                        |                        |                |
|                                          | :math:`f = f_0 + \beta (y - y_c)`      |                        |                |
|                                          |                                        |                        |                |
|                                          | when using beta plane Coriolis         |                        |                |
|                                          |                                        |                        |                |
|                                          | type. :math:`y` is measured from       |                        |                |
|                                          |                                        |                        |                |
|                                          | the southern domain boundary,          |                        |                |
|                                          |                                        |                        |                |
|                                          | :math:`y_c` is the domain center       |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.coriolis_beta**                 | beta-plane constant for                | Real number            | 0.0            |
|                                          |                                        |                        |                |
|                                          | Coriolis param                         |                        |                |
|                                          |                                        |                        |                |
|                                          | :math:`f = f_0 + \beta (y - y_c)`      |                        |                |
|                                          |                                        |                        |                |
|                                          | when using beta plane Coriolis         |                        |                |
|                                          |                                        |                        |                |
|                                          | type. :math:`y` is measured from       |                        |                |
|                                          |                                        |                        |                |
|                                          | the southern domain boundary,          |                        |                |
|                                          |                                        |                        |                |
|                                          | :math:`y_c` is the domain center       |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.horizontal_mixing_type**        | Horizontal mixing type.                | ``analytic`` /         | ``analytic``   |
|                                          |                                        |                        |                |
|                                          | ``analytic`` means function is         | ``constant``           |                |
|                                          |                                        |                        |                |
|                                          | specified in ``prob.cpp``.             | / ``scaled_to_grid``   |                |
|                                          |                                        |                        |                |
|                                          | ``scaled_to_grid`` scales harmonic     |                        |                |
|                                          |                                        |                        |                |
|                                          | viscosity/diffusivity by the grid      |                        |                |
|                                          |                                        |                        |                |
|                                          | cell area. Equivalent to ``DIFF_GRID`` |                        |                |
|                                          |                                        |                        |                |
|                                          | and ``VISC_GRID`` in ROMS.             |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.scaled_to_grid_amr_scaling**    | AMR scaling behavior for               | ``none`` /             | ``none``       |
|                                          |                                        |                        |                |
|                                          | ``scaled_to_grid`` and ``constant``    | ``linear``             |                |
|                                          |                                        |                        |                |
|                                          | coefficients on refined levels.        |                        |                |
|                                          |                                        |                        |                |
|                                          | ``linear`` decreases coefficients in   |                        |                |
|                                          |                                        |                        |                |
|                                          | proportion to the horizontal           |                        |                |
|                                          |                                        |                        |                |
|                                          | refinement ratio.                      |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.visc2**                         | Constant horizontal viscosity,         | Real number            | 0.0            |
|                                          |                                        |                        |                |
|                                          | everywhere. Needed when                |                        |                |
|                                          |                                        |                        |                |
|                                          | ``horizontal_mixing_type`` is          |                        |                |
|                                          |                                        |                        |                |
|                                          | ``constant`` or ``scaled_to_grid``     |                        |                |
|                                          |                                        |                        |                |
|                                          | (in this case, it is the maximum       |                        |                |
|                                          |                                        |                        |                |
|                                          | viscosity over the domain).            |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.tnu2_salt**                     | Constant horizontal diffusivity,       | Real number            | 0.0            |
|                                          |                                        |                        |                |
|                                          | everywhere for salt.                   |                        |                |
|                                          |                                        |                        |                |
|                                          | Needed when                            |                        |                |
|                                          |                                        |                        |                |
|                                          | ``horizontal_mixing_type`` is          |                        |                |
|                                          |                                        |                        |                |
|                                          | ``constant`` or ``scaled_to_grid``     |                        |                |
|                                          |                                        |                        |                |
|                                          | (in this case, it is the maximum       |                        |                |
|                                          |                                        |                        |                |
|                                          | salt diffusivity over the domain).     |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.tnu2_temp**                     | Constant horizontal diffusivity,       | Real number            | 0.0            |
|                                          |                                        |                        |                |
|                                          | everywhere for temperature.            |                        |                |
|                                          |                                        |                        |                |
|                                          | Needed when                            |                        |                |
|                                          |                                        |                        |                |
|                                          | ``horizontal_mixing_type``             |                        |                |
|                                          |                                        |                        |                |
|                                          | is ``constant`` or ``scaled_to_grid``  |                        |                |
|                                          |                                        |                        |                |
|                                          | (in this case, it is the maximum       |                        |                |
|                                          |                                        |                        |                |
|                                          | temperature diffusivity over the       |                        |                |
|                                          |                                        |                        |                |
|                                          | domain).                               |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.nscalar**                       | Number of passive (dye) scalars        | Integer >= 0           | 0              |
|                                          |                                        |                        |                |
|                                          | in addition to temperature and         |                        |                |
|                                          |                                        |                        |                |
|                                          | salinity. Does not count biology       |                        |                |
|                                          |                                        |                        |                |
|                                          | tracers.                               |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.tnu2_scalar**                   | Constant horizontal diffusivity,       | Real number            | 0.0            |
|                                          |                                        |                        |                |
|                                          | everywhere for passive scalar.         |                        |                |
|                                          |                                        |                        |                |
|                                          | Needed when                            |                        |                |
|                                          |                                        |                        |                |
|                                          | ``horizontal_mixing_type``             |                        |                |
|                                          |                                        |                        |                |
|                                          | is ``constant`` or ``scaled_to_grid``  |                        |                |
|                                          |                                        |                        |                |
|                                          | (in this case, it is the maximum       |                        |                |
|                                          |                                        |                        |                |
|                                          | scalar diffusivity over the domain).   |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.tnu2_{var}**                    | Constant horizontal diffusivity for    | Real number            | value of       |
|                                          |                                        |                        |                |
|                                          | one named tracer, overriding           |                        | ``tnu2_scalar``|
|                                          |                                        |                        |                |
|                                          | ``tnu2_scalar``. ``{var}`` is a tracer |                        |                |
|                                          |                                        |                        |                |
|                                          | name: ``tracer``, ``tracer_1``,        |                        |                |
|                                          |                                        |                        |                |
|                                          | ``NO3``, and so on. ``tnu2_temp`` and  |                        |                |
|                                          |                                        |                        |                |
|                                          | ``tnu2_salt`` above are this same form |                        |                |
|                                          |                                        |                        |                |
|                                          | for the first two components.          |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.harmonic_mixing_type**          | Whether harmonic mixing (tracers)      | ``s`` /                | ``s``          |
|                                          |                                        |                        |                |
|                                          | is calculated along s- or geopotential | ``geopotential``       |                |
|                                          |                                        |                        |                |
|                                          | surfaces.                              |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.vertical_mixing_type**          | Vertical mixing type. ``analytic``     | ``analytic`` /         | ``analytic``   |
|                                          |                                        |                        |                |
|                                          | function is specified in               | ``GLS``                |                |
|                                          |                                        |                        |                |
|                                          | ``prob.cpp``.                          |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.gls_stability_type**            | Stability function to use for GLS      | ``Canuto_A`` /         | ``Canuto_A``   |
|                                          |                                        |                        |                |
|                                          |                                        | ``Canuto_B`` /         |                |
|                                          |                                        |                        |                |
|                                          |                                        | ``Galperin``           |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.Akv_bak**                       | Minimum/initial value of Akv           | Real number            | 5.0e-6         |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.Akt_bak**                       | Minimum/initial value of Akt, applied  | Real number            | 1.0e-6         |
|                                          |                                        |                        |                |
|                                          | to temperature and salinity.           |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.Akt_bak_temp**                  | Minimum/initial value of Akt for       | Real number            | value of       |
|                                          |                                        |                        |                |
|                                          | temperature, overriding ``Akt_bak``.   |                        | ``Akt_bak``    |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.Akt_bak_salt**                  | Minimum/initial value of Akt for       | Real number            | value of       |
|                                          |                                        |                        |                |
|                                          | salinity, overriding ``Akt_bak``.      |                        | ``Akt_bak``    |
|                                          |                                        |                        |                |
|                                          | Passive tracers mix with this value.   |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+
| **remora.bulk_fluxes**                   | Whether to use bulk fluxes             | true / false           | false          |
|                                          |                                        |                        |                |
|                                          | parametrization                        |                        |                |
+------------------------------------------+----------------------------------------+------------------------+----------------+

Passive scalars
---------------

REMORA always includes temperature and salinity as the first two conserved
state variables. Passive scalars begin at component ``Tracer_comp = 2``, and the
biology tracers, if any, follow them:

.. code:: text

   temp  salt  | tracer  tracer_1  ...  | NO3  NH4  ...
   0     1     | Tracer_comp = 2        | Bio_comp = 2 + nscalar

The total number of conserved components is

.. math::

   ncons = 2 + \texttt{remora.nscalar} + n_{bio}

so ``remora.nscalar`` counts only the passive scalars — not temperature,
salinity, or biology tracers.

For example:

- ``remora.nscalar = 1`` gives ``temp``, ``salt``, and ``tracer``.
- ``remora.nscalar = 2`` gives ``temp``, ``salt``, ``tracer``, and ``tracer_1``.

Additional passive scalars continue with the names ``tracer_2``,
``tracer_3``, and so on.

``remora.nscalar`` defaults to 0, so dye is opt-in: a run carries temperature and
salinity only until an input asks for more. A tracer nothing asked for is still
advected and diffused every step, and still has to be accounted for in plotfiles and
in boundary and climatology files, so carrying one by default made runs that never
mentioned dye behave as though they had.

Dye and biology coexist freely, so ``remora.nscalar = 2`` with a Fennel model carrying
11 tracers gives 15 components in the order shown above. See :ref:`sec:Fennel`.

With ``remora.ic_type = netcdf`` each dye's initial field is read from the initial
file under its own name -- ``tracer``, ``tracer_1``, and so on -- stored with the
same ``(ocean_time, s_rho, eta_rho, xi_rho)`` layout ROMS uses for ``temp`` and
``salt``. The read is the last thing ``init_data_from_netcdf`` does: dye follows
``remora.ic_type`` along with the physical fields rather than having an
initialization-source flag of its own the way biology does.

Unlike ``temp`` and ``salt``, these fields are optional. A dye whose variable the
initial file does not carry starts at zero, and the run prints the names it did not
find, so an initial file written before the run carried dye still works unchanged and
an idealized file need only carry the dyes it actually wants to specify. The same read
runs on ``remora.nc_init_file_hires`` when initial data is given on a high-resolution
level.

A dye's other file-driven inputs follow the same naming convention -- ``tracer_west``
and friends for boundary data, ``tracer`` for a climatology field, and
``tracer_NudgeCoef`` for a spatially varying nudging coefficient -- but unlike the
initial field they are not optional. A boundary condition or nudging flag that calls
for data the file does not have stops the run and names what is missing; only
``tracer_NudgeCoef`` falls back, to the constant derived from ``remora.tnudg``. See
:ref:`sec:bc-per-tracer` and :ref:`sec:clim-per-tracer`.

Horizontal and vertical mixing coefficients follow different rules, as they do in ROMS.

*Horizontal* diffusivity is per tracer. ``remora.tnu2_scalar`` sets the value every
tracer takes by default -- biology tracers included -- and ``remora.tnu2_{var}``
overrides a single tracer by name:

.. code:: python

   remora.tnu2_scalar  = 5.0      # every dye and biology tracer
   remora.tnu2_NO3     = 50.0     # except NO3

``remora.tnu2_temp`` and ``remora.tnu2_salt`` are the pre-existing spelling of the same
per-name form for the first two components, and they still default to zero rather than
to ``tnu2_scalar``.

*Vertical* diffusivity is not. ROMS computes and stores ``Akt`` for the active tracers
alone, and every passive tracer -- dye and biology alike -- mixes with the salinity
coefficient (``ltrc=MIN(itrc,NAT)`` in ``pre_step3d`` and ``step3d_t``). REMORA does the
same: ``vec_Akt`` has two components, and ``akt_comp()`` in ``REMORA_IndexDefines.H``
performs the mapping.

.. code:: python

   remora.Akt_bak      = 1.0e-6   # temperature and salinity
   remora.Akt_bak_temp = 1.0e-5   # except temperature
   remora.Akt_bak_salt = 1.0e-6   # salinity, and every passive tracer with it

There is deliberately no ``remora.Akt_bak_{var}`` for a passive tracer; supplying one
is an error rather than a silent no-op, since it would name a coefficient that no part
of the vertical diffusion reads.

Note that ``tnu2`` is only consulted when ``remora.horizontal_mixing_type`` is
``constant`` or ``scaled_to_grid``; under the default ``analytic`` the problem's
``init_analytic_hmix`` hook sets the coefficients itself. The same is true of ``Akt``
and ``init_analytic_vmix``; the hooks in ``Source/Prob/`` set the two active-tracer
components directly and ignore ``Akt_bak``, but this can be changed in other problems.

The current analytic *initial condition* hooks seed the first dye tracer only, so with
more than one dye the rest start at zero. Again, this can be changed by the user in
other problems as needed.

The ``remora.sum_interval`` diagnostic reports a volume-weighted sum for every
tracer, labelled by name. Both plotfile writers check every tracer they are about
to write for NaN and inf, and name the offending tracer if one is found; a tracer
absent from ``remora.plot_vars_3d`` is not checked. Both writers also accept
``stflux_{var}`` in ``remora.plot_vars_2d`` for any tracer, so the surface flux array
can be inspected per tracer. Only the temperature and salinity entries are ever filled,
by the bulk-flux or coupling paths; the rest read zero, because Fennel's air-sea gas
exchange acts on the tracer directly rather than through that array.

Refinement accepts any tracer name as its ``field_name``, so a biology tracer can
drive AMR:

.. code:: python

   remora.refinement_indicators = hi_no3
   remora.hi_no3.max_level      = 1
   remora.hi_no3.field_name     = NO3
   remora.hi_no3.value_greater  = 15.0

River sources also work per tracer; see :ref:`sec:rivers`.

Scaled-to-grid horizontal mixing
--------------------------------

If ``remora.horizontal_mixing_type = "scaled_to_grid"``, REMORA follows the ROMS-style approach of
scaling horizontal harmonic mixing coefficients by the grid cell area.

.. math::

   G(i,j) = \sqrt{\frac{1}{pm(i,j)\,pn(i,j)}}

.. math::

   \nu(i,j) = \nu_0 \frac{G(i,j)}{\max(G)}, \qquad
   \kappa_n(i,j) = \kappa_{0,n} \frac{G(i,j)}{\max(G)}

where :math:`\nu_0` is ``remora.visc2`` and :math:`\kappa_{0,n}` are the tracer diffusivities
(``remora.tnu2_temp``, ``remora.tnu2_salt``, ``remora.tnu2_scalar``). This ensures the *maximum*
coefficient over the normalization region equals the user-specified value, while varying spatially
with grid size. Note that if the largest cell area occurs over land, then the maximum over *wet*
cells (and thus what you see after applying ``mask_rho`` in post-processing) may be smaller than
the user-specified value.

- The normalization ``max(G)`` is computed as a global maximum over the level-0 grid (rho points at
  the surface, i.e. ``k=0``) and does not use land/sea masks. Equivalently, it uses the maximum grid-cell
  area :math:`A(i,j) = 1/(pm\,pn)` via :math:`\max(G)=\sqrt{\max(A)}`.

- AMR refinement scaling: if ``remora.scaled_to_grid_amr_scaling = "linear"``, then on AMR level
  :math:`\ell` the coefficients are additionally scaled by the cumulative horizontal refinement ratio.

  .. math::

     \frac{1}{\prod_{m=0}^{\ell-1} \sqrt{r_x(m)\,r_y(m)}}

  For example, with a refinement ratio of ``5 5 1``, level 1 coefficients are reduced by a factor of 5 relative to level 0.

- Ghost cells for the coefficient fields are filled using the same periodic/foextrap boundary fill
  used elsewhere in REMORA. This is done so stencil-based operations (e.g., the psi-point averaging for
  ``visc2_p``) have valid neighbor values near domain boundaries and coarse/fine interfaces. Output
  (plotfiles / NetCDF) writes only the valid region.

When this option is active, the spatially varying coefficients are automatically included in plot output
as 2D (vertically homogeneous) fields:

- ``visc2`` (horizontal viscosity at rho points)
- ``diff2_temp``, ``diff2_salt``, ``diff2_tracer`` (horizontal diffusivities at rho points)
- additional passive scalars appear as ``diff2_tracer_1``, ``diff2_tracer_2``, and so on

Geopotential rotated harmonic tracer diffusion
----------------------------------------------

Harmonic tracer diffusion can be rotated along geopotential (constant-:math:`z`)
surfaces when ``remora.harmonic_mixing_type = "geopotential"``. This formulation
reduces spurious diapycnal mixing over steeply sloping bathymetry, where terrain-following
:math:`s`-levels intersect isopycnal surfaces. This approach corresponds
to the ROMS ``MIX_GEO_TS`` option (``Nonlinear/t3dmix2_geo.h``), although full algorithmic
details are not documented in the ROMS implementation.

Let :math:`z_r(i,j,k)` denote the geopotential (rho-point) vertical coordinate. Local
surface slopes are defined using metric-weighted discrete differences:

.. math::

   \begin{aligned}
   S_x &\equiv dZdx \approx c_x(i,j)\,(z_r(i,j,k)-z_r(i-1,j,k)) \\
   S_y &\equiv dZde \approx c_y(i,j)\,(z_r(i,j,k)-z_r(i,j-1,k))
   \end{aligned}

where :math:`c_x` and :math:`c_y` are C-grid metric factors that include inverse grid spacing
and land–sea masking. These are constructed as face-centered averages:

.. math::

   c_x(i,j) = \tfrac{1}{2}\left(pm(i,j) + pm(i-1,j)\right)\, msku(i,j),

   c_y(i,j) = \tfrac{1}{2}\left(pn(i,j) + pn(i,j-1)\right)\, mskv(i,j).

The rotated diffusion operator can be interpreted in flux-form as:

.. math::

   F_x = -K_h\,H\left(\partial_x T - S_x\,\partial_z T\right),
   \qquad
   F_y = -K_h\,H\left(\partial_y T - S_y\,\partial_z T\right),

where :math:`K_h` is the harmonic diffusivity and :math:`H` is the face-averaged
vertical cell thickness (``Hz``).

In practice, all gradients are computed using finite differences:

- ``dTdx``, ``dTde``: centered horizontal differences on cell faces
- ``dTdz``: vertical differences along rho columns
- ``dZdx``, ``dZde``: metric-weighted slope fields

Diffusivity is interpolated to cell faces and combined with face-averaged vertical thicknesses
prior to flux construction.

The x-face flux (``FX``) is computed as a face-centered diffusivity–thickness product multiplied
by a slope-corrected horizontal tracer gradient:

.. math::

   FX_{i+1/2,j,k} =
   K_{i+1/2,j,k}\,H_{i+1/2,j,k}
   \left[
      dTdx_{i+1/2,j,k}
      - \mathcal{R}_x(S_x, \partial_z T)
   \right],

where :math:`\mathcal{R}_x` denotes a slope-dependent estimate of how much vertical stratification
contaminates the horizontal gradient. This term uses a sign-aware (min/max) stencil that selects
locally appropriate vertical neighbor averages of :math:`\partial_z T` based on the sign of the slope.

Specifically, the reconstruction is given by:

.. math::

   \mathcal{R}_x =
   \tfrac{1}{2}\Big(
   \min(S_x,0)\,(dTdz_{\text{down}}+dTdz_{\text{up}+1})
   +
   \max(S_x,0)\,(dTdz_{\text{down}+1}+dTdz_{\text{up}})
   \Big).

The y-face flux (``FE``) is constructed analogously using :math:`S_y` and the corresponding
y-direction stencil.

A separate vertical coupling term (``FS``) accounts for cross-directional slope–gradient
interactions between horizontal and vertical derivatives. It is constructed using similar
sign-dependent decompositions (min/max splitting) that select locally consistent horizontal
and vertical neighbor contributions.

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

Bulk-flux atmospheric inputs use per-variable source selectors when not received
from the driver in a coupled simulation. The source
selectors accept ``constant``, ``analytic``, or ``netcdf``. ``lwrad_type`` and
``eminusp_type`` also accept ``computed``, which uses REMORA's internal longwave
or evaporation-minus-rain diagnostic path.

The source selector parameters are:

+---------------------------+-----------------------------------------------+-----------------------------+-------------+
| Parameter                 | Definition                                    | Acceptable values           | Default     |
+===========================+===============================================+=============================+=============+
| **remora.uwind_type**     | Source for u-direction wind                   | constant, analytic, netcdf  | analytic    |
+---------------------------+-----------------------------------------------+-----------------------------+-------------+
| **remora.vwind_type**     | Source for v-direction wind                   | constant, analytic, netcdf  | analytic    |
+---------------------------+-----------------------------------------------+-----------------------------+-------------+
| **remora.tair_type**      | Source for air temperature ``Tair``           | constant, analytic, netcdf  | constant    |
+---------------------------+-----------------------------------------------+-----------------------------+-------------+
| **remora.qair_type**      | Source for air humidity ``qair``              | constant, analytic, netcdf  | constant    |
+---------------------------+-----------------------------------------------+-----------------------------+-------------+
| **remora.pair_type**      | Source for air pressure ``Pair``              | constant, analytic, netcdf  | constant    |
+---------------------------+-----------------------------------------------+-----------------------------+-------------+
| **remora.swrad_type**     | Source for shortwave radiation ``swrad``      | constant, analytic, netcdf  | constant    |
+---------------------------+-----------------------------------------------+-----------------------------+-------------+
| **remora.lwrad_type**     | Source for external longwave radiation        | computed, constant,         | computed    |
|                           |                                               | analytic, netcdf            |             |
+---------------------------+-----------------------------------------------+-----------------------------+-------------+
| **remora.cloud_type**     | Source for cloud fraction                     | constant, analytic, netcdf  | constant    |
+---------------------------+-----------------------------------------------+-----------------------------+-------------+
| **remora.rain_type**      | Source for precipitation rate                 | constant, analytic, netcdf  | constant    |
+---------------------------+-----------------------------------------------+-----------------------------+-------------+
| **remora.eminusp_type**   | Source for prescribed ``EminusP``             | computed, constant,         | computed    |
|                           |                                               | analytic, netcdf            |             |
+---------------------------+-----------------------------------------------+-----------------------------+-------------+

The legacy paired wind selector **remora.wind_type** is still accepted, with
accepted values ``analytic`` or ``netcdf``. When ``remora.uwind_type`` or
``remora.vwind_type`` is not set, ``remora.wind_type`` sets the corresponding
wind component source.

Constant-valued bulk inputs are set with:

+--------------------------------------+----------------------------------------+-------------------+----------------+
| Parameter                            | Definition                             | Acceptable values | Default        |
+======================================+========================================+===================+================+
| **remora.uwind**                     | Constant u-direction wind [m/s         | Real number       | 0.0            |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.vwind**                     | Constant v-direction wind [m/s]        | Real number       | 0.0            |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.air_temperature**           | Constant air temperature [C]           | Real number       | 23.567         |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.air_humidity**              | Constant humidity [fraction or kg/kg]  | Real number       | 0.776          |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.air_pressure**              | Constant air pressure [hPa]            | Real number       | 1013.48        |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.surface_radiation_flux**    | Constant shortwave radiation [W/m^2]   | Real number       | 0.0            |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.longwave_radiation_flux**   | Constant external longwave [W/m^2]     | Real number       | 0.0            |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.cloud**                     | Constant cloud cover fraction          | 0 to 1            | 0.0            |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.rain**                      | Constant precipitation rate [kg/m^2/s] | Real number       | 0.0            |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.EminusP**                   | Constant prescribed E-P [kg/m^2/s]     | Real number       | 0.0            |
+--------------------------------------+----------------------------------------+-------------------+----------------+

Additional bulk-flux controls:

+--------------------------------------+----------------------------------------+-------------------+----------------+
| Parameter                            | Definition                             | Acceptable values | Default        |
+======================================+========================================+===================+================+
| **remora.longwave_down**             | Treat external longwave as downward    | true / false      | false          |
|                                      | radiation and compute net longwave     |                   |                |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.longwave_is_net**           | Interpret external longwave as net     | true / false      | false          |
|                                      |                                        |                   |                |
|                                      | longwave and use it as-is              |                   |                |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.longwave_netcdf_varname**   | Name of the NetCDF longwave variable   | String            | ``lwrad``      |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.qair_is_percent**           | Convert NetCDF qair from 0-100 percent | true / false      | false          |
|                                      | to 0-1 fraction                        |                   |                |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.eminusp_correct_ssh**       | Adjust sea surface height for active   | true / false      | false          |
|                                      | E-P source                             |                   |                |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.blk_ZQ**                    | Height [m] of humidity measurements    | Real number       | 10.0           |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.blk_ZT**                    | Height [m] of temperature measurements | Real number       | 10.0           |
+--------------------------------------+----------------------------------------+-------------------+----------------+
| **remora.blk_ZW**                    | Height [m] of wind measurements        | Real number       | 10.0           |
+--------------------------------------+----------------------------------------+-------------------+----------------+

When a source selector is ``constant``, spatially uniform values are set
from values specified in the inputs file.

When a source selector is ``analytic``, REMORA calls the problem-defined
``init_analytic_bulk_flux`` hook and passes pointers only for the fields that
were selected as analytic.

When a source selector is ``netcdf``, the field is
read from **remora.nc_frc_file** using variable names ``Uwind``, ``Vwind``,
``Tair``, ``qair``, ``Pair``, ``swrad``, ``rain``, ``cloud``, ``EminusP``,
and the configured longwave variable name. Multiple files may be specified,
in which case in which case REMORA concatenates their time axes in the
order listed in the inputs file. Forcing variables must have the same dimensions
as the level 0 (coarsest) model grid. Time interpolation is performed automatically
based on the simulation time at level 0, and fields are interpoalted to finer
AMR levels as needed.

The legacy boolean inputs **Tair_from_netcdf**, **qair_from_netcdf**,
**Pair_from_netcdf**, **srflx_from_netcdf**, **longwave_down_from_netcdf**,
**rain_from_netcdf**, **cloud_from_netcdf**, and **EminusP_from_netcdf** are
still accepted for compatibility; setting one to true maps the corresponding
``*_type`` selector to ``netcdf``. The legacy **longwave_netcdf_is_net** input
is also accepted as an alias for ``longwave_is_net``; when true and no
``lwrad_type`` is provided, it selects ``netcdf`` to preserve existing inputs.

The **qair_is_percent** flag should be set to true if the relative humidity
in the NetCDF file is stored as a percentage (0-100) rather than as a
fraction (0-1).

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
| **remora.do_{var}_clim_nudg**         | Whether to nudge tracer     | true / false | false                     |
|                                       |                             |              |                           |
|                                       | ``{var}`` to climatology.   |              |                           |
|                                       |                             |              |                           |
|                                       | ``{var}`` is any tracer     |              |                           |
|                                       |                             |              |                           |
|                                       | name: ``temp``, ``salt``,   |              |                           |
|                                       |                             |              |                           |
|                                       | ``NO3``, ...                |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+
| **remora.nc_clim_his_file**           | NetCDF file name(s) for     | string or    | must be set if one of     |
|                                       |                             |              |                           |
|                                       | climatology data            | list         | ``do_*_clim_nudg``        |
|                                       |                             |              |                           |
|                                       |                             | of strings   | flags is true             |
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
| **remora.clim_{var}_time_varname**    | name of time variable       | string       | ``ocean_time``            |
|                                       |                             |              |                           |
|                                       | for the climatology of      |              |                           |
|                                       |                             |              |                           |
|                                       | tracer ``{var}``            |              |                           |
+---------------------------------------+-----------------------------+--------------+---------------------------+

.. note::

   For AMR runs, climatology fields are read on level 0 and temporally
   interpolated there. REMORA then interpolates those fields to finer AMR
   levels for nudging updates. ``remora.nc_clim_his_file`` may be a single file
   or a space-separated list of files. When multiple files are provided, they
   must be listed in time series order.

.. _sec:clim-per-tracer:

Climatology nudging for individual tracers
------------------------------------------

Every cell-centered tracer can be nudged toward its own climatology, not just
temperature and salinity. The flag is keyed by the tracer's own name, so a run
carrying nitrate uses

.. code-block:: python

   remora.do_NO3_clim_nudg = true

exactly as temperature uses ``remora.do_temp_clim_nudg``. Each tracer turned on
this way needs two things in the input files:

- A 3D field named for the tracer in ``remora.nc_clim_his_file``: ``temp``,
  ``salt``, ``NO3``, and so on. This is the ROMS convention, and it is the same
  name the variable carries elsewhere. If the field is missing, REMORA aborts and
  names it rather than failing deep inside the NetCDF reader.
- Optionally, a spatially varying nudging coefficient named
  ``{var}_NudgeCoef`` in ``remora.nc_clim_coeff_file`` (in inverse days), again
  following ROMS: ``temp_NudgeCoef``, ``salt_NudgeCoef``, ``NO3_NudgeCoef``.

The coefficient is optional in a way the climatology field is not. A tracer with
no ``{var}_NudgeCoef`` in the file falls back to the constant coefficient derived
from ``remora.tnudg``, and REMORA prints which tracers fell back. This means a
nudging coefficient file written for temperature and salinity alone still works
when biology tracers are nudged: those tracers simply use the uniform timescale.

Each tracer may also take its own climatology time axis via
``remora.clim_{var}_time_varname``, which defaults to ``ocean_time``.

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
| **remora.nc_river_file**    | NetCDF file(s) for river         | string or    | must be set if                    |
|                             |                                  |              |                                   |
|                             | sources                          | list         | ``do_rivers``                     |
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
|                             | scalar sources. Default for      |              | if ``do_rivers``                  |
|                             |                                  |              |                                   |
|                             | every dye tracer.                |              |                                   |
+-----------------------------+----------------------------------+--------------+-----------------------------------+
| **remora.do_rivers_{var}**  | Whether rivers are a source of   | true / false | see below                         |
|                             |                                  |              |                                   |
|                             | tracer ``{var}``. Overrides the  |              |                                   |
|                             |                                  |              |                                   |
|                             | defaults above.                  |              |                                   |
+-----------------------------+----------------------------------+--------------+-----------------------------------+

.. note::

   ``remora.nc_river_file`` may be either a single file or a space-separated
   list of files. If multiple river files are provided, they must be listed in
   time series order.

.. _sec:rivers:

River input for individual tracers
----------------------------------

Any cell-centered tracer can take river input, keyed by the tracer's own name, so
temperature and salinity keep ``remora.do_rivers_temp`` and
``remora.do_rivers_salt`` and a biology tracer uses

.. code:: python

   remora.do_rivers_NO3 = true

The defaults are: temperature and salinity on, every passive (dye) scalar following
``remora.do_rivers_scalar``, and biology tracers off. Biology defaults off because a
river concentration for a biogeochemical tracer has to be a deliberate choice, not
something inherited from a dye setting. All of these are ignored unless
``remora.do_rivers`` is true.

Each enabled tracer needs a field named for it in ``remora.nc_river_file``, with
dimensions ``(river_time, s_rho, river)``, following the ROMS convention:
``river_temp``, ``river_salt``, ``river_tracer``, ``river_NO3``, and so on. If the
field is missing REMORA aborts and names it, rather than failing inside the NetCDF
reader.

Tracer fields are concentrations, so a tracer given without an ``s_rho`` dimension,
as ``(river_time, river)``, is used unchanged at every vertical level, and REMORA
prints a warning naming the field. ``river_Vshape`` applies only to
``river_transport``, where it distributes the total transport over the vertical as
ROMS does with ``Qsrc = Qbar * Qshape``; it is not applied to tracers.

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
