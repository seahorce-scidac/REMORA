.. role:: cpp(code)
  :language: c++

.. _sec:Plotfiles:

*********
Plotfiles
*********
.. toctree::
   :maxdepth: 1

Controlling PlotFile Generation
===============================

"Plotfiles" can be written very efficiently in parallel in a native AMReX format or NetCDF (via PnetCDF).

The following options in the inputs file control the generation of plotfiles.

.. _list-of-parameters-9:

List of Parameters
------------------

+----------------------------------------+-----------------------------------+-----------------------+------------+
| Parameter                              | Definition                        | Acceptable            | Default    |
|                                        |                                   |                       |            |
|                                        |                                   | Values                |            |
+========================================+===================================+=======================+============+
| **remora.plotfile_type**               | AMReX or NETCDF plotfiles         | ``amrex`` or          | ``amrex``  |
|                                        |                                   |                       |            |
|                                        |                                   | ``netcdf``            |            |
+----------------------------------------+-----------------------------------+-----------------------+------------+
| **remora.file_min_digits**             | Minimum number of digits          | Integer >= 0          | 5          |
|                                        |                                   |                       |            |
|                                        | in iteration number appended      |                       |            |
|                                        |                                   |                       |            |
|                                        | to plotfile, checkpoint,          |                       |            |
|                                        |                                   |                       |            |
|                                        | or chunked history files          |                       |            |
+----------------------------------------+-----------------------------------+-----------------------+------------+
| **remora.write_history_file**          | do we write                       | false or true         | true       |
|                                        | netcdf files at                   |                       |            |
|                                        |                                   |                       |            |
|                                        | each timestep or one file         |                       |            |
|                                        |                                   |                       |            |
|                                        | for                               |                       |            |
|                                        | all timesteps?                    |                       |            |
+----------------------------------------+-----------------------------------+-----------------------+------------+
| **remora.chunk_history_file**          | do we divide netcdf history       | false or true         | false      |
|                                        |                                   |                       |            |
|                                        | files so that each file contains  |                       |            |
|                                        |                                   |                       |            |
|                                        | only a certain number of time     |                       |            |
|                                        |                                   |                       |            |
|                                        | steps?                            |                       |            |
+----------------------------------------+-----------------------------------+-----------------------+------------+
| **remora.steps_per_history_file**      | Maximum number of steps per       | integer               | -1         |
|                                        |                                   |                       |            |
|                                        | netcdf history file. If <=0,      |                       |            |
|                                        |                                   |                       |            |
|                                        | calculate automatically such      |                       |            |
|                                        |                                   |                       |            |
|                                        | that each file is less than 2GB   |                       |            |
+----------------------------------------+-----------------------------------+-----------------------+------------+
| **remora.plot_file**                   | prefix for                        | String                | “plt”      |
|                                        | plotfiles                         |                       |            |
+----------------------------------------+-----------------------------------+-----------------------+------------+
| **remora.plot_int**                    | how often (by                     | Integer               | -1         |
|                                        | level-0 time                      | :math:`> 0`           |            |
|                                        |                                   |                       |            |
|                                        | steps) to write                   |                       |            |
|                                        | plot files                        |                       |            |
+----------------------------------------+-----------------------------------+-----------------------+------------+
| **remora.plot_int_time**               | how often (in simulation time     | Real                  | -1.0       |
|                                        |                                   | :math:`> 0`           |            |
|                                        |                                   |                       |            |
|                                        | seconds) to write                 |                       |            |
|                                        | plot files                        |                       |            |
+----------------------------------------+-----------------------------------+-----------------------+------------+
| **remora.plot_vars_3d**                | name of                           | list of names         | None       |
|                                        | 3D variables to                   |                       |            |
|                                        | include in                        | (see table below)     |            |
|                                        |                                   |                       |            |
|                                        | plotfiles. Not                    |                       |            |
|                                        | used for netCDF                   |                       |            |
+----------------------------------------+-----------------------------------+-----------------------+------------+
| **remora.plot_vars_2d**                | name of                           | list of names         | None       |
|                                        | 2D variables to                   |                       |            |
|                                        | include in                        | (see table below)     |            |
|                                        |                                   |                       |            |
|                                        | plotfiles. Not                    |                       |            |
|                                        | used for netCDF                   |                       |            |
+----------------------------------------+-----------------------------------+-----------------------+------------+
| **remora.plot_staggered_vels**         | whether to output velocities      | false or true         | false      |
|                                        |                                   |                       |            |
|                                        | on cell faces. They will be in    |                       |            |
|                                        |                                   |                       |            |
|                                        | the UFace/VFace/WFace             |                       |            |
|                                        |                                   |                       |            |
|                                        | multifab. Not used for netCDF     |                       |            |
+----------------------------------------+-----------------------------------+-----------------------+------------+
| **remora.expand_plotvars_to_unif_rr**  | whether to expand a multilevel    | false or true         | false      |
|                                        |                                   |                       |            |
|                                        | plotfile to have a uniform        |                       |            |
|                                        |                                   |                       |            |
|                                        | refinement ratio on levels > 0.   |                       |            |
|                                        |                                   |                       |            |
|                                        | This is necessary for amrvis.     |                       |            |
+----------------------------------------+-----------------------------------+-----------------------+------------+

.. _notes-5:

Notes
-----

-  The NetCDF option is only available if REMORA has been built with USE_PNETCDF enabled.

-  The write_history_file option is only available if **plotfile_type = netcdf**

-  Depending on your PnetCDF build, the code may be unable to write files larger than 2 GB. If the code
   crashes when writing a NetCDF history file (or a single time step, if you have a particularly large grid),
   consider building with MPICH v4.2.2 or instead outputting a native AMReX plotfile instead.

-  Velocity components are defined on faces within the REMORA code, but are averaged onto
   cell centers when written in amrex/native plotfiles. They are not averaged when writing
   NetCDF files.

-  File prefixes can include directories.

-  If both ``remora.plot_int`` and ``remora.plot_int_time`` have been set, plotfile output will occur
  ``plot_int`` steps or ``plot_int_time`` simulation seconds after the last plotfile, whichever happens first.

3D Plotfile Field Options
--------------------------

+--------------------------------+---------------------------+
| Field                          | Definition                |
|                                |                           |
+================================+===========================+
| **salt**                       | salinity                  |
+--------------------------------+---------------------------+
| **temp**                       | temperature               |
+--------------------------------+---------------------------+
| **scalar**                     | passive scalar            |
+--------------------------------+---------------------------+
| **x_velocity**                 | velocity in x-direction   |
+--------------------------------+---------------------------+
| **y_velocity**                 | velocity in y-direction   |
+--------------------------------+---------------------------+
| **z_velocity**                 | velocity in z-direction   |
+--------------------------------+---------------------------+
| **vorticity**                  | vorticity                 |
+--------------------------------+---------------------------+

2D Plotfile Field Options
--------------------------

+--------------------------------+---------------------------+
| Field                          | Definition                |
|                                |                           |
+================================+===========================+
| **zeta**                       |                           |
+--------------------------------+---------------------------+
| **h**                          |                           |
+--------------------------------+---------------------------+
| **ubar**                       |                           |
+--------------------------------+---------------------------+
| **sustr**                      |                           |
+--------------------------------+---------------------------+
| **bustr**                      |                           |
+--------------------------------+---------------------------+
| **vbar**                       |                           |
+--------------------------------+---------------------------+
| **svstr**                      |                           |
+--------------------------------+---------------------------+
| **bvstr**                      |                           |
+--------------------------------+---------------------------+

.. _examples-of-usage-8:

Examples of Usage
-----------------

-  **remora.plotfile_type** = *amrex*

-  **remora.plot_file** = *out/plt_run*

-  **remora.plot_int** = 10

   means that native plot files (actually directories) starting with the prefix
   “*plt_run*” will be generated every 10 level-0 time steps in the directory
   `out`. If using
   amrex format, that directory names will be *plt_run00000*, *plt_run00010*,
   *plt_run00020*, etc. If using NetCDF format, the names will have ".nc" appended.

   AMReX plotfiles will contain data at all of the refinement levels. NetCDF files
   will not be output if there is more than one level.
