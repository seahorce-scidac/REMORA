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

+--------------------------------+------------------+-----------------------+------------+
| Parameter                      | Definition       | Acceptable            | Default    |
|                                |                  | Values                |            |
+================================+==================+=======================+============+
| **remora.plotfile_type**       | AMReX or NETCDF  | "amrex" or            | "amrex"    |
|                                |                  | "netcdf / "NetCDF"    |            |
+--------------------------------+------------------+-----------------------+------------+
| **remora.write_history_file**  | do we write      | false or true         | false      |
|                                | netcdf files at  |                       |            |
|                                | each timestep    |                       |            |
|                                | or one file for  |                       |            |
|                                | all timesteps?   |                       |            |
+--------------------------------+------------------+-----------------------+------------+
| **remora.plot_file**           | prefix for       | String                | “plt”      |
|                                | plotfiles        |                       |            |
+--------------------------------+------------------+-----------------------+------------+
| **remora.plot_int**            | how often (by    | Integer               | -1         |
|                                | level-0 time     | :math:`> 0`           |            |
|                                | steps) to write  |                       |            |
|                                | plot files       |                       |            |
+--------------------------------+------------------+-----------------------+------------+
| **remora.plot_vars**           | name of          | list of names         | None       |
|                                | variables to     |                       |            |
|                                | include in       |                       |            |
|                                | plotfiles        |                       |            |
+--------------------------------+------------------+-----------------------+------------+

.. _notes-5:

Notes
-----

-  The NeTCDF option is only available if REMORA has been built with USE_NETCDF enabled.

-  The write_history_file option is only available if **plotfile_type = netcdf**

-  If  **plotfile_type = netcdf** and **write_history_file = false**, the frequency
   will be determined by **plot_int**

-  Velocity components are defined on faces within the REMORA code, but are averaged onto
   cell centers when written in amrex/native plotfiles. They are not averaged when writing
   NetCDF files.

-  File prefixes can include directories.

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
