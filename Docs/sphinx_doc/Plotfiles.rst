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

"Plotfiles" can be written very efficiently in parallel in a native AMReX format
or in HDF5.  They can also be written in NetCDF.

The following options in the inputs file control the generation of plotfiles.
Note that plotfiles can be written at two different frequencies; the names,
frequency and content of the two streams are controlled separately.

.. _list-of-parameters-9:

List of Parameters
------------------

+------------------------------+------------------+-----------------------+------------+
| Parameter                    | Definition       | Acceptable            | Default    |
|                              |                  | Values                |            |
+==============================+==================+=======================+============+
| **remora.plotfile_type**     | AMReX or NETCDF  | "amrex" or            | "amrex"    |
|                              |                  | "netcdf / "NetCDF"    |            |
+------------------------------+------------------+-----------------------+------------+
| **remora.plot_file_1**       | prefix for       | String                | “*plt_1_*” |
|                              | plotfiles        |                       |            |
|                              | at first freq.   |                       |            |
+------------------------------+------------------+-----------------------+------------+
| **remora.plot_file_2**       | prefix for       | String                | “*plt_2_*” |
|                              | plotfiles        |                       |            |
|                              | at seoncd freq.  |                       |            |
+------------------------------+------------------+-----------------------+------------+
| **remora.plot_int_1**        | how often (by    | Integer               | -1         |
|                              | level-0 time     | :math:`> 0`           |            |
|                              | steps) to write  |                       |            |
|                              | plot files       |                       |            |
|                              | at first freq.   |                       |            |
+------------------------------+------------------+-----------------------+------------+
| **remora.plot_int_2**        | how often (by    | Integer               | -1         |
|                              | level-0 time     | :math:`> 0`           |            |
|                              | steps) to write  |                       |            |
|                              | plot files       |                       |            |
|                              | at seoncd freq.  |                       |            |
+------------------------------+------------------+-----------------------+------------+
| **remora.plot_vars_1**       | name of          | list of names         | None       |
|                              | variables to     |                       |            |
|                              | include in       |                       |            |
|                              | plotfiles        |                       |            |
|                              | at first freq.   |                       |            |
+------------------------------+------------------+-----------------------+------------+
| **remora.plot_vars_2**       | name of          | list of names         | None       |
|                              | variables to     |                       |            |
|                              | include in       |                       |            |
|                              | plotfiles        |                       |            |
|                              | at seoncd freq.  |                       |            |
+------------------------------+------------------+-----------------------+------------+

.. _notes-5:

Notes
-----

-  The NeTCDF option is only available if REMORA has been built with USE_NETCDF enabled.

-  Velocity components are defined on faces within the REMORA code, but are averaged onto
   cell centers when written in amrex/native plotfiles.

-  File prefixes can include directories.

.. _examples-of-usage-8:

Examples of Usage
-----------------

-  **remora.plotfile_type** = *amrex*

-  **remora.plot_file_1** = *out/plt_run*

-  **remora.plot_int_1** = 10

   means that native plot files (actually directories) starting with the prefix
   “*plt_run*” will be generated every 10 level-0 time steps in the directory
   `out`. If using
   amrex format, that directory names will be *plt_run00000*, *plt_run00010*,
   *plt_run00020*, etc. If using NetCDF format, the names will have ".nc" appended.

..   In addition, while the amrex plotfiles will contain data at all of the refinement
   levels,  NetCDF files are separated by level.
