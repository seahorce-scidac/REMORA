.. role:: cpp(code)
  :language: c++

Diffusive Physics
=================

.. _list-of-parameters-12:

List of Parameters
------------------

+----------------------------------+--------------------+---------------------+-------------+
| Parameter                        | Definition         | Acceptable          | Default     |
|                                  |                    | Values              |             |
+==================================+====================+=====================+=============+
| **romsx.spatial_order**          |                    |  2 / 3 / 4 / 5 / 6  | 2           |
+----------------------------------+--------------------+---------------------+-------------+

Forcing Terms
=============

.. _list-of-parameters-19:

List of Parameters
------------------

+----------------------------------+-----------------------------+-------------------+-------------+
| Parameter                        | Definition                  | Acceptable        | Default     |
|                                  |                             | Values            |             |
+==================================+=============================+===================+=============+
| **romsx.use_coriolis**           | Include Coriolis terms.     | true / false      | false       |
|                                  | Coriolis parameter :math`f` |                   |             |
|                                  | is hard-coded for Upwelling |                   |             |
|                                  | problem in                  |                   |             |
|                                  | ``ROMSX::Advance()``        |                   |             |
+----------------------------------+-----------------------------+-------------------+-------------+


Initialization
==============

ROMSX can be initialzed in different ways. These are listed below:

- Custom initialization:
    Several problems under **Exec** are initialized in a custom manner. The state and velocity components are specific to the problem. These problems are meant for demonstration and do not include any terrain or map scale factors.
- Initialization using a NetCDF file:
    Problems in ROMSX can be initialized using a NetCDF file containing the mesoscale data. The state and velocity components of the ROMSX domain are ingested from the mesocale data.
- Initialization using an ``input_sounding`` file:
    Problems in ROMSX can be initialized using an ``input_sounding`` file containing the vertical profile. This file has the same format as used by ``ideal.exe`` executable in WRF. Using this option for initialization, running ``ideal.exe`` and reading from the resulting ``wrfinput_d01`` file are not needed. This option is used for initializing ROMSX domain to a horizontally homogeneous mesoscale state and does not include terrain or map scale factors.

List of Parameters
------------------

+-----------------------------+-------------------+--------------------+------------+
| Parameter                   | Definition        | Acceptable         | Default    |
|                             |                   | Values             |            |
+=============================+===================+====================+============+
| **romsx.init_type**         | Initialization    | “custom”,          | “*custom*” |
|                             | type              | “ideal”,           |            |
|                             |                   | “real”,            |            |
|                             |                   | "input_sounding"   |            |
+-----------------------------+-------------------+--------------------+------------+
| **romsx.nc_init_file**      | NetCDF file with  |  String            | NONE       |
|                             | initial mesocale  |                    |            |
|                             | data              |                    |            |
+-----------------------------+-------------------+--------------------+------------+
| **romsx.nc_bdy_file**       | NetCDF file with  |  String            | NONE       |
|                             | mesocale data at  |                    |            |
|                             | lateral boundaries|                    |            |
+-----------------------------+-------------------+--------------------+------------+

Notes
-----------------

If **romsx.init_type = ideal**”, the problem is initialized with mesoscale data contained in a NetCDF file, provided via ``romsx.nc_init_file``. The mesoscale data are horizontally homogeneous, i.e., there is variation only in vertical direction.

If **romsx.init_type = real**”, the problem is initialized with mesoscale data contained in a NetCDF file, provided via ``romsx.nc_init_file``. The mesoscale data are realistic with variation in all three directions.  In addition, the lateral boundary conditions must be supplied in a NetCDF files specified by **romsx.nc_bdy_file = wrfbdy_d01**”

If **romsx.init_type = custom**” or **romsx.init_type = input_sounding**”, ``romsx.nc_init_file`` and ``romsx.nc_bdy_file`` do not need to be set.

