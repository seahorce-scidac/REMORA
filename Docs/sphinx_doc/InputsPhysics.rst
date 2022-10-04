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
| **romsx.alpha_T**                  | Diffusion coeff.   | Real                | 0.0         |
|                                  | for temperature    |                     |             |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.alpha_C**                  | Diffusion coeff.   | Real                | 0.0         |
|                                  | for scalar         |                     |             |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.rho0_trans**               | Reference density  | Real                | 1.0         |
|                                  | to compute const.  |                     |             |
|                                  | rho*Alpha          |                     |             |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.les_type**                 | Using an LES       | "None",             | "None"      |
|                                  | model, and if so,  | "Smagorinsky",      |             |
|                                  | which type?        | "Deardorff"         |             |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.molec_diff_type**          | Using molecular    | "None",             | "None"      |
|                                  | viscosity and      | "Constant", or      |             |
|                                  | diffusivity?       | "ConstantAlpha"     |             |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.dynamicViscosity**         | Viscous coeff. if  | Real                | 0.0         |
|                                  | DNS                |                     |             |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.Cs**                       | Constant           | Real                | 0.0         |
|                                  | Smagorinsky coeff. |                     |             |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.Pr_t**                     | Turbulent Prandtl  | Real                | 1.0         |
|                                  | Number             |                     |             |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.Sc_t**                     | Turbulent Schmidt  | Real                | 1.0         |
|                                  | Number             |                     |             |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.spatial_order**            |                    |  2 / 3 / 4 / 5 / 6  | 2           |
+----------------------------------+--------------------+---------------------+-------------+

Note: in the equations for the evolution of momentum, potential temperature and advected scalars, the
diffusion coefficients are written as :math:`\mu`, :math:`\rho \alpha_T` and :math:`\rho \alpha_C`, respectively.

If we set ``romsx.molec_diff_type`` to ``Constant``, then

- ``romsx.dynamicViscosity`` is used as the value of :math:`\mu` in the momentum equation, and

- ``romsx.alpha_T`` is multiplied by ``romsx.rho0_trans`` to form the coefficient for potential temperature, and

- ``romsx.alpha_C`` is multiplied by ``romsx.rho0_trans`` to form the coefficient for an advected scalar.

If we set ``romsx.molec_diff_type`` to ``ConstantAlpha``, then

- the dynamic viscosity in the momentum equation is assumed to have the form :math:`\mu = \rho \alpha_M`
  where :math:`\alpha_M` is a momentum diffusivity constant with units of kinematic viscosity, calculated as
  ``romsx.dynamicViscosity`` divided by ``romsx.rho0_trans``;
  this diffusivity is multiplied by the current density :math:`\rho` to form the coefficient in the momentum equation; and

- ``romsx.alpha_T`` is multiplied by the current density :math:`\rho` to form the coefficient for potential temperature, and

- ``romsx.alpha_C`` is multiplied by the current density :math:`\rho` to form the coefficient for an advected scalar.


PBL Scheme
==========

.. _list-of-parameters-13:

List of Parameters
------------------

+----------------------------------+--------------------+---------------------+-------------+
| Parameter                        | Definition         | Acceptable          | Default     |
|                                  |                    | Values              |             |
+==================================+====================+=====================+=============+
| **romsx.pbl_type**                 | Name of PBL Scheme | "None", "MYNN2.5"   | "None"      |
|                                  | to be used         |                     |             |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.pbl_A1**                   | MYNN Constant A1   | Real                | 1.18        |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.pbl_A2**                   | MYNN Constant A2   | Real                | 0.665       |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.pbl_B1**                   | MYNN Constant B1   | Real                | 24.0        |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.pbl_B2**                   | MYNN Constant B2   | Real                | 15.0        |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.pbl_C1**                   | MYNN Constant C1   | Real                | 0.137       |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.pbl_C2**                   | MYNN Constant C1   | Real                | 0.75        |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.pbl_C3**                   | MYNN Constant C3   | Real                | 0.352       |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.pbl_C4**                   | MYNN Constant C4   | Real                | 0.0         |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.pbl_C5**                   | MYNN Constant C5   | Real                | 0.2         |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.advect_QKE**               | Include advection  | bool                | 1           |
|                                  | terms in QKE eqn   |                     |             |
+----------------------------------+--------------------+---------------------+-------------+
| **romsx.diffuse_QKE_3D**           | Include horizontal | bool                | 0           |
|                                  | turb. diffusion    |                     |             |
|                                  | terms in QKE eqn.  |                     |             |
+----------------------------------+--------------------+---------------------+-------------+

Note that the MYNN2.5 scheme must be used in conjunction with a MOST boundary condition
at the surface (Zlo) boundary.

If the PBL scheme is activated, it determines the turbulent diffusivity in the vertical
direction. If an LES model is also specified, it determines only the horizontal turbulent
diffusivity.

Right now, the QKE equation is solved if and only if the MYNN2.5 PBL model is selected. In that
transport equation, it is optional to advect QKE, and to apply LES diffusive transport for QKE
in the horizontal directions (the veritcal component is always computed as part of the PBL
scheme).

Forcing Terms
=============

.. _list-of-parameters-14:

List of Parameters
------------------

+----------------------------------+-------------------+-------------------+-------------+
| Parameter                        | Definition        | Acceptable        | Default     |
|                                  |                   | Values            |             |
+==================================+===================+===================+=============+
| **romsx.abl_driver_type**          | Type of external  | None,             | None        |
|                                  | forcing term      | PressureGradient  |             |
|                                  |                   | GeostrophicWind   |             |
+----------------------------------+-------------------+-------------------+-------------+
| **romsx.abl_pressure_grad**        | Pressure gradient | 3 Reals           | (0.,0.,0.)  |
|                                  | forcing term      |                   |             |
|                                  | (only if          |                   |             |
|                                  | abl.driver_type = |                   |             |
|                                  | PressureGradient) |                   |             |
+----------------------------------+-------------------+-------------------+-------------+
| **romsx.abl_geo_wind**             | Geostrophic       | 3 Reals           | (0.,0.,0.)  |
|                                  | forcing term      |                   |             |
|                                  | (only if          |                   |             |
|                                  | abl.driver_type = |                   |             |
|                                  | GeostrophicWind)  |                   |             |
+----------------------------------+-------------------+-------------------+-------------+
| **romsx.use_gravity**              | Include gravity   | true / false      | false       |
|                                  | in momentum       |                   |             |
|                                  | update?  If true, |                   |             |
|                                  | there is buoyancy |                   |             |
+----------------------------------+-------------------+-------------------+-------------+
| **romsx.use_coriolis**             | Include Coriolis  | true / false      | false       |
|                                  | forcing           |                   |             |
+----------------------------------+-------------------+-------------------+-------------+
| **romsx.use_rayleigh_damping**     | Include explicit  | true / false      | false       |
|                                  | Rayleigh damping  |                   |             |
+----------------------------------+-------------------+-------------------+-------------+


Initialization
==============

ROMSX can be initialzed in different ways. These are listed below:

- Custom initialization:
    Several problems under **Exec** are initialized in a custom manner. The state and velocity components are specific to the problem. These problems are meant for demonstration and do not include any terrain or map scale factors.
- Initialization using a NetCDF file:
    Problems in ROMSX can be initialized using a NetCDF file containing the mesoscale data. The state and velocity components of the ROMSX domain are ingested from the mesocale data. This is a more realistic problem with real atmospheric data used for initialization. The typical filename used for initialization is ``wrfinput_d01``, which is the outcome of running ``ideal.exe`` or ``real.exe`` of the WPS/WRF system.  These problems are run with both terrain and map scale factors.
- Initialization using an ``input_sounding`` file:
    Problems in ROMSX can be initialized using an ``input_sounding`` file containing the vertical profile. This file has the same format as used by ``ideal.exe`` executable in WRF. Using this option for initialization, running ``ideal.exe`` and reading from the resulting ``wrfinput_d01`` file are not needed. This option is used for initializing ROMSX domain to a horizontally homogeneous mesoscale state and does not include terrain or map scale factors.

List of Parameters
------------------

+-----------------------------+-------------------+--------------------+------------+
| Parameter                   | Definition        | Acceptable         | Default    |
|                             |                   | Values             |            |
+=============================+===================+====================+============+
| **romsx.init_type**           | Initialization    | “custom”,          | “*custom*” |
|                             | type              | “ideal”,           |            |
|                             |                   | “real”,            |            |
|                             |                   |"input_sounding"    |            |
+-----------------------------+-------------------+--------------------+------------+
| **romsx.nc_init_file**        | NetCDF file with  |  String            | NONE       |
|                             | initial mesocale  |                    |            |
|                             | data              |                    |            |
+-----------------------------+-------------------+--------------------+------------+
| **romsx.nc_bdy_file**         | NetCDF file with  |  String            | NONE       |
|                             | mesocale data at  |                    |            |
|                             | lateral boundaries|                    |            |
+-----------------------------+-------------------+--------------------+------------+

Notes
-----------------

If **romsx.init_type = ideal**”, the problem is initialized with mesoscale data contained in a NetCDF file, provided via ``romsx.nc_init_file``. The mesoscale data are horizontally homogeneous, i.e., there is variation only in vertical direction.

If **romsx.init_type = real**”, the problem is initialized with mesoscale data contained in a NetCDF file, provided via ``romsx.nc_init_file``. The mesoscale data are realistic with variation in all three directions.  In addition, the lateral boundary conditions must be supplied in a NetCDF files specified by **romsx.nc_bdy_file = wrfbdy_d01**”

If **romsx.init_type = custom**” or **romsx.init_type = input_sounding**”, ``romsx.nc_init_file`` and ``romsx.nc_bdy_file`` do not need to be set.

Terrain Smoothing
=================

Currently, ROMSX has 3 methods of controlling the terrain-fitted coordinates:

- Basic Terain Following (BTF):
    The influence of the terrain decreases linearly with height.
- Smoothed Terrain Following (STF):
    Small-scale terrain structures are progressively smoothed out of the coordinate system as height increases.
- Sullivan Terrain Following (name TBD):
    The influence of the terrain decreases with the cube of height.

List of Parameters
------------------

+-----------------------------+-------------------+--------------------+------------+
| Parameter                   | Definition        | Acceptable         | Default    |
|                             |                   | Values             |            |
+=============================+===================+====================+============+
| **romsx.terrain_smoothing**   | specify terrain   | 0,                 | 0          |
|                             | following         | 1,                 |            |
|                             |                   | 2                  |            |
+-----------------------------+-------------------+--------------------+------------+


Examples of Usage
-----------------

-  **romsx.terrain_smoothing**  = 0
    BTF is used when generating the terrain following coordinate.

-  **romsx.terrain_smoothing**  = 1
    STF is used when generating the terrain following coordinate.

-  **romsx.terrain_smoothing**  = 2
    Sullivan TF is used when generating the terrain following coordinate.
