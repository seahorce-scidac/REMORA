
 .. role:: cpp(code)
    :language: c++

.. _sec:domainBCs:

Physical/Domain Boundary Conditions
===================================

There are two primary types of physical/domain boundary conditions: those which rely only on the
data in the valid regions, and those which rely on externally specified values.

REMORA allows users to specify types of boundary condition with keywords in the inputs file.
The option ``remora.boundary_per_variable`` controls whether conditions are specified on a
per-side (false) or per-variable (true) basis. Conditions are only set on the x- and y-faces.
Boundaries in the z-direction always correspond to the sea floor and surface. Bottom drag is specified with
``remora.rdrag`` (see :ref:`Physics Parameters<list-of-parameters-15>`), and surface wind
stress is specified by the function ``init_analytic_smflux`` in ``prob.cpp``.

REMORA also implements the ROMS bulk fluxes parameterization with optional evaporation-precipitation. The
parametrization can be turned on with ``remora.bulk_fluxes``. Further options can be found in the
:ref:`parameter list <list-of-parameters-bulk-fluxes>`.


Boundary per side
-----------------

To set boundary conditions per domain side for all variables, set
``remora.boundary_per_variable = false``. This is the default behavior.

The information for each face is preceded by
``remora.bc.xlo.type``, ``remora.bc.xhi.type``, ``remora.bc.ylo.type``, or ``remora.bc.yhi.type``. Spelling of the type matters; capitalization does not. The
options for boundary conditions per side are listed below.

- periodic
- inflow
- outflow
- slipwall
- noslipwall
- symmetry
- clamped

Each of these types of physical boundary condition has a mapping to a mathematical boundary condition
for each type; this is summarized in the table below, along with the corresponding ROMS boundary conditions.
If periodic is selected, it must be used for both low and high faces in a direction. The ``remora.is_periodic``
flag must match as described in :ref:`ProblemGeometry <geometry-parameters>`.

For example, setting

::

    remora.bc.xlo.type = "Inflow"
    remora.bc.xhi.type = "Outflow"
    remora.bc.ylo.type = "periodic" #optional
    remora.bc.yhi.type = "periodic" #optional

    remora.is_periodic = 0 1 0

would define a problem with inflow in the low-\ :math:`x` direction,
outflow in the high-\ :math:`x` direction, periodic in the :math:`y`-direction,
and slip wall on the low and high :math:`y`-faces, and
Note that no keyword is needed for a periodic boundary, here only the
specification in ``remora.is_periodic`` is needed.

Boundary per variable
---------------------

To set different boundary conditions for different variables on each side of the domain, set ``remora.boundary_per_variable = true``.
If true, conditions must be set for all of the following variables.

- temperature: ``remora.bc.temp.type``
- salinity: ``remora.bc.salt.type``
- additional scalars: ``remora.bc.scalar.type``, or one entry per tracer (see below)
- 3D u-velocity: ``remora.bc.u.type``
- 3D v-velocity: ``remora.bc.v.type``
- 2D u-velocity: ``remora.bc.ubar.type``
- 2D v-velocity: ``remora.bc.vbar.type``
- sea surface height: ``remora.bc.zeta.type``
- turbulent kinetic energy: ``remora.bc.tke.type``

They must be set to a list of four conditions in the ROMS order West, South, East, North,
mapping onto low-\ :math:`x`, low-\ :math:`y`, high-\ :math:`x`, high-\ :math:`y`. South is the
low-\ :math:`y` face, the one ``remora.bc.ylo.type`` names and ``*_south`` boundary data fills.
The options are

- periodic
- inflow
- outflow
- slipwall
- noslipwall
- symmetry
- clamped
- chapman
- flather
- orlanski_rad
- orlanski_rad_nudg

The corresponding ROMS conditions are indicated in the tables below.

For example, setting

::

    #                        West      South    East      North
    remora.bc.temp.type   =  periodic  clamped  periodic  clamped
    remora.bc.salt.type   =  periodic  clamped  periodic  clamped
    remora.bc.scalar.type =  periodic  clamped  periodic  clamped
    remora.bc.u.type      =  periodic  clamped  periodic  clamped
    remora.bc.v.type      =  periodic  clamped  periodic  clamped
    remora.bc.w.type      =  periodic  clamped  periodic  clamped
    remora.bc.ubar.type   =  periodic  flather  periodic  flather
    remora.bc.vbar.type   =  periodic  flather  periodic  flather
    remora.bc.zeta.type   =  periodic  chapman  periodic  chapman
    remora.bc.tke.type    =  periodic  outflow  periodic  outflow

    remora.is_periodic = 1 0 0

will define a problem that is periodic on the Western and Eastern sides. Temperature, salinity, passive scalar,
3D u-velocity, and 3D v-velocity will be clamped to values given in a NetCDF file specified by
:ref:`remora.nc_bdry_file<icbc-parameters>`. The 2D momentum and zeta BCs are calculated from the
Chapman/Flather conditions, with nudging towards values given in the boundary NetCDF file.

Each velocity entry sets that component's own condition. On the West and East sides, ``u`` and
``ubar`` are the normal components and ``v`` and ``vbar`` the tangential ones; on the South and North
sides it is the reverse. So ``slipwall``, ``symmetry`` or ``flather`` takes its normal-velocity form
on a side only where the normal component's entry names it.

.. _sec:tidal-boundaries:

Tidal boundary forcing
----------------------

Fixed tidal harmonics can be added to the boundary elevation and depth-averaged normal
velocity on Flather boundaries. Configure those boundaries as usual, then set::

    remora.tide_type = netcdf
    remora.nc_tide_file = boundary_tides.nc
    remora.nc_bdry_file = subtidal_boundary.nc

``remora.tide_type`` defaults to ``none``. Reading tides requires a PnetCDF build.
The subtidal boundary file remains required and must cover the simulation time range.
For tides-only forcing, supply zeros at two times spanning the run. Prepare subtidal boundary
values with their tidal signal removed before adding these harmonics.

Every constituent in the tidal file participates in the sum. Amplitudes and phases are fixed
in time and may vary along the boundary. Names are labels; periods determine frequencies.
Changing the constituent set requires a different file, with no constituent-selection input.

Reconstruction uses the same requested time as the subtidal interpolation, including the
current practice of using the start of the level's baroclinic step throughout its fast steps.

Where tides are applied
~~~~~~~~~~~~~~~~~~~~~~~

Tides are applied only on a side whose normal barotropic velocity is ``flather``: ``ubar`` on
the West and East, ``vbar`` on the South and North. On such a side the tide is added to two
boundary values:

- **Normal barotropic velocity.** Flather relaxes the model's normal ``ubar`` or ``vbar`` toward
  the subtidal value plus the tide.
- **Elevation.** The tide is added to the file's ``zeta`` on that side. Flather uses this
  elevation, and so does zeta's own condition there, which must read the boundary file:

  - ``chapman`` radiates the model's elevation and does not use the file value;
  - ``clamped`` holds the elevation at the subtidal value plus the tide;
  - ``orlanski_rad_nudg`` nudges toward the subtidal value plus the tide.

The tangential barotropic velocity keeps its own condition and receives no tide. With
``flather`` it radiates; with ``clamped`` it is held at the subtidal file value.

Every other side gets no tides, and the run prints a warning naming each such side at startup.
This includes a ``clamped`` normal velocity, radiation conditions, and walls. Periodic sides
have no boundary and print nothing.

Conventions
~~~~~~~~~~~

The harmonic sum follows `ROMS set_tides.F
<https://github.com/myroms/roms/blob/develop/ROMS/Nonlinear/set_tides.F>`_: the time argument,
the phase lag in degrees, and the elevation averaged onto the boundary face. REMORA departs from
ROMS in three ways:

- **Sides.** ROMS adds tides wherever the boundary data is read, clamped sides included.
- **Velocity harmonics.** ROMS builds tidal velocities from tidal-ellipse parameters for both
  components. REMORA reads amplitudes and phases of the normal velocity only.
- **Adding to the file values.** ROMS replaces the boundary values with the tide unless it is
  built with ``ADD_FSOBC``/``ADD_M2OBC``. REMORA always adds the tide to the subtidal values.

The sum is

.. math::

    q_{\mathrm{tide}}(s,t) = \sum_c A_c(s)
        \cos\left[\frac{2\pi(t-t_0)}{T_c}
                   - \phi_c(s)\frac{\pi}{180}\right].

Here :math:`q` is elevation or normal barotropic velocity, :math:`s` is position along the
boundary, and :math:`t_0` is the phase-reference epoch. The external reference is
:math:`q_{\mathrm{subtidal}}+q_{\mathrm{tide}}`.

.. list-table:: Units and signs
   :header-rows: 1
   :widths: 25 20 55

   * - Quantity
     - File units
     - Convention
   * - Period
     - ``s``
     - Positive period in seconds.
   * - Elevation amplitude
     - ``m``
     - Nonnegative amplitude. Positive elevation is upward from the subtidal reference.
   * - Velocity amplitude
     - ``m s-1`` or ``meter second-1``
     - Nonnegative amplitude of depth-averaged velocity, rather than transport.
   * - Phase
     - ``degrees``
     - Positive phase is a lag. One constituent's positive maximum occurs at
       :math:`t=t_0+T\phi/360`, modulo its period. Radians are rejected.
   * - ``ubar``
     - ``m s-1``
     - Positive in the model x direction on both west and east boundaries.
   * - ``vbar``
     - ``m s-1``
     - Positive in the model y direction on both south and north boundaries.

Phases are normally in the range 0 to 360 degrees; any finite phase is accepted.
A negative velocity at the epoch is represented by phase, not a negative amplitude.
Elevation and velocity have independent phases. Rotate geographic east/north currents into
model directions before writing the file. Neither velocity uses an outward-normal sign.

``tide_epoch`` is a scalar time with units such as ``days since 2020-01-01 00:00:00``.
Its value is in days, as for REMORA's other NetCDF times; zero selects the date in the units.
The epoch is placed on the model clock using ``remora.time_ref`` and does not reset at restart.
If supplied, ``calendar`` must agree with the model calendar. The reference date uses UTC.

File layout
~~~~~~~~~~~

Use 64-bit-offset NetCDF classic (CDF-2). The dimensions are ``tide`` for constituents,
``eta_rho = n_cell_y + 2`` and ``xi_rho = n_cell_x + 2`` on level 0, and ``name_strlen``
for fixed-length character labels. There is no time-record dimension.

.. list-table:: Tidal variables
   :header-rows: 1
   :widths: 55 45

   * - Variable
     - Dimensions
   * - ``tide_name``
     - ``(tide, name_strlen)``, characters
   * - ``tide_period``
     - ``(tide)``
   * - ``tide_epoch``
     - Scalar
   * - ``tide_zeta_amp_west``, ``tide_zeta_phase_west``; similarly ``east``
     - ``(tide, eta_rho)``
   * - ``tide_ubar_amp_west``, ``tide_ubar_phase_west``; similarly ``east``
     - ``(tide, eta_rho)``
   * - ``tide_zeta_amp_south``, ``tide_zeta_phase_south``; similarly ``north``
     - ``(tide, xi_rho)``
   * - ``tide_vbar_amp_south``, ``tide_vbar_phase_south``; similarly ``north``
     - ``(tide, xi_rho)``

Only sides whose normal barotropic velocity uses Flather require amplitude/phase pairs.
Each such side requires both elevation and normal velocity. Every numeric variable must have
one of the units above exactly, except ``tide_epoch`` which uses date-bearing time units. All pairs share
the constituent order of ``tide_period``. Values must be finite and contain no missing-data
entries; write zero amplitudes where a constituent makes no contribution.

Arrays run in increasing model-coordinate order. West/east entries correspond to level-0
indices ``j = -1, 0, ..., n_cell_y``; south/north entries correspond to
``i = -1, 0, ..., n_cell_x``. The first and last entries are boundary-strip endpoints.
Normal velocities are at their native boundary faces. Elevation is the reference at that
same face; when preparing it from full-grid rho-point harmonics, average the two adjacent
harmonics onto the face, as ROMS does. Rotate and average sine/cosine coefficients rather than
phase angles across 0/360 degrees. REMORA uses those coefficients to interpolate to refined levels.

``Tools/make_boundary_tides.py`` provides ``write_boundary_tides`` for prepared boundary
amplitude/phase arrays. Running it writes a synthetic example::

    python3 Tools/make_boundary_tides.py --nx 128 --ny 64 \
        --epoch "2020-01-01 00:00:00" --output boundary_tides.nc

The script requires NumPy and netCDF4. Replace its demonstration harmonics with values for
the modeled region before using the file for a physical simulation.


.. _sec:bc-per-tracer:

Boundary conditions for individual tracers
------------------------------------------

Every cell-centered tracer can be configured on its own, using the same options and the same
file-driven data path as temperature and salinity. In per-variable mode, address a tracer by its
own name rather than by the shared ``scalar`` keyword:

::

    #                              West      South    East      North
    remora.bc.NO3.type          =  periodic  clamped  periodic  clamped
    remora.bc.oxygen.type       =  periodic  outflow  periodic  outflow

The tracer names are the same ones used for plotfile output: ``temp`` and ``salt``, then the
dye tracers ``tracer``, ``tracer_1``, ... when ``remora.nscalar`` is positive, then the biology
tracer names (``NO3``, ``NH4``, ``chlorophyll``, ...) when :ref:`remora.biology_model<sec:Biology>`
selects a model.

``remora.bc.scalar.type`` still works and still applies to every tracer beyond salt. A per-tracer
entry takes precedence where one is given, so the two can be mixed: set ``scalar`` for the common
case and override only the tracers that differ.

Tracers configured with ``clamped`` or ``orlanski_rad_nudg`` read their boundary values from the
NetCDF boundary file, following the same naming convention as temperature and salinity: a tracer
named ``NO3`` reads ``NO3_west``, ``NO3_east``, ``NO3_south``, and ``NO3_north``. Only the sides
that actually need data are read, so the file need not contain anything for a tracer left on a
local condition such as ``outflow`` or ``slipwall``. If a tracer is given a file-driven condition
and the file does not carry the matching variables, the run stops at setup and names them. The
time axis for a tracer's boundary data can be named with ``remora.bdy_<name>_time_varname`` (for
example ``remora.bdy_NO3_time_varname``), defaulting to ``remora.bdy_time_varname`` as for the
other variables.

This holds in per-side mode as well: a side keyword covers every variable at once, tracers
included, so a per-side ``clamped`` or ``orlanski_rad_nudg`` drives every tracer the run carries
from file, not temperature and salinity alone. A ROMS boundary file that carries temperature and
salinity but nothing for the extra tracers will stop the run at setup, naming the variables it
lacks. Either give the file the missing fields, drop the extra tracers (``remora.nscalar = 0``),
or switch to per-variable mode and state each tracer's condition outright -- ``outflow`` for a
tracer with no boundary data of its own. ``remora.nscalar`` defaults to 0, so this only comes
up for a run that asks for dye or turns on a biology model.

Nudging for ``orlanski_rad_nudg`` uses the tracer timescale ``remora.tnudg`` for every tracer,
matching the ROMS default of a single ``Tnudg`` shared across tracers. A tracer picks up a
spatially varying coefficient instead when climatology nudging is turned on for it and
``remora.nc_clim_coeff_file`` carries a ``{var}_NudgeCoef`` field for that tracer; otherwise it
keeps the constant derived from ``remora.tnudg``. See :ref:`sec:clim-per-tracer`.

.. _sec:bc-options:

Boundary condition options
--------------------------

Boundary types for per-side or per-variable specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scroll table to the right if not all visible. Rightmost column is sea surface height.

+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| Type        | ROMS      | Normal             | Tangential          | Normal             | Tangential          | T, S, etc.         | sea surface        |
|             |           |                    |                     |                    |                     |                    |                    |
|             | name      | vel (3D)           | vel (3D)            | vel (2D)           | vel (2D)            |                    | height             |
+=============+===========+====================+=====================+====================+=====================+====================+====================+
| periodic    | Per       | periodic           | periodic            | periodic           | periodic            | periodic           | periodic           |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| inflow      | Cla       | ext_dir            | ext_dir             | ext_dir            | ext_dir             | ext_dir            | ext_dir            |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| outflow     | Gra       | foextrap           | foextrap            | foextrap           | foextrap            | foextrap           | foextrap           |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| slipwall    | Clo       | ext_dir            | foextrap            | ext_dir            | foextrap            | ext_dir/           | ext_dir/           |
|             |           |                    |                     |                    |                     |                    |                    |
|             |           | (set to 0)         |                     | (set to 0)         |                     | foextrap           | foextrap           |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| noslipwall  | N/A       | ext_dir            | ext_dir             | ext_dir            | ext_dir             | ext_dir/           | ext_dir/           |
|             |           |                    |                     |                    |                     |                    |                    |
|             |           | (set to 0)         | (set to 0)          | (set to 0)         | (set to 0)          | foextrap           | foextrap           |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| symmetry    | N/A       | reflect_odd        | reflect_even        | reflect_odd        | reflect_even        | reflect_even       | reflect_even       |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| clamped*    | Cla       | clamped            | clamped             | clamped            | clamped             | clamped            | clamped            |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+

Boundary types for per-variable specification ONLY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scroll table to the right if not all visible. Rightmost column is sea surface height.

+--------------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| Type               | ROMS      | Normal             | Tangential          | Normal             | Tangential          | T, S, etc.         | sea surface        |
|                    |           |                    |                     |                    |                     |                    |                    |
|                    | name      | vel (3D)           | vel (3D)            | vel (2D)           | vel (2D)            |                    | height             |
+====================+===========+====================+=====================+====================+=====================+====================+====================+
| chapman*           | Che       | N/A                | N/A                 | N/A                | N/A                 | N/A                | chapman            |
+--------------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| flather*           | Fla       | N/A                | N/A                 | flather            | flather             | N/A                | N/A                |
+--------------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| orlanski_rad       | Rad       | orlanski           | orlanski            |  orlanski          | orlanski            | orlanski           | orlanski           |
+--------------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| orlanski_rad_nudg* | RadNud    | orlanski           | orlanski            | orlanski           | orlanski            | orlanski           | orlanski           |
|                    |           |                    |                     |                    |                     |                    |                    |
|                    |           | w/nudging          | w/nudging           | w/nudging          | w/nudging           | w/nudging          | w/nudging          |
+--------------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+

The asterisks (*) indicate conditions that require the specification of a :ref:`boundary file <icbc-parameters>`.

Flather relaxes the normal barotropic velocity toward the boundary file's elevation on that side, so a
side whose normal ``ubar`` or ``vbar`` is ``flather`` must also give ``zeta`` a condition that reads the
file there (``clamped``, ``chapman`` or ``orlanski_rad_nudg``). The run stops at setup otherwise.

Here ``ext_dir``, ``foextrap``, and ``reflect_even`` refer to AMReX keywords.   The ``ext_dir`` type
refers to an "external Dirichlet" boundary, which means the values must be specified by the user, unless
marked as *set to 0* in the table above.
The ``foextrap`` type refers to "first order extrapolation" which sets all the ghost values to the
same value in the last valid cell/face.  (AMReX also has a ``hoextrap``, or "higher order extrapolation"
option, which does a linear extrapolation from the two nearest valid values.)

As an example,

::

    remora.bc.xlo.type                =   "Inflow"
    remora.bc.xlo.velocity            =   1. 0.9  0.
    remora.bc.xlo.temp                =   15.
    remora.bc.xlo.salt                =   35.
    remora.bc.xlo.tracer              =   2.

sets the boundary condition type at the low x face to be an inflow with xlo.type = “Inflow”.
The velocity components take their value from ``velocity``; each tracer takes one keyed by its
own name. There is no keyword standing in for "the remaining tracers": a side covers all of
them, so every tracer the run carries -- ``temp``, ``salt``, each dye, each biology tracer --
needs its own entry on an inflow face, and a value always says which tracer it is for.
``ubar``, ``vbar``, ``zeta`` and ``tke`` have no input of their own and are given zero.

Because ``ubar`` and ``vbar`` are among those, an inflow face carries no barotropic normal
velocity, and the barotropic correction holds the net volume flux through the face near zero
whatever ``velocity`` says. Tracers still enter such a face by horizontal diffusion, and the
inflow values above are what they diffuse from, but ``inflow`` cannot presently be used to drive
a prescribed advective transport through the boundary.

In per-variable mode the prefix already names the variable, so ``value`` sits under it --
``remora.bc.temp.value``, ``remora.bc.NO3.value``, ``remora.bc.u.velocity`` -- and applies to
whichever of that variable's sides are inflow. A tracer that takes its condition from the shared
``remora.bc.scalar.type`` reads ``remora.bc.scalar.value`` along with it; give a tracer its own
``.type`` to give it its own value.

``value`` is read only under a variable prefix, and there is no ``scalar`` keyword for a tracer
inflow value at all.

An inflow face missing a tracer value stops the run at setup, naming the input it wants, as a
missing ``velocity`` does. There is no default value.

We note that ``noslipwall`` allows for non-zero tangential velocities to be specified, such as

::

    remora.is_periodic = 1 0 0

    remora.bc.ylo.type = "NoSlipWall"
    remora.bc.yhi.type = "NoSlipWall"

    remora.bc.ylo.velocity    = 0.0 0.0 0.0
    remora.bc.yhi.velocity    = 2.0 0.0 0.0


It is important to note that external Dirichlet boundary data should be specified
as the value on the face of the cell bounding the domain, even for cell-centered
state data.

.. _sec:nudging-options:

Nudging options
---------------

When using ``orlanski_rad_nudg``, the nudging strength is specified by input parameters. When climatology nudging is turned on, timescales will instead be read from the file specified in ``remora.nc_clim_coeff_file``. ``remora.obcfac`` will still be used.

+-------------------+-------------------------+-------------------+---------------+
| Parameter         | Definition              | Acceptable Values | Default       |
+===================+=========================+===================+===============+
| **remora.tnudg**  | Nudging timescale for   | Positive real     | 0.0           |
|                   |                         |                   |               |
|                   | tracers in days         |                   |               |
+-------------------+-------------------------+-------------------+---------------+
| **remora.m3nudg** | Nudging timescale for   | Positive real     | 0.0           |
|                   |                         |                   |               |
|                   | 3D momentum in days     |                   |               |
+-------------------+-------------------------+-------------------+---------------+
| **remora.obcfac** | Ratio between inflow    | Positive real     | 0.0           |
|                   |                         |                   |               |
|                   | and outflow             |                   |               |
|                   |                         |                   |               |
|                   | boundary conditions     |                   |               |
+-------------------+-------------------------+-------------------+---------------+


