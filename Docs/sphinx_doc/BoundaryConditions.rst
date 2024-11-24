
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
speed is specified by the function ``init_custom_smflux`` in ``prob.cpp``.


Boundary per side
-----------------

To set boundary conditions per domain side for all variables, set
``remora.boundary_per_variable = false``. This is the default behavior.

The information for each face is preceded by
``bc.xlo.type``, ``bc.xhi.type``, ``bc.ylo.type``, or ``bc.yhi.type``. Spelling of the type matters; capitalization does not. The
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
If periodic is selected, it must be used for both low and high faces in a direction. The ``geometry.is_periodic``
flag must match as described in :ref:`ProblemGeometry`<geometry-parameters>`.

For example, setting

::

    bc.xlo.type = "Inflow"
    bc.xhi.type = "Outflow"
    bc.ylo.type = "periodic" #optional
    bc.yhi.type = "periodic" #optional

    geometry.is_periodic = 0 1 0

would define a problem with inflow in the low-\ :math:`x` direction,
outflow in the high-\ :math:`x` direction, periodic in the :math:`y`-direction,
and slip wall on the low and high :math:`y`-faces, and
Note that no keyword is needed for a periodic boundary, here only the
specification in ``geometry.is_periodic`` is needed.

Boundary per variable
---------------------

To set different boundary conditions for different variables on each side of the domain, set ``remora.boundary_per_variable = true``.
If true, conditions must be set for all of the following variables.

- temperature: ``bc.temp.type``
- salinity: ``bc.salt.type``
- passive scalar: ``bc.scalar.type``
- 3D u-velocity: ``bc.u.type``
- 3D v-velocity: ``bc.v.type``
- 2D u-velocity: ``bc.ubar.type``
- 2D v-velocity: ``bc.vbar.type``
- sea surface height: ``bc.zeta.type``
- turbulent kinetic energy: ``bc.tke.type``

They must be set to a list of four conditions in the order West, South, East, North. The options are

- periodic
- inflow
- outflow
- slipwall
- noslipwall
- symmetry
- clamped
- chapman
- flather
- orlanski
- orlanski_nudg

The corresponding ROMS conditions are indicated in the tables below.

For example, setting

::

    #                 West      South    East      North
    bc.temp.type   =  periodic  clamped  periodic  clamped
    bc.salt.type   =  periodic  clamped  periodic  clamped
    bc.scalar.type =  periodic  clamped  periodic  clamped
    bc.u.type      =  periodic  clamped  periodic  clamped
    bc.v.type      =  periodic  clamped  periodic  clamped
    bc.w.type      =  periodic  clamped  periodic  clamped
    bc.ubar.type   =  periodic  flather  periodic  flather
    bc.vbar.type   =  periodic  flather  periodic  flather
    bc.zeta.type   =  periodic  chapman  periodic  chapman
    bc.tke.type    =  periodic  outflow  periodic  outflow

    geometry.is_periodic = 1 0 0

will define a problem that is periodic on the Western and Eastern sides. Temperature, salinity, passive scalar,
3D u-velocity, and 3D v-velocity will be clamped to values given in a NetCDF file specified by
:ref:`remora.nc_bdry_file_0`<icbc-parameters>`. The 2D momentum and zeta BCs are calculated from the
Chapman/Flather conditions, with nudging towards values given in the boundary NetCDF file.

.. _sec:bc-options:

Boundary condition options
--------------------------

Boundary types for per-side or per-variable specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| Type        | ROMS name | Normal vel (3D)    | Tangential vel (3D) | Normal vel (2D)    | Tangential vel (2D) | T, S, etc.         | sea surface height |
+=============+===========+====================+=====================+====================+=====================+====================+====================+
| periodic    | Per       | periodic           | periodic            | periodic           | periodic            | periodic           | periodic           |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| inflow      | Cla       | ext_dir            | ext_dir             | ext_dir            | ext_dir             | ext_dir            | ext_dir            |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| outflow     | Gra       | foextrap           | foextrap            | foextrap           | foextrap            | foextrap           | foextrap           |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| slipwall    | Clo       | ext_dir (set to 0) | foextrap            | ext_dir (set to 0) | foextrap            | ext_dir/foextrap   | ext_dir/foextrap   |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| noslipwall  | N/A       | ext_dir (set to 0) | ext_dir (set to 0)  | ext_dir (set to 0) | ext_dir (set to 0)  | ext_dir/foextrap   | ext_dir/foextrap   |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| symmetry    | N/A       | reflect_odd        | reflect_even        | reflect_odd        | reflect_even        | reflect_even       | reflect_even       |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| clamped*    | Cla       | clamped            | clamped             | clamped            | clamped             | clamped            | clamped            |
+-------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+

Boundary types for per-variable specification ONLY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+----------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| Type           | ROMS name | Normal vel (3D)    | Tangential vel (3D) | Normal vel (2D)    | Tangential vel (2D) | T, S, etc.         | sea surface height |
+================+===========+====================+=====================+====================+=====================+====================+====================+
| chapman*       | Che       | N/A                | N/A                 | N/A                | N/A                 | N/A                | chapman            |
+----------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| flather*       | Fla       | N/A                | N/A                 | flather            | flather             | N/A                | N/A                |
+----------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| orlanski       | Rad       | orlanski           | orlanski            | N/A                | N/A                 | orlanski           | N/A                |
+----------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+
| orlanski_nudg* | RadNud    | orlanski w/nudging | orlanski w/nudging  | N/A                | N/A                 | orlanski w/nudging | N/A                |
+----------------+-----------+--------------------+---------------------+--------------------+---------------------+--------------------+--------------------+

The asterisks (*) indicate conditions that require the specification of a :ref:`boundary file`<icbc-parameters>`.

Here ``ext_dir``, ``foextrap``, and ``reflect_even`` refer to AMReX keywords.   The ``ext_dir`` type
refers to an "external Dirichlet" boundary, which means the values must be specified by the user, unless
marked as *set to 0* in the table above.
The ``foextrap`` type refers to "first order extrapolation" which sets all the ghost values to the
same value in the last valid cell/face.  (AMReX also has a ``hoextrap``, or "higher order extrapolation"
option, which does a linear extrapolation from the two nearest valid values.)

As an example,

::

    bc.xlo.type                =   "Inflow"
    bc.xlo.velocity            =   1. 0.9  0.
    bc.xlo.temp                =   15.
    bc.xlo.scalar              =   2.

sets the boundary condition type at the low x face to be an inflow with xlo.type = “Inflow”.

We note that ``noslipwall`` allows for non-zero tangential velocities to be specified, such as

::

    geometry.is_periodic = 1 0 0

    bc.ylo.type = "NoSlipWall"
    bc.yhi.type = "NoSlipWall"

    bc.ylo.velocity    = 0.0 0.0 0.0
    bc.yhi.velocity    = 2.0 0.0 0.0


It is important to note that external Dirichlet boundary data should be specified
as the value on the face of the cell bounding the domain, even for cell-centered
state data.

.. _sec:nudging-options:

Nudging options
---------------

When using ``orlanski_nudg``, the nudging strength is specified by input parameters. Climatology nudging
has not yet been implemented.

+-------------------+-------------------------+-------------------+---------------+
| Parameter         | Definition              | Acceptable Values | Default       |
+===================+=========================+===================+===============+
| **remora.tnudg**  | Nudging timescale for   | Positive real     | 0.0           |
|                   | tracers in days         |                   |               |
+-------------------+-------------------------+-------------------+---------------+
| **remora.m3nudg** | Nudging timescale for   | Positive real     | 0.0           |
|                   | 3D momentum in days     |                   |               |
+-------------------+-------------------------+-------------------+---------------+
| **remora.obcfac** | Ratio between inflow    | Positive real     | 0.0           |
|                   | and outflow             |                   |               |
|                   | boundary conditions     |                   |               |
+-------------------+-------------------------+-------------------+---------------+


