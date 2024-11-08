
 .. role:: cpp(code)
    :language: c++

.. _sec:domainBCs:

Physical/Domain Boundary Conditions
===================================

There are two primary types of physical/domain boundary conditions: those which rely only on the
data in the valid regions, and those which rely on externally specified values.

REMORA allows users to specify types of boundary condition with keywords in the inputs file.
The option ``remora.boundary_per_variable`` controls whether conditions are specified on a
per-side or per-variable basis. The options are:

- inflow

- outflow

- slipwall

- noslipwall

- symmetry




Boundary per side
-----------------

To set boundary conditions per domain side for all variables, set
``remora.boundary_per_variable = false``. This is the default behavior.

The information for each face is preceded by
``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, or ``zhi``.

+----------------------------------+-----------------+-------------------+-------------+
| Parameter                        | Definition      | Acceptable        | Default     |
|                                  |                 | Values            |             |
+==================================+=================+===================+=============+
| **remora.boundary_per_variable** | physical        | [Real Real -Real] | must be set |
|                                  | location of low |                   |             |
|                                  | corner of the   |                   |             |
|                                  | domain          |                   |             |
+----------------------------------+-----------------+-------------------+-------------+
| **geometry.prob_hi**             | physical        | [Real Real 0]     | must be set |
|                                  | location of     |                   |             |
|                                  | high corner of  |                   |             |
|                                  | the domain      |                   |             |
+----------------------------------+-----------------+-------------------+-------------+
| **geometry.is_periodic**         | is the domain   | 0 if false, 1     | 0 0 0       |
|                                  | periodic in     | if true.          |             |
|                                  | this direction  | Z-component must  |             |
|                                  |                 | be zero           |             |
+----------------------------------+-----------------+-------------------+-------------+

The information for each face is preceded by
``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, or ``zhi``.

Currently available type of boundary conditions are
``inflow``, ``outflow``, ``slipwall``, ``noslipwall``, or ``symmetry``
(Spelling of the type matters; capitalization does not.)

For example, setting

::

    xlo.type = "Inflow"
    xhi.type = "Outflow"
    zlo.type = "SlipWall"
    zhi.type = "SlipWall"

    geometry.is_periodic = 0 1 0

would define a problem with inflow in the low-\ :math:`x` direction,
outflow in the high-\ :math:`x` direction, periodic in the :math:`y`-direction,
and slip wall on the low and high :math:`y`-faces, and
Note that no keyword is needed for a periodic boundary, here only the
specification in ``geometry.is_periodic`` is needed.

Each of these types of physical boundary condition has a mapping to a mathematical boundary condition
for each type; this is summarized in the table below.

.. _sec:dirichlet:

Dirichlet Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

REMORA provides the ability to specify constant Dirichlet BCs in the inputs file. We use the following options
preceded by
``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, and ``zhi``:

+---------------+--------------------+---------------------+--------------------+---------------------+------------------+--------------------+
| Type          | Normal vel (3D)    | Tangential vel (3D) | Normal vel (2D)    | Tangential vel (2D) | T, S, etc.       | sea surface height |
+===============+====================+=====================+====================+=====================+==================+====================+
| inflow        | ext_dir            | ext_dir             | ext_dir            | ext_dir             | ext_dir          | ext_dir            |
+---------------+--------------------+---------------------+--------------------+---------------------+------------------+--------------------+
| outflow       | foextrap           | foextrap            | foextrap           | foextrap            | foextrap         | foextrap           |
+---------------+--------------------+---------------------+--------------------+---------------------+------------------+--------------------+
| slipwall      | ext_dir            | foextrap            | ext_dir            | foextrap            | ext_dir/foextrap | ext_dir/foextrap   |
+---------------+--------------------+---------------------+--------------------+---------------------+------------------+--------------------+
| noslipwall    | ext_dir            | ext_dir             | ext_dir            | ext_dir             | ext_dir/foextrap | ext_dir/foextrap   |
+---------------+--------------------+---------------------+--------------------+---------------------+------------------+--------------------+
| symmetry      | reflect_odd        | reflect_even        | reflect_odd        | reflect_even        | reflect_even     | reflect_even       |
+---------------+--------------------+---------------------+--------------------+---------------------+------------------+--------------------+
| clamped       | clamped            | clamped             | clamped            | clamped             | clamped          | clamped            |
+---------------+--------------------+---------------------+--------------------+---------------------+------------------+--------------------+
| chapman       | N/A                | N/A                 | N/A                | N/A                 | N/A              | chapman            |
+---------------+--------------------+---------------------+--------------------+---------------------+------------------+--------------------+
| flather       | N/A                | N/A                 | flather            | flather             | N/A              | N/A                |
+---------------+--------------------+---------------------+--------------------+---------------------+------------------+--------------------+
| orlanski      | orlanski           | orlanski            | N/A                | N/A                 | orlanski         | N/A                |
+---------------+--------------------+---------------------+--------------------+---------------------+------------------+--------------------+
| orlanski_nudg | orlanski           | orlanski            | N/A                | N/A                 | orlanski         | N/A                |
+---------------+--------------------+---------------------+--------------------+---------------------+------------------+--------------------+

Periodic boundary conditions are specified by the ``geometry.is_periodic`` flag as described in :ref:`Problem Geometry`<geometry-parameters>`.

Boundaries in the z-direction always correspond to the sea floor and surface, so boundary type selection for ``zlo`` and ``zhi``  will not
affect REMORA behavior. Bottom drag is specified with ``remora.rdrag`` (see :ref:`Physics Parameters<list-of-parameters-15>`), and surface wind speed is specified by the function ``init_custom_smflux`` in ``prob.cpp``.

Here ``ext_dir``, ``foextrap``, and ``reflect_even`` refer to AMReX keywords.   The ``ext_dir`` type
refers to an "external Dirichlet" boundary, which means the values must be specified by the user.
The ``foextrap`` type refers to "first order extrapolation" which sets all the ghost values to the
same value in the last valid cell/face.  (AMReX also has a ``hoextrap``, or "higher order extrapolation"
option, which does a linear extrapolation from the two nearest valid values.)

Note that ``outflow`` must be selected for ``xlo``, ``xhi``, ``ylo``, and ``yhi`` when reading boundary data
from a NetCDF file.

As an example,

::

    xlo.type                =   "Inflow"
    xlo.velocity            =   1. 0.9  0.
    xlo.temp                =   15.
    xlo.scalar              =   2.

sets the boundary condition type at the low x face to be an inflow with xlo.type = “Inflow”.

We note that ``noslipwall`` allows for non-zero tangential velocities to be specified, such as

::

    geometry.is_periodic = 1 0 0

    ylo.type = "NoSlipWall"
    yhi.type = "NoSlipWall"

    ylo.velocity    = 0.0 0.0 0.0
    yhi.velocity    = 2.0 0.0 0.0


It is important to note that external Dirichlet boundary data should be specified
as the value on the face of the cell bounding the domain, even for cell-centered
state data.

