.. _sec:Biology:

*******
Biology
*******

REMORA can carry an ecosystem model's tracers in its conserved state and apply
the model's source and sink terms each three-dimensional step, after the GLS
corrector. Two models are available, described on their own pages:

* :ref:`sec:Fennel`, the ROMS Fennel nitrogen-cycle model with optional
  phosphate, carbon, oxygen and river pools, selected at run time.
* :ref:`sec:ECB`, the ROMS Estuarine Carbon Biogeochemistry model, with one
  fixed option set.

This page covers what the two share: selecting a model, the state layout,
initial conditions, output, boundary conditions and restarts.

Selecting a model
=================

.. code:: shell

   remora.biology_model = fennel     # or ecb; none (the default) or off disables biology

Model parameters use the model's own prefix, ``remora.fennel.*`` or
``remora.ecb.*``; the settings on this page have no model prefix.

State layout
============

The conserved state holds temperature, salinity, the passive dye tracers, then
the biology tracers:

.. code:: text

   temp salt [tracer tracer_1 ...] [biology tracers]

The number of biology tracers is set by the model and its options. It is
independent of ``remora.nscalar``, which counts dye tracers only and defaults
to 0, so a biology run carries no dye unless it asks for some.

.. warning::

   ``remora.nscalar`` counts dye *in addition to* the biology tracers. Setting
   it to the number of biology tracers allocates that many extra dye tracers.
   REMORA prints the full component list whenever a run carries both.

Every interface refers to biology tracers by name: initial conditions, plot
variables, boundary conditions and the model kernels. The order on each model's
page only sets the layout of plotfiles and checkpoints.

Initial conditions
==================

``remora.biology_ic_type`` selects the source of the biology initial
condition, independently of ``remora.ic_type``:

.. list-table::
   :header-rows: 1

   * - Value
     - Biology initial condition
   * - ``follow`` (default)
     - NetCDF when ``remora.ic_type = netcdf``, analytic otherwise
   * - ``constant``
     - Each tracer set to ``remora.biology_initial_<name>``
   * - ``analytic``
     - The problem's analytic biology profile
   * - ``netcdf``
     - Read from the NetCDF initial file; needs ``remora.ic_type = netcdf``

Separating the two lets the physical fields come from a NetCDF file that lacks
some biology tracers while biology starts from a constant or analytic state.
Temperature, salinity and dyes always come from the physical initialization.
With ``remora.hires_init_level``, biology is initialized on the
high-resolution level and averaged down with the physical fields.

Constant
--------

Each biology tracer starts at a uniform value set by name, in the style of the
constant forcing values such as ``remora.air_temperature``:

.. code:: shell

   remora.biology_ic_type     = constant
   remora.biology_initial_NO3 = 10.0
   remora.biology_initial_NH4 = 0.1
   # ... one line per tracer

There are no defaults: every tracer the configured model declares must be
given, and the error names any that are missing. Ghost cells are set too.
``Tests/test_files/Upwelling_ECB/Upwelling_ECB.i`` sets all 18 ECB tracers.

NetCDF
------

REMORA reads one variable per biology tracer, named as the tracer. A missing
variable is an error.

Analytic
--------

The problem supplies the profile through
``ProblemBase::init_analytic_biology``. The standard ``Problem`` in
``Exec/REMORA_Prob.cpp`` provides the Fennel profile from ROMS
``ana_biology.h`` for ``BioToy`` and ``Upwelling``, and aborts for any other
problem rather than starting from zeroed tracers; see :ref:`sec:Fennel`.

To add a profile, write ``Source/Prob/REMORA_InitAnalyticBiology_<Problem>.H``
and give it a branch in ``Problem::init_analytic_biology``. The hook receives
``const REMORABiology::Biology&`` and ``REMORABiology::BiologyInitContext&``.
A profile finds each field with
``biology.declaration().component_index(name, context.component_start)``, which
returns -1 for a tracer the run does not carry, and writes only the biology
components. NetCDF and constant initialization do not call the hook.

Output
======

The plot-variable group ``biology`` expands to all of the model's tracers,
excluding dye; each model also accepts its own name, ``fennel`` or ``ecb``, as
an alias. Tracers can also be requested by name:

.. code:: shell

   remora.plot_vars_3d = temp salt biology

Only tracers the run carries are available; REMORA warns about and omits any
other name. AMReX and NetCDF plotfiles use the same names; NetCDF output
requires a PnetCDF build (see :ref:`sec:Plotfiles`). Process diagnostics such
as primary production or air-sea fluxes are not output.

Boundary conditions
===================

Biology tracers offer the same boundary conditions as temperature and
salinity, including the file-driven ``clamped`` and ``orlanski_rad_nudg``.
With per-variable boundary conditions each tracer is addressed by name, for
example ``remora.bc.NO3.type``, and a tracer without its own entry falls back
to ``remora.bc.scalar.type``. A tracer clamped or nudged to file data reads
``<name>_west``, ``<name>_east``, ``<name>_south`` and ``<name>_north`` from the
boundary file, for the sides that need it. A per-side keyword applies its
condition to every tracer, so a per-side ``clamped`` run needs ``NO3_west`` and
the rest in the boundary file; per-variable mode lets a biology tracer keep a
local condition such as ``outflow`` while temperature and salinity are clamped.
See :ref:`sec:bc-per-tracer`.

Restarts
========

Restart with the same model and options as the run that wrote the checkpoint;
changing them changes the tracer layout. Every checkpoint ``Header`` ends with
a ``cons_names`` line listing the cell-centered components in order, and
restart aborts, printing both lists, when it differs from the run's, including
for a changed ``remora.nscalar``. Checkpoints written before this line existed
are checked by component count only. A separate ``BiologyHeader`` records the
model name, which restart also checks; Fennel and no-biology runs accept
checkpoints without it, ECB runs do not. Changing rate parameters across a
restart is allowed.

Adding a model
==============

REMORA holds one ``REMORABiology::Biology`` object, created by
``make_biology`` in ``Source/Biology/REMORA_Biology.cpp``. A model implements
``configure``, which reads its parameters; ``declaration``, which fixes its
tracer names and plot alias before the state is allocated; and ``advance``.
``advance`` receives the old concentrations, the new state, and views of the
geometry and forcing (``BiologyAdvanceContext`` in ``REMORA_Biology.H``), and
adds its concentration change times layer thickness to the new state,
preserving the increments already there. Dispatch is virtual on the host;
device kernels capture plain views and parameters. Adding a model takes its
implementation, its build entries and one line in ``make_biology``; the solver
never tests model names. Shared chemistry lives in
``REMORA_BiologyCarbonate.H`` (surface pCO2) and
``REMORA_BiologyGasExchange.H`` (O2 saturation, CO2 solubility, Schmidt
numbers and transfer velocity).
