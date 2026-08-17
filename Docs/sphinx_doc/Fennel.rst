.. role:: cpp(code)
  :language: c++

.. _sec:Fennel:

**************
Fennel Biology
**************

The Fennel biology module adds biological tracers to REMORA's conserved scalar
state and applies biological source and sink terms during the three-dimensional
time step. The biology update is called after the GLS corrector when a biology
model is enabled.

Turning Fennel On
=================

Enable the module in the inputs file with:

.. code:: shell

   remora.biology_model = fennel

The values ``none`` and ``off`` disable biology. The biology tracer count is set by
the active Fennel tracer set and is independent of ``remora.nscalar``, which counts
passive (dye) scalars only. A run may carry both: the state is laid out as
temperature, salinity, the passive scalars, then the biology tracers. With biology
active ``remora.nscalar`` defaults to 0, so a Fennel run carries no dye unless it
asks for some.

.. warning::

   ``remora.nscalar`` previously had to equal the biology tracer count and was a
   consistency check. It now counts dye scalars *in addition to* biology, so an
   input file carrying ``remora.nscalar = 11`` alongside Fennel will now allocate
   eleven dye tracers on top of the biology ones. Remove the setting to get the
   previous state layout. REMORA prints the full component list whenever a run
   carries both dye and biology.

A minimal nitrogen-only configuration is:

.. code:: shell

   remora.biology_model = fennel
   remora.fennel.po4 = false
   remora.fennel.carbon = false
   remora.fennel.oxygen = false
   remora.fennel.odu = false
   remora.fennel.denitrification = false

Runtime Options
===============

All Fennel options use the ``remora.fennel`` prefix.

+------------------------------------+---------+-----------------------------------------------+
| Option                             | Default | Effect                                        |
+====================================+=========+===============================================+
| ``remora.fennel.po4``              | false   | Adds phosphate and phosphate limitation and   |
|                                    |         | recycling terms.                              |
+------------------------------------+---------+-----------------------------------------------+
| ``remora.fennel.carbon``           | false   | Adds carbon detritus, TIC, alkalinity, carbon |
|                                    |         | cycling, and surface CO2 exchange.            |
+------------------------------------+---------+-----------------------------------------------+
| ``remora.fennel.oxygen``           | false   | Adds dissolved oxygen, oxygen source/sink     |
|                                    |         | terms, and surface O2 exchange.               |
+------------------------------------+---------+-----------------------------------------------+
| ``remora.fennel.odu``              | false   | Adds the optional ``ODU`` tracer. The current |
|                                    |         | port carries this tracer without source/sink  |
|                                    |         | terms.                                        |
+------------------------------------+---------+-----------------------------------------------+
| ``remora.fennel.denitrification``  | false   | Changes bottom nitrogen remineralization when |
|                                    |         | ``bio_sediment`` is active. It does not add a |
|                                    |         | tracer.                                       |
+------------------------------------+---------+-----------------------------------------------+
| ``remora.fennel.bio_sediment``     | true    | Returns sinking particulate material at the   |
|                                    |         | bottom.                                       |
+------------------------------------+---------+-----------------------------------------------+

The core nitrogen tracers are always present, in this order:

1. ``NO3``
2. ``NH4``
3. ``chlorophyll``
4. ``phytoplankton``
5. ``zooplankton``
6. ``LdetritusN``
7. ``SdetritusN``

When ``po4`` is true, ``PO4`` is appended after ``SdetritusN``.
When ``carbon`` is true, the following tracers are appended after optional
``PO4``: ``LdetritusC``, ``SdetritusC``, ``TIC``, and ``alkalinity``.
When ``oxygen`` is true, ``oxygen`` is appended after the carbon tracers, if
present. When ``odu`` is true, ``ODU`` is appended last.

The tracer count is additive:

.. code:: text

   nbio  = 7 + (po4 ? 1 : 0) + (carbon ? 4 : 0) + (oxygen ? 1 : 0) + (odu ? 1 : 0)
   ncons = Tracer_comp + remora.nscalar + nbio

For example, nitrogen-only Fennel uses 7 biology tracers, PO4-only Fennel uses
8, carbon+oxygen+ODU uses 13, and all four optional tracer groups use 14.

The oxygen and carbon options use surface gas-exchange terms and therefore
require surface forcing fields for the transfer velocity: 10-m winds
(``Uwind``/``Vwind``) when ``remora.bulk_fluxes`` is true, and surface stress
(``sustr``/``svstr``) otherwise. The base nitrogen-only configuration does not
require those fields for biology.

Biology Parameters
==================

The following parameters can be overridden in the inputs file as
``remora.fennel.<name> = <value>``.

+-----------------+-------------+---------------------------------------------+
| Parameter       | Default     | Meaning                                     |
+=================+=============+=============================================+
| ``BioIter``     | 1           | Iterations for the nonlinear biology update |
+-----------------+-------------+---------------------------------------------+
| ``AttSW``       | 0.04        | Seawater light attenuation [1/m]            |
+-----------------+-------------+---------------------------------------------+
| ``AttChl``      | 0.02486     | Chlorophyll light attenuation               |
+-----------------+-------------+---------------------------------------------+
| ``PARfrac``     | 0.43        | Shortwave fraction available for biology    |
+-----------------+-------------+---------------------------------------------+
| ``Vp0``         | 1.0         | Phytoplankton growth tuning parameter       |
+-----------------+-------------+---------------------------------------------+
| ``I_thNH4``     | 0.0095      | Radiation threshold for nitrification       |
+-----------------+-------------+---------------------------------------------+
| ``D_p5NH4``     | 0.1         | Half-saturation radiation for nitrification |
+-----------------+-------------+---------------------------------------------+
| ``NitriR``      | 0.05        | Nitrification rate [day-1]                  |
+-----------------+-------------+---------------------------------------------+
| ``K_NO3``       | 2.0         | NO3 uptake inverse half-saturation          |
+-----------------+-------------+---------------------------------------------+
| ``K_NH4``       | 2.0         | NH4 uptake inverse half-saturation          |
+-----------------+-------------+---------------------------------------------+
| ``K_PO4``       | 32.0        | PO4 uptake inverse half-saturation          |
+-----------------+-------------+---------------------------------------------+
| ``K_Phy``       | 1.0         | Zooplankton ingestion half-saturation       |
+-----------------+-------------+---------------------------------------------+
| ``Chl2C_m``     | 0.0535      | Maximum chlorophyll-to-carbon ratio         |
+-----------------+-------------+---------------------------------------------+
| ``ChlMin``      | 0.01        | Chlorophyll minimum threshold               |
+-----------------+-------------+---------------------------------------------+
| ``PhyCN``       | 6.625       | Phytoplankton carbon-to-nitrogen ratio      |
+-----------------+-------------+---------------------------------------------+
| ``R_P2N``       | 0.0625      | Phytoplankton phosphate-to-nitrogen ratio   |
+-----------------+-------------+---------------------------------------------+
| ``PhyIP``       | 1.5         | Phytoplankton NH4 inhibition parameter      |
+-----------------+-------------+---------------------------------------------+
| ``PhyIS``       | 0.025       | Phytoplankton P-I curve initial slope       |
+-----------------+-------------+---------------------------------------------+
| ``PhyMin``      | 0.01        | Phytoplankton minimum threshold             |
+-----------------+-------------+---------------------------------------------+
| ``PhyMR``       | 0.15        | Phytoplankton mortality rate [day-1]        |
+-----------------+-------------+---------------------------------------------+
| ``ZooAE_N``     | 0.75        | Zooplankton nitrogen assimilation fraction  |
+-----------------+-------------+---------------------------------------------+
| ``ZooCN``       | 6.625       | Zooplankton carbon-to-nitrogen ratio        |
+-----------------+-------------+---------------------------------------------+
| ``ZooBM``       | 0.1         | Zooplankton basal metabolism [day-1]        |
+-----------------+-------------+---------------------------------------------+
| ``ZooER``       | 0.1         | Zooplankton excretion rate [day-1]          |
+-----------------+-------------+---------------------------------------------+
| ``ZooGR``       | 0.6         | Zooplankton maximum growth rate [day-1]     |
+-----------------+-------------+---------------------------------------------+
| ``ZooMin``      | 0.01        | Zooplankton minimum threshold               |
+-----------------+-------------+---------------------------------------------+
| ``ZooMR``       | 0.025       | Zooplankton mortality rate [day-1]          |
+-----------------+-------------+---------------------------------------------+
| ``LDeRRN``      | 0.01        | Large nitrogen detritus remineralization    |
+-----------------+-------------+---------------------------------------------+
| ``LDeRRC``      | 0.01        | Large carbon detritus remineralization      |
+-----------------+-------------+---------------------------------------------+
| ``CoagR``       | 0.005       | Coagulation rate [day-1]                    |
+-----------------+-------------+---------------------------------------------+
| ``SDeRRN``      | 0.1         | Small nitrogen detritus remineralization    |
+-----------------+-------------+---------------------------------------------+
| ``SDeRRC``      | 0.03        | Small carbon detritus remineralization      |
+-----------------+-------------+---------------------------------------------+
| ``wPhy``        | 0.1         | Phytoplankton sinking velocity [m/day]      |
+-----------------+-------------+---------------------------------------------+
| ``wLDet``       | 1.0         | Large detritus sinking velocity [m/day]     |
+-----------------+-------------+---------------------------------------------+
| ``wSDet``       | 0.1         | Small detritus sinking velocity [m/day]     |
+-----------------+-------------+---------------------------------------------+
| ``pCO2air``     | 370.0       | Atmospheric CO2 partial pressure [ppmv]     |
+-----------------+-------------+---------------------------------------------+

Initial Data and Restarts
=========================

The source of the biology initial condition is selected independently of
``remora.ic_type``:

+-----------------------------+-------------+-------------------------------------------------------+
| Option                      | Default     | Description                                           |
+=============================+=============+=======================================================+
| ``remora.biology_ic_type``  | ``follow``  | ``follow``: use NetCDF when ``remora.ic_type`` is     |
|                             |             | ``netcdf``, analytic otherwise. ``analytic``: always  |
|                             |             | use the problem's analytic biology profile.           |
|                             |             | ``netcdf``: always read from the NetCDF initial file. |
+-----------------------------+-------------+-------------------------------------------------------+

Note this option has no ``remora.fennel`` prefix; it is a biology-wide setting
rather than a Fennel parameter.

The two are decoupled because a NetCDF initial file supplies only the tracers it
happens to contain. Setting ``remora.biology_ic_type = analytic`` keeps physical
fields coming from NetCDF while biology comes from the problem's analytic
profile, which is useful when the initial file lacks a tracer that an enabled
option requires.

For NetCDF biology initialization, REMORA reads one variable per active tracer,
using the tracer names listed under Runtime Options above. A missing variable is
an error. For analytic biology initialization, the problem must
provide an analytic biology routine; problems that enable biology without one
abort during initialization rather than starting from uninitialized tracers.

``BioToy`` and ``Upwelling`` ship with such a routine, both using the profile from
ROMS ``ana_biology.h``: nitrate follows a cubic in a temperature-derived silicate
proxy, and the remaining tracers start uniform. To add one for another problem, write
a ``Source/Prob/REMORA_InitAnalyticBiology_<Problem>.H`` and give it a branch in
``Problem::init_analytic_biology`` in ``Exec/REMORA_Prob.cpp``. Because that profile
is a function of temperature alone, it carries over unchanged to any problem that
stratifies temperature.

When initial data is specified on a high-resolution level and averaged down
(``remora.hires_init_level``), biology is initialized on that level before the
average-down, using the same ``remora.biology_ic_type`` selection.

Checkpoints contain the active conserved state, including Fennel tracers. Restart
with the same Fennel component options used to create the checkpoint. Changing
``po4``, ``carbon``, ``oxygen``, or ``odu`` across a restart changes the scalar
layout and is not supported unless the checkpoint was produced with a compatible
layout.

Boundary conditions for the Fennel tracers offer the same options as temperature
and salinity, including the file-driven ``clamped`` and ``orlanski_rad_nudg``
conditions. When per-variable boundary conditions are enabled, each tracer can be
addressed by its own name -- ``remora.bc.NO3.type``, ``remora.bc.oxygen.type``,
and so on -- and a tracer with no entry of its own falls back to
``remora.bc.scalar.type``. A tracer clamped or nudged to file data reads
``<name>_west``, ``<name>_east``, ``<name>_south``, and ``<name>_north`` from the
boundary NetCDF file, following the ROMS naming convention, and only for the sides
that need it. Driving a Fennel tracer from file requires per-variable mode
(``remora.boundary_per_variable = true``); a per-side keyword applies the
file-driven condition to temperature and salinity only. See
:ref:`sec:bc-per-tracer` for details.

Validation Against ROMS
=======================

REMORA can optionally build the original ROMS Fennel kernel alongside the native
C++ implementation and select between them at run time, so the two can be
compared without rebuilding. This is validation infrastructure; production runs
need none of it, and both features are compiled out by default.

Build flags:

+-------------------------------------------+--------------------------------------------+
| GNUmake                                   | CMake                                      |
+===========================================+============================================+
| ``USE_FENNEL_FORT=TRUE``                  | ``-DREMORA_ENABLE_FENNEL_FORT=ON``         |
+-------------------------------------------+--------------------------------------------+
| ``USE_BIOLOGY_DIAG=TRUE``                 | ``-DREMORA_ENABLE_BIOLOGY_DIAG=ON``        |
+-------------------------------------------+--------------------------------------------+

The first builds the ROMS reference kernel; the second compiles in the
per-block diagnostic output used to compare the two. They are independent: the
native kernel carries the diagnostic instrumentation whether or not the
reference kernel is built. ``REMORA_ENABLE_FENNEL_FORT`` must be given on the
initial CMake configure line, because it determines whether the Fortran language
is enabled for the project.

Runtime options:

+---------------------------------------+---------+-----------------------------------------------+
| Option                                | Default | Description                                   |
+=======================================+=========+===============================================+
| ``remora.use_biology_cpp_answer``     | 1       | ``1`` uses the native C++ kernel, ``0`` the   |
|                                       |         | ROMS reference kernel. Requires the reference |
|                                       |         | kernel to be compiled in.                     |
+---------------------------------------+---------+-----------------------------------------------+
| ``remora.biology_debug``              | 0       | ``0`` off, ``1`` diagnostics for one target   |
|                                       |         | column, ``2`` for every column. Any value     |
|                                       |         | above zero requires the diagnostics to be     |
|                                       |         | compiled in.                                  |
+---------------------------------------+---------+-----------------------------------------------+
| ``remora.biology_debug_i``            | 2       | Target column index for ``biology_debug = 1`` |
+---------------------------------------+---------+-----------------------------------------------+
| ``remora.biology_debug_j``            | 2       | Target column index for ``biology_debug = 1`` |
+---------------------------------------+---------+-----------------------------------------------+

Requesting something that was not compiled in aborts with a message naming the
flag to rebuild with. Neither option degrades silently: selecting the reference
kernel without it present does not fall back to the native kernel, and enabling
diagnostics without them present does not produce an empty log.

Because ROMS selects its biology options at compile time while REMORA selects
them at run time, the reference kernel is built for one option set at a time.
That set is pinned in ``Source/Biology/Fortran/remora_bio_cppdefs.h`` and is
checked against the runtime options on every call, so a mismatched build stops
rather than producing an invalid comparison.

Diagnostic output volume grows as tags times levels times active tracers, so
``biology_debug = 2`` is only practical on very small grids.

See ``Source/Biology/Fortran/README.md`` for the full build-and-run procedure
and ``Source/Biology/Fortran/validation_lanes.md`` for current validation
coverage.

Output
======

Fennel tracers are output by requesting ``fennel`` in
``remora.plot_vars_3d``, which expands to all active biology tracers. Individual
active tracer names may also be requested. Do not use ``scalar`` for Fennel
output; ``scalar`` is the passive tracer name used when biology is disabled.

Nitrogen-only output example:

.. code:: shell

   remora.plot_vars_3d = temp salt fennel

PO4 output example:

.. code:: shell

   remora.fennel.po4 = true
   remora.plot_vars_3d = temp salt NO3 NH4 chlorophyll phytoplankton zooplankton LdetritusN SdetritusN PO4

Carbon and oxygen output example:

.. code:: shell

   remora.fennel.carbon = true
   remora.fennel.oxygen = true
   remora.plot_vars_3d = temp salt NO3 NH4 chlorophyll phytoplankton zooplankton LdetritusN SdetritusN LdetritusC SdetritusC TIC alkalinity oxygen

Only active tracer names are available. If an inactive or unknown name is listed
in ``remora.plot_vars_3d``, REMORA warns and omits it. Native AMReX plotfiles
and NetCDF plotfiles use the same active tracer names. NetCDF output requires
building REMORA with PnetCDF enabled; see :ref:`sec:Plotfiles`.

Biology diagnostic variables such as primary productivity, nitrification flux,
air-sea CO2 flux, air-sea O2 flux, and denitrification flux are not currently
implemented as separate output variables.
