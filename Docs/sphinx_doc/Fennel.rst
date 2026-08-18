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
temperature, salinity, the passive scalars, then the biology tracers.
``remora.nscalar`` defaults to 0, so a Fennel run carries no dye unless it asks for
some.

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
| ``remora.fennel.river_don``        | false   | Adds a river dissolved organic nitrogen pool  |
|                                    |         | that remineralizes to ``NH4``, plus a carbon  |
|                                    |         | counterpart when ``carbon`` is on. Neither    |
|                                    |         | sinks.                                        |
+------------------------------------+---------+-----------------------------------------------+
| ``remora.fennel.talk_nonconserv``  | false   | Makes alkalinity prognostic instead of a      |
|                                    |         | function of salinity. Requires ``carbon``.    |
+------------------------------------+---------+-----------------------------------------------+
| ``remora.fennel.pco2air_type``     |constant | Source of the atmospheric CO2 partial         |
|                                    |         | pressure: ``constant``, ``data``, or          |
|                                    |         | ``secular``. See below.                       |
+------------------------------------+---------+-----------------------------------------------+
| ``remora.fennel.co2_schmidt``      | ``w92`` | Schmidt-number relation for the CO2 transfer  |
|                                    |         | velocity: ``wanninkhof1992`` (alias ``w92``)  |
|                                    |         | or ``wanninkhof2014`` (alias ``rw14``).       |
+------------------------------------+---------+-----------------------------------------------+
| ``remora.fennel.oxygen_schmidt``   | ``w92`` | Same for O2, with the additional choice       |
|                                    |         | ``ocmip`` (Keeling et al. 1998).              |
+------------------------------------+---------+-----------------------------------------------+

The core nitrogen tracers are always present, in this order:

1. ``NO3``
2. ``NH4``
3. ``chlorophyll``
4. ``phytoplankton``
5. ``zooplankton``
6. ``LdetritusN``
7. ``SdetritusN``

When ``river_don`` is true, ``RdetritusN`` is appended after ``SdetritusN``.
When ``po4`` is true, ``PO4`` follows it.
When ``carbon`` is true, the following tracers are appended after optional
``PO4``: ``LdetritusC``, ``SdetritusC``, ``TIC``, and ``alkalinity``, then
``RdetritusC`` if ``river_don`` is also true.
When ``oxygen`` is true, ``oxygen`` is appended after the carbon tracers, if
present. When ``odu`` is true, ``ODU`` is appended last. This is the order ROMS
uses in ``initialize_biology``, so a NetCDF file written for either code reads
into the same tracers.

The tracer count is additive:

.. code:: text

   nbio  = 7 + (po4 ? 1 : 0) + (carbon ? 4 : 0) + (oxygen ? 1 : 0) + (odu ? 1 : 0)
           + (river_don ? 1 : 0) + (river_don && carbon ? 1 : 0)
   ncons = Tracer_comp + remora.nscalar + nbio

For example, nitrogen-only Fennel uses 7 biology tracers, PO4-only Fennel uses
8, carbon+oxygen+ODU uses 13, and every optional tracer group at once uses 16.

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
| ``RDeRRN``      | 0.03        | River nitrogen detritus remineralization    |
|                 |             | [day-1]. Used only with ``river_don``.      |
+-----------------+-------------+---------------------------------------------+
| ``RDeRRC``      | 0.03        | River carbon detritus remineralization      |
|                 |             | [day-1]. Used only with ``river_don``.      |
+-----------------+-------------+---------------------------------------------+
| ``wPhy``        | 0.1         | Phytoplankton sinking velocity [m/day]      |
+-----------------+-------------+---------------------------------------------+
| ``wLDet``       | 1.0         | Large detritus sinking velocity [m/day]     |
+-----------------+-------------+---------------------------------------------+
| ``wSDet``       | 0.1         | Small detritus sinking velocity [m/day]     |
+-----------------+-------------+---------------------------------------------+
| ``pCO2air``     | 370.0       | Atmospheric CO2 partial pressure [ppmv].    |
|                 |             | Used when ``pco2air_type = constant``.      |
+-----------------+-------------+---------------------------------------------+

Alkalinity
==========

By default alkalinity is diagnostic: at every biology iteration it is
overwritten with :math:`587.05 + 50.56\,S` following Brewer et al. (1986), so
it is a relabelled salinity field and carries no memory of the biology.

With ``remora.fennel.talk_nonconserv = true`` it becomes prognostic. The
overwrite is dropped and alkalinity instead responds to nitrate uptake
(:math:`+`), ammonium uptake (:math:`-`), nitrification (:math:`-2` per mole of
:math:`NH_4`), zooplankton metabolism and excretion (:math:`+`),
remineralization (:math:`+`), and the sediment return (:math:`+`, except under
``denitrification``). Alkalinity feeds the pCO2 solver, so this changes the
air-sea CO2 flux and hence TIC; it changes nothing else.

Atmospheric CO2
===============

``remora.fennel.pco2air_type`` selects where the atmospheric CO2 partial
pressure used by the surface gas exchange comes from:

``constant``
    ``remora.fennel.pCO2air``, unchanging. The default.

``data``
    The annual climatology of Laurent et al. (2017),
    :math:`380.464 + 9.321 \sin(2\pi\, \mathrm{yday}/365.25 + 1.068)`.

``secular``
    A linear trend plus three harmonics, referenced to 1951.

The two time-dependent forms need a calendar. REMORA dates the model clock with
``remora_caldate`` in ``Source/Utils/REMORA_DateClock.H``, a port of ROMS
``caldate``: ``remora.time_ref`` gives the reference date and selects the
calendar, and model time is elapsed time from there. See :ref:`calendar`.

Either of these places a run on 1 January 2020, the first by naming the date and
the second by offsetting the default 0001-01-01 epoch:

.. code:: shell

   remora.fennel.pco2air_type = secular
   remora.time_ref            = 20200101

.. code:: shell

   remora.fennel.pco2air_type = secular
   remora.start_time          = 63713433600   # 2020-01-01, with time_ref = 0

Getting this right matters for ``secular``. Its fit is anchored to 1951, and
extrapolating far back from there gives a negative partial pressure -- at the
0001-01-01 epoch it returns about -2648 ppmv, which would draw CO2 out of the
ocean for the whole run. REMORA aborts on a non-positive value rather than
running with it. ``data`` is an annual climatology with no year dependence, so it
is unaffected, though it still needs the calendar to know the season.

Gas transfer
============

``remora.fennel.co2_schmidt`` and ``remora.fennel.oxygen_schmidt`` select the
Schmidt-number relation for the transfer velocity. Each choice also carries its
own leading rate coefficient -- 0.31 for Wanninkhof (1992), 0.251 for
Wanninkhof (2014) -- so the two move together and are not settable separately.
``ocmip``, available for oxygen only, is Keeling et al. (1998) and pairs with
the 1992 rate coefficient, as it does in ROMS.

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
that need it. A per-side keyword applies the file-driven condition to every
tracer alike, Fennel tracers included, so a per-side ``clamped`` run needs
``NO3_west`` and the rest in the boundary file; per-variable mode is what lets a
biology tracer keep a local condition such as ``outflow`` while temperature and
salinity are clamped. See :ref:`sec:bc-per-tracer` for details.

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
