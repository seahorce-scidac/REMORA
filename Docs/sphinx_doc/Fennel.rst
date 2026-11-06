.. role:: cpp(code)
  :language: c++

.. _sec:Fennel:

**************
Fennel Biology
**************

``remora.biology_model = fennel`` selects a C++ port of the ROMS Fennel
nitrogen-cycle model, with optional phosphate, carbon, oxygen and river pools
selected at run time. Settings shared with other models, such as initial
conditions, output groups and restarts, are on :ref:`sec:Biology`.

With every option at its default, Fennel carries the seven nitrogen tracers:

.. code:: shell

   remora.biology_model = fennel

Options
=======

All Fennel settings use the ``remora.fennel`` prefix.

+------------------------------------+-----------+---------------------------------------------+
| Option                             | Default   | Effect                                      |
+====================================+===========+=============================================+
| ``remora.fennel.po4``              | false     | Adds phosphate, with phosphate limitation   |
|                                    |           | and recycling.                              |
+------------------------------------+-----------+---------------------------------------------+
| ``remora.fennel.carbon``           | false     | Adds carbon detritus, TIC and alkalinity,   |
|                                    |           | carbon cycling, and surface CO2 exchange.   |
+------------------------------------+-----------+---------------------------------------------+
| ``remora.fennel.oxygen``           | false     | Adds dissolved oxygen, its source and sink  |
|                                    |           | terms, and surface O2 exchange.             |
+------------------------------------+-----------+---------------------------------------------+
| ``remora.fennel.odu``              | false     | Adds the ``ODU`` tracer, which has no       |
|                                    |           | source or sink terms in this port.          |
+------------------------------------+-----------+---------------------------------------------+
| ``remora.fennel.river_don``        | false     | Adds a river dissolved organic nitrogen     |
|                                    |           | pool that remineralizes to ``NH4``, plus a  |
|                                    |           | carbon counterpart with ``carbon``. Neither |
|                                    |           | sinks.                                      |
+------------------------------------+-----------+---------------------------------------------+
| ``remora.fennel.bio_sediment``     | true      | Returns sinking particulate material at the |
|                                    |           | bottom.                                     |
+------------------------------------+-----------+---------------------------------------------+
| ``remora.fennel.denitrification``  | false     | Changes the bottom nitrogen return under    |
|                                    |           | ``bio_sediment``. Adds no tracer.           |
+------------------------------------+-----------+---------------------------------------------+
| ``remora.fennel.talk_nonconserv``  | false     | Makes alkalinity prognostic; see            |
|                                    |           | `Alkalinity`_. Requires ``carbon``.         |
+------------------------------------+-----------+---------------------------------------------+
| ``remora.fennel.pco2air_type``     | constant  | Atmospheric CO2: ``constant``, ``data`` or  |
|                                    |           | ``secular``; see `Atmospheric CO2`_.        |
+------------------------------------+-----------+---------------------------------------------+
| ``remora.fennel.co2_schmidt``      | ``w92``   | CO2 gas transfer relation:                  |
|                                    |           | ``wanninkhof1992`` (``w92``) or             |
|                                    |           | ``wanninkhof2014`` (``rw14``).              |
+------------------------------------+-----------+---------------------------------------------+
| ``remora.fennel.oxygen_schmidt``   | ``w92``   | O2 gas transfer relation: as for CO2, or    |
|                                    |           | ``ocmip`` (Keeling et al. 1998).            |
+------------------------------------+-----------+---------------------------------------------+

The oxygen and carbon options need surface forcing for the gas transfer
velocity: 10 m winds when ``remora.bulk_fluxes`` is true, surface stress
otherwise.

Tracers
=======

The seven nitrogen tracers are always present, followed by the optional ones
in this order, each only when its option is on:

.. code:: text

   NO3 NH4 chlorophyll phytoplankton zooplankton LdetritusN SdetritusN
   RdetritusN                              (river_don)
   PO4                                     (po4)
   LdetritusC SdetritusC TIC alkalinity    (carbon)
   RdetritusC                              (carbon and river_don)
   oxygen                                  (oxygen)
   ODU                                     (odu)

This is the ROMS ``initialize_biology`` order, so a NetCDF file written for
either code reads into the same tracers. Nitrogen-only Fennel has 7 tracers
and every option at once 16.

Parameters
==========

Each parameter is set as ``remora.fennel.<name>``. Defaults follow ROMS; units
are in ``REMORA_FennelBiology.H``. ``K_Phy`` must be positive and ``AttSW`` and
``AttChl`` non-negative, and with ``carbon`` a constant ``pCO2air`` must be
positive.

+-------------+-------------+---------------------------------------------+
| Parameter   | Default     | Meaning                                     |
+=============+=============+=============================================+
| ``BioIter`` | 1           | Iterations of the nonlinear update; 0 or    |
|             |             | less disables the source terms              |
+-------------+-------------+---------------------------------------------+
| ``AttSW``   | 0.04        | Seawater light attenuation [1/m]            |
+-------------+-------------+---------------------------------------------+
| ``AttChl``  | 0.02486     | Chlorophyll light attenuation               |
+-------------+-------------+---------------------------------------------+
| ``PARfrac`` | 0.43        | Shortwave fraction available for biology    |
+-------------+-------------+---------------------------------------------+
| ``Vp0``     | 1.0         | Phytoplankton growth tuning parameter       |
+-------------+-------------+---------------------------------------------+
| ``I_thNH4`` | 0.0095      | Radiation threshold for nitrification       |
+-------------+-------------+---------------------------------------------+
| ``D_p5NH4`` | 0.1         | Half-saturation radiation for nitrification |
+-------------+-------------+---------------------------------------------+
| ``NitriR``  | 0.05        | Nitrification rate [day-1]                  |
+-------------+-------------+---------------------------------------------+
| ``K_NO3``   | 2.0         | NO3 uptake inverse half-saturation          |
+-------------+-------------+---------------------------------------------+
| ``K_NH4``   | 2.0         | NH4 uptake inverse half-saturation          |
+-------------+-------------+---------------------------------------------+
| ``K_PO4``   | 32.0        | PO4 uptake inverse half-saturation          |
+-------------+-------------+---------------------------------------------+
| ``K_Phy``   | 1.0         | Zooplankton ingestion half-saturation       |
+-------------+-------------+---------------------------------------------+
| ``Chl2C_m`` | 0.0535      | Maximum chlorophyll-to-carbon ratio         |
+-------------+-------------+---------------------------------------------+
| ``ChlMin``  | 0.01        | Chlorophyll minimum threshold               |
+-------------+-------------+---------------------------------------------+
| ``PhyCN``   | 6.625       | Phytoplankton carbon-to-nitrogen ratio      |
+-------------+-------------+---------------------------------------------+
| ``R_P2N``   | 0.0625      | Phytoplankton phosphate-to-nitrogen ratio   |
+-------------+-------------+---------------------------------------------+
| ``PhyIP``   | 1.5         | Phytoplankton NH4 inhibition parameter      |
+-------------+-------------+---------------------------------------------+
| ``PhyIS``   | 0.025       | Phytoplankton P-I curve initial slope       |
+-------------+-------------+---------------------------------------------+
| ``PhyMin``  | 0.01        | Phytoplankton minimum threshold             |
+-------------+-------------+---------------------------------------------+
| ``PhyMR``   | 0.15        | Phytoplankton mortality rate [day-1]        |
+-------------+-------------+---------------------------------------------+
| ``ZooAE_N`` | 0.75        | Zooplankton nitrogen assimilation fraction  |
+-------------+-------------+---------------------------------------------+
| ``ZooCN``   | 6.625       | Zooplankton carbon-to-nitrogen ratio        |
+-------------+-------------+---------------------------------------------+
| ``ZooBM``   | 0.1         | Zooplankton basal metabolism [day-1]        |
+-------------+-------------+---------------------------------------------+
| ``ZooER``   | 0.1         | Zooplankton excretion rate [day-1]          |
+-------------+-------------+---------------------------------------------+
| ``ZooGR``   | 0.6         | Zooplankton maximum growth rate [day-1]     |
+-------------+-------------+---------------------------------------------+
| ``ZooMin``  | 0.01        | Zooplankton minimum threshold               |
+-------------+-------------+---------------------------------------------+
| ``ZooMR``   | 0.025       | Zooplankton mortality rate [day-1]          |
+-------------+-------------+---------------------------------------------+
| ``LDeRRN``  | 0.01        | Large nitrogen detritus remineralization    |
+-------------+-------------+---------------------------------------------+
| ``LDeRRC``  | 0.01        | Large carbon detritus remineralization      |
+-------------+-------------+---------------------------------------------+
| ``CoagR``   | 0.005       | Coagulation rate [day-1]                    |
+-------------+-------------+---------------------------------------------+
| ``SDeRRN``  | 0.1         | Small nitrogen detritus remineralization    |
+-------------+-------------+---------------------------------------------+
| ``SDeRRC``  | 0.03        | Small carbon detritus remineralization      |
+-------------+-------------+---------------------------------------------+
| ``RDeRRN``  | 0.03        | River nitrogen detritus remineralization    |
|             |             | [day-1], with ``river_don``                 |
+-------------+-------------+---------------------------------------------+
| ``RDeRRC``  | 0.03        | River carbon detritus remineralization      |
|             |             | [day-1], with ``river_don``                 |
+-------------+-------------+---------------------------------------------+
| ``wPhy``    | 0.1         | Phytoplankton sinking velocity [m/day]      |
+-------------+-------------+---------------------------------------------+
| ``wLDet``   | 1.0         | Large detritus sinking velocity [m/day]     |
+-------------+-------------+---------------------------------------------+
| ``wSDet``   | 0.1         | Small detritus sinking velocity [m/day]     |
+-------------+-------------+---------------------------------------------+
| ``pCO2air`` | 370.0       | Atmospheric CO2 partial pressure [ppmv],    |
|             |             | with ``pco2air_type = constant``            |
+-------------+-------------+---------------------------------------------+

Alkalinity
==========

By default alkalinity is diagnostic: each biology iteration overwrites it with
:math:`587.05 + 50.56\,S` (Brewer et al. 1986), so it carries no memory of the
biology.

With ``talk_nonconserv = true`` it is prognostic, responding to nitrate uptake
(:math:`+`), ammonium uptake (:math:`-`), nitrification (:math:`-2` per mole
of :math:`NH_4`), zooplankton metabolism and excretion (:math:`+`),
remineralization (:math:`+`) and the sediment return (:math:`+`, except under
``denitrification``). Alkalinity enters only the pCO2 calculation, so this
changes the air-sea CO2 flux and hence TIC, and nothing else.

Atmospheric CO2
===============

``pco2air_type`` sets the atmospheric CO2 partial pressure for the surface gas
exchange:

``constant``
    ``remora.fennel.pCO2air``. The default.

``data``
    The annual climatology of Laurent et al. (2017),
    :math:`380.464 + 9.321 \sin(2\pi\, \mathrm{yday}/365.25 + 1.068)`.

``secular``
    A linear trend plus three harmonics, referenced to 1951.

Both time-dependent forms date the model clock with ``remora_caldate``, a port
of ROMS ``caldate``: ``remora.time_ref`` gives the reference date and calendar
(see :ref:`calendar`). Either of these places a run on 1 January 2020:

.. code:: shell

   remora.time_ref   = 20200101

   # or, keeping the default 0001-01-01 epoch:
   remora.start_time = 63713433600

This matters for ``secular``: extrapolated far from 1951 its fit goes negative,
to about -2648 ppmv at year 1, which would draw CO2 out of the ocean. REMORA
aborts on a non-positive value. ``data`` has no year dependence but still needs
the calendar for the season.

Gas transfer
============

``co2_schmidt`` and ``oxygen_schmidt`` choose the Schmidt-number relation for
the transfer velocity. Each carries its own leading rate coefficient, 0.31 for
Wanninkhof (1992) and 0.251 for Wanninkhof (2014), so the two are not set
separately. ``ocmip`` pairs with the 1992 coefficient, as in ROMS.

Analytic initial condition
==========================

The standard ``Problem`` gives ``BioToy`` and ``Upwelling`` the Fennel profile
from ROMS ``ana_biology.h``
(``Source/Prob/REMORA_InitAnalyticBiology_BioToy.H``): nitrate follows a cubic
in a temperature-derived silicate proxy, and the other tracers start uniform,
including any optional pools the run carries. The profile depends on
temperature alone, so another problem that stratifies temperature can use it by
adding its name to the BioToy branch of ``Problem::init_analytic_biology``.
It aborts for any model but Fennel. See :ref:`sec:Biology` for the hook and
the other initial-condition choices.

Validation against ROMS
=======================

REMORA can build the original ROMS Fennel kernel alongside the C++ port and
select between them at run time. This is validation infrastructure, compiled
out by default.

+-------------------------------+--------------------------------------+-----------------------------------+
| GNUmake                       | CMake                                | Builds                            |
+===============================+======================================+===================================+
| ``USE_FENNEL_FORT=TRUE``      | ``-DREMORA_ENABLE_FENNEL_FORT=ON``   | The ROMS reference kernel         |
+-------------------------------+--------------------------------------+-----------------------------------+
| ``USE_BIOLOGY_DIAG=TRUE``     | ``-DREMORA_ENABLE_BIOLOGY_DIAG=ON``  | Per-block diagnostic output       |
+-------------------------------+--------------------------------------+-----------------------------------+

The two are independent. ``REMORA_ENABLE_FENNEL_FORT`` must be on the first
CMake configure line, since it enables the Fortran language.

+-----------------------------------+---------+-----------------------------------------------+
| Option                            | Default | Effect                                        |
+===================================+=========+===============================================+
| ``remora.use_biology_cpp_answer`` | 1       | ``1`` runs the C++ kernel, ``0`` the ROMS     |
|                                   |         | kernel, which must be compiled in.            |
+-----------------------------------+---------+-----------------------------------------------+
| ``remora.biology_debug``          | 0       | ``1`` writes diagnostics for one column,      |
|                                   |         | ``2`` for every column; needs the             |
|                                   |         | diagnostics compiled in.                      |
+-----------------------------------+---------+-----------------------------------------------+
| ``remora.biology_debug_i``,       | 2, 2    | The column for ``biology_debug = 1``          |
| ``remora.biology_debug_j``        |         |                                               |
+-----------------------------------+---------+-----------------------------------------------+

Requesting something not compiled in aborts with the flag to rebuild with,
rather than falling back. ROMS fixes its options at compile time, so the
reference kernel is built for the one option set pinned in
``Source/Biology/Fortran/remora_bio_cppdefs.h``, which is checked against the
runtime options on every call. ``biology_debug = 2`` output grows with tags,
levels and tracers, so it suits only very small grids. See
``Source/Biology/Fortran/README.md`` for the full procedure.

Regression coverage
===================

``Upwelling_Fennel`` (41x80x16, ten steps, analytic physics and biology)
compares against a gold file with a deliberately wide option set: ``carbon``,
``oxygen``, ``river_don`` and ``talk_nonconserv`` on, ``pco2air_type =
secular``, both Schmidt relations ``rw14``, and ``denitrification`` off.
Where an option has a default and an alternative, the test generally takes the
alternative, since the Fortran-bridge comparison already checks the defaults
against ROMS. ``denitrification`` stays off because the alkalinity increment
from the sediment return exists only in the branch ROMS takes without it.

The test sets ``remora.time_ref = 20200301``, just after a leap day, so the
secular trend depends on the calendar's leap-day branch; moving the run one
day changes the gold file by about 3e-8 relative, well above the 1e-11
tolerance.

``Upwelling_Fennel_tiled`` reruns the case with small horizontal CPU tiles
against the same gold file. ``Upwelling_Fennel_hires_init`` and
``Upwelling_Fennel_hires_init_ic`` check the analytic profile and its
average-down from a high-resolution level. ``Upwelling_Fennel_restart`` checks
that a restart reproduces the run exactly and that restart rejects a
checkpoint with a different model or component layout (see
:ref:`sec:Biology`). The ``BiologyInterface`` unit test checks the tracer
layout and name lookup for every option combination, the default parameters,
constant initialization, the analytic profile and the disabled update. Left to the Fortran-bridge comparison: ``po4``, ``odu``,
``denitrification``, ``bio_sediment = false``, ``BioIter > 1``, the ``data``
form of ``pco2air_type``, the ``w92`` and ``ocmip`` relations, and NetCDF
biology initialization.
