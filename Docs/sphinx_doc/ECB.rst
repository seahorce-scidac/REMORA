.. _sec:ECB:

***********
ECB Biology
***********

``remora.biology_model = ecb`` selects a C++ port of the ROMS Estuarine Carbon
Biogeochemistry (ECB) model. Settings shared with other models, such as
initial conditions, output groups and restarts, are on :ref:`sec:Biology`.
``Tests/test_files/Upwelling_ECB/Upwelling_ECB.i`` is a complete example; its
parameters exercise the equations and are not a calibrated setup.

Reference and option set
========================

The port follows ``ecb_tile`` and ``pCO2_water`` in
``ROMS/Nonlinear/Biology/ecb.h``, with parameters and indices from
``ecb_mod.h``, at `ROMS commit 31cabd318c
<https://github.com/myroms/roms/tree/31cabd318c9379429c365c308b9e0f2ee166a6fc>`_,
without local changes. ROMS's MIT license notice is kept with the kernel.

ECB has one fixed option set, the ROMS options:

* ``CARBON``, ``OXYGEN``, ``BIO_SEDIMENT`` and ``DENITRIFICATION``, which ECB
  requires;
* ``PO4`` and ``TALK_NONCONSERV``;
* ``RW14_OXYGEN_SC`` and ``RW14_CO2_SC``, the Wanninkhof (2014) gas transfer
  relations;
* constant atmospheric pCO2, the OCMIP carbonate solver, land masking, and the
  default monotone sinking reconstruction.

The Fennel option names (``remora.ecb.carbon``, ``po4``, ``co2_schmidt`` and
so on) are rejected rather than ignored. ROMS's ``WET_DRY`` branch, other
carbonate solvers, pCO2 trends, sediment-proxy optics, net sulfate reduction
and diagnostics are not ported, and ``remora.use_biology_cpp_answer`` does not
select a Fortran ECB.

Tracers
=======

.. list-table::
   :header-rows: 1

   * - Tracers, in state order
     - Units
   * - ``NO3``, ``NH4``
     - mmol N/m3
   * - ``chlorophyll``
     - mg chlorophyll/m3
   * - ``phytoplankton``, ``zooplankton``, ``LdetritusN``, ``SdetritusN``, ``DON``
     - mmol N/m3
   * - ``PO4``
     - mmol P/m3
   * - ``LdetritusC``, ``SdetritusC``, ``TIC``
     - mmol C/m3
   * - ``alkalinity``
     - milliequivalents/m3
   * - ``DOC``
     - mmol C/m3
   * - ``oxygen``
     - mmol O2/m3
   * - ``Rdon``
     - mmol N/m3
   * - ``Rdoc``
     - mmol C/m3
   * - ``calcium``
     - mmol Ca/m3

``DON`` and ``DOC`` are the semilabile pools, and ``Rdon`` and ``Rdoc`` the
refractory ones, unlike Fennel's river pools. In this ROMS version ``Rdon``,
``Rdoc`` and calcium have no biological sources or sinks: they are only
transported, and ``Rdoc`` adds to light attenuation. ECB has no built-in
analytic profile, so start it from constant or NetCDF values, or supply a
profile (see :ref:`sec:Biology`).

Parameters
==========

Each of the 42 parameters is set as ``remora.ecb.<name>``, with ROMS names and
units. Each defaults to the value ROMS ``External/bio_ECB.in`` sets at the
pinned commit, not the Fennel value of the same name; ``REMORA_ECBBiology.H``
documents each one.

.. code:: text

   AttSW AttChl Chl2C_m ChlMin CoagR D_p5NH4 I_thNH4 K_NH4
   K_NO3 K_PO4 K_Phy LDeRRN LDeRRC NitriR PARfrac PhyCN R_P2N
   PhyIS PhyMin PhyMR SDeRRN SDeRRC RDeRRN RDeRRC Vp0
   wLDet wPhy wSDet ZooAE_N ZooBM ZooCN ZooER ZooGR ZooMin ZooMR
   ElDON EsDON deltN deltC gammaC slCexc pCO2air

``BioIter`` (default 1) sets the iterations of the update; 0 or less disables
it. ``K_NO3``, ``K_NH4`` and ``K_PO4`` are inverse half-saturations
[m3/mmol], sinking velocities are in m/day, and rates in day-1 unless stated
otherwise. ``slCexc``, ``K_Phy`` and ``pCO2air`` must be positive and ``AttSW`` and
``AttChl`` non-negative.

Forcing
=======

ECB reads shortwave radiation, both bottom-stress components, and the 10 m
winds under ``remora.bulk_fluxes`` or the surface stress otherwise; the
shortwave needs ``bulk_fluxes`` or ``atm2ocn_flux_mode``. REMORA's shortwave
is in W/m\ :sup:`2`, so surface PAR is ``PARfrac * srflx``; ROMS also
multiplies by ``rho0 * Cp`` because its ``srflx`` is kinematic. A column must
have at least two layers.

Validation
==========

The tests parallel Fennel's. ``Upwelling_ECB`` (8x8x8, three steps, analytic
physics and constant biology) compares against a gold file.
``Seamount_ECB_tiled`` runs ECB on Seamount, whose biology varies in both x and
y, untiled and with small horizontal CPU tiles, and requires the two runs to
agree bitwise. Upwelling is uniform in x, so it could not show a kernel that
mixes up columns across x tile edges. The ``BiologyInterface`` unit test checks
the tracer layout and name lookup, the default parameters, constant
initialization and the disabled update. During the port the kernel matched the
reference ROMS Fortran exactly on single-column cases; that comparison is not
kept as a test.

Not yet validated: a full ROMS and REMORA simulation compared side by side,
NetCDF initialization, river, boundary and climatology forcing, mesh
refinement, and GPU execution.
