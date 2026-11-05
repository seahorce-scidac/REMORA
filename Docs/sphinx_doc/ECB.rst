.. _sec:ECB:

ECB Biology
===========

``remora.biology_model = ecb`` selects a native C++ port of the Estuarine
Carbon Biogeochemistry model through the common ``Biology`` interface.
ECB owns its parameters, tracer inventory and column kernel. Analytic profiles
belong to the application, following the reviewed Stage A API. Temperature,
salinity and passive dyes precede its biology block as for Fennel.

Scientific reference
--------------------

The port follows clean `ROMS trunk (develop), commit 31cabd318c9379429c365c308b9e0f2ee166a6fc
<https://github.com/myroms/roms/tree/31cabd318c9379429c365c308b9e0f2ee166a6fc>`_,
retrieved on 3 November 2026. The source routines are ``ecb_tile`` and
``pCO2_water`` in ``ROMS/Nonlinear/Biology/ecb.h``; parameters and indices
come from ``ecb_mod.h``. No local ROMS scientific patches are incorporated.
ROMS's MIT license notice is retained with the translated kernel.

This first port implements one explicit option set:

* Carbon, oxygen, instantaneous sediment recycling and denitrification,
  which the reference ECB code requires implicitly.
* Phosphate limitation (``PO4``) and nonconservative alkalinity
  (``TALK_NONCONSERV``).
* Wanninkhof 2014 oxygen and CO2 Schmidt numbers and their paired gas
  transfer coefficient (``RW14_OXYGEN_SC``, ``RW14_CO2_SC``).
* Constant atmospheric pCO2, OCMIP carbonate chemistry with its 30 bisections,
  zero phosphate/silicate chemistry arguments as in the reference call,
  and the default monotone sinking boundary reconstruction.
* Land masking, with the common physical wet mask applied to increments.
  ROMS's separate ``WET_DRY`` branch is not implemented.

The code retains reference mathematics, including the **maximum** with 4.0
in the temperature growth expression, the extraction/final TIC bounds,
and the ordered resuspension/burial updates. Refractory organic nitrogen,
refractory organic carbon and calcium have no biological source/sink in this
reference; they remain transported tracers. The source does not establish
that all 18 pools participate in every process.

Alternative carbonate solvers, atmospheric pCO2 trends, sediment-proxy optics,
net sulfate reduction, diagnostics and other CPP variants require separate
ports and validation. There are no runtime ECB switches for these variants.
The Fennel backend switch does not select a Fortran implementation of ECB.

Configuration and initialization
--------------------------------

``Tests/test_files/Upwelling_ECB/Upwelling_ECB.i`` is a complete synthetic
example using analytic physics, a uniform ECB initial condition, one passive
dye and bulk forcing. Its numerical parameters exercise the equations; they
are not a calibrated Chesapeake Bay setup.

The 42 rate, optical and stoichiometric parameters are read from
``remora.ecb``, with names and units following ROMS. Each defaults to the value
ROMS ``External/bio_ECB.in`` sets at the pinned commit, not to the Fennel value
of the same name; ``REMORA_ECBBiology.H`` documents each one. ``BioIter``
defaults to 1; a nonpositive value disables the source update. ``AttSW``,
``slCexc`` and ``K_Phy`` must be positive and ``AttChl`` non-negative. Input half-saturation coefficients ``K_NO3``, ``K_NH4`` and
``K_PO4`` are the inverse values used by the ROMS kernel (m3/mmol).
``wPhy``, ``wSDet`` and ``wLDet`` are in m/day; rates are in day^-1 unless
specified otherwise in the reference parameter descriptions.

The input names are:

.. code:: text

   AttSW AttChl Chl2C_m ChlMin CoagR D_p5NH4 I_thNH4 K_NH4
   K_NO3 K_PO4 K_Phy LDeRRN LDeRRC NitriR PARfrac PhyCN R_P2N
   PhyIS PhyMin PhyMR SDeRRN SDeRRC RDeRRN RDeRRC Vp0
   wLDet wPhy wSDet ZooAE_N ZooBM ZooCN ZooER ZooGR ZooMin ZooMR
   ElDON EsDON deltN deltC gammaC slCexc pCO2air

For a uniform start, set ``remora.biology_ic_type = constant`` and give each of
the 18 tracers below by name as ``remora.biology_initial_<name>``, for example
``remora.biology_initial_TIC = 2100.0``; ``Tests/test_files/Upwelling_ECB``
sets all 18. Every tracer must be given, and ghost cells are set too.
Temperature, salinity and dyes are initialized by the physical application and
left unchanged. ECB has no built-in analytic profile; an application can
implement the common ``init_analytic_biology`` hook for a spatial one and use
``remora.biology_ic_type = analytic``. NetCDF initialization uses the names
below and requires every active tracer; that path still needs an ECB execution
validation case.

Every interface refers to ECB tracers by name: the constant initial values,
NetCDF initialization, plot variables and the kernel itself, which looks up
each tracer's state component from its name. The order below is only the
order of the components in plotfiles and checkpoints.

.. list-table:: ECB tracer names and units
   :header-rows: 1

   * - Names
     - Units
   * - NO3, NH4
     - mmol N/m3
   * - chlorophyll
     - mg chlorophyll/m3
   * - phytoplankton, zooplankton, LdetritusN, SdetritusN, DON
     - mmol N/m3
   * - PO4
     - mmol P/m3
   * - LdetritusC, SdetritusC, TIC
     - mmol C/m3
   * - alkalinity
     - milliequivalents/m3
   * - DOC
     - mmol C/m3
   * - oxygen
     - mmol O2/m3
   * - Rdon
     - mmol N/m3
   * - Rdoc, calcium
     - mmol C/m3 and mmol Ca/m3, respectively

``DON``/``DOC`` are semilabile pools; ``Rdon``/``Rdoc`` are refractory pools.
These meanings differ from Fennel's river-detritus fields.
``biology`` and the ECB-declared alias ``ecb`` expand to all 18 fields.

Update and restart
------------------

The source update stays at the existing biology call in ``advance_3d``.
Old concentrations drive scratch chemistry; only ``mask * Hz * (Bio-Bio_old)``
is added to the incoming weighted state. The selected surface wind/stress and
both bottom-stress components are supplied through explicit common views.
The solver must allocate shortwave forcing through ``bulk_fluxes`` or
``atm2ocn_flux_mode``. REMORA's shortwave flux is in W/m\ :sup:`2`, so surface
PAR is ``PARfrac * srflx``. ROMS ``ecb.h`` also multiplies by ``rho0 * Cp``
because its ``srflx`` is kinematic; REMORA does not.

``configure`` rejects ``AttSW <= 0``, ``AttChl < 0``, ``slCexc <= 0`` and
``K_Phy <= 0``, each of which divides by zero in the kernel.
A column must contain at least two layers. Like Fennel, the kernel works on
whole columns, which REMORA's grids and CPU tiles always provide; the answer
does not depend on the horizontal tiling.

Every checkpoint ``Header`` ends with a ``cons_names`` line listing the
cell-centered components in order, and restart rejects a checkpoint whose list
differs from the run's, even when the counts agree. This check is not specific
to biology: it also catches a changed ``remora.nscalar``. Checkpoints written
before the line existed are checked by count only. A separate versioned
``BiologyHeader`` records the biology model name, which restart also checks;
it does not record the ROMS version the port follows. Changing rate parameters alone is permitted. Legacy checkpoints without
``BiologyHeader`` remain accepted for Fennel and none; ECB requires it.

Validation and remaining work
-----------------------------

``BiologyInterface`` compares native increments with 24 independent clean-ROMS
column fixtures, generated by the reference Fortran code with bounds checking.
It covers day/night, low nutrients/oxygen, gas exchange, resuspension, sinking,
masks, iteration counts, layer counts, dye offsets, immutable old state and
preservation of incoming weighted state. The generator, source hashes, option
set and fixture format are in ``Tests/reference/ECB``. Fortran and NumPy are
needed only to regenerate fixtures, not to build or run the native tests.

``Upwelling_ECB`` exercises analytic initialization, solver advancement,
transport, output and restart using the existing regression machinery.
These tests establish a native second ecosystem and focused reference parity.
They do not complete the v4 plan's first Stage D: matched full ROMS/REMORA
simulation output, NetCDF/river/boundary forcing, AMR/regridding, MPI and GPU
execution still need validation before optional module infrastructure (Stage B).
No oyster or harmful-algal-bloom model is part of this change.
