.. _sec:ECB:

ECB Biology
===========

``remora.biology_model = ecb`` selects a native C++ port of the Estuarine
Carbon Biogeochemistry model through the common ``Biology`` interface.
ECB owns its parameters, ordered inventory and column kernel. Analytic profiles
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

All 42 rate/optical/stoichiometric parameters used by this variant are required
under ``remora.ecb``. Their names and units follow ROMS. There are no inherited
Fennel defaults. ``BioIter`` defaults to 1; a nonpositive value disables the
source update. Input half-saturation coefficients ``K_NO3``, ``K_NH4`` and
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

Analytic initialization uses 18 explicitly supplied concentrations in
``remora.biology_initial``, in the following fixed order. Temperature, salinity
and dyes are initialized by the physical application and left unchanged.
The standard ``Problem`` implements this generic uniform profile. A custom
application can instead implement the common ``init_analytic_biology``
hook for a spatial profile. Missing uniform values fail clearly. Existing
NetCDF initialization uses the names below and requires every active tracer;
that path still needs an ECB execution validation case.

.. list-table:: Ordered ECB fields and NetCDF/output names
   :header-rows: 1

   * - Names, in order
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
A column must contain at least two layers, and each grid box must span its
entire vertical column. Horizontal decomposition is allowed; vertically
split grid boxes fail explicitly. CPU tiling does not split this kernel.

New checkpoints contain a versioned ``BiologyHeader`` sidecar recording model
variant identity and every conserved component name in order. Restart rejects
an identity or ordering mismatch even when component counts agree. Changing
rate parameters alone is permitted. The existing checkpoint data formats are
reused. Legacy checkpoints without the sidecar remain accepted for Fennel and
none, with their original count-only identity limits; ECB requires the sidecar.

Validation and remaining work
-----------------------------

``BiologyInterface`` compares native increments with 20 independent clean-ROMS
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
