.. role:: cpp(code)
  :language: c++

.. _sec:SeaIce:

*******
Sea Ice
*******

REMORA includes a port of the native sea-ice model of myROMS
(``ROMS/Nonlinear/SeaIce``: Budgell 2005; Hedstrom; Durski and Kurapov; rewritten
by Arango). It has Mellor and Kantha (1989) thermodynamics with snow and melt
ponds, elastic-viscous-plastic (EVP) rheology (Hunke and Dukowicz 1997), and
Smolarkiewicz (MPDATA) or upwind advection of the ice tracers. The ice is a 2D
model coupled to the ocean surface: it reads the bulk-flux atmosphere and the
upper ocean, and sets the surface stress, heat flux and salt flux the ocean
uses.

The port reproduces myROMS, quirks included (see :ref:`sec:SeaIceQuirks`). On
the LAKE_ICE test of myROMS it matches the first ice step to round-off in every
ice field (see :ref:`sec:SeaIceValidation`).

Turning Sea Ice On
==================

.. code:: shell

   remora.ice_model = mellor_kantha     # none (default) or mellor_kantha
   remora.bulk_fluxes = true            # required
   remora.bulk_flux_scheme = coare35    # recommended: the bulk flux of current myROMS

The ice model requires bulk fluxes (the air-ice stress and the ice surface heat
budget come out of the bulk flux), salinity, and a single level
(``amr.max_level = 0``). Regridding is not supported.

``remora.bulk_flux_scheme`` selects which ROMS bulk flux REMORA follows. The
default, ``legacy``, is REMORA's original port. ``coare35`` follows myROMS: COARE
3.5 roughness lengths, the 2024 wind-stress normalization (Wilkin and
Parela-Roman), the buoyancy flux with ``Qstar`` in the denominator, and negative
rain read as snow. The ice model reuses the bulk flux's gustiness and humidities,
so use ``coare35`` to reproduce myROMS. ROMS uses :math:`g = 9.81`; set
``remora.g = 9.81`` to match it.

Options
=======

The CPP options of the ROMS ice model are runtime switches under the ``ice``
prefix.

+----------------------------+---------+------------------------------------------------------------+
| Option                     | Default | Effect (ROMS option)                                       |
+============================+=========+============================================================+
| ``ice.do_thermo``          | true    | Thermodynamics (``ICE_THERMO``, ``ICE_MK``)                |
+----------------------------+---------+------------------------------------------------------------+
| ``ice.do_momentum``        | true    | EVP momentum (``ICE_MOMENTUM``, ``ICE_EVP``)               |
+----------------------------+---------+------------------------------------------------------------+
| ``ice.do_advect``          | true    | Advection of the ice tracers (``ICE_ADVECT``,              |
|                            |         | ``ICE_SMOLAR``)                                            |
+----------------------------+---------+------------------------------------------------------------+
| ``ice.advect_upwind``      | false   | First-order upwind instead of MPDATA (``ICE_UPWIND``)      |
+----------------------------+---------+------------------------------------------------------------+
| ``ice.mom_bulk``           | false   | Quadratic ice-water drag :math:`C_{d,io}|\Delta u|`        |
|                            |         | instead of the log layer (``ICE_MOM_BULK``)                |
+----------------------------+---------+------------------------------------------------------------+
| ``ice.convsnow``           | false   | Convert flooded snow to ice (``ICE_CONVSNOW``)             |
+----------------------------+---------+------------------------------------------------------------+
| ``ice.albedo_type``        | ec92    | ``ec92`` (``ICE_ALB_EC92``), ``default`` (ROMS's default   |
|                            |         | branch) or ``none``. ``csim`` is refused: ROMS's           |
|                            |         | ``ALBEDO_CSIM`` branch does not compile                    |
+----------------------------+---------+------------------------------------------------------------+
| ``ice.copy_old_to_new``    | false   | Not in ROMS: carry the ice tracers to the new time level   |
|                            |         | when advection is off (see :ref:`sec:SeaIceQuirks`)        |
+----------------------------+---------+------------------------------------------------------------+
| ``ice.bc_type``            | closed  | Four entries (xlo ylo xhi yhi): ``closed`` or ``gradient``.|
|                            | x4      | Periodic sides follow ``remora.is_periodic``               |
+----------------------------+---------+------------------------------------------------------------+
| ``ice.gamma2``             | 1       | Slipperiness of closed ice boundaries (ROMS ``GAMMA2``)    |
+----------------------------+---------+------------------------------------------------------------+

The physical parameters are those of ROMS's ``ice.in``, with its defaults:
``ice.nEVP`` (60), ``ice.AirRho`` (1.267), ``ice.IceRho`` (900),
``ice.SnowDryRho`` (330), ``ice.SnowWetRho`` (450), ``ice.Cd_ai`` (3e-3),
``ice.Cd_io`` (7.5e-3), ``ice.Astrength`` (20), ``ice.Pstar`` (2.5e4),
``ice.zetaMin`` (0), ``ice.zetaMax`` (1e14), ``ice.ellip_sq`` (4), ``ice.min_ai``
(0), ``ice.max_ai`` (1), ``ice.min_hi`` (0), ``ice.max_hmelt`` (0.1),
``ice.stressAng`` (0, degrees; read but unused, as in ROMS), ``ice.ice_emiss``,
``ice.spec_heat_air``, ``ice.trans_coeff`` and ``ice.sublimation``.

Time Stepping
=============

There is one ice step per baroclinic step, called from ``setup_step`` after the
bulk flux and the surface flux assembly, which is where ROMS's ``main3d`` calls
``seaice``. In order:

#. The ice albedo, just before the bulk flux, which then computes the air-ice
   stress, the ice surface heat budget and the shortwave left for the ocean.
#. ``ice_spdiw``: the ocean velocity averaged over the top 10 m, and the speed of
   the ice relative to it.
#. ``ice_vbc``: the ice-water transfer coefficient and the ocean surface stress,
   the ice-water stress over the ice-covered fraction and the air-sea stress over
   the rest.
#. With momentum, ``nEVP`` EVP substeps of ``ice_evp``, ``ice_evp_sig`` and
   ``ice_elastic``.
#. With advection, ``ice_enthalpy``, ``ice_advect`` and ``ice_limit``.
#. With thermodynamics, ``ice_mk`` and ``ice_limit``. ``ice_mk`` overwrites the
   ocean surface heat flux and adds to its salt flux.

As in ROMS, every ice state variable has two time levels, indexed by six level
indices (``liold``/``linew`` for the ice tracers, ``liuol``/``liunw`` for the ice
velocity, ``lieol``/``lienw`` for the EVP substeps), swapped in ROMS's order.

Initial State, Input and Output
===============================

With ``remora.ic_type = netcdf``, the ice state is read from the initial file
under its ROMS names (``Aice``, ``ice_thickness``, ``meltpond_thickness``,
``snow_thickness``, ``ice_age``, ``ice_Sxx``, ``ice_Sxy``, ``ice_Syy``, ``Tice``,
``ice_enthalpy``, ``ice_age_thickness``, ``Uice``, ``Vice``, ``Uevp``, ``Vevp``).
Variables the file lacks start at zero. As in ROMS, the ice fields (such as the
molecular sublayer temperature and salinity) are not read, so a ROMS restart
file does not continue exactly. With analytic initial conditions, a problem can
supply an initial state through ``Problem::init_analytic_seaice``; otherwise the
ice starts at zero, as in ROMS.

Checkpoints hold the full ice state, its time-level indices and the ice fields,
so a restart starts from exactly the ice state that was saved.

The NetCDF history file writes every ice state variable and the ice fields ROMS
can write, under their ROMS names, along with ``shflux`` and ``stflx_salt``, the
surface heat and salt fluxes the ocean actually takes. In AMReX plotfiles, list
ice variables in ``remora.plot_vars_2d`` by name, or use ``seaice`` for the whole
state. State variables are written at time level ``liunw``, as ROMS writes them.

Limitations
===========

* Single level only.
* Ice lateral boundaries: closed, gradient and periodic. ROMS's clamped,
  radiation, nudging, mixed and nested types are not ported, nor are boundary or
  climatology files for the ice.
* Not ported: frazil ice (ROMS never calls ``ice_frazil``), ``ICE_SHOREFAST``,
  ``ICE_KEEL``, ``AICLM_NUDGING``, ``MICLM_NUDGING``, ``FASTICE_CLIMATOLOGY``,
  ``OUTFLOW_MASK``, ``ICE_BIO`` and ``ICESHELF`` (most do not compile in ROMS),
  and ``ICE_STRENGTH_QUAD`` (ROMS never sets its gridded strength, so it gives
  zero strength).
* The mixed-layer depth of ``ice_spdiw`` is 10 m. ROMS uses the KPP boundary layer
  depth where deeper; REMORA has no KPP.
* REMORA's Stefan-Boltzmann constant (5.670374419e-8) is more precise than
  ROMS's (5.67e-8), which changes the longwave fluxes by 6.6e-5 relative.

.. _sec:SeaIceQuirks:

Behavior Kept from ROMS
=======================

The port keeps myROMS's behavior, including these quirks, so that it can be
validated against ROMS:

* The ice-ocean shear and the ocean stress use the ice velocity from before this
  step's EVP swap, which from the second step on is the velocity of two steps
  back.
* Without advection, nothing copies the ice tracers from the old time level to
  the new one, so the thermodynamics works on the state of two steps back and a
  run alternates between two independent states. ``ice.copy_old_to_new`` carries
  them forward instead.
* The default albedo branch compares a Celsius temperature as if it were Kelvin,
  so it always gives the dry albedo.
* The ice thermodynamics overwrites the ocean surface heat flux.
* ``ice_vibc`` sets a closed western ghost row from itself rather than from the
  interior, and masks a gradient northern row with the interior row's mask;
  ROMS's ``ice_bc2d`` has no gradient branch, so gradient sides of the rho-point
  fields keep their boundary values.

One ROMS behavior is not kept: ROMS negates the whole enthalpy array inside each
tile's advection call, which races under OpenMP. REMORA negates it once.

.. _sec:SeaIceValidation:

Validation
==========

The port was validated against myROMS (commit 31cabd31) on its LAKE_ICE test
(myroms/roms_test ``lake_ice``), restarted from a day-60 spin-up in which a
quarter of the lake is ice covered, with GLS mixing in both models, and ROMS
writing double-precision output (``OUT_DOUBLE``). On the first ice step, where
the ocean states agree, every ice state variable, ice field, rate and ocean
flux matches ROMS to round-off (at most 1e-13 relative) when REMORA uses ROMS's
Stefan-Boltzmann constant; with REMORA's own constant the differences are at
most 7e-5. This holds for the stock configuration (upwind advection, quadratic
drag), with MPDATA, with the log-layer drag, with thermodynamics and momentum
without advection, and with thermodynamics alone. ROMS's ice age depends on its
tiling, because ``ice_limit`` leaves ``Iage`` out of its MPI exchange; REMORA's
results are bitwise identical across ranks and box layouts, and match a
single-tile ROMS run.

After the first step the two models drift apart through their ocean components,
which differ in ways that have nothing to do with the ice.

``Exec/LakeIce`` holds the inputs of that comparison, and the ``SeaIceBox``
regression test runs a small analytic case.
