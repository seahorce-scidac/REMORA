/*
** remora_bio_cppdefs.h
**
** CPP configuration for the REMORA Fennel Fortran bridge (Path A).
**
** ROMS selects biology options at compile time; REMORA selects them at run
** time through remora.fennel.*.  The bridge can therefore only be built for
** ONE ROMS option set at a time.  The set chosen here must match the runtime
** configuration of the case being validated, and the isohelper validates that
** at every call (see fennel_check_config in REMORA_fennel_roms.F) so a
** mismatch aborts loudly instead of producing a silently invalid comparison.
**
** Default set below matches Exec/BioToy/inputs:
**     remora.fennel.carbon          = 1
**     remora.fennel.denitrification = 1
**     remora.fennel.bio_sediment    = 1
**     remora.fennel.po4             = 0   (default)
**     remora.fennel.oxygen          = 0   (default)
**     remora.fennel.odu             = 0   (default)
**
** To validate a different option set, change the defines here and rebuild.
** Only this bridge directory recompiles; the native C++ path is unaffected.
*/

#ifndef REMORA_BIO_CPPDEFS_H
#define REMORA_BIO_CPPDEFS_H

/* Land/sea masking: REMORA always carries rmask, and the native kernel
** multiplies the final increment by it, so MASKING must stay on. */
#define MASKING

/* Fennel option set under validation. */
#define CARBON
#define DENITRIFICATION
#define BIO_SEDIMENT

/* Deliberately NOT defined.
#define OXYGEN
#define PO4
#define ODU
*/

/*
** BULK_FLUXES selects how the gas-exchange transfer velocity is computed:
** from Uwind/Vwind when defined, from sustr/svstr when not.  Both branches
** exist natively in REMORA_Biology.cpp and are chosen at run time by
** remora.bulk_fluxes, so this define must agree with the case being
** validated.  fennel_check_config enforces that on every call.
**
** Defined here because Exec/BioToy/inputs sets remora.bulk_fluxes = true.
*/
#define BULK_FLUXES

/* Not supported by the native port; leaving these off keeps the two paths
** on the same code path.  Enabling any of them requires a native
** implementation and a separate validation lane.
**   RIVER_DON, TALK_NONCONSERV, pCO2_RZ, PCO2AIR_DATA, PCO2AIR_SECULAR,
**   OCMIP_OXYGEN_SC, RW14_OXYGEN_SC, RW14_CO2_SC, SEDBIO_COUP
*/

/* Diagnostics and wet/dry are out of scope for the parity campaign. */
/* #define DIAGNOSTICS_BIO */
/* #define WET_DRY */

/* Explicit-shape dummy arguments: the bridge passes contiguous packed
** buffers with explicit bounds rather than ROMS pointer slices. */
#undef ASSUMED_SHAPE

/* Serial bridge; no ROMS distributed-memory or profiling scaffolding. */
#undef DISTRIBUTE
#undef PROFILE

#endif /* REMORA_BIO_CPPDEFS_H */
