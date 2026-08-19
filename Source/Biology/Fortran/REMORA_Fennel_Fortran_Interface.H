#ifndef REMORA_FENNEL_FORTRAN_INTERFACE_H
#define REMORA_FENNEL_FORTRAN_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * C interface to the ROMS Fennel kernel (Path A oracle).
 *
 * Wraps biology_tile from the tracked copy of ROMS fennel.h. See
 * Source/Biology/Fortran/REMORA_fennel_roms.F for the isohelper and
 * Source/Biology/Fortran/tag_map.md for the diagnostic contract.
 *
 * @note All arrays are contiguous packed buffers owned by the caller, in
 *       Fortran (column-major, i fastest) order. Nothing is aliased directly
 *       from an AMReX FArrayBox: REMORA FABs carry ghost cells and a
 *       0-based k index, while ROMS wants 1-based k over the valid region,
 *       so the pack/unpack is explicit and lives in REMORA_FennelBridge.cpp.
 *
 * @note Storage bounds (ilo..ihi, jlo..jhi, nz) and active bounds
 *       (istr..iend, jstr..jend) are passed separately, matching the ROMS
 *       memory-extent versus compute-extent distinction.
 *
 * @note i and j are REMORA's own 0-based indices; the Fortran side sets its
 *       array lower bounds to match, so no horizontal index translation
 *       happens anywhere. Only k is shifted (REMORA k -> ROMS k+1), by the
 *       packer.
 *
 * Array shapes, in Fortran notation:
 *   rmask, srflx, uwind, vwind, pH   (ilo:ihi, jlo:jhi)
 *   sustr, svstr                     (ilo:ihi+1, jlo:jhi+1)
 *   Hz, z_r                          (ilo:ihi, jlo:jhi, 1:nz)
 *   z_w                              (ilo:ihi, jlo:jhi, 0:nz)
 *   t_nstp, t_nnew                   (ilo:ihi, jlo:jhi, 1:nz, 1:nbt+2)
 *
 * The tracer dimension is ROMS-ordered: index 1 is temperature, 2 is
 * salinity, 3.. are the biology tracers in Fennel order. That equals the
 * REMORA component index plus one throughout.
 *
 * @param nbt        biology tracer count; validated against the compiled
 *                   ROMS option set, which aborts the run on mismatch
 * @param dt_sec     full baroclinic step in seconds (NOT divided by BioIter;
 *                   the ROMS kernel does that division itself)
 * @param rho0       reference density, for the ROMS kinematic-srflx
 *                   conversion PARsur = PARfrac*srflx*rho0*Cp
 * @param cp_heat    seawater specific heat, same conversion
 * @param pH         inout for interface fidelity; write-only in practice,
 *                   since fennel.h fixes DoNewton = 0 and the bisection
 *                   brackets are constants
 * @param dbg_level  0 off, 1 target column, 2 all columns
 */
void fennel_bridge_advance_c (
    int ilo, int ihi, int jlo, int jhi, int nz, int nbt,
    int istr, int iend, int jstr, int jend,
    int use_po4, int use_carbon, int use_oxygen, int use_odu,
    int use_denit, int use_biosed, int use_bulk,
    int use_rivdon, int use_talknc,
    int pco2air_kind, int co2_sc_kind, int o2_sc_kind,
    int bioiter, double dt_sec, double rho0, double cp_heat,
    double p_AttSW, double p_AttChl, double p_PARfrac,
    double p_Vp0, double p_I_thNH4, double p_D_p5NH4,
    double p_NitriR, double p_K_NO3, double p_K_NH4,
    double p_K_PO4, double p_K_Phy, double p_Chl2C_m,
    double p_ChlMin, double p_PhyCN, double p_R_P2N,
    double p_PhyIP, double p_PhyIS, double p_PhyMin,
    double p_PhyMR, double p_ZooAE_N, double p_ZooCN,
    double p_ZooBM, double p_ZooER, double p_ZooGR,
    double p_ZooMin, double p_ZooMR, double p_LDeRRN,
    double p_LDeRRC, double p_CoagR, double p_SDeRRN,
    double p_SDeRRC, double p_RDeRRN, double p_RDeRRC,
    double p_wPhy, double p_wLDet,
    double p_wSDet, double p_pCO2air,
    const double* rmask,
    const double* Hz,
    const double* z_r,
    const double* z_w,
    const double* srflx,
    const double* sustr,
    const double* svstr,
    const double* uwind,
    const double* vwind,
    double* pH,
    const double* t_nstp,
    double* t_nnew,
    int dbg_level, int dbg_i, int dbg_j);

#ifdef __cplusplus
}
#endif

#endif // REMORA_FENNEL_FORTRAN_INTERFACE_H
