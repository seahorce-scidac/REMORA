/**
 * \file REMORA_FennelBridge.cpp
 *
 * Path A of the Fennel parity campaign: pack REMORA state into ROMS-shaped
 * buffers, call the unmodified ROMS kernel through the isohelper, and unpack
 * the result. Selected at run time by remora.use_biology_cpp_answer = 0.
 *
 * The pack/unpack boundary is deliberately explicit rather than aliasing
 * FArrayBox memory: REMORA FABs carry ghost cells and index k from 0 at the
 * bottom, while ROMS wants a contiguous 1-based k over the valid region. An
 * explicit copy keeps both paths auditable and immune to ghost-layout
 * surprises, at a cost that is irrelevant for validation runs.
 *
 * See Source/Biology/Fortran/tag_map.md for the diagnostic contract.
 */

#include <AMReX_MultiFab.H>

#include <REMORA.H>
#include <REMORA_Biology.H>
#include <REMORA_Constants.H>

#include "REMORA_Fennel_Fortran_Interface.H"

using namespace amrex;

namespace {

// Column-major (Fortran) linear index into a packed buffer covering
// [ilo,ihi] x [jlo,jhi] with nk levels and nc components.
AMREX_FORCE_INLINE
int fidx (int i, int j, int k, int n,
          int ilo, int jlo, int ni, int nj, int nk) noexcept
{
    return ((n * nk + k) * nj + (j - jlo)) * ni + (i - ilo);
}

} // namespace

void
REMORA::advance_biology_fortran (int lev, MultiFab const& mf_cons_old,
                                 MultiFab& mf_cons_new,
                                 int N, Real dt_lev)
{
#ifdef AMREX_USE_GPU
    amrex::ignore_unused(lev, mf_cons_old, mf_cons_new, N, dt_lev);
    amrex::Abort("Fennel Fortran bridge is host-only. Serial-first acceptance "
                 "is the default validation scope; GPU parity is a separate "
                 "lane and is not opened by this campaign.");
#else

    const auto parms = fennel_params;

    const bool need_gas = parms.oxygen || parms.carbon;
    const bool do_bulk_flux = solverChoice.bulk_fluxes;

    if (vec_Hz[lev] == nullptr || vec_z_w[lev] == nullptr ||
        vec_z_r[lev] == nullptr || vec_mskr[lev] == nullptr ||
        vec_srflx[lev] == nullptr ||
        (need_gas && do_bulk_flux &&
         (vec_uwind[lev] == nullptr || vec_vwind[lev] == nullptr)) ||
        (need_gas && !do_bulk_flux &&
         (vec_sustr[lev] == nullptr || vec_svstr[lev] == nullptr))) {
        amrex::Abort("Fennel Fortran bridge requires Hz, z_r, z_w, rmask, "
                     "srflx, and, when carbon or oxygen is active, either "
                     "uwind/vwind (bulk_fluxes) or sustr/svstr");
    }

    // ROMS tracer index equals REMORA component index plus one throughout:
    // itemp=1/Temp_comp=0, isalt=2/Salt_comp=1, iNO3_=3/Tracer_comp=2. The
    // Fennel tracer ordering in REMORA_Biology.H::components() matches
    // fennel_mod.h:527-561 exactly, so a single offset covers all of them.
    const int nbt  = REMORABiology::Fennel::num_tracers(parms);
    const int ntrc = nbt + 2;

    // ROMS N(ng) is a level count; REMORA N is the last valid k index.
    const int nz = N + 1;

    // REMORA stores shortwave radiation in W/m2. ROMS fennel.h expects
    // kinematic srflx and converts it back with rho0*Cp when forming PARsur.
    const Real srflx_to_roms_kinematic =
        Real(1.0) / (solverChoice.rho0 * Cp);

    for (MFIter mfi(mf_cons_new, false); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();

        const int ilo = bx.smallEnd(0);
        const int ihi = bx.bigEnd(0);
        const int jlo = bx.smallEnd(1);
        const int jhi = bx.bigEnd(1);

        const int ni = ihi - ilo + 1;
        const int nj = jhi - jlo + 1;

        // Stress buffers are one cell wider in each direction because
        // fennel.h averages sustr(i+1,j) and svstr(i,j+1) to cell centres.
        const int nis = ni + 1;
        const int njs = nj + 1;

        Array4<Real const> const& state_old = mf_cons_old.const_array(mfi);
        Array4<Real const> const& state_new = mf_cons_new.const_array(mfi);
        Array4<Real>       const& state_out = mf_cons_new.array(mfi);
        Array4<Real const> const& Hz_a    = vec_Hz[lev]->const_array(mfi);
        Array4<Real const> const& z_r_a   = vec_z_r[lev]->const_array(mfi);
        Array4<Real const> const& z_w_a   = vec_z_w[lev]->const_array(mfi);
        Array4<Real const> const& srflx_a = vec_srflx[lev]->const_array(mfi);
        Array4<Real const> const& mskr_a  = vec_mskr[lev]->const_array(mfi);

        Array4<Real const> sustr_a;
        Array4<Real const> svstr_a;
        Array4<Real const> uwind_a;
        Array4<Real const> vwind_a;
        if (need_gas) {
            if (do_bulk_flux) {
                uwind_a = vec_uwind[lev]->const_array(mfi);
                vwind_a = vec_vwind[lev]->const_array(mfi);
            } else {
                sustr_a = vec_sustr[lev]->const_array(mfi);
                svstr_a = vec_svstr[lev]->const_array(mfi);
            }
        }

        Vector<double> b_rmask(std::size_t(ni) * nj, 0.0);
        Vector<double> b_srflx(std::size_t(ni) * nj, 0.0);
        Vector<double> b_pH   (std::size_t(ni) * nj, 0.0);
        Vector<double> b_Hz   (std::size_t(ni) * nj * nz, 0.0);
        Vector<double> b_z_r  (std::size_t(ni) * nj * nz, 0.0);
        Vector<double> b_z_w  (std::size_t(ni) * nj * (nz + 1), 0.0);
        Vector<double> b_sustr(std::size_t(nis) * njs, 0.0);
        Vector<double> b_svstr(std::size_t(nis) * njs, 0.0);
        Vector<double> b_uwind(std::size_t(ni) * nj, 0.0);
        Vector<double> b_vwind(std::size_t(ni) * nj, 0.0);
        Vector<double> b_told (std::size_t(ni) * nj * nz * ntrc, 0.0);
        Vector<double> b_tnew (std::size_t(ni) * nj * nz * ntrc, 0.0);

        // ---- pack -------------------------------------------------------
        for (int j = jlo; j <= jhi; ++j) {
            for (int i = ilo; i <= ihi; ++i) {
                const int s = fidx(i, j, 0, 0, ilo, jlo, ni, nj, 1);
                b_rmask[s] = mskr_a(i, j, 0);
                b_srflx[s] = srflx_a(i, j, 0) * srflx_to_roms_kinematic;
            }
        }

        // ROMS k is 1-based over the valid column: ROMS k = REMORA k + 1.
        // This is the only index-base translation in the bridge.
        for (int k = 0; k < nz; ++k) {
            for (int j = jlo; j <= jhi; ++j) {
                for (int i = ilo; i <= ihi; ++i) {
                    const int s = fidx(i, j, k, 0, ilo, jlo, ni, nj, nz);
                    b_Hz [s] = Hz_a (i, j, k);
                    b_z_r[s] = z_r_a(i, j, k);
                }
            }
        }

        // z_w has no shift: ROMS z_w(0:N(ng)) and REMORA z_w k=0..N+1 both
        // index faces from the bottom, and both hold nz+1 values.
        for (int k = 0; k <= nz; ++k) {
            for (int j = jlo; j <= jhi; ++j) {
                for (int i = ilo; i <= ihi; ++i) {
                    b_z_w[fidx(i, j, k, 0, ilo, jlo, ni, nj, nz + 1)] =
                        z_w_a(i, j, k);
                }
            }
        }

        if (need_gas && !do_bulk_flux) {
            for (int j = jlo; j <= jhi + 1; ++j) {
                for (int i = ilo; i <= ihi + 1; ++i) {
                    const int s = fidx(i, j, 0, 0, ilo, jlo, nis, njs, 1);
                    b_sustr[s] = sustr_a(i, j, 0);
                    b_svstr[s] = svstr_a(i, j, 0);
                }
            }
        }

        if (need_gas && do_bulk_flux) {
            for (int j = jlo; j <= jhi; ++j) {
                for (int i = ilo; i <= ihi; ++i) {
                    const int s = fidx(i, j, 0, 0, ilo, jlo, ni, nj, 1);
                    b_uwind[s] = uwind_a(i, j, 0);
                    b_vwind[s] = vwind_a(i, j, 0);
                }
            }
        }

        // Buffer slot n is ROMS tracer n+1, which is REMORA component n.
        for (int n = 0; n < ntrc; ++n) {
            const int comp = n;
            for (int k = 0; k < nz; ++k) {
                for (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {
                        const int s =
                            fidx(i, j, k, n, ilo, jlo, ni, nj, nz);
                        b_told[s] = state_old(i, j, k, comp);
                        b_tnew[s] = state_new(i, j, k, comp);
                    }
                }
            }
        }

        // ---- call -------------------------------------------------------
        fennel_bridge_advance_c(
            ilo, ihi, jlo, jhi, nz, nbt,
            ilo, ihi, jlo, jhi,
            parms.po4 ? 1 : 0,
            parms.carbon ? 1 : 0,
            parms.oxygen ? 1 : 0,
            parms.odu ? 1 : 0,
            parms.denitrification ? 1 : 0,
            parms.bio_sediment ? 1 : 0,
            do_bulk_flux ? 1 : 0,
            parms.BioIter,
            double(dt_lev), double(solverChoice.rho0), double(Cp),
            double(parms.AttSW),   double(parms.AttChl), double(parms.PARfrac),
            double(parms.Vp0),     double(parms.I_thNH4), double(parms.D_p5NH4),
            double(parms.NitriR),  double(parms.K_NO3),  double(parms.K_NH4),
            double(parms.K_PO4),   double(parms.K_Phy),  double(parms.Chl2C_m),
            double(parms.ChlMin),  double(parms.PhyCN),  double(parms.R_P2N),
            double(parms.PhyIP),   double(parms.PhyIS),  double(parms.PhyMin),
            double(parms.PhyMR),   double(parms.ZooAE_N), double(parms.ZooCN),
            double(parms.ZooBM),   double(parms.ZooER),  double(parms.ZooGR),
            double(parms.ZooMin),  double(parms.ZooMR),  double(parms.LDeRRN),
            double(parms.LDeRRC),  double(parms.CoagR),  double(parms.SDeRRN),
            double(parms.SDeRRC),  double(parms.wPhy),   double(parms.wLDet),
            double(parms.wSDet),   double(parms.pCO2air),
            b_rmask.data(), b_Hz.data(), b_z_r.data(), b_z_w.data(),
            b_srflx.data(), b_sustr.data(), b_svstr.data(),
            b_uwind.data(), b_vwind.data(), b_pH.data(),
            b_told.data(), b_tnew.data(),
            biology_debug, biology_debug_i, biology_debug_j);

        // ---- unpack -----------------------------------------------------
        // biology_tile updates only the biology tracers, so only those are
        // written back; temperature and salinity in cons_new stay untouched.
        for (int n = Tracer_comp; n < ntrc; ++n) {
            const int comp = n;
            for (int k = 0; k < nz; ++k) {
                for (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {
                        state_out(i, j, k, comp) =
                            Real(b_tnew[fidx(i, j, k, n, ilo, jlo,
                                             ni, nj, nz)]);
                    }
                }
            }
        }
    }
#endif // AMREX_USE_GPU
}
