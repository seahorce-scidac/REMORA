#include <REMORA.H>

using namespace amrex;
using namespace REMORASeaIce;

#ifdef REMORA_USE_NETCDF
#include "REMORA_NCFile.H"
#endif

/**
 * Allocate the ice state and fields of a level and zero them, as ROMS's
 * initialize_ice does (IniVal = 0 for every state variable and field).
 *
 * @param[in] lev level to operate on
 * @param[in] ba  BoxArray of the (3D) level
 * @param[in] dm  DistributionMapping of the level
 */
void
REMORA::init_seaice (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

    const IntVect ng(ice_ngrow, ice_ngrow, 0);

    vec_ice_si_r[lev] = std::make_unique<MultiFab>(ba2d, dm, 2*IceR::NumVars, ng);
    vec_ice_si_u[lev] = std::make_unique<MultiFab>(convert(ba2d,IntVect(1,0,0)), dm, 2*IceU::NumVars, ng);
    vec_ice_si_v[lev] = std::make_unique<MultiFab>(convert(ba2d,IntVect(0,1,0)), dm, 2*IceU::NumVars, ng);
    vec_ice_fi_r[lev] = std::make_unique<MultiFab>(ba2d, dm, IceFR::NumVars, ng);
    vec_ice_fi_u[lev] = std::make_unique<MultiFab>(convert(ba2d,IntVect(1,0,0)), dm, IceFU::NumVars, ng);
    vec_ice_fi_v[lev] = std::make_unique<MultiFab>(convert(ba2d,IntVect(0,1,0)), dm, IceFU::NumVars, ng);

    vec_ice_si_r[lev]->setVal(zero);
    vec_ice_si_u[lev]->setVal(zero);
    vec_ice_si_v[lev]->setVal(zero);
    vec_ice_fi_r[lev]->setVal(zero);
    vec_ice_fi_u[lev]->setVal(zero);
    vec_ice_fi_v[lev]->setVal(zero);

    ice_lvl = IceLevels{};
}

#ifdef REMORA_USE_NETCDF
/**
 * Read the ice state from the initial file into time level 1 (index 0), as
 * ROMS's get_state does: each state variable the file carries is read, the
 * rest stay zero. The ice fields (T0mk, S0mk, IOmt, ...) are not read, also as
 * in ROMS, so they start at zero even when restarting from a ROMS file.
 *
 * Uevp and Vevp are read at u and v points. ROMS reads them as rho-point
 * fields, which fails on its own restart files; the values do not matter
 * since the first EVP substep overwrites them from Uice and Vice.
 *
 * @param[in] lev level to operate on
 */
void
REMORA::init_seaice_from_netcdf (int lev)
{
    const std::string& fname = nc_init_file[lev][0];
    if (num_boxes_at_level[lev] != 1) {
        amrex::Abort("init_seaice_from_netcdf: only a single initial-data box is supported");
    }
    const Box& domain = boxes_at_level[lev][0];

    auto read_into = [&] (MultiFab& mf, int comp, const std::string& nc_name)
    {
        if (!QueryNetCDFHasVars(fname, {nc_name})) {
            return false;
        }
        FArrayBox nc_fab;
        Vector<FArrayBox*> fabs {&nc_fab};
        Vector<std::string> names {nc_name};
        Vector<NC_Data_Dims_Type> dims {NC_Data_Dims_Type::Time_SN_WE};
        BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, names, dims, fabs);

        for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
            // Copy on intersection, valid and ghost cells alike: the file carries
            // ROMS's boundary rows, which sit in REMORA's first ghost row.
            mf[mfi].template copy<RunOn::Device>(nc_fab, 0, comp, 1);
        }
        return true;
    };

    Vector<std::string> found;
    Vector<std::string> missing;
    auto note = [&] (bool ok, const std::string& nm) { (ok ? found : missing).push_back(nm); };

    const auto& rnames = r_state_nc_names();
    for (int n = 0; n < IceR::NumVars; ++n) {
        note(read_into(*vec_ice_si_r[lev], rcomp(n, 0), rnames[n]), rnames[n]);
    }
    const auto& unames = u_state_nc_names();
    const auto& vnames = v_state_nc_names();
    for (int n = 0; n < IceU::NumVars; ++n) {
        note(read_into(*vec_ice_si_u[lev], ucomp(n, 0), unames[n]), unames[n]);
        note(read_into(*vec_ice_si_v[lev], ucomp(n, 0), vnames[n]), vnames[n]);
    }

    amrex::Print() << "Sea ice initial state from " << fname << ":";
    for (const auto& nm : found) { amrex::Print() << " " << nm; }
    amrex::Print() << std::endl;
    if (!missing.empty()) {
        amrex::Print() << "  not in the file, starting at zero:";
        for (const auto& nm : missing) { amrex::Print() << " " << nm; }
        amrex::Print() << std::endl;
    }

    // Masked, as ROMS's nf_fread2d does under MASKING
    const int nr = vec_ice_si_r[lev]->nComp();
    const int nu = vec_ice_si_u[lev]->nComp();
    for (MFIter mfi(*vec_ice_si_r[lev], false); mfi.isValid(); ++mfi) {
        Array4<Real> const& sr = vec_ice_si_r[lev]->array(mfi);
        Array4<Real> const& su = vec_ice_si_u[lev]->array(mfi);
        Array4<Real> const& sv = vec_ice_si_v[lev]->array(mfi);
        Array4<Real const> const& mskr = vec_mskr[lev]->const_array(mfi);
        Array4<Real const> const& msku = vec_msku[lev]->const_array(mfi);
        Array4<Real const> const& mskv = vec_mskv[lev]->const_array(mfi);
        const Box rbx = mfi.growntilebox(IntVect(1,1,0));
        const Box ubx = mfi.grownnodaltilebox(0, IntVect(1,1,0));
        const Box vbx = mfi.grownnodaltilebox(1, IntVect(1,1,0));
        ParallelFor(rbx, nr, [=] AMREX_GPU_DEVICE (int i, int j, int, int n)
        {
            sr(i,j,0,n) *= mskr(i,j,0);
        });
        ParallelFor(ubx, nu, [=] AMREX_GPU_DEVICE (int i, int j, int, int n)
        {
            su(i,j,0,n) *= msku(i,j,0);
        });
        ParallelFor(vbx, nu, [=] AMREX_GPU_DEVICE (int i, int j, int, int n)
        {
            sv(i,j,0,n) *= mskv(i,j,0);
        });
    }

    ice_fill_state(lev);
}
#endif
