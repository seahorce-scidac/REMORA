#include <REMORA.H>

using namespace amrex;
using namespace REMORASeaIce;

/**
 * @param[in]  lev     level
 * @param[in]  name    output name (ROMS spelling, e.g. Aice, ice_thickness, Uice)
 * @param[out] mf      MultiFab holding it
 * @param[out] comp    its component
 * @param[out] stagger 0 for rho points, 1 for u points, 2 for v points
 */
bool
REMORA::ice_output_source (int lev, const std::string& name, MultiFab*& mf,
                           int& comp, int& stagger) const
{
    if (!has_seaice()) { return false; }

    const int lout = ice_lvl.liunw;

    const auto& rnames = r_state_nc_names();
    for (int n = 0; n < IceR::NumVars; ++n) {
        if (name == rnames[n]) {
            mf = vec_ice_si_r[lev].get(); comp = rcomp(n, lout); stagger = 0;
            return true;
        }
    }
    const auto& unames = u_state_nc_names();
    const auto& vnames = v_state_nc_names();
    for (int n = 0; n < IceU::NumVars; ++n) {
        if (name == unames[n]) {
            mf = vec_ice_si_u[lev].get(); comp = ucomp(n, lout); stagger = 1;
            return true;
        }
        if (name == vnames[n]) {
            mf = vec_ice_si_v[lev].get(); comp = ucomp(n, lout); stagger = 2;
            return true;
        }
    }
    const auto& fnames = r_field_nc_names();
    for (int n = 0; n < IceFR::NumVars; ++n) {
        if (!fnames[n].empty() && name == fnames[n]) {
            mf = vec_ice_fi_r[lev].get(); comp = n; stagger = 0;
            return true;
        }
    }
    return false;
}

/**
 * @param[in]    name      output name (ROMS spelling)
 * @param[inout] dst       one MultiFab per level, on the variable's staggering
 * @param[in]    dst_comp  component of dst to fill
 */
bool
REMORA::copy_ice_output (const std::string& name, Vector<MultiFab>& dst, int dst_comp) const
{
    MultiFab* mf = nullptr;
    int comp = 0;
    int stagger = 0;
    if (!ice_output_source(0, name, mf, comp, stagger)) { return false; }
    for (int lev = 0; lev <= finest_level; ++lev) {
        ice_output_source(lev, name, mf, comp, stagger);
        MultiFab::Copy(dst[lev], *mf, comp, dst_comp, 1, 0);
    }
    return true;
}

Vector<std::string>
REMORA::ice_output_names () const
{
    Vector<std::string> names;
    if (!has_seaice()) { return names; }
    for (const auto& nm : r_state_nc_names()) { names.push_back(nm); }
    for (const auto& nm : u_state_nc_names()) { names.push_back(nm); }
    for (const auto& nm : v_state_nc_names()) { names.push_back(nm); }
    for (const auto& nm : r_field_nc_names()) {
        if (!nm.empty()) { names.push_back(nm); }
    }
    return names;
}
