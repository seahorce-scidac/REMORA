#include <REMORA_SeaIce.H>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_String.H>

#include <cmath>

using namespace amrex;

namespace REMORASeaIce {

IceModel parse_ice_model (const std::string& name)
{
    const std::string lc = amrex::toLower(name);
    if (lc == "none") {
        return IceModel::none;
    } else if (lc == "mellor_kantha" || lc == "mk") {
        return IceModel::mellor_kantha;
    }
    amrex::Abort("remora.ice_model = " + name + " is not recognized; use none or mellor_kantha");
    return IceModel::none;
}

std::string ice_model_name (IceModel model)
{
    switch (model) {
        case IceModel::mellor_kantha: return "mellor_kantha";
        default:                      return "none";
    }
}

namespace {

IceAlbedoType parse_albedo_type (const std::string& name)
{
    const std::string lc = amrex::toLower(name);
    if (lc == "none")    { return IceAlbedoType::none; }
    if (lc == "ec92")    { return IceAlbedoType::ec92; }
    if (lc == "csim")    { return IceAlbedoType::csim; }
    if (lc == "default") { return IceAlbedoType::default_mk; }
    amrex::Abort("ice.albedo_type = " + name + " is not recognized; use none, ec92, csim or default");
    return IceAlbedoType::none;
}

IceBCType parse_bc_type (const std::string& name)
{
    const std::string lc = amrex::toLower(name);
    if (lc == "closed")   { return IceBCType::closed; }
    if (lc == "gradient") { return IceBCType::gradient; }
    amrex::Abort("ice.bc_type entry " + name + " is not recognized; use closed or gradient"
                 " (periodic sides follow remora.is_periodic)");
    return IceBCType::closed;
}

} // namespace

void SeaIceParams::init_params ()
{
    ParmParse pp("ice");

    pp.queryAdd("do_thermo", do_thermo);
    pp.queryAdd("do_momentum", do_momentum);
    pp.queryAdd("do_advect", do_advect);
    pp.queryAdd("advect_upwind", advect_upwind);
    pp.queryAdd("mom_bulk", mom_bulk);
    pp.queryAdd("convsnow", convsnow);
    pp.queryAdd("copy_old_to_new", copy_old_to_new);

    std::string albedo_string = "ec92";
    pp.queryAdd("albedo_type", albedo_string);
    albedo_type = parse_albedo_type(albedo_string);
    if (albedo_type == IceAlbedoType::csim) {
        // ice_albedo.F's ALBEDO_CSIM branch uses fh where only fhi is computed, so ROMS
        // cannot compile it; there is no ROMS answer to reproduce.
        amrex::Abort("ice.albedo_type = csim is not available: the ROMS ALBEDO_CSIM branch"
                     " does not compile (it reads an undefined fh)");
    }

    pp.queryAdd("nEVP", nEVP);
    pp.queryAdd("AirRho", AirRho);
    pp.queryAdd("IceRho", IceRho);
    pp.queryAdd("SnowDryRho", SnowDryRho);
    pp.queryAdd("SnowWetRho", SnowWetRho);
    pp.queryAdd("Cd_ai", Cd_ai);
    pp.queryAdd("Cd_io", Cd_io);
    pp.queryAdd("Astrength", Astrength);
    pp.queryAdd("Pstar", Pstar);
    pp.queryAdd("zetaMin", zetaMin);
    pp.queryAdd("zetaMax", zetaMax);
    pp.queryAdd("ellip_sq", ellip_sq);
    pp.queryAdd("min_ai", min_ai);
    pp.queryAdd("max_ai", max_ai);
    pp.queryAdd("min_hi", min_hi);
    pp.queryAdd("max_hmelt", max_hmelt);

    // Degrees in the input, radians internally, as in ice_inp.h. ROMS reads it and never
    // uses it; neither does this port.
    Real stressAng_deg = stressAng * Real(180.0) / Real(M_PI);
    pp.queryAdd("stressAng", stressAng_deg);
    stressAng = stressAng_deg * Real(M_PI) / Real(180.0);

    pp.queryAdd("ice_emiss", ice_emiss);
    pp.queryAdd("spec_heat_air", spec_heat_air);
    pp.queryAdd("trans_coeff", trans_coeff);
    pp.queryAdd("sublimation", sublimation);
    pp.queryAdd("gamma2", gamma2);

    if (nEVP < 1) {
        amrex::Abort("ice.nEVP must be at least 1");
    }

    // Optional: closed on every side unless given
    Vector<std::string> bc_strings(4, "closed");
    pp.queryarr("bc_type", bc_strings);
    if (bc_strings.size() != 4) {
        amrex::Abort("ice.bc_type needs four entries: xlo ylo xhi yhi");
    }
    for (int n = 0; n < 4; ++n) {
        bc_type[n] = parse_bc_type(bc_strings[n]);
    }
}

const Vector<std::string>& r_state_nc_names ()
{
    static const Vector<std::string> names = {
        "Aice", "ice_thickness", "meltpond_thickness", "snow_thickness", "ice_age",
        "ice_Sxx", "ice_Sxy", "ice_Syy", "Tice", "ice_enthalpy", "ice_age_thickness"
    };
    AMREX_ASSERT(names.size() == IceR::NumVars);
    return names;
}

const Vector<std::string>& u_state_nc_names ()
{
    static const Vector<std::string> names = { "Uice", "Uevp" };
    return names;
}

const Vector<std::string>& v_state_nc_names ()
{
    static const Vector<std::string> names = { "Vice", "Vevp" };
    return names;
}

const Vector<std::string>& r_field_nc_names ()
{
    // varinfo.yaml names of the Fi entries ROMS can write; the forcing block has none.
    static const Vector<std::string> names = {
        "",                           // Bvis
        "",                           // Hsse
        "io_friction_vel",            // IOfv
        "io_massflux",                // IOmf
        "io_transfer_coef",           // IOmt
        "",                           // IOvs
        "ice_sst",                    // Isst
        "",                           // Pgrd
        "",                           // Pice
        "ice_snow_heat_conductivity", // Qcon
        "ice_snow_heat_flux",         // Qrhs
        "",                           // Svis
        "under_ice_salt",             // S0mk
        "under_ice_temp",             // T0mk
        "ice_divergence",             // Wdiv
        "ai_melt_freeze",             // W_ai
        "ao_melt_freeze",             // W_ao
        "frazil_freeze",              // W_fr
        "io_melt_freeze",             // W_io
        "melt_freeze_runoff",         // W_ro
        "", "", "", "", "", "", ""    // forcing block
    };
    AMREX_ASSERT(names.size() == IceFR::NumVars);
    return names;
}

} // namespace REMORASeaIce
