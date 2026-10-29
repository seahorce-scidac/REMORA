#ifndef REMORA_SEAICE_H_
#define REMORA_SEAICE_H_

#include <array>
#include <string>

#include <AMReX_Box.H>
#include <AMReX_Extension.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

/**
 * Sea-ice model, ported from the native myROMS ice model (ROMS/Nonlinear/SeaIce,
 * Budgell 2005; Hedstrom; Durski and Kurapov; rewritten by Arango).
 *
 * The state keeps ROMS's layout: every prognostic variable holds two time
 * levels, and the six level indices of mod_stepping (liold/linew for the ice
 * tracers, liuol/liunw for the ice velocity, lieol/lienw for the EVP substeps)
 * are carried over unchanged, including the order in which seaice.F swaps them.
 * The port reproduces ROMS's results, so the quirks that follow from that
 * order are kept on purpose; see the comments where they bite.
 */
namespace REMORASeaIce {

enum class IceModel {
    none,
    mellor_kantha
};

/** Ice/snow albedo formulation (ROMS ice_albedo.F branches). */
enum class IceAlbedoType {
    none,       //!< ICE_ALBEDO off: albedo_ice stays 0
    ec92,       //!< ICE_ALB_EC92, Ebert and Curry (1992)
    csim,       //!< ALBEDO_CSIM
    default_mk  //!< the default branch of ice_albedo.F
};

/** Ice lateral boundary condition types (the ROMS LBC keywords that are ported). Periodic
 *  sides are not a type here: they follow the geometry, and override the side's type. */
enum class IceBCType {
    closed,     //!< Clo
    gradient    //!< Gra
};

/**
 * Rho-point ice state variables. Each is stored at two time levels in the
 * rho-point state MultiFab, component `lvl*NumVars + var`.
 */
struct IceR {
    enum {
        Aice = 0,  //!< ice concentration
        Hice,      //!< ice thickness (mass per area / IceRho)
        Hmel,      //!< surface meltwater thickness
        Hsno,      //!< snow thickness
        Iage,      //!< ice age
        ISxx,      //!< internal ice xx-stress
        ISxy,      //!< internal ice xy-stress
        ISyy,      //!< internal ice yy-stress
        Tice,      //!< ice interior temperature
        Enth,      //!< ice/brine enthalpy
        Hage,      //!< thickness weighted by ice age
        NumVars
    };
};

/** U-point ice state variables (V-point ones mirror them). */
struct IceU {
    enum {
        Uice = 0,  //!< ice U-velocity
        Uevp,      //!< ice U-velocity within the EVP subcycle
        NumVars
    };
};

/**
 * Rho-point ice fields, a single time level. The first block is ROMS's
 * ICE%Fi array (rho-located entries); the second holds what bulk_flux hands
 * the ice model through FORCES in ROMS.
 */
struct IceFR {
    enum {
        Bvis = 0,  //!< ice bulk viscosity
        Hsse,      //!< sea surface elevation
        IOfv,      //!< ice-ocean friction velocity
        IOmf,      //!< ice-ocean mass flux
        IOmt,      //!< ice-ocean momentum transfer coefficient
        IOvs,      //!< ice-ocean velocity shear
        Isst,      //!< ice/snow surface temperature
        Pgrd,      //!< gridded ice strength (never set in ROMS)
        Pice,      //!< ice pressure or strength
        Qcon,      //!< ice/snow heat conductivity
        Qrhs,      //!< RHS heat flux over ice/snow
        Svis,      //!< ice shear viscosity
        S0mk,      //!< molecular sublayer salinity
        T0mk,      //!< molecular sublayer temperature
        Wdiv,      //!< ice divergence rate
        W_ai,      //!< melt/freeze rate at air/ice
        W_ao,      //!< melt/freeze rate at air/ocean
        W_fr,      //!< ice accretion rate by frazil
        W_io,      //!< melt/freeze rate at ice/ocean
        W_ro,      //!< melt/freeze runoff rate
        // Forcing from the bulk flux (ROMS FORCES arrays)
        Qnet_ai,   //!< net heat flux, ice to atmosphere [W/m2]
        Qnet_ao,   //!< net heat flux, ocean to atmosphere [W/m2]
        snow,      //!< snowfall [m/s]
        srflx_ice, //!< shortwave radiation reaching the ice [W/m2]
        albedo_ice,//!< ice/snow albedo
        sustr_ai,  //!< air-ice U-stress [N/m2], rho points
        svstr_ai,  //!< air-ice V-stress [N/m2], rho points
        NumVars
    };
};

/** U-point ice fields (V-point ones mirror them). */
struct IceFU {
    enum {
        AIus = 0,  //!< air-ice U-stress
        Uavg,      //!< mixed-layer averaged ocean U-velocity
        sustr_ao,  //!< air-ocean U-stress (the bulk-flux sustr)
        NumVars
    };
};

/** Ghost cells on every ice MultiFab: MPDATA reaches two cells out, plus one. */
constexpr int ice_ngrow = 3;

/** The two time levels of a state variable, as components of its MultiFab. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rcomp (int var, int lvl) noexcept { return lvl * IceR::NumVars + var; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int ucomp (int var, int lvl) noexcept { return lvl * IceU::NumVars + var; }

/** The valid faces of a nodal box minus those on a non-periodic domain edge (ROMS IstrP..Iend). */
inline amrex::Box interior_faces (const amrex::Box& nodal_bx, int dir, const amrex::Geometry& geom)
{
    amrex::Box b = nodal_bx;
    if (!geom.isPeriodic(dir)) {
        const amrex::Box& dom = geom.Domain();
        b.setSmall(dir, std::max(b.smallEnd(dir), dom.smallEnd(dir)+1));
        b.setBig  (dir, std::min(b.bigEnd(dir),   dom.bigEnd(dir)));
    }
    return b;
}

/**
 * Time-level indices, 0-based versions of ROMS's mod_stepping ones. ROMS
 * starts all six at 1, so all six start at 0 here.
 */
struct IceLevels {
    int liold = 0;
    int linew = 0;
    int liuol = 0;
    int liunw = 0;
    int lieol = 0;
    int lienw = 0;
};

/** Runtime parameters: ROMS's ice.in plus the CPP options of the ice model. */
struct SeaIceParams {
    IceModel model = IceModel::none;

    // CPP options of the ice model, as runtime switches
    bool do_thermo = true;        //!< ICE_THERMO / ICE_MK
    bool do_momentum = true;      //!< ICE_MOMENTUM / ICE_EVP
    bool do_advect = true;        //!< ICE_ADVECT / ICE_SMOLAR
    bool advect_upwind = false;   //!< ICE_UPWIND: no antidiffusive MPDATA pass
    bool mom_bulk = false;        //!< ICE_MOM_BULK: quadratic ice-water drag
    bool convsnow = false;        //!< ICE_CONVSNOW
    /** Not in ROMS: copy the ice tracers from liold to linew at the start of a step
     *  when advection is off. ROMS does not, so its thermodynamics then works on the
     *  state of two steps back; this matches the patched ROMS used to validate
     *  thermodynamics alone. Off by default. */
    bool copy_old_to_new = false;
    IceAlbedoType albedo_type = IceAlbedoType::ec92;

    // ice.in
    int nEVP = 60;
    amrex::Real AirRho = amrex::Real(1.267);
    amrex::Real IceRho = amrex::Real(900.0);
    amrex::Real SnowDryRho = amrex::Real(330.0);
    amrex::Real SnowWetRho = amrex::Real(450.0);
    amrex::Real Cd_ai = amrex::Real(3.0e-3);
    amrex::Real Cd_io = amrex::Real(7.5e-3);
    amrex::Real Astrength = amrex::Real(20.0);
    amrex::Real Pstar = amrex::Real(2.5e4);
    amrex::Real zetaMin = amrex::Real(0.0);
    amrex::Real zetaMax = amrex::Real(1.0e14);
    amrex::Real ellip_sq = amrex::Real(4.0);
    amrex::Real min_ai = amrex::Real(0.0);
    amrex::Real max_ai = amrex::Real(1.0);
    amrex::Real min_hi = amrex::Real(0.0);
    amrex::Real max_hmelt = amrex::Real(0.1);
    amrex::Real stressAng = amrex::Real(0.0);   //!< radians; read but unused, as in ROMS
    amrex::Real ice_emiss = amrex::Real(0.97);
    amrex::Real spec_heat_air = amrex::Real(1004.0);
    amrex::Real trans_coeff = amrex::Real(1.75e-3);
    amrex::Real sublimation = amrex::Real(2.834e6);
    amrex::Real gamma2 = amrex::Real(1.0);      //!< slipperiness for closed ice BCs (ROMS GAMMA2)

    /** Lateral boundary type per side (xlo, ylo, xhi, yhi), one set for every variable. */
    std::array<IceBCType, 4> bc_type {IceBCType::closed, IceBCType::closed,
                                      IceBCType::closed, IceBCType::closed};

    /** Read the `ice.*` parameters. */
    void init_params ();
};

/** Parse remora.ice_model. */
IceModel parse_ice_model (const std::string& name);

/** The input name of an ice model. */
std::string ice_model_name (IceModel model);

AMREX_FORCE_INLINE bool has_seaice (IceModel model) { return model != IceModel::none; }

/** ROMS NetCDF names (varinfo.yaml) of the state variables. */
const amrex::Vector<std::string>& r_state_nc_names ();
const amrex::Vector<std::string>& u_state_nc_names ();
const amrex::Vector<std::string>& v_state_nc_names ();

/** ROMS NetCDF names of the rho-point ice fields that ROMS writes; empty where it writes none. */
const amrex::Vector<std::string>& r_field_nc_names ();

} // namespace REMORASeaIce

#endif
