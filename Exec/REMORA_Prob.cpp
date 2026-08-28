#include "REMORA_Prob.H"
#include "REMORA_DepthStretchTransform.H"

using namespace amrex;

#ifdef ERF_REMORA_FORCE_PROBINIT_LINK
// Force archive extraction of this TU when REMORA is linked as a static library
// inside a parent coupled executable and amrex_probinit is weak.
void remora_probinit_link_anchor_func () noexcept {}
#endif

std::unique_ptr<ProblemBase>
amrex_probinit (const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem(const amrex::Real* /*problo*/, const amrex::Real* /*probhi*/)
{}

/**
 * \brief Initializes bathymetry h and surface height Zeta
 */
void Problem::init_analytic_bathymetry (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_h)
{
    ParmParse pp("remora");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    auto geomdata = geom.data();
    bool EWPeriodic = geomdata.isPeriodic(0);
    bool NSPeriodic = geomdata.isPeriodic(1);

    if (my_prob_name_ci == "advection") {
#include "Prob/REMORA_InitAnalyticBathymetry_Advection.H"

    } else if ( (my_prob_name_ci == "blankproblem")   ||
                (my_prob_name_ci == "dogbone")        ||
                (my_prob_name_ci == "idealminigrid")  ||
                (my_prob_name_ci == "idealminiriv") ) {
// No initialization of bathymetry occurs with these prob_name's

    } else if (my_prob_name_ci == "coupletoerf") {
#include "Prob/REMORA_InitAnalyticBathymetry_CoupleToERF.H"

    } else if (my_prob_name_ci == "boundarylayer") {
#include "Prob/REMORA_InitAnalyticBathymetry_BoundaryLayer.H"

    } else if (my_prob_name_ci == "channeltest") {
#include "Prob/REMORA_InitAnalyticBathymetry_ChannelTest.H"

    } else if (my_prob_name_ci == "dogboneanalytic") {
#include "Prob/REMORA_InitAnalyticBathymetry_DogboneAnalytic.H"

    } else if (my_prob_name_ci == "doublegyre") {
        mf_h.setVal(Real(500.0));

    } else if (my_prob_name_ci == "doublyperiodic") {
#include "Prob/REMORA_InitAnalyticBathymetry_DoublyPeriodic.H"

    } else if (my_prob_name_ci == "particleadvectionflat") {
#include "Prob/REMORA_InitAnalyticBathymetry_ParticleAdvectionFlat.H"

    } else if (my_prob_name_ci == "seamount") {
#include "Prob/REMORA_InitAnalyticBathymetry_Seamount.H"

    } else if (my_prob_name_ci == "upwelling") {
#include "Prob/REMORA_InitAnalyticBathymetry_Upwelling.H"

    } else if (my_prob_name_ci == "upwellingcoupling") {
#include "Prob/REMORA_InitAnalyticBathymetry_Upwelling_Coupling.H"

    } else if (my_prob_name_ci == "upwellingml") {
#include "Prob/REMORA_InitAnalyticBathymetry_Upwelling_ML.H"
    }

    amrex::Gpu::streamSynchronize();
}

/**
 * \brief Initializes custom sea surface height
 */
void Problem::init_analytic_zeta (
        int lev, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        MultiFab& mf_zeta)
{
    ParmParse pp("remora");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    // This is currently the only case that does anything but set mf_zeta to 0
    if (my_prob_name_ci == "dogboneanalytic") {
#include "Prob/REMORA_InitAnalyticZeta_DogboneAnalytic.H"

    } else {
        mf_zeta.setVal(zero);
    }

    Gpu::streamSynchronize();
}

void Problem::init_analytic_prob(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_cons,
        amrex::MultiFab& mf_xvel,
        amrex::MultiFab& mf_yvel)
{
    ParmParse pp("remora");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    // Initialize to zero to ensure any otherwise uninitialized tracers are 0
    mf_cons.setVal(zero);

    // This is currently the only case that does anything but set mf_zeta to 0
    if (my_prob_name_ci == "advection") {
#include "Prob/REMORA_InitAnalyticProb_Advection.H"

    } else if (my_prob_name_ci == "boundarylayer") {
#include "Prob/REMORA_InitAnalyticProb_BoundaryLayer.H"

    } else if (my_prob_name_ci == "channeltest") {
#include "Prob/REMORA_InitAnalyticProb_ChannelTest.H"

    } else if (my_prob_name_ci == "coupletoerf") {
#include "Prob/REMORA_InitAnalyticProb_CoupleToERF.H"

    } else if (my_prob_name_ci == "dogboneanalytic") {
#include "Prob/REMORA_InitAnalyticProb_DogboneAnalytic.H"

    } else if (my_prob_name_ci == "doublegyre") {
#include "Prob/REMORA_InitAnalyticProb_DoubleGyre.H"

    } else if (my_prob_name_ci == "doublyperiodic") {
#include "Prob/REMORA_InitAnalyticProb_DoublyPeriodic.H"

    } else if (my_prob_name_ci == "particleadvectionflat") {
#include "Prob/REMORA_InitAnalyticProb_ParticleAdvectionFlat.H"

    } else if (my_prob_name_ci == "seamount") {
#include "Prob/REMORA_InitAnalyticProb_Seamount.H"

    } else if (my_prob_name_ci == "upwelling") {
#include "Prob/REMORA_InitAnalyticProb_Upwelling.H"

    } else if (my_prob_name_ci == "upwellingcoupling") {
#include "Prob/REMORA_InitAnalyticProb_Upwelling_Coupling.H"

    } else if (my_prob_name_ci == "upwellingml") {
#include "Prob/REMORA_InitAnalyticProb_Upwelling_ML.H"

    }

    Gpu::streamSynchronize();
}

void Problem::init_analytic_vmix(
        int lev,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        MultiFab& mf_Akv, MultiFab& mf_Akt)
{
    ParmParse pp("remora");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if ( (my_prob_name_ci == "dogbone") ||
         (my_prob_name_ci == "dogboneanalytic") ) {
        mf_Akv.setVal(m_solverChoice.Akv_bak);
        // Akt_bak has one value per active tracer, so leave the per-component values REMORA
        // already set rather than flattening temperature and salinity to a single number.

    } else if (my_prob_name_ci == "advection") {
#include "Prob/REMORA_InitAnalyticVMix_Advection.H"

    } else if (my_prob_name_ci == "boundarylayer") {
#include "Prob/REMORA_InitAnalyticVMix_BoundaryLayer.H"

    } else if (my_prob_name_ci == "channeltest") {
#include "Prob/REMORA_InitAnalyticVMix_ChannelTest.H"

    } else if (my_prob_name_ci == "doublegyre") {
#include "Prob/REMORA_InitAnalyticVMix_DoubleGyre.H"

    } else if (my_prob_name_ci == "doublyperiodic") {
#include "Prob/REMORA_InitAnalyticVMix_DoublyPeriodic.H"

    } else if (my_prob_name_ci == "idealminigrid") {
#include "Prob/REMORA_InitAnalyticVMix_IdealMiniGrid.H"

    } else if (my_prob_name_ci == "idealrivgrid") {
#include "Prob/REMORA_InitAnalyticVMix_IdealRivGrid.H"

    } else if (my_prob_name_ci == "particleadvectionflat") {
#include "Prob/REMORA_InitAnalyticVMix_ParticleAdvectionFlat.H"

    } else if (my_prob_name_ci == "seamount") {
#include "Prob/REMORA_InitAnalyticVMix_Seamount.H"

    } else if ( (my_prob_name_ci == "upwelling") ||
                (my_prob_name_ci == "upwellingcoupling") ) {
#include "Prob/REMORA_InitAnalyticVMix_Upwelling.H"

    } else if (my_prob_name_ci == "upwellingml") {
#include "Prob/REMORA_InitAnalyticVMix_Upwelling_ML.H"
    }

    Gpu::streamSynchronize();
}

void Problem::init_analytic_hmix(
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        MultiFab& mf_visc2_p,
        MultiFab& mf_visc2_r,
        MultiFab& mf_diff2)
{
    ParmParse pp("remora");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if ( (my_prob_name_ci == "dogbone") ||
         (my_prob_name_ci == "dogboneanalytic") ) {
        mf_visc2_p.setVal(zero);
        mf_visc2_r.setVal(zero);
        mf_diff2.setVal(zero);

    } else if (my_prob_name_ci == "advection") {
#include "Prob/REMORA_InitAnalyticHMix_Advection.H"

    } else if (my_prob_name_ci == "boundarylayer") {
#include "Prob/REMORA_InitAnalyticHMix_BoundaryLayer.H"

    } else if (my_prob_name_ci == "channeltest") {
#include "Prob/REMORA_InitAnalyticHMix_ChannelTest.H"

    } else if (my_prob_name_ci == "doublegyre") {
#include "Prob/REMORA_InitAnalyticHMix_DoubleGyre.H"

    } else if (my_prob_name_ci == "doublyperiodic") {
#include "Prob/REMORA_InitAnalyticHMix_DoublyPeriodic.H"

    } else if (my_prob_name_ci == "idealminigrid") {
#include "Prob/REMORA_InitAnalyticHMix_IdealMiniGrid.H"

    } else if (my_prob_name_ci == "idealrivgrid") {
#include "Prob/REMORA_InitAnalyticHMix_IdealRivGrid.H"

    } else if (my_prob_name_ci == "particleadvectionflat") {
#include "Prob/REMORA_InitAnalyticHMix_ParticleAdvectionFlat.H"

    } else if (my_prob_name_ci == "seamount") {
#include "Prob/REMORA_InitAnalyticHMix_Seamount.H"

    } else if ( (my_prob_name_ci == "upwelling") ||
                (my_prob_name_ci == "upwellingcoupling") ) {
#include "Prob/REMORA_InitAnalyticHMix_Upwelling.H"

    } else if (my_prob_name_ci == "upwellingml") {
#include "Prob/REMORA_InitAnalyticHMix_Upwelling_ML.H"

    }

    Gpu::streamSynchronize();
}

void Problem::init_analytic_smflux(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        MultiFab& mf_sustr, MultiFab& mf_svstr)
{
    ParmParse pp("remora");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    auto geomdata = geom.data();
    bool EWPeriodic = geomdata.isPeriodic(0);
    bool NSPeriodic = geomdata.isPeriodic(1);

    if ((my_prob_name_ci == "advection")       ||
        (my_prob_name_ci == "dogbone")         ||
        (my_prob_name_ci == "dogboneanalytic") ||
        (my_prob_name_ci == "idealminigrid")   ||
        (my_prob_name_ci == "idealrivgrid")    ||
        (my_prob_name_ci == "seamount")        ||
        (my_prob_name_ci == "upwellingcoupling") ) {
        mf_sustr.setVal(zero);
        mf_svstr.setVal(zero);

    } else if (my_prob_name_ci == "channeltest") {
#include "Prob/REMORA_InitAnalyticSMFlux_ChannelTest.H"

    } else if (my_prob_name_ci == "doublegyre") {
#include "Prob/REMORA_InitAnalyticSMFlux_DoubleGyre.H"

    } else if (my_prob_name_ci == "doublyperiodic") {
#include "Prob/REMORA_InitAnalyticSMFlux_DoublyPeriodic.H"

    } else if (my_prob_name_ci == "particleadvectionflat") {
#include "Prob/REMORA_InitAnalyticSMFlux_ParticleAdvectionFlat.H"
    } else if (my_prob_name_ci == "upwelling") {
#include "Prob/REMORA_InitAnalyticSMFlux_Upwelling.H"

    } else if (my_prob_name_ci == "upwellingml") {
#include "Prob/REMORA_InitAnalyticSMFlux_UpwellingML.H"

    }

    Gpu::streamSynchronize();
}

void Problem::init_analytic_grid_scale (
        int /*lev*/, const amrex::Geometry& geom,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& mf_pm, amrex::MultiFab& mf_pn)
{
    ParmParse pp("remora");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if (my_prob_name_ci == "boundarylayer") {
#include "Prob/REMORA_InitAnalyticGridScale_BoundaryLayer.H"
    }

    Gpu::streamSynchronize();
}

void Problem::init_analytic_masks(
        int lev,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        MultiFab& mf_mskr)
{
    ParmParse pp("remora");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if (my_prob_name_ci == "dogboneanalytic") {
#include "Prob/REMORA_InitAnalyticMasks_DogboneAnalytic.H"
    }

    Gpu::streamSynchronize();
}

void Problem::init_analytic_coriolis (
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& /*mf_fcor*/)
{
}

void Problem::init_analytic_surface_var (
        int lev,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        amrex::MultiFab& mf_Uwind,
        amrex::MultiFab& mf_Vwind,
        amrex::MultiFab& /*mf_Tair*/,
        amrex::MultiFab& /*mf_qair*/,
        amrex::MultiFab& /*mf_Pair*/,
        amrex::MultiFab& /*mf_srflx*/,
        amrex::MultiFab& /*mf_longwave_down*/,
        amrex::MultiFab& /*mf_rain*/,
        amrex::MultiFab& /*mf_cloud*/,
        amrex::MultiFab& /*mf_EminusP*/)
{
    ParmParse pp("remora");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if ( (my_prob_name_ci == "coupletoerf")   ||
         (my_prob_name_ci == "dogbone")     ) {
        mf_Uwind.setVal(zero);
        mf_Vwind.setVal(zero);

    } else if (my_prob_name_ci == "boundarylayer")  {
#include "Prob/REMORA_InitAnalyticSurfaceVar_BoundaryLayer.H"
    }

    Gpu::streamSynchronize();
}

void Problem::init_analytic_biology (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORABiology::FennelParameters const& m_fennel_params,
        REMORA const& remora,
        amrex::MultiFab& mf_cons)
{
    ParmParse pp("remora");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);
    const auto bio_comp = REMORABiology::Fennel::components(m_fennel_params, remora.bio_comp_start());
    // Upwelling reuses BioToy's profile, which is ROMS ana_biology.h: it is a function of
    // temperature alone, so it carries over to any problem that stratifies temperature.
    if ((my_prob_name_ci == "biotoy") || (my_prob_name_ci == "upwelling")) {
#include "Prob/REMORA_InitAnalyticBiology_BioToy.H"

    } else {
        // Falling through would leave the biology tracers at whatever cons_new was zeroed
        // to, and a Fennel run starting from all-zero tracers looks plausible rather than
        // broken. The base-class hook in REMORA_prob_common.H errors for the same reason;
        // this override has to do so too, or that guard is dead.
        amrex::Error("remora.prob_name = " + my_prob_name + " has no analytic biology initial"
                     " condition. Add one to Source/Prob/ and a branch for it in"
                     " Problem::init_analytic_biology, or set remora.biology_ic_type = netcdf.");
    }
}

