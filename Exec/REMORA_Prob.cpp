#include "REMORA_Prob.H"

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
                (my_prob_name_ci == "coupletoerf")    ||
                (my_prob_name_ci == "dogbone")        ||
                (my_prob_name_ci == "idealminigrid")  ||
                (my_prob_name_ci == "idealminiriv") ) {
// No initialization of bathymetry occurs with these prob_name's

    } else if (my_prob_name_ci == "boundarylayer") {
#include "Prob/REMORA_InitAnalyticBathymetry_BoundaryLayer.H"

    } else if (my_prob_name_ci == "channeltest") {
#include "Prob/REMORA_InitAnalyticBathymetry_ChannelTest.H"

    } else if (my_prob_name_ci == "dogboneanalytic") {
#include "Prob/REMORA_InitAnalyticBathymetry_DogboneAnalytic.H"

    } else if (my_prob_name_ci == "doublegyre") {
        mf_h.setVal(500.0_rt);

    } else if (my_prob_name_ci == "doublyperiodic") {
#include "Prob/REMORA_InitAnalyticBathymetry_DoublyPeriodic.H"

    } else if (my_prob_name_ci == "particles") {
#include "Prob/REMORA_InitAnalyticBathymetry_ParticlesOverSeamount.H"

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
        mf_zeta.setVal(0.0_rt);
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

    } else if (my_prob_name_ci == "particlesoverseamount") {
#include "Prob/REMORA_InitAnalyticProb_ParticlesOverSeamount.H"

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
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& m_solverChoice,
        REMORA const& /*remora*/,
        MultiFab& mf_Akv, MultiFab& mf_Akt)
{
    ParmParse pp("remora");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if ( (my_prob_name_ci == "dogbone") ||
         (my_prob_name_ci == "dogboneanalytic") ) {
        mf_Akv.setVal(m_solverChoice.Akv_bak);
        mf_Akt.setVal(m_solverChoice.Akt_bak);
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
        mf_visc2_p.setVal(0.0_rt);
        mf_visc2_r.setVal(0.0_rt);
        mf_diff2.setVal(0.0_rt);

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

    } else if (my_prob_name_ci == "particlesoverseamount") {
#include "Prob/REMORA_InitAnalyticHMix_ParticlesOverSeamount.H"

    } else if (my_prob_name_ci == "seamount") {
#include "Prob/REMORA_InitAnalyticHMix_Seamount.H"

    } else if (my_prob_name_ci == "upwelling") {
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
        (my_prob_name_ci == "seamount") ) {
        mf_sustr.setVal(0.0);
        mf_svstr.setVal(0.0);

    } else if (my_prob_name_ci == "channeltest") {
#include "Prob/REMORA_InitAnalyticSMFlux_ChannelTest.H"

    } else if (my_prob_name_ci == "doublegyre") {
#include "Prob/REMORA_InitAnalyticSMFlux_DoubleGyre.H"

    } else if (my_prob_name_ci == "doublyperiodic") {
#include "Prob/REMORA_InitAnalyticSMFlux_DoublyPeriodic.H"

    } else if (my_prob_name_ci == "particlesoverseamount") {
#include "Prob/REMORA_InitAnalyticSMFlux_ParticlesOverSeamount.H"

    } else if (my_prob_name_ci == "upwelling") {
#include "Prob/REMORA_InitAnalyticSMFlux_Upwelling.H"

    } else if (my_prob_name_ci == "upwellingcoupling") {
#include "Prob/REMORA_InitAnalyticSMFlux_UpwellingCoupling.H"

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

void Problem::init_analytic_wind(
        int lev,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        amrex::MultiFab& mf_Uwind, amrex::MultiFab& mf_Vwind)
{
    ParmParse pp("remora");
    std::string my_prob_name; pp.get("prob_name",my_prob_name);
    std::string my_prob_name_ci = amrex::toLower(my_prob_name);

    if ( (my_prob_name_ci == "coupletoerf")   ||
         (my_prob_name_ci == "dogbone")     ) {
        mf_Uwind.setVal(0.0_rt);
        mf_Vwind.setVal(0.0_rt);

    } else if (my_prob_name_ci == "boundarylayer")  {
#include "Prob/REMORA_InitAnalyticWind_BoundaryLayer.H"
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
