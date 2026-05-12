#include "prob.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "REMORA_IndexDefines.H"
#include "REMORA_DepthStretchTransform.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem(const amrex::Real* /*problo*/, const amrex::Real* /*probhi*/)
{
    // Parse params
    ParmParse pp("remora.prob");

    // Add problem-specific parameters below, e.g.

    /*
    pp.query("u_0", parms.u_0);
    pp.query("v_0", parms.v_0);
    */
}


/**
 * \brief Initializes custom sea surface height
 */
void Problem::init_analytic_zeta (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        MultiFab& mf_zeta)
{
    mf_zeta.setVal(0.0_rt);
}

void Problem::init_analytic_prob(
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& m_solverChoice,
        REMORA const& /*remora*/,
        amrex::MultiFab& mf_cons,
        amrex::MultiFab& mf_xvel,
        amrex::MultiFab& mf_yvel,
        amrex::MultiFab& mf_zvel)
{
    auto T0 = m_solverChoice.T0;
    auto S0 = m_solverChoice.S0;

    mf_cons.setVal(T0,Temp_comp,1);
    mf_cons.setVal(S0,Salt_comp,1);
    mf_cons.setVal(0.0_rt,Tracer_comp,1);
    mf_xvel.setVal(0.0_rt);
    mf_yvel.setVal(0.0_rt);
    mf_zvel.setVal(0.0_rt);
}

void Problem::init_analytic_wind(
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        MultiFab& mf_Uwind, MultiFab& mf_Vwind)
{
    mf_Uwind.setVal(0.0_rt);
    mf_Vwind.setVal(0.0_rt);
}
