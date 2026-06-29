#include "prob.H"
#include "REMORA_prob_common.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "REMORA_IndexDefines.H"
#include "REMORA_DepthStretchTransform.H"

using namespace amrex;

void Problem::init_analytic_vmix(
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& m_solverChoice,
        REMORA const& /*remora*/,
        MultiFab& mf_Akv, MultiFab& mf_Akt)
{
    mf_Akv.setVal(m_solverChoice.Akv_bak);
    mf_Akt.setVal(m_solverChoice.Akt_bak);
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
    mf_visc2_p.setVal(0.0_rt);
    mf_visc2_r.setVal(0.0_rt);
    mf_diff2.setVal(0.0_rt);
}

void Problem::init_analytic_smflux(
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        MultiFab& mf_sustr, MultiFab& mf_svstr)
{
    mf_svstr.setVal(0.0_rt);
    mf_sustr.setVal(0.0_rt);
}
