#ifndef _PROB_H_
#define _PROB_H_

#include "AMReX_REAL.H"
#include "REMORA_prob_common.H"

struct ProbParm {
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem(const amrex::Real* /*problo*/, const amrex::Real* /*probhi*/);

    void init_analytic_bathymetry (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_h) override;

    void init_analytic_grid_scale (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_pm, amrex::MultiFab& mf_pn) override;

    void init_analytic_zeta (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& mf_zeta) override;

    void init_analytic_prob(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_cons,
        amrex::MultiFab& mf_xvel,
        amrex::MultiFab& mf_yvel) override;

    void init_analytic_vmix(
        int lev,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        amrex::MultiFab& mf_Akv, amrex::MultiFab& mf_Akt) override;

    void init_analytic_hmix(
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& mf_visc2_p,
        amrex::MultiFab& mf_visc2_r,
        amrex::MultiFab& mf_diff2) override;

    void init_analytic_wind(
        int lev,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        amrex::MultiFab& mf_Uwind, amrex::MultiFab& mf_Vwind) override;

private:
    [[maybe_unused]] ProbParm parms;
};

#endif

