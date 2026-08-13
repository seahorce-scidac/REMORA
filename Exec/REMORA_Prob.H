#ifndef REMORA_PROB_H_
#define REMORA_PROB_H_

#include <string>

#include "AMReX_REAL.H"
#include "REMORA_prob_common.H"

class Problem : public ProblemBase
{
public:
    Problem (const amrex::Real* problo, const amrex::Real* probhi);

    void init_analytic_bathymetry (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_h) override;

    void init_analytic_prob(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_cons,
        amrex::MultiFab& mf_xvel,
        amrex::MultiFab& mf_yvel) override;

     void init_analytic_zeta (
        int lev, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        amrex::MultiFab& mf_zeta) override;

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

    void init_analytic_masks(
        int lev,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        amrex::MultiFab& mf_mskr) override;

    void init_analytic_smflux(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_sustr,
        amrex::MultiFab& mf_svstr) override;

    void init_analytic_grid_scale (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_pm, amrex::MultiFab& mf_pn) override;

    void init_analytic_coriolis (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_fcor) override;

    void init_analytic_surface_var (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_Uwind,
        amrex::MultiFab& mf_Vwind,
        amrex::MultiFab& mf_Tair,
        amrex::MultiFab& mf_qair,
        amrex::MultiFab& mf_Pair,
        amrex::MultiFab& mf_srflx,
        amrex::MultiFab& mf_longwave_down,
        amrex::MultiFab& mf_rain,
        amrex::MultiFab& mf_cloud,
        amrex::MultiFab& mf_EminusP) override;
};

#endif
