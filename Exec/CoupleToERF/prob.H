#ifndef _PROB_H_
#define _PROB_H_

#include "AMReX_REAL.H"
#include "REMORA_prob_common.H"

struct ProbParm {
    // Add problem-specific parameters here
}; // namespace ProbParm


// DELETE ANY UNUSED FUNCTIONS
class Problem : public ProblemBase
{
public:
    Problem(const amrex::Real* /*problo*/, const amrex::Real* /*probhi*/);

    void init_analytic_zeta (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_zeta) override;

    void init_analytic_prob(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_cons,
        amrex::MultiFab& mf_xvel,
        amrex::MultiFab& mf_yvel) override;

    void init_analytic_wind(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_Uwind, amrex::MultiFab& mf_Vwind) override;
};

#endif

