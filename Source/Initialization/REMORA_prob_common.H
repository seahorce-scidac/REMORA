#ifndef _REMORA_PROBCOMMON_H_
#define _REMORA_PROBCOMMON_H_

#include <REMORA.H>
#include <REMORA_DataStruct.H>

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

// Forward declaration to resolve circular dependency
class REMORA;

/**
 * Class to hold problem-specific routines
 */

class ProblemBase
{
public:
    /**
     * Virtual destructor to avoid data leakage with derived class
     */
    virtual ~ProblemBase () = default;

    virtual void
    init_analytic_prob (
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& /*mf_cons*/,
        amrex::MultiFab& /*mf_xvel*/,
        amrex::MultiFab& /*mf_yvel*/,
        amrex::MultiFab& /*mf_zvel*/) {
          amrex::Error("Initialization set to analytic but init_analytic_prob not defined in prob.cpp");
    }

    virtual void
    init_analytic_bathymetry (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& /*mf_h*/) {
          amrex::Error("Doing analytic bathymetry, but not defined in prob.cpp");
    }

    virtual void
    init_analytic_grid_scale (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& /*mf_pm*/, amrex::MultiFab& /*mf_pn*/) {
          amrex::Error("Doing analytic grid scale (pm, pn), but not defined in prob.cpp");
    }

    virtual void
    init_analytic_zeta (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& /*mf_zeta*/) {
          amrex::Error("Doing analytic sea surface height (zeta), but not defined in prob.cpp");
    }

    virtual void
    init_analytic_coriolis (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& /*mf_fcor*/) {
          amrex::Error("Doing analytic Coriolis forcing, but not defined in prob.cpp");
    }

    virtual void
    init_analytic_vmix (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& /*mf_Akv*/,
        amrex::MultiFab& /*mf_Akt*/) {
          amrex::Error("Doing analytic vertical mixing, but not defined in prob.cpp");
    }

    virtual void
    init_analytic_hmix (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& /*mf_visc2_p*/,
        amrex::MultiFab& /*mf_visc2_r*/,
        amrex::MultiFab& /*mf_diff2*/) {
          amrex::Error("Doing analytic horizontal mixing, but not defined in prob.cpp");
    }

    virtual void
    init_analytic_smflux (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& /*mf_sustr*/,
        amrex::MultiFab& /*mf_svstr*/) {
          amrex::Error("Doing analytic surface momentum flux, but not defined in prob.cpp");
    }

    virtual void
    init_analytic_wind (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& /*mf_uwind*/,
        amrex::MultiFab& /*mf_vwind*/) {
          amrex::Error("Doing analytic wind, but not defined in prob.cpp");
    }
};

/**
 * Function to init the physical bounds of the domain
 * and instantiate a Problem derived from ProblemBase
*/
extern std::unique_ptr<ProblemBase> amrex_probinit (
  const amrex_real* problo,
  const amrex_real* probhi) AMREX_ATTRIBUTE_WEAK;

#endif
