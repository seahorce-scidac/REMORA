#include "prob.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "REMORA_IndexDefines.H"
#include "REMORA_DepthStretchTransform.H"

using namespace amrex;

// DELETE ANY UNUSED FUNCTIONS
/**
 * \brief Initializes bathymetry h and surface height Zeta
 */
void Problem::init_analytic_bathymetry (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_h)
{}

void Problem::init_analytic_grid_scale (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_pm, amrex::MultiFab& mf_pn)
{}

/**
 * \brief Initializes custom sea surface height
 */
void Problem::init_analytic_zeta (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        MultiFab& mf_zeta)
{}

void Problem::init_analytic_prob(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_cons,
        amrex::MultiFab& mf_xvel,
        amrex::MultiFab& mf_yvel)
{}

void Problem::init_analytic_vmix(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        MultiFab& mf_Akv, MultiFab& mf_Akt)
{}

void Problem::init_analytic_hmix(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        MultiFab& mf_visc2_p,
        MultiFab& mf_visc2_r,
        MultiFab& mf_diff2)
{}

void Problem::init_analytic_wind(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        MultiFab& mf_Uwind, MultiFab& mf_Vwind)
{}

void Problem::init_analytic_smflux(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        MultiFab& mf_sustr, MultiFab& mf_svstr)
{}

void Problem::init_analytic_coriolis (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_fcor)
{}
