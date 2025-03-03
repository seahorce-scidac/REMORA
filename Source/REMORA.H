#ifndef REMORA_H_
#define REMORA_H_

#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>
#include <AMReX_InterpFaceRegister.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_YAFluxRegister.H>
#include <AMReX_ErrorList.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <REMORA_Math.H>
#include <REMORA_PhysBCFunct.H>
#include <REMORA_FillPatcher.H>

#include <REMORA_IndexDefines.H>
#include <REMORA_TimeInterpolatedData.H>
#include <REMORA_DataStruct.H>
#include <REMORA_Derive.H>

#ifdef REMORA_USE_PARTICLES
#include "REMORA_ParticleData.H"
#endif

#ifdef REMORA_USE_NETCDF
#include "REMORA_NCInterface.H"
#endif

#ifdef REMORA_USE_MOAB
#include "REMORA_MOAB.H"
#endif

#include <iostream>

#ifdef AMREX_LAZY
#include <AMReX_Lazy.H>
#endif

namespace InterpType {
    enum {
        PCInterp = 0,
        NodeBilinear,
        CellConservativeLinear,
        CellBilinear,
        CellQuadratic,
        CellConservativeProtected,
        CellConservativeQuartic
    };
}

class REMORA
    : public amrex::AmrCore
{
public:

    ////////////////
    // public member functions

    /**
     * constructor - reads in parameters from inputs file
     *             - sizes multilevel arrays and data structures
     */
    REMORA ();
    virtual ~REMORA();

    /** Advance solution to final time */
    void Evolve ();

    /** Tag cells for refinement */
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    /** Initialize multilevel data */
    void InitData ();

    /** Init (NOT restart or regrid) */
    void init_only (int lev, amrex::Real time);

    /** Restart */
    void restart ();

    /** Called after every level 0 timestep */
    void post_timestep (int nstep, amrex::Real time, amrex::Real dt_lev);

    // Diagnostics

    /** Integrate conserved quantities for diagnostics */
    void sum_integrated_quantities(amrex::Real time);

    /* Perform the volume-weighted sum */
    amrex::Real
    volWgtSumMF(int lev,
      const amrex::MultiFab& mf, int comp, bool local, bool finemask);

    /* Decide if it is time to take an action */
    bool is_it_time_for_action(int nstep, amrex::Real time, amrex::Real dt,
                               int action_interval, amrex::Real action_per);

    /** Make a new level using provided BoxArray and DistributionMapping and
     * fill with interpolated coarse level data.
     * Overrides the pure virtual function in AmrCore
     */
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                     const amrex::DistributionMapping& dm) override;

    /** Remake an existing level using provided BoxArray and DistributionMapping and
     * fill with existing fine and coarse data.
     * Overrides the pure virtual function in AmrCore
     */
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                  const amrex::DistributionMapping& dm) override;

    /** Delete level data
     * Overrides the pure virtual function in AmrCore
     */
    virtual void ClearLevel (int lev) override;

    /** Make a new level from scratch using provided BoxArray and DistributionMapping.
     * Only used during initialization.
     * Overrides the pure virtual function in AmrCore
     */
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm) override;

    /** Set pm and pn arrays on level lev. Only works if using Custom initialization */
    void set_pm_pn (int lev);

    /** Set zeta components to be equal to time-averaged Zt_avg1 */
    void set_zeta_to_Ztavg (int lev);

    /** Set psi-point mask to be consistent with rho-point mask */
    void update_mskp (int lev);

    /** compute dt from CFL considerations */
    amrex::Real estTimeStep (int lev) const;

    /** Interface for advancing the data at one level by one "slow" timestep */
    void remora_advance(int level,
                        amrex::MultiFab& cons_old,  amrex::MultiFab& cons_new,
                        amrex::MultiFab& xvel_old,  amrex::MultiFab& yvel_old,  amrex::MultiFab& zvel_old,
                        amrex::MultiFab& xvel_new,  amrex::MultiFab& yvel_new,  amrex::MultiFab& zvel_new,
                        amrex::MultiFab& source,
                        const amrex::Geometry fine_geom,
                        const amrex::Real dt, const amrex::Real time);

    /** Make mask to zero out covered cells*/
    amrex::MultiFab& build_fine_mask(int lev);

    /** write plotfile */
    void WritePlotFile ();

    void WriteMultiLevelPlotfileWithBathymetry (const std::string &plotfilename,
                                                int nlevels,
                                                const amrex::Vector<const amrex::MultiFab*> &mf,
                                                const amrex::Vector<const amrex::MultiFab*> &mf_nd,
                                                const amrex::Vector<std::string> &varnames,
                                                amrex::Real time,
                                                const amrex::Vector<int> &level_steps,
                                                const std::string &versionName = "HyperCLaw-V1.1",
                                                const std::string &levelPrefix = "Level_",
                                                const std::string &mfPrefix = "Cell",
                                                const amrex::Vector<std::string>& extra_dirs = amrex::Vector<std::string>()) const;


    void WriteGenericPlotfileHeaderWithBathymetry (std::ostream &HeaderFile,
                                                   int nlevels,
                                                   const amrex::Vector<amrex::BoxArray> &bArray,
                                                   const amrex::Vector<std::string> &varnames,
                                                   amrex::Real time,
                                                   const amrex::Vector<int> &level_steps,
                                                   const std::string &versionName,
                                                   const std::string &levelPrefix,
                                                   const std::string &mfPrefix) const;

    std::string pp_prefix {"remora"};

    amrex::Vector<amrex::MultiFab*> cons_old;
    amrex::Vector<amrex::MultiFab*> xvel_old;
    amrex::Vector<amrex::MultiFab*> yvel_old;
    amrex::Vector<amrex::MultiFab*> zvel_old;

    amrex::Vector<amrex::MultiFab*> cons_new;
    amrex::Vector<amrex::MultiFab*> xvel_new;
    amrex::Vector<amrex::MultiFab*> yvel_new;
    amrex::Vector<amrex::MultiFab*> zvel_new;

    // Program state data is represented by vectors of pointers to AMReX Multifabs.
    // There is one pointer per level

    /** Bathymetry data (2D, positive valued, h in ROMS) **/
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_hOfTheConfusingName;

    /** Width of cells in the vertical (z-) direction (3D, Hz in ROMS) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Hz;
    /** u-volume flux (3D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Huon;
    /** v-volume flux (3D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Hvom;
    /** u velocity RHS (3D, includes horizontal and vertical advection) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_ru;
    /** v velocity RHS (3D, includes horizontal and vertical advection) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rv;
    /** u velocity RHS (2D, includes horizontal and vertical advection) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_ru2d;
    /** v velocity RHS (2D, includes horizontal and vertical advection) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rv2d;
    /** u velocity RHS, integrated, including advection and bottom/surface stresses (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rufrc;
    /** v velocity RHS, integrated, including advection and bottom/surface stresses (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rvfrc;
    /** Vertical viscosity coefficient (3D), set in .in file */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Akv;
    /** Vertical diffusion coefficient (3D), set in .in file */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Akt;
    /** Harmonic viscosity defined on the psi points (corners of horizontal grid cells) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_visc2_p;
    /** Harmonic viscosity defined on the rho points (centers) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_visc2_r;
    /** Harmonic diffusivity for temperature / salinity */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_diff2;

    /** x coordinates at rho points (cell centers) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_x_r;

    /** y coordinates at rho points (cell centers) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_y_r;

    /** z coordinates at rho points (cell centers) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_z_r;

    /** z coordinates at w points (faces between z-cells) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_z_w;

    /** Scaled vertical coordinate (range [0,1]) that transforms to z, defined at rho points (cell centers) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_s_r;

    /** Scaled vertical coordinate (range [0,1]) that transforms to z, defined at w-points (cell faces) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_s_w;

    /** z coordinates at psi points (cell nodes) **/
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_z_phys_nd;

    /** Average of the free surface, zeta (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Zt_avg1;

    /** Surface stress in the u direction */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_sustr;

    /** Surface stress in the v direction */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_svstr;

    /** Linear drag coefficient [m/s], defined at rho points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rdrag;
    /** Quadratic drag coefficient [unitless], defined at rho points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rdrag2;
    /** Bottom roughness length [m], defined at rho points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_ZoBot;

    /** Bottom stress in the u direction */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_bustr;

    /** Bottom stress in the v direction */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_bvstr;

    /** time average of barotropic x velocity flux (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_DU_avg1;
    /** correct time average of barotropic x velocity flux for coupling (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_DU_avg2;
    /** time average of barotropic y velocity flux */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_DV_avg1;
    /** correct time average of barotropic y velocity flux for coupling (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_DV_avg2;
    /** barotropic x velocity for the RHS (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rubar;
    /** barotropic y velocity for the RHS (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rvbar;
    /** free surface height for the RHS (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rzeta;
    /** barotropic x velocity (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_ubar;
    /** barotropic y velocity (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_vbar;
    /** free surface height (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_zeta;

    /** land/sea mask at cell centers (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_mskr;

    /** land/sea mask at x-faces (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_msku;

    /** land/sea mask at y-faces (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_mskv;

    /** land/sea mask at cell corners (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_mskp;

    /** map factor (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_pm;

    /** map factor (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_pn;

    /** coriolis factor (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_fcor;

    /** x_grid on rho points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_xr;

    /** y_grid on rho points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_yr;

    /** x_grid on u-points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_xu;

    /** y_grid on u-points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_yu;

    /** x_grid on v-points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_xv;

    /** y_grid on v-points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_yv;

    /** x_grid on psi-points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_xp;

    /** y_grid on psi-points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_yp;

    /** saltstore, tempstore, etc*/
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_sstore;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rhoS;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rhoA;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_bvf;

    amrex::Vector<amrex::Real> vec_weight1;
    amrex::Vector<amrex::Real> vec_weight2;

    // GLS
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_tke;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_gls;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Lscale;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Akk;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Akp;

    /** advance a single level for a single time step */
    void Advance (int lev, amrex::Real time, amrex::Real dt_lev, int iteration, int ncycle);

    /** Set everything up for a step on a level */
    void setup_step(int lev, amrex::Real time, amrex::Real dt_lev);

    /** 3D advance on a single level */
    void advance_3d_ml (int lev, amrex::Real dt_lev);

    /** 2D advance, one predictor/corrector step */
    void advance_2d_onestep (int lev, amrex::Real dt_lev, amrex::Real dtfast_lev, int my_iif, int nfast_counter);

    /** Perform a 2D predictor (predictor_2d_step=True) or corrector (predictor_2d_step=False) step */
    void advance_2d (int lev,
                     amrex::MultiFab const* mf_rhoS,
                     amrex::MultiFab const* mf_rhoA,
                     amrex::MultiFab      * mf_ru,
                     amrex::MultiFab      * mf_rv,
                     amrex::MultiFab      * mf_rufrc,
                     amrex::MultiFab      * mf_rvfrc,
                     amrex::MultiFab      * mf_Zt_avg1,
                     std::unique_ptr<amrex::MultiFab>& mf_DU_avg1,
                     std::unique_ptr<amrex::MultiFab>& mf_DU_avg2,
                     std::unique_ptr<amrex::MultiFab>& mf_DV_avg1,
                     std::unique_ptr<amrex::MultiFab>& mf_DV_avg2,
                     std::unique_ptr<amrex::MultiFab>& mf_rubar,
                     std::unique_ptr<amrex::MultiFab>& mf_rvbar,
                     std::unique_ptr<amrex::MultiFab>& mf_rzeta,
                     std::unique_ptr<amrex::MultiFab>& mf_ubar,
                     std::unique_ptr<amrex::MultiFab>& mf_vbar,
                     amrex::MultiFab      * mf_zeta,
                     amrex::MultiFab const* mf_h,
                     amrex::MultiFab const* mf_pm,
                     amrex::MultiFab const* mf_pn,
                     amrex::MultiFab const* mf_fcor,
                     amrex::MultiFab const* mf_visc2_p,
                     amrex::MultiFab const* mf_visc2_r,
                     amrex::MultiFab const* mf_mskr,
                     amrex::MultiFab const* mf_msku,
                     amrex::MultiFab const* mf_mskv,
                     amrex::MultiFab const* mf_mskp,
                     const amrex::Real dtfast_lev,
                     const bool predictor_2d_step,
                     const bool first_2d_step, int my_iif,
                     int & next_indx1);

    /** Advance the 3D variables */
    void advance_3d (int lev,
                     amrex::MultiFab& mf_cons,
                     amrex::MultiFab& mf_u        , amrex::MultiFab& mf_v,
                     amrex::MultiFab* mf_sstore,
                     amrex::MultiFab* mf_ru       , amrex::MultiFab* mf_rv,
                     std::unique_ptr<amrex::MultiFab>& mf_DU_avg1,
                     std::unique_ptr<amrex::MultiFab>& mf_DU_avg2,
                     std::unique_ptr<amrex::MultiFab>& mf_DV_avg1,
                     std::unique_ptr<amrex::MultiFab>& mf_DV_avg2,
                     std::unique_ptr<amrex::MultiFab>& mf_ubar,
                     std::unique_ptr<amrex::MultiFab>& mf_vbar,
                     std::unique_ptr<amrex::MultiFab>& mf_Akv,
                     std::unique_ptr<amrex::MultiFab>& mf_Akt,
                     std::unique_ptr<amrex::MultiFab>& mf_Hz,
                     std::unique_ptr<amrex::MultiFab>& mf_Huon,
                     std::unique_ptr<amrex::MultiFab>& mf_Hvom,
                     std::unique_ptr<amrex::MultiFab>& mf_z_w,
                     amrex::MultiFab const* mf_h,
                     amrex::MultiFab const* mf_pm,
                     amrex::MultiFab const* mf_pn,
                     amrex::MultiFab const* mf_msku,
                     amrex::MultiFab const* mf_mskv,
                     const int N,
                     const amrex::Real dt_lev);

    /** Wrapper function for prestep */
    void prestep (int lev,
                  amrex::MultiFab& mf_uold, amrex::MultiFab& mf_vold,
                  amrex::MultiFab& mf_u, amrex::MultiFab& mf_v,
                        amrex::MultiFab* mf_ru,
                        amrex::MultiFab* mf_rv,
                  amrex::MultiFab& S_old,
                  amrex::MultiFab& S_new,
                  amrex::MultiFab& mf_W, amrex::MultiFab& mf_DC,
                  /* MF mf_FC? */
                  const amrex::MultiFab* mf_z_r,
                  const amrex::MultiFab* mf_z_w,
                  const amrex::MultiFab* mf_h,
                  const amrex::MultiFab* mf_pm,
                  const amrex::MultiFab* mf_pn,
                  const amrex::MultiFab* mf_sustr,
                  const amrex::MultiFab* mf_svstr,
                  const amrex::MultiFab* mf_bustr,
                  const amrex::MultiFab* mf_bvstr,
                  const amrex::MultiFab* mf_msku,
                  const amrex::MultiFab* mf_mskv,
                  const int iic, const int nfirst,
                  const int nnew, int nstp, int nrhs,
                  int N, const amrex::Real dt_lev);

    /** Prestep calculations for the tracers */
    void prestep_t_advection (const amrex::Box& tbx,
                        const amrex::Box& gbx,
                        const amrex::Array4<amrex::Real      >& tempold,
                        const amrex::Array4<amrex::Real      >& tempcache,
                        const amrex::Array4<amrex::Real      >& Hz,
                        const amrex::Array4<amrex::Real      >& Huon,
                        const amrex::Array4<amrex::Real      >& Hvom,
                        const amrex::Array4<amrex::Real      >& W,
                        const amrex::Array4<amrex::Real      >& DC,
                        const amrex::Array4<amrex::Real      >& FC,
                        const amrex::Array4<amrex::Real      >& sstore,
                        const amrex::Array4<amrex::Real const>& z_w,
                        const amrex::Array4<amrex::Real const>& h,
                        const amrex::Array4<amrex::Real const>& pm,
                        const amrex::Array4<amrex::Real const>& pn,
                        const amrex::Array4<amrex::Real const>& msku,
                        const amrex::Array4<amrex::Real const>& mskv,
                        int iic, int ntfirst, int nrhs, int N,
                        const amrex::Real dt_lev);

    void rhs_t_3d (const amrex::Box& bx,
                   const amrex::Box& gbx,
                   const amrex::Array4<amrex::Real      >& t,
                   const amrex::Array4<amrex::Real const>& tempstore,
                   const amrex::Array4<amrex::Real const>& Huon,
                   const amrex::Array4<amrex::Real const>& Hvom,
                   const amrex::Array4<amrex::Real const>& Hz,
                   const amrex::Array4<amrex::Real const>& pn,
                   const amrex::Array4<amrex::Real const>& pm,
                   const amrex::Array4<amrex::Real const>& W,
                   const amrex::Array4<amrex::Real      >& FC,
                   const amrex::Array4<amrex::Real const>& msku,
                   const amrex::Array4<amrex::Real const>& mskv,
                   int nrhs, int nnew, int N, const amrex::Real dt_lev);

    /** Calculation of the RHS */
    void rhs_uv_3d (const amrex::Box& xbx,
                    const amrex::Box& ybx,
                    const amrex::Array4<amrex::Real const>& uold,
                    const amrex::Array4<amrex::Real const>& vold,
                    const amrex::Array4<amrex::Real      >& ru,
                    const amrex::Array4<amrex::Real      >& rv,
                    const amrex::Array4<amrex::Real      >& rufrc,
                    const amrex::Array4<amrex::Real      >& rvfrc,
                    const amrex::Array4<amrex::Real const>& sustr,
                    const amrex::Array4<amrex::Real const>& svstr,
                    const amrex::Array4<amrex::Real const>& bustr,
                    const amrex::Array4<amrex::Real const>& bvstr,
                    const amrex::Array4<amrex::Real const>& Huon,
                    const amrex::Array4<amrex::Real const>& Hvom,
                    const amrex::Array4<amrex::Real const>& pm,
                    const amrex::Array4<amrex::Real const>& pn,
                    const amrex::Array4<amrex::Real const>& W,
                    const amrex::Array4<amrex::Real      >& FC,
                    int nrhs, int N);

    void rhs_uv_2d (const amrex::Box& xbx,
                    const amrex::Box& ybx,
                    const amrex::Array4<amrex::Real const>& uold,
                    const amrex::Array4<amrex::Real const>& vold,
                    const amrex::Array4<amrex::Real      >& ru,
                    const amrex::Array4<amrex::Real      >& rv,
                    const amrex::Array4<amrex::Real const>& Duon,
                    const amrex::Array4<amrex::Real const>& Dvom,
                    const int nrhs);

    void rho_eos (const amrex::Box& bx,
                  const amrex::Array4<amrex::Real const>& state,
                  const amrex::Array4<amrex::Real      >& rho,
                  const amrex::Array4<amrex::Real      >& rhoA,
                  const amrex::Array4<amrex::Real      >& rhoS,
                  const amrex::Array4<amrex::Real      >& bvf,
                  const amrex::Array4<amrex::Real const>& Hz,
                  const amrex::Array4<amrex::Real const>& z_w,
                  const amrex::Array4<amrex::Real const>& z_r,
                  const amrex::Array4<amrex::Real const>& h,
                  const amrex::Array4<amrex::Real const>& mskr,
                  const int N);

    void lin_eos (const amrex::Box& bx,
                  const amrex::Array4<amrex::Real const>& state,
                  const amrex::Array4<amrex::Real      >& rho,
                  const amrex::Array4<amrex::Real      >& rhoA,
                  const amrex::Array4<amrex::Real      >& rhoS,
                  const amrex::Array4<amrex::Real      >& bvf,
                  const amrex::Array4<amrex::Real const>& Hz,
                  const amrex::Array4<amrex::Real const>& z_w,
                  const amrex::Array4<amrex::Real const>& z_r,
                  const amrex::Array4<amrex::Real const>& h,
                  const amrex::Array4<amrex::Real const>& mskr,
                  const int N);

    void nonlin_eos (const amrex::Box& bx,
                  const amrex::Array4<amrex::Real const>& state,
                  const amrex::Array4<amrex::Real      >& rho,
                  const amrex::Array4<amrex::Real      >& rhoA,
                  const amrex::Array4<amrex::Real      >& rhoS,
                  const amrex::Array4<amrex::Real      >& bvf,
                  const amrex::Array4<amrex::Real const>& Hz,
                  const amrex::Array4<amrex::Real const>& z_w,
                  const amrex::Array4<amrex::Real const>& z_r,
                  const amrex::Array4<amrex::Real const>& h,
                  const amrex::Array4<amrex::Real const>& mskr,
                  const int N);


    void prsgrd (const amrex::Box& bx,
                 const amrex::Box& gbx,
                 const amrex::Box& utbx,
                 const amrex::Box& vtbx,
                 const amrex::Array4<amrex::Real      >& ru,
                 const amrex::Array4<amrex::Real      >& rv,
                 const amrex::Array4<amrex::Real const>& pn,
                 const amrex::Array4<amrex::Real const>& pm,
                 const amrex::Array4<amrex::Real const>& rho,
                 const amrex::Array4<amrex::Real      >& FC,
                 const amrex::Array4<amrex::Real const>& Hz,
                 const amrex::Array4<amrex::Real const>& z_r,
                 const amrex::Array4<amrex::Real const>& z_w,
                 const amrex::Array4<amrex::Real const>& msku,
                 const amrex::Array4<amrex::Real const>& mskv,
                 const int nrhs, const int N);

    /** Update velocities or tracers with diffusion/viscosity
     * as the last part of the prestep */
    void prestep_diffusion (const amrex::Box& bx,
                        const amrex::Box& gbx,
                        const int ioff, const int joff,
                        const amrex::Array4<amrex::Real      >& vel,
                        const amrex::Array4<amrex::Real const>& vel_old,
                        const amrex::Array4<amrex::Real      >& rvel,
                        const amrex::Array4<amrex::Real const>& Hz,
                        const amrex::Array4<amrex::Real const>& Akv,
                        const amrex::Array4<amrex::Real      >& DC,
                        const amrex::Array4<amrex::Real      >& FC,
                        const amrex::Array4<amrex::Real const>& sstr,
                        const amrex::Array4<amrex::Real const>& bstr,
                        const amrex::Array4<amrex::Real const>& z_r,
                        const amrex::Array4<amrex::Real const>& pm,
                        const amrex::Array4<amrex::Real const>& pn,
                        const int iic, const int ntfirst, const int nnew, int nstp, int nrhs, int N,
                        const amrex::Real lambda, const amrex::Real dt_lev);

    void vert_visc_3d (const amrex::Box& bx,
                       const int ioff, const int joff,
                       const amrex::Array4<amrex::Real      >& phi,
                       const amrex::Array4<amrex::Real const>& Hz,
                       const amrex::Array4<amrex::Real      >& Hzk,
                       const amrex::Array4<amrex::Real      >& AK,
                       const amrex::Array4<amrex::Real      >& Akv,
                       const amrex::Array4<amrex::Real      >& BC,
                       const amrex::Array4<amrex::Real      >& DC,
                       const amrex::Array4<amrex::Real      >& FC,
                       const amrex::Array4<amrex::Real      >& CF,
                       const int nnew, const int N,
                       const amrex::Real dt_lev);

    void update_massflux_3d (const amrex::Box& bx,
                             const int ioff, const int joff,
                             const amrex::Array4<amrex::Real      >& phi,
                             const amrex::Array4<amrex::Real      >& phibar,
                             const amrex::Array4<amrex::Real      >& Hphi,
                             const amrex::Array4<amrex::Real const>& Hz,
                             const amrex::Array4<amrex::Real const>& pm_or_pn,
                             const amrex::Array4<amrex::Real const>& Dphi1,
                             const amrex::Array4<amrex::Real const>& Dphi2,
                             const amrex::Array4<amrex::Real      >& DC,
                             const amrex::Array4<amrex::Real      >& FC,
                             const amrex::Array4<amrex::Real const>& msk,
                             const int nnew);

    void vert_mean_3d (const amrex::Box& bx,
                       const int ioff, const int joff,
                       const amrex::Array4<amrex::Real      >& phi,
                       const amrex::Array4<amrex::Real const>& Hz,
                       const amrex::Array4<amrex::Real const>& Dphi_avg1,
                       const amrex::Array4<amrex::Real      >& DC,
                       const amrex::Array4<amrex::Real      >& CF,
                       const amrex::Array4<amrex::Real const>& dxlen,
                       const amrex::Array4<amrex::Real const>& msk,
                       const int nnew, const int N);

    /** Harmonic viscosity */
    void uv3dmix  (const amrex::Box& xbx,
                   const amrex::Box& ybx,
                   const amrex::Array4<amrex::Real      >& u,
                   const amrex::Array4<amrex::Real      >& v,
                   const amrex::Array4<amrex::Real const>& uold,
                   const amrex::Array4<amrex::Real const>& vold,
                   const amrex::Array4<amrex::Real      >& rufrc,
                   const amrex::Array4<amrex::Real      >& rvfrc,
                   const amrex::Array4<amrex::Real const>& visc2_p,
                   const amrex::Array4<amrex::Real const>& visc2_r,
                   const amrex::Array4<amrex::Real const>& Hz,
                   const amrex::Array4<amrex::Real const>& pm,
                   const amrex::Array4<amrex::Real const>& pn,
                   const amrex::Array4<amrex::Real const>& mskp,
                   int nrhs, int nnew,
                   const amrex::Real dt_lev);

    /** Harmonic diffusivity for tracers */
    void t3dmix  (const amrex::Box& bx,
                  const amrex::Array4<amrex::Real      >& state,
                  const amrex::Array4<amrex::Real      >& state_rhs,
                  const amrex::Array4<amrex::Real const>& diff2,
                  const amrex::Array4<amrex::Real const>& Hz,
                  const amrex::Array4<amrex::Real const>& pm,
                  const amrex::Array4<amrex::Real const>& pn,
                  const amrex::Array4<amrex::Real const>& msku,
                  const amrex::Array4<amrex::Real const>& mskv,
                  const amrex::Real dt_lev,
                  const int ncomp);

    void coriolis (const amrex::Box& xbx,
                   const amrex::Box& ybx,
                   const amrex::Array4<amrex::Real const>& uold,
                   const amrex::Array4<amrex::Real const>& vold,
                   const amrex::Array4<amrex::Real      >& ru,
                   const amrex::Array4<amrex::Real      >& rv,
                   const amrex::Array4<amrex::Real const>& Hz,
                   const amrex::Array4<amrex::Real const>& fomn,
                   int nrhs, int nr);

    void set_2darrays (int lev);

    void set_zeta_average (int lev);

    void set_zeta (int lev);

    void set_bathymetry (int lev);

    void set_coriolis (int lev);

    void stretch_transform (int lev);

    void init_set_vmix (int lev);
    void set_analytical_vmix (int lev);
    void init_gls_vmix (int lev, SolverChoice solver_choice);
    void gls_prestep (int lev, amrex::MultiFab* mf_gls, amrex::MultiFab* mf_tke,
                      amrex::MultiFab& mf_W,
                      amrex::MultiFab* mf_msku, amrex::MultiFab* mf_mskv,
                      const int nstp, const int nnew, const int iic, const int ntfirst,
                      const int N, const amrex::Real dt_lev);
    void gls_corrector (int lev, amrex::MultiFab* mf_gls, amrex::MultiFab* mf_tke,
                        amrex::MultiFab& mf_W, amrex::MultiFab* mf_Akv, amrex::MultiFab* mf_Akt,
                        amrex::MultiFab* mf_Akk, amrex::MultiFab* mf_Akp,
                        amrex::MultiFab* mf_mskr,
                        amrex::MultiFab* mf_msku, amrex::MultiFab* mf_mskv,
                        const int nstp, const int nnew,
                        const int N, const amrex::Real dt_lev);

    void scale_rhs_vars ();
    void scale_rhs_vars_inv ();

    void set_smflux (int lev, amrex::Real time);

    void set_hmixcoef (int lev);

    void set_drag (int lev);

    void set_weights (int lev);

    // Fill the physical boundary conditions for cell-centered velocity (diagnostic only)
    void FillBdyCCVels (int lev, amrex::MultiFab& mf_cc_vel);

    // Fill a new MultiFab by copying in phi from valid region and filling ghost cells
    void FillPatch (int lev, amrex::Real time,
                    amrex::MultiFab& mf_to_be_filled,
                    amrex::Vector<amrex::MultiFab*> const& mfs,
                    const int bccomp,
                    const int bdy_var_type = BdyVars::null,
                    const int icomp=0,
                    const bool fill_all=true,
                    const bool fill_set=true,
                    const int  n_not_fill=0,
                    const int icomp_calc=0,
                    const amrex::Real dt = amrex::Real(0.0),
                    const amrex::MultiFab& mf_calc = amrex::MultiFab());

    void FillPatchNoBC (int lev, amrex::Real time,
                    amrex::MultiFab& mf_to_be_filled,
                    amrex::Vector<amrex::MultiFab*> const& mfs,
                    const int bdy_var_type = BdyVars::null,
                    const int icomp=0,
                    const bool fill_all=true,
                    const bool fill_set=true);

   void fill_from_bdyfiles (amrex::MultiFab& mf_to_fill,
                            const amrex::MultiFab& mf_mask,
                            const amrex::Real time,
                            const int bccomp,
                            const int bdy_var_type,
                            const int icomp_to_fill,
                            const int icomp_calc = 0,
                            const amrex::MultiFab& mf_calc = amrex::MultiFab(),
                            const amrex::Real = amrex::Real(0.0));

    void init_beta_plane_coriolis (int lev);

    void init_data_from_netcdf (int lev);
    void init_bdry_from_netcdf ();
    void init_masks_from_netcdf (int lev);
    void init_bathymetry_from_netcdf (int lev);
    void init_zeta_from_netcdf (int lev);
    void init_coriolis_from_netcdf (int lev);

    void mask_arrays_for_write (int lev, amrex::Real fill_value);


#ifdef REMORA_USE_NETCDF
    static int total_nc_plot_file_step;
#endif

private:

    ///////////////////////////
    // private member functions
    ///////////////////////////

    /** read in some parameters from inputs file */
    void ReadParameters();

    /** set covered coarse cells to be the average of overlying fine cells */
    void AverageDown ();

    void init_bcs();

    void init_custom(int lev);

    void init_stuff (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    void init_masks (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    void resize_stuff (int lev);

    /** more flexible version of AverageDown() that lets you average down across multiple levels */
    void AverageDownTo (int crse_lev);

    void Construct_REMORAFillPatchers (int lev);

    void Define_REMORAFillPatchers (int lev);

    /** fill an entire multifab by interpolating from the coarser level
     * this comes into play when a new level of refinement appears
     */
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab* mf_fine, amrex::MultiFab* mf_crse,
          const int  icomp = 0,
          const bool fill_all = true);

    /** utility to copy in data from old and/or new state into another multifab */
    TimeInterpolatedData GetDataAtTime (int lev, amrex::Real time);

    /** advance a level by dt,  includes a recursive call for finer levels */
    void timeStep (int lev, amrex::Real time, int iteration);

    /** advance all levels by dt, loops over finer levels */
    void timeStepML (amrex::Real time, int iteration);

    /** a wrapper for estTimeStep() */
    void ComputeDt ();

    /** get plotfile name */
    std::string PlotFileName (int lev) const;

    // set which variables and derived quantities go into plotfiles
    void setPlotVariables (const std::string& pp_plot_var_names);

    // append variables to plot
    void appendPlotVariables (const std::string& pp_plot_var_names);

#ifdef REMORA_USE_NETCDF
    //! Write plotfile using NETCDF
    void WriteNCPlotFile       (int istep);
    void WriteNCPlotFile_which (int lev, int which_subdomain,
                                bool write_header, ncutils::NCFile& ncf,
                                bool is_history);

    //! Write MultiFab in NetCDF format
    void WriteNCMultiFab (const amrex::FabArray<amrex::FArrayBox> &fab,
                          const std::string& name,
                          bool set_ghost = false) const;

    //! Read MultiFab in NetCDF format
    void ReadNCMultiFab (amrex::FabArray<amrex::FArrayBox> &fab,
                         const std::string &name,
                         int coordinatorProc = amrex::ParallelDescriptor::IOProcessorNumber(),
                         int allow_empty_mf = 0);

    // Vectors (over time) of Vector (over variables) of FArrayBoxs for holding the data read from the wrfbdy NetCDF file
    amrex::Vector<amrex::Vector<amrex::FArrayBox>> bdy_data_xlo;
    amrex::Vector<amrex::Vector<amrex::FArrayBox>> bdy_data_xhi;
    amrex::Vector<amrex::Vector<amrex::FArrayBox>> bdy_data_ylo;
    amrex::Vector<amrex::Vector<amrex::FArrayBox>> bdy_data_yhi;

    amrex::Real start_bdy_time;
    amrex::Real bdy_time_interval;

    static bool write_history_file;

    int bdy_width{0};
    int bdy_set_width{0};

#endif // REMORA_USE_NETCDF

#ifdef REMORA_USE_MOAB
    void InitMOABMesh();
#endif
    // Fillpatcher classes for coarse-fine boundaries
    int cf_width{0};
    int cf_set_width{0};
    amrex::Vector<REMORAFillPatcher> FPr_c;
    amrex::Vector<REMORAFillPatcher> FPr_u;
    amrex::Vector<REMORAFillPatcher> FPr_v;
    amrex::Vector<REMORAFillPatcher> FPr_w;

    // Fillpatcher classes for 2d velocity variables
    amrex::Vector<REMORAFillPatcher> FPr_ubar;
    amrex::Vector<REMORAFillPatcher> FPr_vbar;

    // write checkpoint file to disk
    void WriteCheckpointFile ();

    // read checkpoint file from disk
    void ReadCheckpointFile ();

    // Read the file passed to amr.restart and use it as an initial condition for
    // the current simulation. Supports a different number of components and
    // ghost cells.
    void InitializeFromFile ();

    // Initialize the new-time data at a level from the initial_data MultiFab
    void InitializeLevelFromData (int lev, const amrex::MultiFab& initial_data);

    // utility to skip to next line in Header
    static void GotoNextLine (std::istream& is);

    // Single level functions called by advance()
    void post_update (amrex::MultiFab& state_mf, const amrex::Real time, const amrex::Geometry& geom);
    void fill_rhs (amrex::MultiFab& rhs_mf, const amrex::MultiFab& state_mf, const amrex::Real time, const amrex::Geometry& geom);

    ////////////////
    // private data members

    amrex::Vector<int> num_boxes_at_level;                   // how many boxes specified at each level by tagging criteria
    amrex::Vector<int> num_files_at_level;                   // how many wrfinput files specified at each level
    amrex::Vector<amrex::Vector<amrex::Box>> boxes_at_level; //      the boxes specified at each level by tagging criteria

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // whether to set boundary conditions by variable rather than just by side
    bool set_bcs_by_var;

    amrex::Vector<std::unique_ptr<REMORAPhysBCFunct>> physbcs;
    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap".

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mapfac_m;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mapfac_u;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mapfac_v;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> sst;

    // array of flux registers
    amrex::Vector<amrex::YAFluxRegister*> advflux_reg;

    // A BCRec is essentially a 2*DIM integer array storing the boundary
    // condition type at each lo/hi walls in each direction. We have one BCRec
    // for each component of the cell-centered variables and each velocity component.
    amrex::Vector           <amrex::BCRec> domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d;

    // We store these so that we can print them out in the job_info file
    amrex::Array<std::string,2*AMREX_SPACEDIM> domain_bc_type;

    // These hold the Dirichlet values at walls which need them ...
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NCONS+8> m_bc_extdir_vals;

    // These are the "physical" boundary condition types (e.g. "inflow")
    amrex::GpuArray<amrex::GpuArray<REMORA_BC, AMREX_SPACEDIM*2>,BCVars::NumTypes> phys_bc_type;

    int last_plot_file_step;

    int last_check_file_step;
    int plot_file_on_restart = 1;

    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    // Time of the start of the simulation, in seconds
    amrex::Real start_time = 0.0;

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";

    // Time step controls
    static amrex::Real cfl;
    static amrex::Real init_shrink;
    static amrex::Real change_max;

    // Fixed dt for level 0 timesteps (only used if positive)
    static amrex::Real fixed_dt;
    static amrex::Real fixed_fast_dt;
    static int fixed_ndtfast_ratio;
    int nfast;

    // Whether to substep fine levels in time
    int do_substep = 0;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // plotfile prefix and frequency
    std::string plot_file_name {"plt_"};
    int plot_int = -1;

    // Checkpoint type, prefix and frequency
    std::string check_file {"chk"};
    int check_int = -1;

    // Whether to chunk netcdf history file
    bool chunk_history_file = false;

    // Number of time steps per netcdf history file.
    // -1 is the default and means code will automatically compute
    // number of steps per file based on grid size to keep files < 2 GB
    int steps_per_history_file = -1;

    amrex::Vector<std::string> plot_var_names;
    const amrex::Vector<std::string>     cons_names {"temp", "salt", "scalar"};

    // Note that the order of variable names here must match the order in PlotFile.cpp
    const amrex::Vector<std::string> derived_names {"vorticity"};


    // algorithm choices
    static SolverChoice solverChoice;

#ifdef REMORA_USE_PARTICLES
    // Particle container with all particle species
    ParticleData particleData;

    // variables and functions for tracers particles
    bool  m_use_tracer_particles; /*!< tracer particles that advect with flow */
    bool  m_use_hydro_particles;  /*!< tracer particles that fall with gravity */

    /*! Read tracer and hydro particles parameters */
    void readTracersParams();

    /*! Initialize tracer and hydro particles */
    void initializeTracers ( amrex::ParGDBBase*,
                             const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& );

    /*! Evolve tracers and hydro particles */
    void evolveTracers( int lev,
                        amrex::Real dt,
                        amrex::Vector<amrex::MultiFab const*>& flow_vel,
                        const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& );

#endif

    static int verbose;

    // Diagnostic output interval
    static int sum_interval;
    static amrex::Real sum_per;

    // Minimum number of digits in plotfile name or chunked history file
    static int file_min_digits;

    // Native or NetCDF
    static PlotfileType plotfile_type;

    // NetCDF initialization files for solution and grid
    static amrex::Vector<amrex::Vector<std::string>> nc_init_file;
    static amrex::Vector<amrex::Vector<std::string>> nc_grid_file;

    // NetCDF file of time series of boundary data
    static std::string nc_bdry_file;

    // Counter for which time index we are writing to in the netcdf history file
    int history_count = 0;

    void refinement_criteria_setup();

    //
    // Holds info for dynamically generated tagging criteria
    //
    static amrex::Vector<amrex::AMRErrorTag> ref_tags;

    //
    // Build a mask that zeroes out values on a coarse level underlying
    //     grids on the next finest level
    //
    amrex::MultiFab fine_mask;

    AMREX_FORCE_INLINE
    int
    ComputeGhostCells(const int& spatial_order) {
        int nGhostCells = 0;

        switch (spatial_order) {
            case 2:
                nGhostCells = 2; // We need this many to compute the eddy viscosity in the ghost cells
                break;
            case 3:
                nGhostCells = 2;
                break;
            case 4:
                nGhostCells = 2;
                break;
            case 5:
                nGhostCells = 3;
                break;
            case 6:
                nGhostCells = 3;
                break;
            default:
                amrex::Error("Must specify spatial order to be 2,3,4,5 or 6");
        }

        return nGhostCells;
    }

    AMREX_FORCE_INLINE
    amrex::YAFluxRegister* getAdvFluxReg (int lev)
    {
        return advflux_reg[lev];
    }

    AMREX_FORCE_INLINE
    std::ostream&
    DataLog (int i)
    {
        return *datalog[i];
    }

    AMREX_FORCE_INLINE
    int
    NumDataLogs () noexcept
    {
        return datalog.size();
    }

    static amrex::Real startCPUTime;
    static amrex::Real previousCPUTimeUsed;

    amrex::Real
    getCPUTime() const
    {
        int numCores = amrex::ParallelDescriptor::NProcs();
#ifdef _OPENMP
        numCores = numCores * omp_get_max_threads();
#endif

        amrex::Real T =
            numCores * (amrex::Real(amrex::ParallelDescriptor::second()) - startCPUTime) +
            previousCPUTimeUsed;

        return T;
    }

    void setRecordDataInfo (int i, const std::string& filename)
    {
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            datalog[i] = std::make_unique<std::fstream>();
            datalog[i]->open(filename.c_str(),std::ios::out|std::ios::app);
            if (!datalog[i]->good()) {
                amrex::FileOpenFailed(filename);
            }
        }
        amrex::ParallelDescriptor::Barrier("REMORA::setRecordDataInfo");
    }

    amrex::Vector<std::unique_ptr<std::fstream> > datalog;
    amrex::Vector<std::string> datalogname;

    //! The filename of the ith datalog file.
    const std::string DataLogName (int i) const noexcept { return datalogname[i]; }

public:
    void writeJobInfo (const std::string& dir) const;
    static void writeBuildInfo (std::ostream& os);

    static void print_banner(MPI_Comm /*comm*/, std::ostream& /*out*/);
    static void print_usage(MPI_Comm /*comm*/, std::ostream& /*out*/);
    static void print_error(MPI_Comm /*comm*/, const std::string& msg);
    static void print_summary(std::ostream&);
    static void print_tpls(std::ostream& /*out*/);
};

#endif
