#ifndef ROMSX_H_
#define ROMSX_H_

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
#include <AMReX_FluxRegister.H>
#include <AMReX_ErrorList.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <ROMSX_Math.H>
#include <ROMSX_PhysBCFunct.H>
#include <IndexDefines.H>
#include <TimeInterpolatedData.H>
#include <DataStruct.H>
#include <Derive.H>

#ifdef ROMSX_USE_NETCDF
#include "NCWpsFile.H"
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

template<typename V, typename T>
bool containerHasElement(const V& iterable, const T& query) {
    return std::find(iterable.begin(), iterable.end(), query) != iterable.end();
}

class ROMSX
    : public amrex::AmrCore
{
public:

    ////////////////
    // public member functions

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel arrays and data structures
    ROMSX ();
    virtual ~ROMSX();

    // Advance solution to final time
    void Evolve ();

    // Tag cells for refinement
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    // Initialize multilevel data
    void InitData ();

    // Init (NOT restart or regrid)
    void init_only (int lev, amrex::Real time);

    // Restart
    void restart ();

    // Called after every level 0 timestep
    void post_timestep (int nstep, amrex::Real time, amrex::Real dt_lev);

    // Diagnostics
    void sum_integrated_quantities(amrex::Real time);

    // Perform the volume-weighted sum
    amrex::Real
    volWgtSumMF(int lev,
      const amrex::MultiFab& mf, int comp, bool local, bool finemask);

    // Decide if it is time to take an action
    bool is_it_time_for_action(int nstep, amrex::Real time, amrex::Real dt,
                               int action_interval, amrex::Real action_per);

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                     const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                  const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm) override;

    // compute dt from CFL considerations
    amrex::Real estTimeStep (int lev) const;

    // Interface for advancing the data at one level by one "slow" timestep
    void romsx_advance(int level,
                       amrex::MultiFab& cons_old,  amrex::MultiFab& cons_new,
                       amrex::MultiFab& xvel_old,  amrex::MultiFab& yvel_old,  amrex::MultiFab& zvel_old,
                       amrex::MultiFab& xvel_new,  amrex::MultiFab& yvel_new,  amrex::MultiFab& zvel_new,
                       amrex::MultiFab& source,
                       const amrex::Geometry fine_geom,
                       const amrex::Real dt, const amrex::Real time
                      );

    amrex::MultiFab& build_fine_mask(int lev);

    void MakeHorizontalAverages();

    // write plotfile to disk
    void WritePlotFile  (int which, amrex::Vector<std::string> plot_var_names);

    void WriteMultiLevelPlotfileWithTerrain (const std::string &plotfilename,
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


    void WriteGenericPlotfileHeaderWithTerrain (std::ostream &HeaderFile,
                                                int nlevels,
                                                const amrex::Vector<amrex::BoxArray> &bArray,
                                                const amrex::Vector<std::string> &varnames,
                                                amrex::Real time,
                                                const amrex::Vector<int> &level_steps,
                                                const std::string &versionName,
                                                const std::string &levelPrefix,
                                                const std::string &mfPrefix) const;

#ifdef ROMSX_USE_NETCDF
    //! Write a timestep to 1D vertical column output for coupling
    void writeToNCColumnFile(const int lev,
                             const std::string& colfile_name, const amrex::Real xloc, const amrex::Real yloc,
                             const amrex::Real time);

    // Read the netcdf wrfinput file once
    void read_from_wrfinput(int lev, int idx,
                            amrex::Vector<amrex::FArrayBox>& NC_xvel_fab, amrex::Vector<amrex::FArrayBox>& NC_yvel_fab,
                            amrex::Vector<amrex::FArrayBox>& NC_zvel_fab, amrex::Vector<amrex::FArrayBox>& NC_rho_fab,
                            amrex::Vector<amrex::FArrayBox>& NC_rhop_fab, amrex::Vector<amrex::FArrayBox>& NC_rhoth_fab,
                            amrex::Vector<amrex::FArrayBox>& NC_MUB_fab,
                            amrex::Vector<amrex::FArrayBox>& NC_MSFU_fab, amrex::Vector<amrex::FArrayBox>& NC_MSFV_fab,
                            amrex::Vector<amrex::FArrayBox>& NC_MSFM_fab,
                            amrex::Vector<amrex::FArrayBox>& NC_SST_fab,
                            amrex::Vector<amrex::FArrayBox>& NC_C1H_fab,  amrex::Vector<amrex::FArrayBox>& NC_C2H_fab,
                            amrex::Vector<amrex::FArrayBox>& NC_RDNW_fab,
                            amrex::Vector<amrex::FArrayBox>& NC_PH_fab, amrex::Vector<amrex::FArrayBox>& NC_PHB_fab,
                            amrex::Vector<amrex::FArrayBox>& NC_ALB_fab, amrex::Vector<amrex::FArrayBox>& NC_PB_fab);

    // Copy from the NC*fabs into the MultiFabs holding the initial data
    void init_state_from_wrfinput(int lev,
                            amrex::FArrayBox& state_vel_fab, amrex::FArrayBox& x_vel_fab,
                            amrex::FArrayBox&     y_vel_fab, amrex::FArrayBox& z_vel_fab,
                            const amrex::Vector<amrex::FArrayBox>& NC_xvel_fab, const amrex::Vector<amrex::FArrayBox>& NC_yvel_fab,
                            const amrex::Vector<amrex::FArrayBox>& NC_zvel_fab, const amrex::Vector<amrex::FArrayBox>& NC_rho_fab,
                            const amrex::Vector<amrex::FArrayBox>& NC_rhoth_fab);

    void init_msfs_from_wrfinput(int lev,
                            amrex::FArrayBox& msfu_fab, amrex::FArrayBox& msfv_fab,
                            amrex::FArrayBox& msfm_fab,
                            const amrex::Vector<amrex::FArrayBox>& NC_MSFU_fab,
                            const amrex::Vector<amrex::FArrayBox>& NC_MSFV_fab,
                            const amrex::Vector<amrex::FArrayBox>& NC_MSFM_fab);

    // Copy from NC_PH* fabs into the MultiFabs holding the terrain data
    void init_terrain_from_wrfinput(int lev,  amrex::FArrayBox& z_phys,
                                    const amrex::Vector<amrex::FArrayBox>& NC_PH_fab,
                                    const amrex::Vector<amrex::FArrayBox>& NC_PHB_fab);

    // Copy from NC_PB fabs into the MultiFab holding the base state
    void init_base_state_from_wrfinput(int lev, const amrex::Box& bx,
                                       amrex::FArrayBox& p_hse_fab,
                                       amrex::FArrayBox& pi_hse_fab,
                                       amrex::FArrayBox& r_hse_fab,
                                       const amrex::FArrayBox& z_phys_nd_fab,
                                       const amrex::FArrayBox& z_phys_cc_fab,
                                       const amrex::Vector<amrex::FArrayBox>& NC_ALB_fab,
                                       const amrex::Vector<amrex::FArrayBox>& NC_PB_fab);
#endif //ROMSX_USE_NETCDF

    std::string pp_prefix {"romsx"};

    amrex::Vector<amrex::Vector<amrex::MultiFab> > vars_new;
    amrex::Vector<amrex::Vector<amrex::MultiFab> > vars_old;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_hOfTheConfusingName;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Hz;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Huon;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Hvom;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_ru;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rv;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rufrc;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rvfrc;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Akv;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_visc3d_r;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_visc2_p;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_visc2_r;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_visc2_salinity;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_visc2_temperature;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_x_r;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_y_r;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_z_r;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_z_w;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_s_r;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Zt_avg1;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_sustr;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_svstr;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rdrag;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_bustr;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_bvstr;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_DU_avg1;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_DU_avg2;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_DV_avg1;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_DV_avg2;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rubar;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rvbar;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rzeta;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_ubar;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_vbar;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_zeta;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_t3;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_s3;

    amrex::Vector<amrex::Real> vec_weight1;
    amrex::Vector<amrex::Real> vec_weight2;

    // advance a single level for a single time step
    void Advance (int lev, amrex::Real time, amrex::Real dt_lev, int iteration, int ncycle);

    void advance_2d (int lev,
                     amrex::MultiFab& mf_u , amrex::MultiFab& mf_v ,
                     std::unique_ptr<amrex::MultiFab>& mf_ru,
                     std::unique_ptr<amrex::MultiFab>& mf_rv,
                     std::unique_ptr<amrex::MultiFab>& mf_rufrc,
                     std::unique_ptr<amrex::MultiFab>& mf_rvfrc,
                     std::unique_ptr<amrex::MultiFab>& mf_Zt_avg1,
                     std::unique_ptr<amrex::MultiFab>& mf_DU_avg1,
                     std::unique_ptr<amrex::MultiFab>& mf_DU_avg2,
                     std::unique_ptr<amrex::MultiFab>& mf_DV_avg1,
                     std::unique_ptr<amrex::MultiFab>& mf_DV_avg2,
                     std::unique_ptr<amrex::MultiFab>& mf_rubar,
                     std::unique_ptr<amrex::MultiFab>& mf_rvbar,
                     std::unique_ptr<amrex::MultiFab>& mf_rzeta,
                     std::unique_ptr<amrex::MultiFab>& mf_ubar,
                     std::unique_ptr<amrex::MultiFab>& mf_vbar,
                     std::unique_ptr<amrex::MultiFab>& mf_zeta,
                     std::unique_ptr<amrex::MultiFab>& mf_h,
                     std::unique_ptr<amrex::MultiFab>& mf_visc2_p,
                     std::unique_ptr<amrex::MultiFab>& mf_visc2_r,
                     const int ncomp, const amrex::Real dt_lev, const amrex::Real dtfast_lev,
                     const bool predictor_2d_step,
                     const bool first_2d_step, int my_iif, int nfast,
                     int & next_indx1);

    void advance_3d (int lev,
                     amrex::MultiFab& mf_u  , amrex::MultiFab& mf_v ,
                     amrex::MultiFab& mf_tempold  , amrex::MultiFab& mf_saltold ,
                     amrex::MultiFab& mf_temp  , amrex::MultiFab& mf_salt ,
                     std::unique_ptr<amrex::MultiFab>& mf_tempstore,
                     std::unique_ptr<amrex::MultiFab>& mf_saltstore,
                     std::unique_ptr<amrex::MultiFab>& mf_ru,
                     std::unique_ptr<amrex::MultiFab>& mf_rv,
                     std::unique_ptr<amrex::MultiFab>& mf_DU_avg1,
                     std::unique_ptr<amrex::MultiFab>& mf_DU_avg2,
                     std::unique_ptr<amrex::MultiFab>& mf_DV_avg1,
                     std::unique_ptr<amrex::MultiFab>& mf_DV_avg2,
                     std::unique_ptr<amrex::MultiFab>& mf_ubar,
                     std::unique_ptr<amrex::MultiFab>& mf_vbar,
                     amrex::MultiFab& mf_AK , amrex::MultiFab& mf_DC,
                     amrex::MultiFab& mf_Hzk,
                     std::unique_ptr<amrex::MultiFab>& mf_Akv,
                     std::unique_ptr<amrex::MultiFab>& mf_Hz,
                     std::unique_ptr<amrex::MultiFab>& mf_Huon,
                     std::unique_ptr<amrex::MultiFab>& mf_Hvom,
                     const int ncomp, const int N,
                     const amrex::Real dt_lev);

    void prestep_t_3d (const amrex::Box& bx,
                        amrex::Array4<amrex::Real> uold,
                        amrex::Array4<amrex::Real> vold,
                        amrex::Array4<amrex::Real> u,
                        amrex::Array4<amrex::Real> v,
                        amrex::Array4<amrex::Real> tempold,
                        amrex::Array4<amrex::Real> saltold,
                        amrex::Array4<amrex::Real> temp,
                        amrex::Array4<amrex::Real> salt,
                        amrex::Array4<amrex::Real> ru,
                        amrex::Array4<amrex::Real> rv,
                        amrex::Array4<amrex::Real> Hz,
                        amrex::Array4<amrex::Real> Akv,
                        amrex::Array4<amrex::Real> on_u,
                        amrex::Array4<amrex::Real> om_v,
                        amrex::Array4<amrex::Real> Huon,
                        amrex::Array4<amrex::Real> Hvom,
                        amrex::Array4<amrex::Real> pm,
                        amrex::Array4<amrex::Real> pn,
                        amrex::Array4<amrex::Real> W,
                        amrex::Array4<amrex::Real> DC,
                        amrex::Array4<amrex::Real> FC,
                        amrex::Array4<amrex::Real> tempstore,
                        amrex::Array4<amrex::Real> saltstore,
                        amrex::Array4<amrex::Real> FE,
                        amrex::Array4<amrex::Real> FX,
                        amrex::Array4<amrex::Real> z_r,
                        int iic, int ntfirst, int nnew, int nstp, int nrhs, int N,
                        const amrex::Real lambda, const amrex::Real dt_lev);

    void prestep_uv_3d (const amrex::Box& bx,
                        amrex::Array4<amrex::Real> uold,
                        amrex::Array4<amrex::Real> vold,
                        amrex::Array4<amrex::Real> u,
                        amrex::Array4<amrex::Real> v,
                        amrex::Array4<amrex::Real> ru,
                        amrex::Array4<amrex::Real> rv,
                        amrex::Array4<amrex::Real> Hz,
                        amrex::Array4<amrex::Real> Akv,
                        amrex::Array4<amrex::Real> on_u,
                        amrex::Array4<amrex::Real> om_v,
                        amrex::Array4<amrex::Real> Huon,
                        amrex::Array4<amrex::Real> Hvom,
                        amrex::Array4<amrex::Real> pm,
                        amrex::Array4<amrex::Real> pn,
                        amrex::Array4<amrex::Real> W,
                        amrex::Array4<amrex::Real> DC,
                        amrex::Array4<amrex::Real> FC,
                        amrex::Array4<amrex::Real> z_r,
                        amrex::Array4<amrex::Real> sustr,
                        amrex::Array4<amrex::Real> svstr,
                        amrex::Array4<amrex::Real> bustr,
                        amrex::Array4<amrex::Real> bvstr,
                        int iic, int ntfirst, int nnew, int nstp, int nrhs, int N,
                        const amrex::Real lambda, const amrex::Real dt_lev);

    void rhs_t_3d (const amrex::Box& bx,
                   amrex::Array4<amrex::Real> told,
                   amrex::Array4<amrex::Real> t,
                   amrex::Array4<amrex::Real> tempstore,
                   amrex::Array4<amrex::Real> Huon,
                   amrex::Array4<amrex::Real> Hvom,
                   amrex::Array4<amrex::Real> oHz,
                   amrex::Array4<amrex::Real> pn,
                   amrex::Array4<amrex::Real> pm,
                   amrex::Array4<amrex::Real> W,
                   amrex::Array4<amrex::Real> FC,
                   int nrhs, int nnew, int N, const amrex::Real dt_lev);

    void rhs_3d (const amrex::Box& bx,
                 amrex::Array4<amrex::Real> uold,
                 amrex::Array4<amrex::Real> vold,
                 amrex::Array4<amrex::Real> ru,
                 amrex::Array4<amrex::Real> rv,
                 amrex::Array4<amrex::Real> rufrc,
                 amrex::Array4<amrex::Real> rvfrc,
                 amrex::Array4<amrex::Real> sustr,
                 amrex::Array4<amrex::Real> svstr,
                 amrex::Array4<amrex::Real> bustr,
                 amrex::Array4<amrex::Real> bvstr,
                 amrex::Array4<amrex::Real> Huon,
                 amrex::Array4<amrex::Real> Hvom,
                 amrex::Array4<amrex::Real> on_u,
                 amrex::Array4<amrex::Real> om_v,
                 amrex::Array4<amrex::Real> om_u,
                 amrex::Array4<amrex::Real> on_v,
                 amrex::Array4<amrex::Real> W,
                 amrex::Array4<amrex::Real> FC,
                 int nrhs, int N);

    void rhs_2d (const amrex::Box& bx,
                 amrex::Array4<amrex::Real> uold,
                 amrex::Array4<amrex::Real> vold,
                 amrex::Array4<amrex::Real> ru,
                 amrex::Array4<amrex::Real> rv,
                 //                 amrex::Array4<amrex::Real> rufrc,
                 //                 amrex::Array4<amrex::Real> rvfrc,
                 //                 amrex::Array4<amrex::Real> sustr,
                 //                 amrex::Array4<amrex::Real> svstr,
                 amrex::Array4<amrex::Real> Huon,
                 amrex::Array4<amrex::Real> Hvom,
                 //                 amrex::Array4<amrex::Real> on_u,
                 //                 amrex::Array4<amrex::Real> om_v,
                 //                 amrex::Array4<amrex::Real> om_u,
                 //                 amrex::Array4<amrex::Real> on_v,
                 //                 amrex::Array4<amrex::Real> W,
                 //                 amrex::Array4<amrex::Real> FC,
                 int nrhs, int N);

    void update_vel_3d (const amrex::Box& bx,
                        const int ioff, const int joff,
                        amrex::Array4<amrex::Real> vel,
                        amrex::Array4<amrex::Real> vel_old,
                        amrex::Array4<amrex::Real> rvel,
                        amrex::Array4<amrex::Real> Hz,
                        amrex::Array4<amrex::Real> Akv,
                        amrex::Array4<amrex::Real> pm,
                        amrex::Array4<amrex::Real> pn,
                        amrex::Array4<amrex::Real> DC,
                        amrex::Array4<amrex::Real> FC,
                        amrex::Array4<amrex::Real> sstr,
                        amrex::Array4<amrex::Real> bstr,
                        amrex::Array4<amrex::Real> z_r,
                        const int iic, const int ntfirst, const int nnew, int nstp, int nrhs, int N,
                        const amrex::Real lambda, const amrex::Real dt_lev);

    void vert_visc_3d (const amrex::Box& bx,
                       const int ioff, const int joff,
                       amrex::Array4<amrex::Real> phi,
                       amrex::Array4<amrex::Real> Hz,
                       amrex::Array4<amrex::Real> Hzk,
                       amrex::Array4<amrex::Real> oHz,
                       amrex::Array4<amrex::Real> AK,
                       amrex::Array4<amrex::Real> Akv,
                       amrex::Array4<amrex::Real> BC,
                       amrex::Array4<amrex::Real> DC,
                       amrex::Array4<amrex::Real> FC,
                       amrex::Array4<amrex::Real> CF,
                       const int nnew, const int N,
                       const amrex::Real dt_lev);

    void set_massflux_3d (const amrex::Box& bx,
                       const int ioff, const int joff,
                       amrex::Array4<amrex::Real> phi,
                       amrex::Array4<amrex::Real> Hphi,
                       amrex::Array4<amrex::Real> Hz,
                       amrex::Array4<amrex::Real> on_u_or_om_v,
                       const int nnew);

    void update_massflux_3d (const amrex::Box& bx,
                       const int ioff, const int joff,
                       amrex::Array4<amrex::Real> phi,
                       amrex::Array4<amrex::Real> Hphi,
                       amrex::Array4<amrex::Real> Hz,
                       amrex::Array4<amrex::Real> on_u_or_om_v,
                       amrex::Array4<amrex::Real> Dphi1,
                       amrex::Array4<amrex::Real> Dphi2,
                       amrex::Array4<amrex::Real> DC,
                       amrex::Array4<amrex::Real> FC,
                       amrex::Array4<amrex::Real> CF,
                       const int nnew);

    void vert_mean_3d (const amrex::Box& bx,
                       const int ioff, const int joff,
                       amrex::Array4<amrex::Real> phi,
                       amrex::Array4<amrex::Real> Hz,
                       amrex::Array4<amrex::Real> Hzk,
                       amrex::Array4<amrex::Real> oHz,
                       amrex::Array4<amrex::Real> Dphi_avg1,
                       amrex::Array4<amrex::Real> Akv,
                       amrex::Array4<amrex::Real> BC,
                       amrex::Array4<amrex::Real> DC,
                       amrex::Array4<amrex::Real> FC,
                       amrex::Array4<amrex::Real> CF,
                       amrex::Array4<amrex::Real> dxlen,
                       const int nnew, const int N,
                       const amrex::Real dt_lev);

    void uv3dmix  (const amrex::Box& bx,
                 amrex::Array4<amrex::Real> u,
                 amrex::Array4<amrex::Real> v,
                 amrex::Array4<amrex::Real> rufrc,
                 amrex::Array4<amrex::Real> rvfrc,
                 amrex::Array4<amrex::Real> visc2_p,
                 amrex::Array4<amrex::Real> visc2_r,
                 amrex::Array4<amrex::Real> Hz,
                 amrex::Array4<amrex::Real> om_r,
                 amrex::Array4<amrex::Real> on_r,
                 amrex::Array4<amrex::Real> om_p,
                 amrex::Array4<amrex::Real> on_p,
                 amrex::Array4<amrex::Real> pm,
                 amrex::Array4<amrex::Real> pn,
                 int nrhs, int nnew,
                 const amrex::Real dt_lev);

    void t3dmix  (const amrex::Box& bx,
                 amrex::Array4<amrex::Real> t,
                 amrex::Array4<amrex::Real> visc2,
                 amrex::Array4<amrex::Real> Hz,
                 amrex::Array4<amrex::Real> pm,
                 amrex::Array4<amrex::Real> pn,
                 int nrhs, int nnew,
                 const amrex::Real dt_lev);

#ifdef UV_COR
    void coriolis (const amrex::Box& bx,
                   amrex::Array4<amrex::Real> uold,
                   amrex::Array4<amrex::Real> vold,
                   amrex::Array4<amrex::Real> ru,
                   amrex::Array4<amrex::Real> rv,
                   amrex::Array4<amrex::Real> Hz,
                   amrex::Array4<amrex::Real> fomn,
                   int nrhs, int nr);
#endif

    void set_2darrays (int lev);

    void set_depth (int lev);

    void set_vmix (int lev);

    void set_smflux (int lev, amrex::Real time);

    void set_hmixcoef (int lev);

    void set_drag (int lev);

    void set_weights (int lev);

private:

    ///////////////////////////
    // private member functions
    ///////////////////////////

    // read in some parameters from inputs file
    void ReadParameters();

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown ();

    void init1DArrays();

    void init_bcs();

    void init_custom(int lev);

#ifdef ROMSX_USE_NETCDF
    void init_from_wrfinput(int lev);
#endif // ROMSX_USE_NETCDF

    // more flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo (int crse_lev);

    // compute a new multifab by copying in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    void FillPatch (int lev, amrex::Real time, const amrex::Vector<amrex::MultiFab*>& mf);

    // compute a new multifab by copying in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    void FillPatch (int lev, amrex::Real time,
                    amrex::Vector<amrex::MultiFab>& mf);

    // compute new multifabs by copying in data from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    // unlike FillPatch, FillIntermediatePatch will use the supplied multifabs instead of fine level data.
    // This is to support filling boundary cells at an intermediate time between old/new times
    // on the fine level when valid data at a specific time is already available (such as
    // at each RK stage when integrating between initial and final times at a given level).
    // NOTE: mfs should always contain {cons, xvel, yvel, zvel} multifab data.
    // if which_var is supplied, then only fill the specified variable in the vector of mfs
    void FillIntermediatePatch (int lev, amrex::Real time,
                                amrex::Vector<std::reference_wrapper<amrex::MultiFab> > mfs,
                                int ng_cons, int ng_vel, bool cons_only, int scomp_cons, int ncomp_cons);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch (int lev, amrex::Real time,
                          amrex::MultiFab& mf, int icomp, int ncomp, int var_idx);

    // fill all multifabs (and all components) in a vector of multifabs corresponding to the
    // grid variables defined in vars_old and vars_new just as FillCoarsePatch.
    void FillCoarsePatchAllVars (int lev, amrex::Real time,
                                 amrex::Vector<amrex::MultiFab>& vmf);

    // utility to copy in data from old and/or new state into another multifab
    TimeInterpolatedData GetDataAtTime (int lev, amrex::Real time);

    // advance a level by dt
    // includes a recursive call for finer levels
    void timeStep (int lev, amrex::Real time, int iteration);

    //! Initialize HSE
    void initHSE();

    //! Initialize Rayleigh damping profiles
    void initRayleigh();

    // a wrapper for estTimeStep()
    void ComputeDt ();

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames (amrex::Vector<std::string> plot_var_names) const;

    // set which variables and derived quantities go into plotfiles
    void setPlotVariables (const std::string& pp_plot_var_names, amrex::Vector<std::string>& plot_var_names);

#ifdef ROMSX_USE_NETCDF
    //! Write plotfile using NETCDF
    void writeNCPlotFile(int lev, int which, const std::string& dir,
                         const amrex::Vector<const amrex::MultiFab*> &mf,
                         const amrex::Vector<std::string> &varnames,
                         const amrex::Vector<int> level_steps, const amrex::Real time) const;

    //! Write checkpointFile using NetCdf
    void WriteNCCheckpointFile () const;

    //! Read checkpointFile for restart
    void ReadNCCheckpointFile ();

    //! Write MultiFab in NetCDF format
    void WriteNCMultiFab (const amrex::FabArray<amrex::FArrayBox> &fab,
                          const std::string& name,
                          bool set_ghost = false) const;

    //! Read MultiFab in NetCDF format
    void ReadNCMultiFab (amrex::FabArray<amrex::FArrayBox> &fab,
                         const std::string &name,
                         int coordinatorProc = amrex::ParallelDescriptor::IOProcessorNumber(),
                         int allow_empty_mf = 0);

    //! Create 1D vertical column output for coupling
    void createNCColumnFile(const int lev,
                            const std::string& colfile_name, const amrex::Real xloc, const amrex::Real yloc);

    // Copy from the NC*fabs into the MultiFabs holding the boundary data
    void init_from_wrfbdy(amrex::Vector<amrex::FArrayBox*> x_vel_lateral,
                          amrex::Vector<amrex::FArrayBox*> y_vel_lateral,
                          amrex::Vector<amrex::FArrayBox*> z_vel_lateral,
                          amrex::Vector<amrex::FArrayBox*> T_lateral);

    // *** *** FArrayBox's for holding the SURFACE data
    // Data read from the wrfinput NetCDF file
    amrex::FArrayBox NC_HGT_fab;    // Terrain height (variable HGT in the wrfinput file)
    amrex::IArrayBox NC_IVGTYP_fab; // Vegetation type (IVGTYP); Discrete numbers;
    amrex::FArrayBox NC_z0_fab;     // Surface Roughness, z0 = z0 (IVGTYP)
    amrex::FArrayBox NC_PSFC_fab;   // Surface pressure

    // TODO: Clarify the relation between SST and TSK
    amrex::FArrayBox NC_SST_fab;    // Sea Surface Temperature; Defined even for land area
    amrex::FArrayBox NC_TSK_fab;    // Surface Skin Temperature; Appears to be same as SST...

    amrex::FArrayBox NC_eta_fab, NC_phb_fab, NC_z_phy_fab;
#endif // ROMSX_USE_NETCDF

#ifdef ROMSX_USE_NETCDF
    // Vectors (over time) of Vector (over variables) of FArrayBoxs for holding the data read from the wrfbdy NetCDF file
    amrex::Vector<amrex::Vector<amrex::FArrayBox>> bdy_data_xlo;
    amrex::Vector<amrex::Vector<amrex::FArrayBox>> bdy_data_xhi;
    amrex::Vector<amrex::Vector<amrex::FArrayBox>> bdy_data_ylo;
    amrex::Vector<amrex::Vector<amrex::FArrayBox>> bdy_data_yhi;

    amrex::Real bdy_time_interval;
#endif // ROMSX_USE_NETCDF

    // write checkpoint file to disk
    void WriteCheckpointFile () const;

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

    amrex::Vector<std::unique_ptr<ROMSXPhysBCFunct>> physbcs;
    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap".

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_phys_nd;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> z_phys_cc;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>   detJ_cc;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mapfac_m;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mapfac_u;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> mapfac_v;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> sst;

    amrex::Vector<amrex::MultiFab> base_state;

    // array of flux registers
    amrex::Vector<amrex::FluxRegister*> flux_registers;

    // A BCRec is essentially a 2*DIM integer array storing the boundary
    // condition type at each lo/hi walls in each direction. We have one BCRec
    // for each component of the cell-centered variables and each velocity component.
    amrex::Vector           <amrex::BCRec> domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d;

    // We store these so that we can print them out in the job_info file
    amrex::Array<std::string,2*AMREX_SPACEDIM> domain_bc_type;

    // These hold the Dirichlet values at walls which need them ...
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_extdir_vals;

    // These are the "physical" boundary condition types (e.g. "inflow")
    amrex::GpuArray<ROMSX_BC, AMREX_SPACEDIM*2> phys_bc_type;

    int last_plot_file_step_1;
    int last_plot_file_step_2;

    int total_plot_file_step_1 = -1;

    int last_check_file_step;
    int plot_file_on_restart = 1;

    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

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

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // plotfile prefix and frequency
    std::string plot_file_1 {"plt_1_"};
    std::string plot_file_2 {"plt_2_"};
    int plot_int_1 = -1;
    int plot_int_2 = -1;

    // Checkpoint type, prefix and frequency
    std::string check_file {"chk"};
    std::string check_type {"native"};
    std::string restart_type {"native"};
    int check_int = -1;

    amrex::Vector<std::string> plot_var_names_1;
    amrex::Vector<std::string> plot_var_names_2;
    const amrex::Vector<std::string> velocity_names {"x_velocity", "y_velocity", "z_velocity"};
    const amrex::Vector<std::string> cons_names { "density", "temp", "rhoKE", "rhoQKE", "rhoadv_0", "omega", "salt"};

    // Note that the order of variable names here must match the order in Derive.cpp
    const amrex::Vector<std::string> derived_names {"KE", "QKE", "scalar", "pres_hse",
                                                    "dens_hse", "pert_pres", "pert_dens", "dpdx", "dpdy",
                                                    "pres_hse_x", "pres_hse_y", "z_phys", "detJ"
                                                   };

    // algorithm choices
    static SolverChoice solverChoice;

    static int verbose;
    static int use_native_mri;
    static int no_substepping;

    // mesh refinement
    static std::string coupling_type;
    static int do_reflux;
    static int do_avg_down;

    // Diagnostic output interval
    static int sum_interval;
    static amrex::Real sum_per;

    // Native or NetCDF
    static std::string plotfile_type;

    // init_type:  "custom", "ideal", "real"
    static std::string init_type;

    // NetCDF initialization (wrfinput) file
    static amrex::Vector<amrex::Vector<std::string>> nc_init_file;

    // NetCDF initialization (wrfbdy) file
    static std::string nc_bdy_file;

    // Text input_sounding file
    static std::string input_sounding_file;

    // 1D CDF output (for ingestion in AMR-Wind)
    static int         output_1d_column;
    static int         column_interval;
    static amrex::Real column_per;
    static amrex::Real column_loc_x;
    static amrex::Real column_loc_y;
    static std::string column_file_name;

    // 2D BndryRegister output (for ingestion in AMR-Wind)
    static int         output_bndry_planes;
    static int         bndry_output_planes_interval;
    static amrex::Real bndry_output_planes_per;
    static amrex::Real bndry_output_planes_start_time;

    // 2D BndryRegister input
    static int          input_bndry_planes;

    static int ng_dens_hse;
    static int ng_pres_hse;

    amrex::Vector<amrex::Real> h_havg_density;
    amrex::Vector<amrex::Real> h_havg_temperature;
    amrex::Vector<amrex::Real> h_havg_pressure;
    amrex::Gpu::DeviceVector<amrex::Real> d_havg_density;
    amrex::Gpu::DeviceVector<amrex::Real> d_havg_temperature;
    amrex::Gpu::DeviceVector<amrex::Real> d_havg_pressure;

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
      int nGhostCells;

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
    amrex::FluxRegister&
    get_flux_reg(int lev)
    {
      return *flux_registers[lev];
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
        numCores * (amrex::ParallelDescriptor::second() - startCPUTime) +
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
        amrex::ParallelDescriptor::Barrier("ROMSX::setRecordDataInfo");
    }

    amrex::Vector<std::unique_ptr<std::fstream> > datalog;
    amrex::Vector<std::string> datalogname;

    //! The filename of the ith datalog file.
    const std::string DataLogName (int i) const noexcept { return datalogname[i]; }

public:
    void writeJobInfo(const std::string& dir) const;
    static void writeBuildInfo(std::ostream& os);
};

#endif
