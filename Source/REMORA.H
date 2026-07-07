#ifndef REMORA_H_
#define REMORA_H_

#include <string>
#include <limits>
#include <array>
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
#include "REMORA_prob_common.H"

#ifdef REMORA_USE_PARTICLES
#include "REMORA_ParticleData.H"
#endif

#ifdef REMORA_USE_NETCDF
#include "REMORA_NCInterface.H"
#include "REMORA_NCTimeSeries.H"
#include "REMORA_NCTimeSeriesRiver.H"
#include "REMORA_NCTimeSeriesBoundary.H"
#endif

#ifdef REMORA_USE_MOAB
#include "REMORA_MOAB.H"
#endif

#include <iostream>

#ifdef AMREX_LAZY
#include <AMReX_Lazy.H>
#endif

#ifndef AMREX_USE_MPI
using amrex::MPI_COMM_WORLD;
using amrex::MPI_Comm;
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

// Forward declarations to resolve circular dependency
class ProblemBase;
#ifdef REMORA_USE_NETCDF
class NCTimeSeries;
class NCTimeSeriesRiver;
class NCTimeSeriesBoundary;
#endif

/*!
 * \brief Class that stores all relevant simulation state data with methods for time stepping
 */
class REMORA
    : public amrex::AmrCore
{
public:

    REMORA ();
    REMORA (const amrex::RealBox& rb, int max_level_in, const amrex::Vector<int>& n_cell_in, int coord, const amrex::Vector<amrex::IntVect>& ref_ratio, const amrex::Array<int,AMREX_SPACEDIM>& is_per, std::string prefix);
    virtual ~REMORA ();

    /** \brief Advance solution to final time */
    void Evolve ();

    /** \brief Write checkpoint and plotfiles at end of simulation */
    void WriteAtFinalTime ();

    /** \brief Write checkpoint and plotfiles at intermediate point of simulation, if needed */
    void WriteAtIntermediateTime (int step, amrex::Real cur_time);

    /** \brief Tag cells for refinement */
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    /** \brief Initialize multilevel data */
    void InitData ();
    amrex::Real EvolveOneStep (amrex::Real time, amrex::Real dt_request);

    /** \brief Init (NOT restart or regrid) */
    void init_only (int lev, amrex::Real time);

    /** \briefRestart */
    void restart ();

    /** \brief Called after every level 0 timestep */
    void post_timestep (int nstep, amrex::Real time, amrex::Real dt_lev);

    /**
     * \brief Extracts SST from the 3D conservative state for the atmospheric driver.
     *
     * Reads Temp_comp at the top water-column cell (k_sfc), converts from
     * Celsius to Kelvin, and copies the result into state[SSTIndex].
     *
     * @param[in,out] state  OCN2ATM slab buffer sized by the driver (one MultiFab
     *                       per ocean-to-atmosphere export layer). state[SSTIndex]
     *                       (index 0) is overwritten with sea-surface temperature
     *                       (SST) sampled from the Temp_comp tracer at the uppermost
     *                       sigma level (k = Nz) and converted from degrees Celsius
     *                       to Kelvin for the atmospheric driver. An empty vector or
     *                       null state[0] is treated as a no-op.
     * @param[in    ] time   Current ocean model time (unused; retained for driver
     *                       interface conformance).
     */
    void PackSurfaceState (amrex::Vector<amrex::MultiFab*>& state, amrex::Real time);

    /**
     * \brief Receives atmospheric states from the driver and applies unit conversions.
     *
     * Fills REMORA's internal forcing MultiFabs from states and records which
     * lanes were successfully updated in driver_atmos_state_from_driver.
     * Unit conversions applied: Pair Pa to mb; Tair K to Celsius.
     *
     * @param[in] states  ATM2OCN forcing slab buffer from the driver (Warner et al.
     *                    2010, Block B state-passing contract), indexed by AtmosState.
     *                    Expected units per lane:
     *                    Uwind/Vwind: 10-m winds [m/s];
     *                    Pair: mean sea-level pressure [Pa, converted to mb];
     *                    Qair: near-surface specific humidity [kg/kg];
     *                    Tair: 2-m air temperature [K, converted to degC];
     *                    Cloud: cloud fraction [0-1];
     *                    Rain: precipitation rate [kg/m^2/s];
     *                    SWrad/LWrad: downwelling shortwave/longwave radiation [W/m^2].
     *                    Missing lanes (null pointer or index out of range) are skipped;
     *                    driver_atmos_state_from_driver tracks populated lanes for the
     *                    bulk-flux parameterization fallback logic.
     * @param[in] time    Current ocean model time (unused; retained for driver
     *                    interface conformance).
     */
    void ApplyAtmosphericStates (const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time);

    // Diagnostics

    /** \brief Integrate conserved quantities for diagnostics */
    void sum_integrated_quantities (amrex::Real time);

    /** \brief Perform the volume-weighted sum */
    amrex::Real
    volWgtSumMF (int lev,
      const amrex::MultiFab& mf, int comp, bool local, bool finemask);

    /** \brief Decide if it is time to take an action */
    bool is_it_time_for_action (int nstep, amrex::Real time, amrex::Real dt,
                               int action_interval, amrex::Real action_per);

    /** \brief Make a new level using provided BoxArray and DistributionMapping and
     * fill with interpolated coarse level data.
     * Overrides the pure virtual function in AmrCore
     */
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                     const amrex::DistributionMapping& dm) override;

    /** \brief Remake an existing level using provided BoxArray and DistributionMapping and
     * fill with existing fine and coarse data.
     * Overrides the pure virtual function in AmrCore
     */
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                  const amrex::DistributionMapping& dm) override;

    /** \brief Delete level data
     * Overrides the pure virtual function in AmrCore
     */
    virtual void ClearLevel (int lev) override;

    /** \brief Make a new level from scratch using provided BoxArray and DistributionMapping.
     * Only used during initialization.
     * Overrides the pure virtual function in AmrCore
     */
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm) override;

    /** \brief Set pm and pn arrays and x/y coords on level lev */
    void set_grid_scale (int lev);

    /** \brief Set x/y coords on level lev based on pm and pn */
    void set_grid_coords_from_grid_scale (int lev);

    /** \brief Set curvilinear derivative terms on level lev based on pm and pn */
    void set_curvilinear_terms_from_grid_scale (int lev);

    /** \brief Set pm/pn by averaging down from higher-resolution grid */
    void set_grid_vars_averaged_down (int lev);

    /** \brief Set zeta components to be equal to time-averaged Zt_avg1 */
    void set_zeta_to_Ztavg (int lev);

    /** \brief Calculate u-, v-, and psi-point masks based on rho-point masks after analytic initialization */
    void calculate_nodal_masks (int lev);

    /** \brief Copy maskr to all z levels */
    void fill_3d_masks (int lev);

    /** \brief Set psi-point mask to be consistent with rho-point mask */
    void update_mskp (int lev);

    /** \brief compute dt from CFL considerations */
    amrex::Real estTimeStep (int lev) const;

    /** \brief Interface for advancing the data at one level by one "slow" timestep */
    void remora_advance(int level,
                        amrex::MultiFab& cons_old,  amrex::MultiFab& cons_new,
                        amrex::MultiFab& xvel_old,  amrex::MultiFab& yvel_old,  amrex::MultiFab& zvel_old,
                        amrex::MultiFab& xvel_new,  amrex::MultiFab& yvel_new,  amrex::MultiFab& zvel_new,
                        amrex::MultiFab& source,
                        const amrex::Geometry fine_geom,
                        const amrex::Real dt, const amrex::Real time);

    /** \brief Make mask to zero out covered cells (for mesh refinement) */
    amrex::MultiFab& build_fine_mask(int lev);

    /** \brief main driver for writing AMReX plotfiles */
    void WritePlotFile (int istep);

    /** \brief write out particular data to an AMReX plotfile */
    void WriteMultiLevelPlotfileWithBathymetry (const std::string &plotfilename,
                                                int nlevels,
                                                const amrex::Vector<const amrex::MultiFab*> &mf,
                                                const amrex::Vector<const amrex::MultiFab*> &mf_nd,
                                                const amrex::Vector<const amrex::MultiFab*> &mf_u,
                                                const amrex::Vector<const amrex::MultiFab*> &mf_v,
                                                const amrex::Vector<const amrex::MultiFab*> &mf_w,
                                                const amrex::Vector<const amrex::MultiFab*> &mf_2d_rho,
                                                const amrex::Vector<const amrex::MultiFab*> &mf_2d_u,
                                                const amrex::Vector<const amrex::MultiFab*> &mf_2d_v,
                                                const amrex::Vector<std::string> &varnames_3d,
                                                const amrex::Vector<std::string> &varnames_2d_rho,
                                                const amrex::Vector<std::string> &varnames_2d_u,
                                                const amrex::Vector<std::string> &varnames_2d_v,
                                                const amrex::Vector<amrex::Geometry>& my_geom,
                                                amrex::Real time,
                                                const amrex::Vector<int> &level_steps,
                                                const amrex::Vector<amrex::IntVect>& rr,
                                                const std::string &versionName = "HyperCLaw-V1.1",
                                                const std::string &levelPrefix = "Level_",
                                                const std::string &mfPrefix = "Cell",
                                                const amrex::Vector<std::string>& extra_dirs = amrex::Vector<std::string>()) const;



    /** \brief write out header data for an AMReX plotfile */
    void WriteGenericPlotfileHeaderWithBathymetry (std::ostream &HeaderFile,
                                                   int nlevels,
                                                   const amrex::Vector<amrex::BoxArray> &bArray,
                                                   const amrex::Vector<std::string> &varnames_3d,
                                                   const amrex::Vector<std::string> &varnames_2d_rho,
                                                   const amrex::Vector<std::string> &varnames_2d_u,
                                                   const amrex::Vector<std::string> &varnames_2d_v,
                                                   const amrex::Vector<amrex::Geometry>& my_geom,
                                                   amrex::Real time,
                                                   const amrex::Vector<int> &level_steps,
                                                   const amrex::Vector<amrex::IntVect>& rr,
                                                   const std::string &versionName,
                                                   const std::string &levelPrefix,
                                                   const std::string &mfPrefix) const;

    /** \brief default prefix for input file parameters */
    std::string pp_prefix {"remora"};

    /** \brief multilevel data container for last step's scalar data: temperature, salinity, passive tracer  */
    amrex::Vector<amrex::MultiFab*> cons_old;
    /** \brief multilevel data container for last step's x velocities (u in ROMS) */
    amrex::Vector<amrex::MultiFab*> xvel_old;
    /** \brief multilevel data container for last step's y velocities (v in ROMS) */
    amrex::Vector<amrex::MultiFab*> yvel_old;
    /** \brief multilevel data container for last step's z velocities (largely unused; W stored separately) */
    amrex::Vector<amrex::MultiFab*> zvel_old;

    /** \brief multilevel data container for current step's scalar data: temperature, salinity, passive tracer  */
    amrex::Vector<amrex::MultiFab*> cons_new;
    /** \brief multilevel data container for current step's x velocities (u in ROMS) */
    amrex::Vector<amrex::MultiFab*> xvel_new;
    /** \brief multilevel data container for current step's y velocities (v in ROMS) */
    amrex::Vector<amrex::MultiFab*> yvel_new;
    /** \brief multilevel data container for current step's z velocities (largely unused; W stored separately) */
    amrex::Vector<amrex::MultiFab*> zvel_new;

    /** \brief multilevel data container for high res initial data: temperature, salinity, passive tracer  */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_cons_full_domain;
    /** \brief multilevel data container for high res initial x velocities (u in ROMS) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_xvel_full_domain;
    /** \brief multilevel data container for high res initial y velocities (v in ROMS) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_yvel_full_domain;
    /** \brief multilevel data container for current step's z velocities (largely unused; W stored separately) */

    // Program state data is represented by vectors of pointers to AMReX Multifabs.
    // There is one pointer per level

    /** \brief Bathymetry data (2D, positive valued, h in ROMS) **/
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_h;

    /** \brief Bathymetry data on the whole domain at each potential level **/
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_h_full_domain;

    /** \brief Width of cells in the vertical (z-) direction (3D, Hz in ROMS) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Hz;
    /** \brief u-volume flux (3D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Huon;
    /** \brief v-volume flux (3D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Hvom;
    /** \brief u velocity RHS (3D, includes horizontal and vertical advection) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_ru;
    /** \brief v velocity RHS (3D, includes horizontal and vertical advection) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rv;
    /** \brief u velocity RHS (2D, includes horizontal and vertical advection) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_ru2d;
    /** \brief v velocity RHS (2D, includes horizontal and vertical advection) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rv2d;
    /** \brief u velocity RHS, integrated, including advection and bottom/surface stresses (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rufrc;
    /** \brief v velocity RHS, integrated, including advection and bottom/surface stresses (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rvfrc;
    /** \brief Vertical viscosity coefficient (3D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Akv;
    /** \brief Vertical diffusion coefficient (3D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Akt;
    /** \brief Harmonic viscosity defined on the psi points (corners of horizontal grid cells) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_visc2_p;
    /** \brief Harmonic viscosity defined on the rho points (centers) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_visc2_r;
    /** \brief Harmonic diffusivity for temperature / salinity */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_diff2;

    /** \brief z coordinates at rho points (cell centers) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_z_r;

    /** \brief z coordinates at w points (faces between z-cells) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_z_w;
    /** \brief Scaled vertical coordinate (range [0,1]) that transforms to z, defined at rho points (cell centers) */
    amrex::Gpu::DeviceVector<amrex::Real> s_r;
    /** \brief Scaled vertical coordinate (range [0,1]) that transforms to z, defined at w-points (cell faces) */
    amrex::Gpu::DeviceVector<amrex::Real> s_w;

    //amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_s_r;
    //amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_s_w;
    //amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Cs_r;
    //amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Cs_w;

    /** \brief Stretching coefficients at rho points */
    amrex::Gpu::DeviceVector<amrex::Real> Cs_r;
    /** \brief Stretching coefficients at w points */
    amrex::Gpu::DeviceVector<amrex::Real> Cs_w;

    /** \brief z coordinates at psi points (cell nodes) **/
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_z_phys_nd;

    /** \brief Average of the free surface, zeta (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Zt_avg1;

    /** \brief Surface stress in the u direction */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_sustr;
    /** \brief Surface stress in the v direction */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_svstr;

    /** \brief Wind in the u direction, defined at rho-points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_uwind;
    /** \brief Wind in the v direction, defined at rho-points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_vwind;
    /** \brief Air temperature [°C], defined at rho-points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Tair;
    /** \brief Specific humidity [kg/kg], defined at rho-points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_qair;
    /** \brief Air pressure [mb], defined at rho-points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Pair;

    /** \brief Shortwave radiation flux [W/m²], defined at rho-points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_srflx;
    /** \brief longwave radiation */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_lrflx;
    /** \brief Downward longwave radiation */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_longwave_down;
    /** \brief latent heat flux */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_lhflx;
    /** \brief sensible heat flux */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_shflx;

    /** \brief Surface tracer flux; working arrays */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_stflx;
    /** \brief Surface tracer flux; input arrays */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_stflux;
    /** \brief Bottom tracer flux; working arrays */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_btflx;
    /** \brief Bottom tracer flux; input arrays */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_btflux;

    /** \brief precipitation rate [kg/m^2/s] */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rain;
    /** \brief evaporation rate [kg/m^2/s] */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_evap;
    /** \brief cloud cover fraction [0-1], defined at rho-points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_cloud;
    /** \brief evaporation minus precipitation [kg/m^2/s], defined at rho-points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_EminusP;

    /** \brief provenance flags for driver-supplied atmospheric forcing lanes */
    std::array<bool, AtmosState::NumTypes> driver_atmos_state_from_driver {};

    /** \brief Linear drag coefficient [m/s], defined at rho points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rdrag;
    /** \brief Quadratic drag coefficient [unitless], defined at rho points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rdrag2;
    /** \brief Bottom roughness length [m], defined at rho points */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_ZoBot;

    /** \brief Bottom stress in the u direction */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_bustr;
    /** \brief Bottom stress in the v direction */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_bvstr;

    /** \brief time average of barotropic x velocity flux (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_DU_avg1;
    /** \brief correct time average of barotropic x velocity flux for coupling (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_DU_avg2;
    /** \brief time average of barotropic y velocity flux */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_DV_avg1;
    /** \brief correct time average of barotropic y velocity flux for coupling (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_DV_avg2;
    /** \brief barotropic x velocity for the RHS (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rubar;
    /** \brief barotropic y velocity for the RHS (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rvbar;
    /** \brief free surface height for the RHS (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rzeta;
    /** \brief barotropic x velocity (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_ubar;
    /** \brief barotropic y velocity (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_vbar;
    /** \brief free surface height (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_zeta;

    /** \brief high resolution initial free surface height (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_zeta_full_domain;

    /** \brief land/sea mask at cell centers (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_mskr;
    /** \brief land/sea mask at x-faces (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_msku;
    /** \brief land/sea mask at y-faces (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_mskv;
    /** \brief land/sea mask at cell corners (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_mskp;
    /** \brief land/sea mask at cell centers, copied to all z levels (3D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_mskr3d;

    /** \brief horizontal scaling factor: 1 / dx (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_pm;
    /** \brief horizontal scaling factor: 1 / dy (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_pn;
    /** \brief horizontal scaling factor: 1 / dx (2D) on whole domain */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_pm_full_domain;
    /** \brief horizontal scaling factor: 1 / dy (2D) on whole domain */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_pn_full_domain;

    /** \brief coriolis factor (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_fcor;

    /** \brief x_grid on rho points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_xr;
    /** \brief y_grid on rho points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_yr;

    /** \brief x_grid on u-points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_xu;
    /** \brief y_grid on u-points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_yu;

    /** \brief x_grid on v-points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_xv;
    /** \brief y_grid on v-points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_yv;

    /** \brief x_grid on psi-points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_xp;
    /** \brief y_grid on psi-points (2D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_yp;

    /** \brief d(1/n)/d(xi) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_dndx;
    /** \brief d(1/m)/d(eta) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_dmde;

    /** \brief additional scratch space for calculations on temp, salt, etc */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_sstore;

    /** \brief density perturbation */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rhoS;
    /** \brief vertically-averaged density */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_rhoA;
    /** \brief Brunt-Vaisala frequency (3D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_bvf;
    /** \brief Thermal expansion coefficient (3D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_alpha;
    /** \brief Saline contraction coefficient (3D) */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_beta;

    /** \brief Weights for calculating avg1 in 2D advance */
    amrex::Vector<amrex::Real> vec_weight1;
    /** \brief Weights for calculating avg2 in 2D advance */
    amrex::Vector<amrex::Real> vec_weight2;

    // GLS
    /** \brief Turbulent kinetic energy */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_tke;
    /** \brief Turbulent generic length scale */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_gls;
    /** \brief Vertical mixing turbulent length scale */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Lscale;
    /** \brief Turbulent kinetic energy vertical diffusion coefficient */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Akk;
    /** \brief Turbulent length scale vertical diffusion coefficient */
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vec_Akp;

    /** \brief Climatology nudging coefficients */
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> vec_nudg_coeff;

    /** \brief advance a single level for a single time step */
    void Advance (int lev, amrex::Real time, amrex::Real dt_lev, int iteration, int ncycle);

    /** \brief Set everything up for a step on a level */
    void setup_step(int lev, amrex::Real time, amrex::Real dt_lev);

    /** \brief 3D advance on a single level */
    void advance_3d_ml (int lev, amrex::Real dt_lev);

    /** \brief 2D advance, one predictor/corrector step */
    void advance_2d_onestep (int lev, amrex::Real dt_lev, amrex::Real dtfast_lev, int my_iif, int nfast_counter);

    /** \brief Perform a 2D predictor (predictor_2d_step=True) or corrector (predictor_2d_step=False) step */
    void advance_2d (int lev,
                     amrex::MultiFab const* mf_rhoS,
                     amrex::MultiFab const* mf_rhoA,
                     amrex::MultiFab      * mf_ru2d,
                     amrex::MultiFab      * mf_rv2d,
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
                     amrex::Real dtfast_lev,
                     bool predictor_2d_step,
                     bool first_2d_step, int my_iif,
                     int & next_indx1);

    /** \brief Advance the 3D variables */
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
                     amrex::MultiFab const* mf_mskr,
                     amrex::MultiFab const* mf_msku,
                     amrex::MultiFab const* mf_mskv,
                     const int N,
                     const amrex::Real dt_lev);

    /** \brief Calculate bulk temperature, salinity, wind fluxes */
    void bulk_fluxes (int lev,
                      amrex::MultiFab* mf_cons,
                      amrex::MultiFab* mf_uwind,
                      amrex::MultiFab* mf_vwind,
                      amrex::MultiFab* mf_Tair,
                      amrex::MultiFab* mf_qair,
                      amrex::MultiFab* mf_Pair,
                      amrex::MultiFab* mf_srflx,
                      amrex::MultiFab* mf_longwave_down,
                      amrex::MultiFab* mf_evap,
                      amrex::MultiFab* mf_sustr,
                      amrex::MultiFab* mf_svstr,
                      amrex::MultiFab* mf_stflux,
                      amrex::MultiFab* mf_lrflx,
                      amrex::MultiFab* mf_lhflx,
                      amrex::MultiFab* mf_shflx,
                      const int N);

    /** \brief Wrapper function for prestep */
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

    /** \brief Prestep advection calculations for the tracers */
    void prestep_t_advection (int lev, const amrex::Box& tbx,
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
                        const amrex::Array4<int         const>& river_pos,
                        const amrex::Array4<amrex::Real const>& river_source,
                        int iic, int ntfirst, int nrhs, int N,
                        const amrex::Real dt_lev);

    /** \brief RHS terms for tracer */
    void rhs_t_3d (int lev,
                   const amrex::Box& bx,
                   const amrex::Array4<amrex::Real      >& t,
                   const amrex::Array4<amrex::Real const>& tempstore,
                   const amrex::Array4<amrex::Real const>& Huon,
                   const amrex::Array4<amrex::Real const>& Hvom,
                   const amrex::Array4<amrex::Real const>& Hz,
                   const amrex::Array4<amrex::Real const>& pn,
                   const amrex::Array4<amrex::Real const>& pm,
                   const amrex::Array4<amrex::Real const>& W,
                   const amrex::Array4<amrex::Real      >& FC,
                   const amrex::Array4<amrex::Real const>& mskr,
                   const amrex::Array4<amrex::Real const>& msku,
                   const amrex::Array4<amrex::Real const>& mskv,
                   const amrex::Array4<int         const>& river_pos,
                   const amrex::Array4<amrex::Real const>& river_source,
                   int nrhs, int nnew, int N, const amrex::Real dt_lev);

    /** \brief RHS terms for 3D momentum */
    void rhs_uv_3d (int lev,
                    const amrex::Box& xbx,
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

    /** \brief RHS terms for 2D momentum */
    void rhs_uv_2d (int lev,
                    const amrex::Box& xbx,
                    const amrex::Box& ybx,
                    const amrex::Array4<amrex::Real const>& uold,
                    const amrex::Array4<amrex::Real const>& vold,
                    const amrex::Array4<amrex::Real      >& ru,
                    const amrex::Array4<amrex::Real      >& rv,
                    const amrex::Array4<amrex::Real const>& Duon,
                    const amrex::Array4<amrex::Real const>& Dvom,
                    const int nrhs);

    /** \brief Wrapper around equation of state calculation */
    void rho_eos (const amrex::Box& bx,
                  const amrex::Array4<amrex::Real const>& state,
                  const amrex::Array4<amrex::Real      >& rho,
                  const amrex::Array4<amrex::Real      >& rhoA,
                  const amrex::Array4<amrex::Real      >& rhoS,
                  const amrex::Array4<amrex::Real      >& bvf,
                  const amrex::Array4<amrex::Real      >& alpha,
                  const amrex::Array4<amrex::Real      >& beta,
                  const amrex::Array4<amrex::Real const>& Hz,
                  const amrex::Array4<amrex::Real const>& z_w,
                  const amrex::Array4<amrex::Real const>& z_r,
                  const amrex::Array4<amrex::Real const>& h,
                  const amrex::Array4<amrex::Real const>& mskr,
                  const int N);

    /** \brief Calculate density and related quantities from linear equation of state */
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

    /** \brief Calculate density and related quantities from nonlinear equation of state */
    void nonlin_eos (const amrex::Box& bx,
                  const amrex::Array4<amrex::Real const>& state,
                  const amrex::Array4<amrex::Real      >& rho,
                  const amrex::Array4<amrex::Real      >& rhoA,
                  const amrex::Array4<amrex::Real      >& rhoS,
                  const amrex::Array4<amrex::Real      >& bvf,
                  const amrex::Array4<amrex::Real      >& alpha,
                  const amrex::Array4<amrex::Real      >& beta,
                  const amrex::Array4<amrex::Real const>& Hz,
                  const amrex::Array4<amrex::Real const>& z_w,
                  const amrex::Array4<amrex::Real const>& z_r,
                  const amrex::Array4<amrex::Real const>& h,
                  const amrex::Array4<amrex::Real const>& mskr,
                  const int N);

    /** \brief Calculate pressure gradient */
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

    /** \brief Update velocities or tracers with diffusion/viscosity
     * as the last part of the prestep */
    void prestep_diffusion (const amrex::Box& bx,
                        const amrex::Box& gbx,
                        const int ioff, const int joff,
                        const amrex::Array4<amrex::Real      >& vel,
                        const amrex::Array4<amrex::Real const>& vel_old,
                        const amrex::Array4<amrex::Real      >& rvel,
                        const amrex::Array4<amrex::Real const>& Hz,
                        const amrex::Array4<amrex::Real const>& Akv,
                        const amrex::Array4<amrex::Real      >& FC,
                        const amrex::Array4<amrex::Real const>& sstr,
                        const amrex::Array4<amrex::Real const>& bstr,
                        const amrex::Array4<amrex::Real const>& z_r,
                        const amrex::Array4<amrex::Real const>& pm,
                        const amrex::Array4<amrex::Real const>& pn,
                        const int iic, const int ntfirst, const int nnew, int nstp, int nrhs, int N,
                        const amrex::Real lambda, const amrex::Real dt_lev);

    /** \brief Calculate effects of vertical viscosity or diffusivity */
    void vert_visc_3d (const amrex::Box& bx,
                       const int ioff, const int joff,
                       const amrex::Array4<amrex::Real      >& phi,
                       const amrex::Array4<amrex::Real const>& Hz,
                       const amrex::Array4<amrex::Real      >& Hzk,
                       const amrex::Array4<amrex::Real      >& AK,
                       const amrex::Array4<amrex::Real const>& Akv,
                       const amrex::Array4<amrex::Real      >& BC,
                       const amrex::Array4<amrex::Real      >& DC,
                       const amrex::Array4<amrex::Real      >& FC,
                       const amrex::Array4<amrex::Real      >& CF,
                       const int nnew, const int N,
                       const amrex::Real dt_lev);

    /** \brief Correct mass flux */
    void update_massflux_3d (int lev, const amrex::Box& bx,
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

    /** \brief Adjust 3D momentum variables based on vertical mean momentum */
    void vert_mean_3d (const amrex::Box& bx,
                       const int ioff, const int joff,
                       const amrex::Array4<amrex::Real      >& phi,
                       const amrex::Array4<amrex::Real const>& Hz,
                       const amrex::Array4<amrex::Real const>& Dphi_avg1,
                       const amrex::Array4<amrex::Real      >& DC,
                       const amrex::Array4<amrex::Real      >& CF,
                       const amrex::Array4<amrex::Real const>& pm_or_pn,
                       const amrex::Array4<amrex::Real const>& msk,
                       const int nnew, const int N);

    /** \brief Harmonic viscosity */
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

    /** \brief Wrapper for harmonic diffusivity for tracers */
    void t3dmix2 (const amrex::Box& bx,
                  const amrex::Array4<amrex::Real      >& state,
                  const amrex::Array4<amrex::Real      >& state_rhs,
                  const amrex::Array4<amrex::Real const>& diff2,
                  const amrex::Array4<amrex::Real const>& Hz,
                  const amrex::Array4<amrex::Real const>& z_r,
                  const amrex::Array4<amrex::Real const>& pm,
                  const amrex::Array4<amrex::Real const>& pn,
                  const amrex::Array4<amrex::Real const>& msku,
                  const amrex::Array4<amrex::Real const>& mskv,
                  const amrex::Real dt_lev,
                  const int ncomp, const int N);

    /** \brief Harmonic diffusivity for tracers along S-coordinate level surfaces */
    void t3dmix2_s(const amrex::Box& bx,
                   const amrex::Array4<amrex::Real      >& state,
                   const amrex::Array4<amrex::Real      >& state_rhs,
                   const amrex::Array4<amrex::Real const>& diff2,
                   const amrex::Array4<amrex::Real const>& Hz,
                   const amrex::Array4<amrex::Real const>& pm,
                   const amrex::Array4<amrex::Real const>& pn,
                   const amrex::Array4<amrex::Real const>& msku,
                   const amrex::Array4<amrex::Real const>& mskv,
                   const amrex::Real dt_lev, const int ncomp);

    /** \brief Harmonic diffusivity for tracers along geopotential surfaces */
    void t3dmix2_geo (const amrex::Box& bx,
                      const amrex::Array4<amrex::Real      >& state,
                      const amrex::Array4<amrex::Real      >& state_rhs,
                      const amrex::Array4<amrex::Real const>& diff2,
                      const amrex::Array4<amrex::Real const>& Hz,
                      const amrex::Array4<amrex::Real const>& z_r,
                      const amrex::Array4<amrex::Real const>& pm,
                      const amrex::Array4<amrex::Real const>& pn,
                      const amrex::Array4<amrex::Real const>& msku,
                      const amrex::Array4<amrex::Real const>& mskv,
                      const amrex::Real dt_lev,
                      const int ncomp, const int N);

    /** \brief Calculate Coriolis terms */
    void coriolis (const amrex::Box& xbx,
                   const amrex::Box& ybx,
                   const amrex::Array4<amrex::Real const>& uold,
                   const amrex::Array4<amrex::Real const>& vold,
                   const amrex::Array4<amrex::Real      >& ru,
                   const amrex::Array4<amrex::Real      >& rv,
                   const amrex::Array4<amrex::Real const>& Hz,
                   const amrex::Array4<amrex::Real const>& fomn,
                   int nrhs, int nr);

    /** \brief Calculate curvilinear advection terms */
    void curvilinear (const amrex::Box& bx,
                   const amrex::Box& xbx,
                   const amrex::Box& ybx,
                   const amrex::Array4<amrex::Real const>& uold,
                   const amrex::Array4<amrex::Real const>& vold,
                   const amrex::Array4<amrex::Real      >& ru,
                   const amrex::Array4<amrex::Real      >& rv,
                   const amrex::Array4<amrex::Real const>& Hz,
                   const amrex::Array4<amrex::Real const>& dndx,
                   const amrex::Array4<amrex::Real const>& dmde,
                   int nrhs, int nr);

    /** \brief Apply climatology nudging */
    void apply_clim_nudg (const amrex::Box& bx,
                   int ioff, int joff,
                   const amrex::Array4<amrex::Real      >& var,
                   const amrex::Array4<amrex::Real const>& var_old,
                   const amrex::Array4<amrex::Real const>& var_clim,
                   const amrex::Array4<amrex::Real const>& clim_coeff,
                   const amrex::Array4<amrex::Real const>& Hz,
                   const amrex::Array4<amrex::Real const>& pm,
                   const amrex::Array4<amrex::Real const>& pn,
                   const amrex::Real dt_lev = amrex::Real(0.0));

    /** \brief Set 2D momentum arrays from 3D momentum */
    void set_2darrays (int lev);

    /** \brief Set Zt_avg1 to zeta */
    void set_zeta_average (int lev);

    /** \brief Initialize zeta from file or analytic */
    void set_zeta (int lev);

    /** \brief Copy over zeta data that has been averaged down from high res */
    void set_zeta_averaged_down (int lev);

    /** \brief Initialize bathymetry from file or analytic */
    void set_bathymetry (int lev);

    /** \brief Copy over bathymetry data that has been averaged down from high resolution
     * input netcdf file */
    void set_bathymetry_averaged_down (int lev);

    /** \brief Initialize Coriolis factor from file or analytic */
    void set_coriolis (int lev);

    /** \brief Initialize land-sea masks from file or analytic */
    void set_masks (int lev);

    /** \brief Calculate vertical stretched coordinates */
    void stretch_transform (int lev);

    /** \brief Initialize vertical mixing coefficients from file or analytic */
    void init_set_vmix (int lev);
    /** \brief Set vertical mixing coefficients from analytic */
    void set_analytic_vmix (int lev);
    /** \brief Initialize GLS variables */
    void init_gls_vmix (int lev, SolverChoice solver_choice);
    /** \brief Prestep for GLS calculation */
    void gls_prestep (int lev, amrex::MultiFab* mf_gls, amrex::MultiFab* mf_tke,
                      amrex::MultiFab& mf_W,
                      amrex::MultiFab* mf_msku, amrex::MultiFab* mf_mskv,
                      const int nstp, const int nnew, const int iic, const int ntfirst,
                      const int N, const amrex::Real dt_lev);
    /** \brief Corrector step for GLS calculation */
    void gls_corrector (int lev, amrex::MultiFab* mf_gls, amrex::MultiFab* mf_tke,
                        amrex::MultiFab& mf_W, amrex::MultiFab* mf_Akv, amrex::MultiFab* mf_Akt,
                        amrex::MultiFab* mf_Akk, amrex::MultiFab* mf_Akp,
                        amrex::MultiFab* mf_mskr,
                        amrex::MultiFab* mf_msku, amrex::MultiFab* mf_mskv,
                        const int nstp, const int nnew,
                        const int N, const amrex::Real dt_lev);

    /** \brief Scale RHS momentum variables by 1/cell area, needed before FillPatch to different levels*/
    void scale_rhs_vars ();
    /** \brief Scale RHS momentum variables by cell area, needed after FillPatch to different levels*/
    void scale_rhs_vars_inv ();

    /** \brief Initialize or calculate surface momentum flux from file or analytic */
    void set_smflux (int lev);

    /** \brief Initialize or calculate wind speed from file or analytic */
    void set_wind (int lev);

    /** \brief Initialize horizontal mixing coefficients */
    void set_hmixcoef (int lev);

    /** \brief Initialize flat bathymetry to value from problo */
    void init_flat_bathymetry (int lev);

    /** \brief Initialize or calculate bottom drag */
    void set_drag (int lev);

    /** \brief Set weights for averaging 3D variables to 2D */
    void set_weights (int lev);

    /** \brief Fill the physical boundary conditions for cell-centered velocity (diagnostic only) */
    void FillBdyCCVels (int lev, amrex::MultiFab& mf_cc_vel);

    /** \brief Fill a new MultiFab by copying in phi from valid region and filling ghost cells */
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

    /** \brief Fill a new MultiFab by copying in phi from valid region and filling ghost cells without applying boundary conditions */
    void FillPatchNoBC (int lev, amrex::Real time,
                    amrex::MultiFab& mf_to_be_filled,
                    amrex::Vector<amrex::MultiFab*> const& mfs,
                    const int bdy_var_type = BdyVars::null,
                    const int icomp=0,
                    const bool fill_all=true,
                    const bool fill_set=true);
    /** \brief Fill boundary data from netcdf file */
    void fill_from_bdyfiles (int lev,
                             amrex::MultiFab& mf_to_fill,
                            const amrex::MultiFab& mf_mask,
                            const amrex::Real time,
                            const int bccomp,
                            const int bdy_var_type,
                            const int icomp_to_fill,
                            const int icomp_calc = 0,
                            const amrex::MultiFab& mf_calc = amrex::MultiFab(),
                            const amrex::Real = amrex::Real(0.0));

    /** \brief fill an entire multifab by interpolating from the coarser level
     *         using the piecewise constant interpolater */
    void FillCoarsePatchPC (int lev, amrex::Real time, amrex::MultiFab* mf_fine,
                amrex::MultiFab* mf_crse,
                const int bccomp,
                const int bdy_var_type = BdyVars::null,
                const int  icomp = 0,
                const bool fill_all = true,
                const int  n_not_fill=0,
                const int icomp_calc=0,
                const amrex::Real dt = amrex::Real(0.0),
                const amrex::MultiFab& mf_calc = amrex::MultiFab());

    /** \brief fill an entire multifab by interpolating from the coarser level */
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab* mf_fine,
                amrex::MultiFab* mf_crse,
                const int bccomp,
                const int bdy_var_type = BdyVars::null,
                const int  icomp = 0,
                const bool fill_all = true,
                const int  n_not_fill=0,
                const int icomp_calc=0,
                const amrex::Real dt = amrex::Real(0.0),
                const amrex::MultiFab& mf_calc = amrex::MultiFab());

    /** \brief fill an entire multifab by interpolating from the coarser level,
     *         explicitly specifying interpolator to use */
    void FillCoarsePatchMap (int lev, amrex::Real time, amrex::MultiFab* mf_fine,
                amrex::MultiFab* mf_crse,
                const int bccomp,
                const int bdy_var_type = BdyVars::null,
                const int  icomp = 0,
                const bool fill_all = true,
                const int  n_not_fill=0,
                const int icomp_calc=0,
                const amrex::Real dt = amrex::Real(0.0),
                const amrex::MultiFab& mf_calc = amrex::MultiFab(),
                      amrex::Interpolater* mapper = nullptr);

    /** \brief initialize and calculate stretch coefficients */
    void init_stretch_coeffs ();
    /** \brief calculate vertical stretch coefficients */
    void calc_stretch_coeffs ();

    /** \brief Calculate Coriolis parameters from beta plane parametrization */
    void init_beta_plane_coriolis (int lev);

    /** \brief Full domain bathymetry data initialization from analytic */
    void init_bathymetry_full_domain_from_analytic ();

    /** \brief Problem initialization from averaged-down high resolution data */
    void set_init_data_averaged_down (int lev);
    /** \brief Problem initialization from NetCDF file */
    void init_data_from_netcdf (int lev);
    /** \brief High resolution roblem initialization from NetCDF file */
    void init_data_full_domain_from_netcdf ();
    /** \brief Boundary data initialization from NetCDF file */
    void init_bdry_from_netcdf (int lev);
    /** \brief Mask data initialization from NetCDF file */
    void init_masks_from_netcdf (int lev);
    /** \brief Bathymetry data initialization from NetCDF file */
    void init_bathymetry_from_netcdf (int lev);
    /** \brief Full domain high-res bathymetry data initialization from NetCDF file */
    void init_bathymetry_full_domain_from_netcdf ();

    /** \brief Grid variable initialization from NetCDF file */
    void init_grid_vars_from_netcdf (int lev);
    /** \brief Full domain high-res grid variable initialization from NetCDF file */
    void init_grid_vars_full_domain_from_netcdf ();

    /** \brief Extrapolate grid metrics to edge of MultiFab */
    void extrapolate_metric_to_physical_boundaries(amrex::MultiFab& mf, const amrex::Geometry& geom);

    /** \brief Sea-surface height data initialization from NetCDF file */
    void init_zeta_from_netcdf (int lev);
    /** \brief Full-domain high res sea-surface height data initialization from NetCDF file */
    void init_zeta_full_domain_from_netcdf ();
    /** \brief Coriolis parameter data initialization from NetCDF file */
    void init_coriolis_from_netcdf (int lev);
    /** \brief Climatology nudging coefficient initialization from NetCDF file */
    void init_clim_nudg_coeff_from_netcdf (int lev);
    /** \brief Wrapper to initialize climatology nudging coefficient */
    void init_clim_nudg_coeff (int lev);
    void init_riv_pos_from_netcdf (int lev);

    /** \brief Allocate multifabs for storing full-domain bathymetry and grid vars data */
    void allocate_bathymetry_grid_vars_full_domain ();
    /** \brief Allocate multifabs for storing full-domain high resolution initial data */
    void allocate_init_full_domain ();

    /** \brief Convert data in a multifab from inverse days to inverse seconds */
    void convert_inv_days_to_inv_s (amrex::MultiFab*);

    /** \brief Mask data arrays before writing output */
    void mask_arrays_for_write (int lev, amrex::Real fill_value, amrex::Real fill_where);

    /** \brief Average down from level lev+1 to lev in mf, including grow cells */
    void average_down_with_grow_cells(int lev, amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mf);

    int xvel_bc () const noexcept { return BCVars::xvel_bc(ncons); }
    int yvel_bc () const noexcept { return BCVars::yvel_bc(ncons); }
    int zvel_bc () const noexcept { return BCVars::zvel_bc(ncons); }
    int ubar_bc () const noexcept { return BCVars::ubar_bc(ncons); }
    int vbar_bc () const noexcept { return BCVars::vbar_bc(ncons); }
    int zeta_bc () const noexcept { return BCVars::zeta_bc(ncons); }
    int tke_bc  () const noexcept { return BCVars::tke_bc(ncons); }
    int foextrap_periodic_bc () const noexcept { return BCVars::foextrap_periodic_bc(ncons); }
    int foextrap_bc () const noexcept { return BCVars::foextrap_bc(ncons); }
    int u2d_simple_bc () const noexcept { return BCVars::u2d_simple_bc(ncons); }
    int v2d_simple_bc () const noexcept { return BCVars::v2d_simple_bc(ncons); }
    int num_bc_vars () const noexcept { return BCVars::NumTypes(ncons); }


#ifdef REMORA_USE_NETCDF
    static int total_nc_plot_file_step;
#endif

private:

    ///////////////////////////
    // private member functions
    ///////////////////////////

    /** \brief read in some parameters from inputs file */
    void ReadParameters();

    /** \brief Build runtime scalar names after nscalar is known */
    void init_scalar_metadata();

    /** \brief set covered coarse cells to be the average of overlying fine cells */
    void AverageDown ();

    /** \brief Read in boundary parameters from input file and set up data structures */
    void init_bcs();

    /** \brief Initialize initial problem data from analytic functions */
    void init_analytic(int lev);
    /** \brief Initialize high resolution initial problem data from analytic functions */
    void init_full_domain_from_analytic();

    /** \brief Initialize high resolution initial sea surface height from analytic functions */
    void init_full_domain_zeta_from_analytic();

    /** \brief Allocate MultiFabs for state and evolution variables */
    void init_stuff (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    /** \brief Allocate MultiFabs for masks */
    void init_masks (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    /** \brief Resize variable containers to accommodate data on levels 0 to max_lev */
    void resize_stuff (int lev);

    /** \brief more flexible version of AverageDown() that lets you average down across multiple levels */
    void AverageDownTo (int crse_lev);

    /** \brief Construct FillPatchers */
    void Construct_REMORAFillPatchers (int lev);

    /** \brief Define FillPatchers */
    void Define_REMORAFillPatchers (int lev);

    /** \brief utility to copy in data from old and/or new state into another multifab */
    TimeInterpolatedData GetDataAtTime (int lev, amrex::Real time);

    /** \brief advance a level by dt,  includes a recursive call for finer levels */
    void timeStep (int lev, amrex::Real time, int iteration);

    /** \brief advance all levels by dt, loops over finer levels */
    void timeStepML (amrex::Real time, int iteration);

    /** \brief a wrapper for estTimeStep() */
    void ComputeDt ();

    /** \brief get plotfile name */
    std::string PlotFileName (int lev) const;

    /* set which variables and derived quantities go into plotfiles */
    void set2DPlotVariables (const std::string& pp_plot_var_names_2d);
    void set3DPlotVariables (const std::string& pp_plot_var_names_3d);

    /* append variables to plot */
    void append2DPlotVariables (const std::string& pp_plot_var_names_2d);
    void append3DPlotVariables (const std::string& pp_plot_var_names_3d);

#ifdef REMORA_USE_NETCDF
    /** \brief Write plotfile using NetCDF (wrapper) */
    void WriteNCPlotFile       (int istep, amrex::MultiFab const* plotMF);

    /** \brief Write a particular NetCDF plotfile */
    void WriteNCPlotFile_which (int lev, int which_subdomain,
                                amrex::MultiFab const* plotMF,
                                bool write_header, ncutils::NCFile& ncf,
                                bool is_history);

    /** \brief Write MultiFab in NetCDF format */
    void WriteNCMultiFab (const amrex::FabArray<amrex::FArrayBox> &fab,
                          const std::string& name,
                          bool set_ghost = false) const;

    /** \brief Read MultiFab in NetCDF format */
    void ReadNCMultiFab (amrex::FabArray<amrex::FArrayBox> &fab,
                         const std::string &name,
                         int coordinatorProc = amrex::ParallelDescriptor::IOProcessorNumber(),
                         int allow_empty_mf = 0);

    /** \brief Start time in the time series of boundary data */
    amrex::Real start_bdy_time;
    /** \brief Interval between boundary data times */
    amrex::Real bdy_time_interval;

    /** \brief Whether to output NetCDF files as a single history file with several time steps */
    static bool write_history_file;

    /** \brief Data container for u-component surface momentum flux read from file */
    std::unique_ptr<NCTimeSeries> sustr_data_from_file;
    /** \brief Data container for v-component surface momentum flux read from file */
    std::unique_ptr<NCTimeSeries> svstr_data_from_file;
    /** \brief Data container for u-direction wind read from file */
    std::unique_ptr<NCTimeSeries> Uwind_data_from_file;
    /** \brief Data container for v-direction wind read from file */
    std::unique_ptr<NCTimeSeries> Vwind_data_from_file;
    /** \brief Data container for air temperature read from file */
    std::unique_ptr<NCTimeSeries> Tair_data_from_file;
    /** \brief Data container for specific humidity read from file */
    std::unique_ptr<NCTimeSeries> qair_data_from_file;
    /** \brief Data container for air pressure read from file */
    std::unique_ptr<NCTimeSeries> Pair_data_from_file;
    /** \brief Data container for shortwave radiation flux read from file */
    std::unique_ptr<NCTimeSeries> srflx_data_from_file;
    /** \brief Data container for downward longwave radiation flux read from file */
    std::unique_ptr<NCTimeSeries> longwave_down_data_from_file;
    /** \brief Data container for precipitation rate read from file */
    std::unique_ptr<NCTimeSeries> rain_data_from_file;
    /** \brief Data container for cloud cover fraction read from file */
    std::unique_ptr<NCTimeSeries> cloud_data_from_file;
    /** \brief Data container for evaporation minus precipitation read from file */
    std::unique_ptr<NCTimeSeries> EminusP_data_from_file;

    /** \brief Data container for ubar climatology data read from file */
    std::unique_ptr<NCTimeSeries> ubar_clim_data_from_file;
    /** \brief Data container for vbar climatology data read from file */
    std::unique_ptr<NCTimeSeries> vbar_clim_data_from_file;
    /** \brief Data container for u-velocity climatology data read from file */
    std::unique_ptr<NCTimeSeries> u_clim_data_from_file;
    /** \brief Data container for v-velocity climatology data read from file */
    std::unique_ptr<NCTimeSeries> v_clim_data_from_file;
    /** \brief Data container for temperature climatology data read from file */
    std::unique_ptr<NCTimeSeries> temp_clim_data_from_file;
    /** \brief Data container for salinity climatology data read from file */
    std::unique_ptr<NCTimeSeries> salt_clim_data_from_file;

    /** \brief Vector of data containers for scalar data in rivers */
    amrex::Vector<std::unique_ptr<NCTimeSeriesRiver>> river_source_cons;
    /** \brief Data container for momentum transport in rivers */
    std::unique_ptr<NCTimeSeriesRiver> river_source_transport;
    /** \brief Data container for vertically integrated momentum transport in rivers */
    std::unique_ptr<NCTimeSeriesRiver> river_source_transportbar;

    /** \brief Vector over BdyVars of boundary series data containers */
    amrex::Vector<amrex::Vector<std::unique_ptr<NCTimeSeriesBoundary>>> boundary_series;

#endif // REMORA_USE_NETCDF
    /** \brief iMultiFab for river positions; contents are indices of rivers */
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> vec_river_position;
    /** \brief Vector over rivers of river direction: 0: u-face; 1: v-face; 2: w-face */
    amrex::Gpu::DeviceVector<int> river_direction;


#ifdef REMORA_USE_MOAB
    /** \brief Initialize connection to MOAB */
    void InitMOABMesh();
#endif
    /** \brief Nudging width at coarse-fine interface */
    int cf_width{0};
    /** \brief Width for fixing values at coarse-fine interface */
    int cf_set_width{0};
    // Fillpatcher classes for coarse-fine boundaries
    /** \brief Vector over levels of FillPatchers for scalars */
    amrex::Vector<REMORAFillPatcher> FPr_c;
    /** \brief Vector over levels of FillPatchers for u (3D) */
    amrex::Vector<REMORAFillPatcher> FPr_u;
    /** \brief Vector over levels of FillPatchers for v (3D) */
    amrex::Vector<REMORAFillPatcher> FPr_v;
    /** \brief Vector over levels of FillPatchers for w */
    amrex::Vector<REMORAFillPatcher> FPr_w;

    // Fillpatcher classes for 2d velocity variables
    /** \brief Vector over levels of FillPatchers for ubar (2D) */
    amrex::Vector<REMORAFillPatcher> FPr_ubar;
    /** \brief Vector over levels of FillPatchers for vbar (2D) */
    amrex::Vector<REMORAFillPatcher> FPr_vbar;

    /** \brief write checkpoint file to disk */
    void WriteCheckpointFile ();

    /** \brief read checkpoint file from disk */
    void ReadCheckpointFile ();

    /** \brief Read the file passed to remora.restart and use it as an initial condition for the current simulation */
    void InitializeFromFile ();

    /** \brief Initialize the new-time data at a level from the initial_data MultiFab */
    void InitializeLevelFromData (int lev, const amrex::MultiFab& initial_data);

    /** \brief utility to skip to next line in Header */
    static void GotoNextLine (std::istream& is);

    ////////////////
    // private data members

    /** \brief Pointer to container of analytical functions for problem definition */
    std::unique_ptr<ProblemBase> prob = nullptr;

    /** \brief how many boxes specified at each level by tagging criteria */
    amrex::Vector<int> num_boxes_at_level;
    /** \brief how many netcdf input files specified at each level */
    amrex::Vector<int> num_files_at_level;
    /** \brief the boxes specified at each level by tagging criteria */
    amrex::Vector<amrex::Vector<amrex::Box>> boxes_at_level;

    /** \brief which step? */
    amrex::Vector<int> istep;
    /** \brief How many substeps on each level? */
    amrex::Vector<int> nsubsteps;
    /** \brief new time at each level */
    amrex::Vector<amrex::Real> t_new;
    /** \brief old time at each level */
    amrex::Vector<amrex::Real> t_old;
    /** \brief time step at each level */
    amrex::Vector<amrex::Real> dt;

    /** \brief whether to set boundary conditions by variable rather than just by side */
    bool set_bcs_by_var;

    /** \brief Vector (over level) of functors to apply physical boundary conditions */
    amrex::Vector<std::unique_ptr<REMORAPhysBCFunct>> physbcs;

    /** \brief array of flux registers for refluxing in multilevel */
    amrex::Vector<std::unique_ptr<amrex::YAFluxRegister>> advflux_reg;

    ////////////////
    // boundary data members

    // A BCRec is essentially a 2*DIM integer array storing the boundary
    // condition type at each lo/hi walls in each direction. We have one BCRec
    // for each component of the cell-centered variables and each velocity component.
    /** \brief vector (over BCVars) of BCRecs */
    amrex::Vector           <amrex::BCRec> domain_bcs_type;
    /** \brief GPU vector (over BCVars) of BCRecs */
    amrex::Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d;

    /** \brief Array of strings describing domain boundary conditions */
    amrex::Array<std::string,2*AMREX_SPACEDIM> domain_bc_type;

    /** \brief Array holding the Dirichlet values at walls which need them */
    amrex::Vector<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2>> m_bc_extdir_vals;

    /** \brief Array holding the "physical" boundary condition types (e.g. "inflow") */
    amrex::Vector<amrex::GpuArray<REMORA_BC, AMREX_SPACEDIM*2>> phys_bc_type;

    /** \brief These are flags that indicate whether we need to read in boundary data from file */
    amrex::GpuArray<amrex::GpuArray<bool, AMREX_SPACEDIM*2>,BdyVars::NumTypes+1> phys_bc_need_data;

    /** \brief Container to connect boundary data being read in boundary condition containers.
     *
     * Needed for phys_bc_need_data to work correctly
     */
    amrex::Vector<int> bdy_index;

    /** \brief Step when we last output a plotfile */
    int last_plot_file_step;
    /** \brief Simulation time when we last output a plotfile */
    amrex::Real last_plot_file_time;

    /** \brief Step when we last output a checkpoint file */
    int last_check_file_step;
    /** \brief Simulation time when we last output a checkpoint file */
    amrex::Real last_check_file_time;
    /** \brief Whether to output a plotfile on restart from checkpoint */
    int plot_file_on_restart = 1;

    ////////////////
    // runtime parameters

    /** \brief maximum number of steps */
    int max_step = std::numeric_limits<int>::max();
    /** \brief Time to stop */
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    /** \brief Time of the start of the simulation, in seconds */
    amrex::Real start_time = 0.0;

    /** \brief If set, restart from this checkpoint file */
    std::string restart_chkfile = "";

    /** \brief Number of passive scalars carried in the state */
    int nscalar = 1;
    /** \brief Number of conserved scalars in the state (temperature + salt + passive scalars) */
    int ncons = Tracer_comp + 1;

    // Time step controls
    /** \brief CFL condition */
    static amrex::Real cfl;
    /** \brief Fraction maximum change in subsequent time steps */
    static amrex::Real change_max;
    /** \brief User specified fixed baroclinic time step */
    static amrex::Real fixed_dt;
    /** \brief User specified fixed barotropic time step */
    static amrex::Real fixed_fast_dt;
    /** \brief User specified, number of barotropic steps per baroclinic step */
    static int fixed_ndtfast_ratio;
    /** \brief Number of fast steps to take */
    int nfast;

    /** \brief Whether to substep fine levels in time */
    int do_substep = 0;

    /** \brief how often each level regrids the higher levels of refinement (after a level advances that many time steps) */
    int regrid_int = 2;

    // I/O controls
    /** \brief Plotfile prefix */
    std::string plot_file_name {"plt_"};
    /** \brief Plotfile output interval in iterations */
    int plot_int = -1;
    /** \brief Plotfile output interval in seconds */
    amrex::Real plot_int_time = amrex::Real(-1.0);
    /** \brief Checkpoint file prefix */
    std::string check_file {"chk"};
    /** \brief Checkpoint output interval in iterations */
    int check_int = -1;
    /** \brief Checkpoint output interval in seconds */
    amrex::Real check_int_time = amrex::Real(-1.0);

    /** \brief Whether to chunk netcdf history file */
    bool chunk_history_file = false;
    /** \brief Number of time steps per netcdf history file.
     *
        -1 is the default and means code will automatically compute
         number of steps per file based on grid size to keep files < 2 GB */
    int steps_per_history_file = -1;
    /** \brief Counter for which time index we are writing to in the netcdf history file */
    int history_count = 0;

    /** \brief Names of 3D variables to output to AMReX plotfile */
    amrex::Vector<std::string> plot_var_names_3d;
    /** \brief Names of 2D variables to output to AMReX plotfile */
    amrex::Vector<std::string> plot_var_names_2d;
    /** \brief Names of scalars for plotfile output */
    amrex::Vector<std::string> cons_names {"temp", "salt", "tracer"};
    // Note that the order of variable names here must match the order in PlotFile.cpp
    /** \brief Names of derived fields for plotfiles */
    const amrex::Vector<std::string> derived_names {"vorticity"};

    /** \brief Container for algorithmic choices */
    static SolverChoice solverChoice;

    /** \brief whether plotfile variables should be expanded to a uniform refinement ratio */
    bool expand_plotvars_to_unif_rr = false;

    /** \brief fill value for masked arrays in amrex plotfiles */
    amrex::Real plotfile_fill_value = amrex::Real(0.0);
    /** \brief fill value for masked arrays in netcdf output */
    amrex::Real netcdf_fill_value = amrex::Real(1.0e37);


#ifdef REMORA_USE_PARTICLES
    /** \brief Particle container with all particle species */
    ParticleData particleData;

    // variables and functions for tracers particles
    /** \brief tracer particles that advect with flow */
    bool  m_use_tracer_particles;
    /** \brief tracer particles that fall with gravity */
    bool  m_use_hydro_particles;

    /** \brief Read tracer and hydro particles parameters */
    void readTracersParams();

    /** \brief Initialize tracer and hydro particles */
    void initializeTracers ( amrex::ParGDBBase*,
                             const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& );

    /** \brief Evolve tracers and hydro particles */
    void evolveTracers( int lev,
                        amrex::Real dt,
                        amrex::Vector<amrex::MultiFab const*>& flow_vel,
                        const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& );

#endif
    /** \brief Verbosity level of output */
    static int verbose;

    /** \brief Diagnostic sum output interval in number of steps */
    static int sum_interval;
    /** \brief Diagnostic sum output interval in time */
    static amrex::Real sum_per;

    /** \brief Minimum number of digits in plotfile name or chunked history file */
    static int file_min_digits;

      /** \brief Whether to write the staggered velocities (not averaged to cell centers) */
      static bool plot_staggered_vels;

      /** \brief Whether to write nodal data (Nu_nd) to plotfiles */
      static bool plot_nodal_data;

     /** \brief Native or NetCDF plotfile output */
     static PlotfileType plotfile_type;

    /** \brief NetCDF initialization file */
    static amrex::Vector<amrex::Vector<std::string>> nc_init_file;
    /** \brief NetCDF grid file */
    static amrex::Vector<amrex::Vector<std::string>> nc_grid_file;
    /** \brief NetCDF boundary data */
    static amrex::Vector<std::string> nc_bdry_file;

    /** \brief Which level the high resolution bathymetry is at */
    int hires_grid_level = -1;
    /** \brief Grid file for high resolution bathymetry */
    std::string nc_grid_file_hires;
    /** \brief Box for the full domain at nc_hires_grid_level */
    amrex::Box nc_hires_grid_box;

    /** \brief Which level the high resolution initialization data is at */
    int hires_init_level = -1;
    /** \brief Init file for high resolution */
    std::string nc_init_file_hires;
    /** \brief Box for the full domain at nc_hires_init_level */
    amrex::Box nc_hires_init_box;

    /** \brief Cumulative refinement ratio between level 0 and level i */
    amrex::Vector<amrex::IntVect> cum_ref_ratios;

    /** \brief NetCDF forcing file(s) */
    amrex::Vector<std::string> nc_frc_file;
    /** \brief NetCDF river file(s) */
    amrex::Vector<std::string> nc_riv_file;

    /** \brief NetCDF climatology history file(s) */
    amrex::Vector<std::string> nc_clim_his_file;
    /** \brief NetCDF climatology coefficient file */
    std::string nc_clim_coeff_file;

    /** \brief Default name of time field for boundary data */
    std::string bdry_time_varname = "ocean_time";
    /** \brief Name of time fields for boundary data */
    amrex::Vector<std::string> bdry_time_name_byvar;

    /** \brief Name of time field for ubar climatology data */
    std::string clim_ubar_time_varname = "ocean_time";
    /** \brief Name of time field for vbar climatology data */
    std::string clim_vbar_time_varname = "ocean_time";
    /** \brief Name of time field for u climatology data */
    std::string clim_u_time_varname = "ocean_time";
    /** \brief Name of time field for v climatology data */
    std::string clim_v_time_varname = "ocean_time";
    /** \brief Name of time field for salinity climatology data */
    std::string clim_salt_time_varname = "ocean_time";
    /** \brief Name of time field for temperature climatology data */
    std::string clim_temp_time_varname = "ocean_time";
    /** \brief Name of time field for river time */
    std::string riv_time_varname = "river_time";
    /** \brief Name of time field for forcing data */
    std::string frc_time_varname = "";

    /** \brief Set refinement criteria */
    void refinement_criteria_setup();

    /** \brief Holds info for dynamically generated tagging criteria */
    static amrex::Vector<amrex::AMRErrorTag> ref_tags;

    /** \brief Mask that zeroes out values on a coarse level underlying grids on the next finest level */
    amrex::MultiFab fine_mask;

    /** \brief Helper function to determine number of ghost cells */
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

    /** \brief Helper function for IO stream */
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

    /** \brief Variable for CPU timing */
    static amrex::Real startCPUTime;
    /** \brief Accumulator variable for CPU time used thusfar */
    static amrex::Real previousCPUTimeUsed;

    /** \brief Get CPU time used */
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

    /** \brief The filename of the ith datalog file. */
    const std::string DataLogName (int i) const noexcept { return datalogname[i]; }

public:
    /** \brief Write job info to stdout */
    void writeJobInfo (const std::string& dir) const;
    /** \brief Write build info to os */
    static void writeBuildInfo (std::ostream& os);

    static void print_banner(MPI_Comm /*comm*/, std::ostream& /*out*/);
    static void print_usage(MPI_Comm /*comm*/, std::ostream& /*out*/);
    static void print_error(MPI_Comm /*comm*/, const std::string& msg);
    static void print_summary(std::ostream&);
    static void print_tpls(std::ostream& /*out*/);

    /** \brief Accessor method for t_old to expose to outside classes */
    amrex::Real get_t_old(int lev) const;

    /** \brief Evaluate stability function psi for wind speed */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real bulk_psiu(amrex::Real ZoL) {
        // Compute stability function, psi
        using namespace amrex;

        Real psiu;
        if (ZoL < 0.0_rt) {
            // Unstable conditions.
            Real x    = std::pow(1.0_rt-15.0_rt*ZoL,0.25_rt);
            Real psik = 2.0_rt*std::log(0.5_rt*(1.0_rt+x))+
                        std::log(0.5_rt*(1.0_rt+x*x))-
                        2.0_rt*std::atan(x)+0.5_rt*PI;
            //  For very unstable conditions, use free-convection (Fairall).
            Real sqrt3 = std::sqrt(3.0_rt);
            Real y = std::pow(1.0_rt-10.15_rt*ZoL,1.0_rt/3.0_rt);
            Real psic = 1.5_rt*std::log(1.0_rt/3.0_rt*(1.0_rt+y+y*y))-
                   sqrt3*std::atan((1.0_rt+2.0_rt*y)/sqrt3)+PI/sqrt3;
            //  Match Kansas and free-convection forms with weighting Fw.
            Real Zol2 = ZoL*ZoL;
            Real Fw = Zol2/(1.0_rt+Zol2);

            psiu =  (1.0_rt-Fw)*psik+Fw*psic;
        } else {
            // Stable conditions
            Real cff=std::min(50.0_rt,0.35_rt*ZoL);
            psiu = -((1.0_rt+ZoL)+0.6667_rt*(ZoL-14.28_rt)/
                     std::exp(cff)+8.525_rt);
        }
        return psiu;
    }

    /** \brief Evaluate stability function psi for moisture and heat */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real bulk_psit(amrex::Real ZoL) {
        // Compute stability function, psi
        using namespace amrex;

        Real psit;
        if (ZoL < 0.0_rt) {
            // Unstable conditions.
            Real x = std::pow(1.0_rt-15.0_rt*ZoL,0.5_rt);
            Real psik = 2.0_rt*std::log(0.5_rt*(1.0_rt+x));
            //  For very unstable conditions, use free-convection (Fairall).
            Real sqrt3=std::sqrt(3.0_rt);
            Real y=std::pow(1.0_rt-34.15_rt*ZoL,1.0_rt/3.0_rt);
            Real psic=1.5_rt*std::log(1.0_rt/3.0_rt*(1.0_rt+y+y*y))-
                sqrt3*std::atan((1.0_rt+2.0_rt*y)/sqrt3)+PI/sqrt3;

            // Match Kansas and free-convection forms with weighting Fw.
            Real ZoL2=ZoL*ZoL;
            Real Fw=ZoL2/(1.0_rt+ZoL2);
            psit = (1.0_rt-Fw)*psik+Fw*psic;
        } else {
            //  Stable conditions.
            Real cff=std::min(50.0_rt,0.35_rt*ZoL);
            psit=-(std::pow(1.0_rt+2.0_rt*ZoL,1.5_rt)+
                    0.6667_rt*(ZoL-14.28_rt)/std::exp(cff)+8.525_rt);
        }
        return psit;
    }
};

#endif
