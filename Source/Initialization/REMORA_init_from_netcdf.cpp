/**
 * \file REMORA_init_from_netcdf.cpp
 */

#include <REMORA.H>
#include <REMORA_Constants.H>
#include <REMORA_prob_common.H>
#include <REMORA_DataStruct.H>

using namespace amrex;

#ifdef REMORA_USE_NETCDF

/** \brief helper function for reading in initial state data from netcdf */
void
read_data_from_netcdf (int /*lev*/, const Box& domain, const std::string& fname,
                       FArrayBox& NC_temp_fab, FArrayBox& NC_salt_fab,
                       FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                       FArrayBox& NC_ubar_fab, FArrayBox& NC_vbar_fab);

/** \brief helper function for reading in land-sea masks from netcdf */
void
read_masks_from_netcdf (int /*lev*/, const Box& domain, const std::string& fname,
                       FArrayBox& NC_mskr_fab, FArrayBox& NC_msku_fab,
                       FArrayBox& NC_mskv_fab);

/** \brief helper function for reading boundary data from netcdf */
Real
read_bdry_from_netcdf (const Box& domain, const std::string& fname,
                       Vector<Vector<FArrayBox>>& bdy_data_xlo,
                       Vector<Vector<FArrayBox>>& bdy_data_xhi,
                       Vector<Vector<FArrayBox>>& bdy_data_ylo,
                       Vector<Vector<FArrayBox>>& bdy_data_yhi,
                       int& width, amrex::Real& start_bdy_time,
                       std::string bdry_time_varname,
                       amrex::GpuArray<amrex::GpuArray<bool, AMREX_SPACEDIM*2>,BdyVars::NumTypes+1>&);

/** \brief helper function to initialize state from netcdf */
void
init_state_from_netcdf (int lev,
                        FArrayBox&  temp_fab, FArrayBox&  salt_fab,
                        FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                        FArrayBox&  ubar_fab, FArrayBox&  vbar_fab,
                        const Vector<FArrayBox>& NC_temp_fab,
                        const Vector<FArrayBox>& NC_salt_fab,
                        const Vector<FArrayBox>& NC_xvel_fab,
                        const Vector<FArrayBox>& NC_yvel_fab,
                        const Vector<FArrayBox>& NC_ubar_fab,
                        const Vector<FArrayBox>& NC_vbar_fab);

/** \brief helper function to read bathymetry from netcdf */
void
read_bathymetry_from_netcdf (int lev, const Box& domain, const std::string& fname,
                             FArrayBox& NC_h_fab,
                             FArrayBox& NC_pm_fab, FArrayBox& NC_pn_fab,
                             FArrayBox& NC_xr_fab, FArrayBox& NC_yr_fab,
                             FArrayBox& NC_xu_fab, FArrayBox& NC_yu_fab,
                             FArrayBox& NC_xv_fab, FArrayBox& NC_yv_fab,
                             FArrayBox& NC_xp_fab, FArrayBox& NC_yp_fab);

/** \brief helper function to read coriolis factor from netcdf */
void
read_coriolis_from_netcdf (int lev, const Box& domain, const std::string& fname, FArrayBox& NC_fcor_fab);

/** \brief helper function to read sea surface height from netcdf */
void
read_zeta_from_netcdf (int lev, const Box& domain, const std::string& fname,
                             FArrayBox& NC_zeta_fab);

/** \brief helper function to read climatology nudging from netcdf */
void
read_clim_nudg_coeff_from_netcdf (int lev, const Box& domain, const std::string& fname,
                                  bool do_m2_clim_nudg,
                                  bool do_m3_clim_nudg,
                                  bool do_temp_clim_nudg,
                                  bool do_salt_clim_nudg,
                                  FArrayBox& NC_M2NC_fab,
                                  FArrayBox& NC_M3NC_fab,
                                  FArrayBox& NC_TempNC_fab,
                                  FArrayBox& NC_SaltNC_fab);

/** \brief helper function to read in vector of data from netcdf */
void read_vec_from_netcdf (int lev, const amrex::Vector<std::string>& fnames, const std::string& field_name, amrex::Vector<int>& vec_dat);

/**
 * @param lev Integer specifying the current level
 */
void
REMORA::init_data_from_netcdf (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_temp_fab ; NC_temp_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_salt_fab ; NC_salt_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_xvel_fab ; NC_xvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_yvel_fab ; NC_yvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_ubar_fab ; NC_ubar_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_vbar_fab ; NC_vbar_fab.resize(num_boxes_at_level[lev]);

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_data_from_netcdf(lev, boxes_at_level[lev][idx], nc_init_file[lev][idx],
                              NC_temp_fab[idx], NC_salt_fab[idx],
                              NC_xvel_fab[idx], NC_yvel_fab[idx],
                              NC_ubar_fab[idx], NC_vbar_fab[idx]);
    }


    MultiFab mf_temp(*cons_new[lev], make_alias, Temp_comp, 1);
    MultiFab mf_salt(*cons_new[lev], make_alias, Salt_comp, 1);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
    // Don't tile this since we are operating on full FABs in this routine
    for ( MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi )
    {
        // Define fabs for holding the initial data
        FArrayBox &temp_fab = mf_temp[mfi];
        FArrayBox &salt_fab = mf_salt[mfi];
        FArrayBox &xvel_fab = (*xvel_new[lev])[mfi];
        FArrayBox &yvel_fab = (*yvel_new[lev])[mfi];
        FArrayBox &ubar_fab = (*vec_ubar[lev])[mfi];
        FArrayBox &vbar_fab = (*vec_vbar[lev])[mfi];

        init_state_from_netcdf(lev, temp_fab, salt_fab,
                               xvel_fab, yvel_fab,
                               ubar_fab, vbar_fab,
                               NC_temp_fab, NC_salt_fab,
                               NC_xvel_fab, NC_yvel_fab,
                               NC_ubar_fab, NC_vbar_fab);
    } // mf
    } // omp
}

/**
 * @param lev Integer specifying the current level
 */
void
REMORA::init_zeta_from_netcdf (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_zeta_fab     ; NC_zeta_fab.resize(num_boxes_at_level[lev]);

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_zeta_from_netcdf(lev,boxes_at_level[lev][idx], nc_init_file[lev][idx],
                                    NC_zeta_fab[idx]);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
        // Don't tile this since we are operating on full FABs in this routine
        for ( MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi )
        {
            FArrayBox &zeta_fab  = (*vec_zeta[lev])[mfi];

            //
            // FArrayBox to FArrayBox copy does "copy on intersection"
            // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
            //

            zeta_fab.template    copy<RunOn::Device>(NC_zeta_fab[idx],0,0,1);
        } // mf
        } // omp
    } // idx

    vec_zeta[lev]->FillBoundary(geom[lev].periodicity());
    (*physbcs[lev])(*vec_zeta[lev],*vec_mskr[lev].get(),0,1,vec_zeta[lev]->nGrowVect(),t_new[lev],BCVars::zeta_bc,0,*vec_zeta[lev],*vec_msku[lev],*vec_mskv[lev]);
//    (*physbcs[lev])(*vec_zeta[lev],*vec_mskr[lev].get(),1,1,vec_zeta[lev]->nGrowVect(),t_new[lev],BCVars::zeta_bc);
//    (*physbcs[lev])(*vec_zeta[lev],*vec_mskr[lev].get(),2,1,vec_zeta[lev]->nGrowVect(),t_new[lev],BCVars::zeta_bc);

    if (solverChoice.boundary_from_netcdf) {
        Real told = t_new[lev];
        fill_from_bdyfiles(lev, *vec_zeta[lev], *vec_mskr[lev], told, BCVars::zeta_bc,BdyVars::zeta,0,0,
                           *vec_zeta[lev]);
    }
    if (lev>0) {
        FillPatch(lev, t_old[lev], *vec_zeta[lev], GetVecOfPtrs(vec_zeta), BCVars::zeta_bc, BdyVars::zeta,
                  0, false,false,0,0,0.0,*vec_zeta[lev]);
    }
//    fill_from_bdyfiles(lev, *vec_zeta[lev], *vec_mskr[lev], told, BCVars::zeta_bc,BdyVars::zeta,1,1);
//    fill_from_bdyfiles(lev, *vec_zeta[lev], *vec_mskr[lev], told, BCVars::zeta_bc,BdyVars::zeta,2,2);
}
/**
 * @param lev Integer specifying the current level
 */
void
REMORA::init_bathymetry_from_netcdf (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_h_fab     ; NC_h_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_pm_fab    ; NC_pm_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_pn_fab    ; NC_pn_fab.resize(num_boxes_at_level[lev]);

    Vector<FArrayBox> NC_xr_fab    ; NC_xr_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_yr_fab    ; NC_yr_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_xu_fab    ; NC_xu_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_yu_fab    ; NC_yu_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_xv_fab    ; NC_xv_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_yv_fab    ; NC_yv_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_xp_fab    ; NC_xp_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_yp_fab    ; NC_yp_fab.resize(num_boxes_at_level[lev]);

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_bathymetry_from_netcdf(lev, boxes_at_level[lev][idx], nc_grid_file[lev][idx],
                                    NC_h_fab[idx],
                                    NC_pm_fab[idx], NC_pn_fab[idx],
                                    NC_xr_fab[idx], NC_yr_fab[idx],
                                    NC_xu_fab[idx], NC_yu_fab[idx],
                                    NC_xv_fab[idx], NC_yv_fab[idx],
                                    NC_xp_fab[idx], NC_yp_fab[idx]);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
        // Don't tile this since we are operating on full FABs in this routine
        for ( MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi )
        {
            FArrayBox &h_fab     = (*vec_h[lev])[mfi];
            FArrayBox &pm_fab    = (*vec_pm[lev])[mfi];
            FArrayBox &pn_fab    = (*vec_pn[lev])[mfi];
            FArrayBox &xr_fab    = (*vec_xr[lev])[mfi];
            FArrayBox &yr_fab    = (*vec_yr[lev])[mfi];
            FArrayBox &xu_fab    = (*vec_xu[lev])[mfi];
            FArrayBox &yu_fab    = (*vec_yu[lev])[mfi];
            FArrayBox &xv_fab    = (*vec_xv[lev])[mfi];
            FArrayBox &yv_fab    = (*vec_yv[lev])[mfi];
            FArrayBox &xp_fab    = (*vec_xp[lev])[mfi];
            FArrayBox &yp_fab    = (*vec_yp[lev])[mfi];

            //
            // FArrayBox to FArrayBox copy does "copy on intersection"
            // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
            //

            // Copy into both components of h
            h_fab.template     copy<RunOn::Device>(NC_h_fab[idx],0,0,1);
            h_fab.template     copy<RunOn::Device>(NC_h_fab[idx],0,1,1);

            pm_fab.template    copy<RunOn::Device>(NC_pm_fab[idx]);
            pn_fab.template    copy<RunOn::Device>(NC_pn_fab[idx]);

            xr_fab.template    copy<RunOn::Device>(NC_xr_fab[idx]);
            yr_fab.template    copy<RunOn::Device>(NC_yr_fab[idx]);
            xu_fab.template    copy<RunOn::Device>(NC_xu_fab[idx]);
            yu_fab.template    copy<RunOn::Device>(NC_yu_fab[idx]);
            xv_fab.template    copy<RunOn::Device>(NC_xv_fab[idx]);
            yv_fab.template    copy<RunOn::Device>(NC_yv_fab[idx]);
            xp_fab.template    copy<RunOn::Device>(NC_xp_fab[idx]);
            yp_fab.template    copy<RunOn::Device>(NC_yp_fab[idx]);
        } // mf
        } // omp
    } // idx

    const double dummy_time = 0.0_rt;
    // Unconditional foextrap will overwrite periodicity, but EnforcePeriodicity will
    // be called on h afterwards
    FillPatch(lev,dummy_time,*vec_h[lev],GetVecOfPtrs(vec_h),
            BCVars::foextrap_periodic_bc,
            BdyVars::null,0,false,true,1);
    FillPatch(lev,dummy_time,*vec_h[lev],GetVecOfPtrs(vec_h),
            BCVars::foextrap_periodic_bc,
            BdyVars::null,1,false,true,1);

    if (lev > 0) {
        FillPatch(lev,dummy_time,*vec_pm[lev],GetVecOfPtrs(vec_pm),
                BCVars::foextrap_periodic_bc,
                BdyVars::null,0,false,true);
        FillPatch(lev,dummy_time,*vec_pn[lev],GetVecOfPtrs(vec_pn),
                BCVars::foextrap_periodic_bc,
                BdyVars::null,0,false,true);
    }


    int ng = vec_pm[lev]->nGrow();

    const auto& dom_lo = amrex::lbound(geom[lev].Domain());
    const auto& dom_hi = amrex::ubound(geom[lev].Domain());

    //
    // We need values of pm and pn outside the domain so we fill
    //    them here with foextrap
    //
    // We first fill interior ghost cells because we will need to extrapolate
    //    from ghost cells inside the domain to ghost cells outside the domain
    //
    vec_pm[lev]->FillBoundary(geom[lev].periodicity());
    vec_pn[lev]->FillBoundary(geom[lev].periodicity());

    vec_xr[lev]->FillBoundary(geom[lev].periodicity());
    vec_yr[lev]->FillBoundary(geom[lev].periodicity());
    vec_xu[lev]->FillBoundary(geom[lev].periodicity());
    vec_yu[lev]->FillBoundary(geom[lev].periodicity());
    vec_xv[lev]->FillBoundary(geom[lev].periodicity());
    vec_yv[lev]->FillBoundary(geom[lev].periodicity());
    vec_xp[lev]->FillBoundary(geom[lev].periodicity());
    vec_yp[lev]->FillBoundary(geom[lev].periodicity());

    for ( MFIter mfi(*vec_pm[lev]); mfi.isValid(); ++mfi )
    {
        Box bx   = mfi.tilebox();

        auto pm_fab = vec_pm[lev]->array(mfi);
        auto pn_fab = vec_pn[lev]->array(mfi);

        Box gbx_lox = adjCellLo(bx,0,ng); gbx_lox.grow(1,ng); gbx_lox.setBig  (0,dom_lo.x-2);
        Box gbx_hix = adjCellHi(bx,0,ng); gbx_hix.grow(1,ng); gbx_hix.setSmall(0,dom_hi.x+2);
        Box gbx_loy = adjCellLo(bx,1,ng); gbx_loy.grow(0,ng); gbx_loy.setBig  (1,dom_lo.y-2);
        Box gbx_hiy = adjCellHi(bx,1,ng); gbx_hiy.grow(0,ng); gbx_hiy.setSmall(1,dom_hi.y+2);

        // if (gbx_lox.ok()) amrex::AllPrint() << "GBX_XLO " << gbx_lox << std::endl;
        // if (gbx_hix.ok()) amrex::AllPrint() << "GBX_XHI " << gbx_hix << std::endl;
        // if (gbx_loy.ok()) amrex::AllPrint() << "GBX_YLO " << gbx_loy << std::endl;
        // if (gbx_hiy.ok()) amrex::AllPrint() << "GBX_YHI " << gbx_hiy << std::endl;

        if (gbx_lox.ok()) {
            ParallelFor(gbx_lox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                pm_fab(i,j,k,0) = pm_fab(dom_lo.x-1,j,k,0);
                pn_fab(i,j,k,0) = pn_fab(dom_lo.x-1,j,k,0);
            });
        }
        if (gbx_hix.ok()) {
            ParallelFor(gbx_hix, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                pm_fab(i,j,k,0) = pm_fab(dom_hi.x+1,j,k,0);
                pn_fab(i,j,k,0) = pn_fab(dom_hi.x+1,j,k,0);
            });
        }
        if (gbx_loy.ok()) {
            ParallelFor(gbx_loy, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                pm_fab(i,j,k,0) = pm_fab(i,dom_lo.y-1,k,0);
                pn_fab(i,j,k,0) = pn_fab(i,dom_lo.y-1,k,0);
            });
        }
        if (gbx_hiy.ok()) {
            ParallelFor(gbx_hiy, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                pm_fab(i,j,k,0) = pm_fab(i,dom_hi.y+1,k,0);
                pn_fab(i,j,k,0) = pn_fab(i,dom_hi.y+1,k,0);
            });
        }
    } // mfi

    vec_h[lev]->FillBoundary(geom[lev].periodicity());
}

/**
 * @param lev Integer specifying the current level
 */
void
REMORA::init_coriolis_from_netcdf (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_fcor_fab     ; NC_fcor_fab.resize(num_boxes_at_level[lev]);

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_coriolis_from_netcdf(lev,boxes_at_level[lev][idx], nc_grid_file[lev][idx],
                                    NC_fcor_fab[idx]);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
        // Don't tile this since we are operating on full FABs in this routine
        for ( MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi )
        {
            FArrayBox &fcor_fab  = (*vec_fcor[lev])[mfi];

            //
            // FArrayBox to FArrayBox copy does "copy on intersection"
            // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
            //

            fcor_fab.template    copy<RunOn::Device>(NC_fcor_fab[idx]);
        } // mf
        } // omp
    } // idx
    vec_fcor[lev]->FillBoundary(geom[lev].periodicity());
}

/**
 * @param lev Integer specifying the current level
 */
void
REMORA::init_masks_from_netcdf (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_mskr_fab     ; NC_mskr_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_msku_fab     ; NC_msku_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_mskv_fab     ; NC_mskv_fab.resize(num_boxes_at_level[lev]);

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_masks_from_netcdf(lev,boxes_at_level[lev][idx], nc_grid_file[lev][idx],
                                    NC_mskr_fab[idx],NC_msku_fab[idx],
                                    NC_mskv_fab[idx]);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
        // Don't tile this since we are operating on full FABs in this routine
        for ( MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi )
        {
            FArrayBox &mskr_fab  = (*vec_mskr[lev])[mfi];
            FArrayBox &msku_fab  = (*vec_msku[lev])[mfi];
            FArrayBox &mskv_fab  = (*vec_mskv[lev])[mfi];

            //
            // FArrayBox to FArrayBox copy does "copy on intersection"
            // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
            //

            mskr_fab.template    copy<RunOn::Device>(NC_mskr_fab[idx]);
            msku_fab.template    copy<RunOn::Device>(NC_msku_fab[idx]);
            mskv_fab.template    copy<RunOn::Device>(NC_mskv_fab[idx]);
        } // mf
        } // omp
    } // idx

    update_mskp(lev);
    vec_mskr[lev]->FillBoundary(geom[lev].periodicity());
    vec_msku[lev]->FillBoundary(geom[lev].periodicity());
    vec_mskv[lev]->FillBoundary(geom[lev].periodicity());
    vec_mskp[lev]->FillBoundary(geom[lev].periodicity());
}

/**
 * @param lev Integer specifying the current level
 */
void
REMORA::init_bdry_from_netcdf (int lev)
{
    if (nc_bdry_file.empty() || nc_bdry_file[0].empty()) {
        amrex::Error("NetCDF boundary file name must be provided via input");
    }

    amrex::Vector<std::string> field_name  = {"u", "v", "temp", "salt", "ubar", "vbar", "zeta"};
    amrex::Vector<IntVect    > index_types = {IntVect(1,0,0), IntVect(0,1,0),
                                             IntVect(0,0,0), IntVect(0,0,0),
                                             IntVect(1,0,0), IntVect(0,1,0),
                                             IntVect(0,0,0)};
    std::vector<bool       > is_2d       = {false, false, false, false, true, true, true};

    amrex::Print() << "DOING INIT AT LEVEL " << lev << std::endl;
    int rx = 1; int ry = 1;
    if (lev > 0) {
        for (int k = lev-1; k >= 0; k--) {
            rx *= ref_ratio[k][0];
            ry *= ref_ratio[k][1];
        }
    }
    for (int ivar = 0; ivar < BdyVars::NumTypes; ivar++) {
        boundary_series[lev].push_back(std::unique_ptr<NCTimeSeriesBoundary>(new NCTimeSeriesBoundary(lev, geom, nc_bdry_file, field_name[ivar],
                                                                bdry_time_name_byvar[ivar],
                                                                index_types[ivar],
                                                                &phys_bc_need_data[ivar], is_2d[ivar], rx, ry)));
        boundary_series[lev][ivar]->Initialize();
    }
}

/**
 * \brief Helper function to initialize state and velocity data in a Fab from a REMORAdataset.
 *
 * @param lev Integer specifying current level
 * @param state_fab FArrayBox object holding the state data we initialize
 * @param temp_fab  FArrayBox object holding the temperature data we initialize
 * @param salt_fab  FArrayBox object holding the salt        data we initialize
 * @param x_vel_fab FArrayBox object holding the x-velocity data we initialize
 * @param y_vel_fab FArrayBox object holding the y-velocity data we initialize
 * @param ubar_fab  FArrayBox object holding the ubar       data we initialize
 * @param vbar_fab  FArrayBox object holding the vbar       data we initialize
 * @param zeta_fab  FArrayBox object holding the zeta       data we initialize
 * @param NC_temp_fab Vector of FArrayBox objects with the REMORA dataset specifying temperature
 * @param NC_salt_fab Vector of FArrayBox objects with the REMORA dataset specifying salinity
 * @param NC_xvel_fab Vector of FArrayBox objects with the REMORA dataset specifying x-velocity
 * @param NC_yvel_fab Vector of FArrayBox objects with the REMORA dataset specifying y-velocity
 * @param NC_ubar_fab Vector of FArrayBox objects with the REMORA dataset specifying ubar
 * @param NC_vbar_fab Vector of FArrayBox objects with the REMORA dataset specifying vbar
 * @param NC_zeta_fab Vector of FArrayBox objects with the REMORA dataset specifying zeta
 */
void
init_state_from_netcdf (int /*lev*/,
                        FArrayBox&  temp_fab, FArrayBox&  salt_fab,
                        FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                        FArrayBox&  ubar_fab, FArrayBox&  vbar_fab,
                        const Vector<FArrayBox>& NC_temp_fab,
                        const Vector<FArrayBox>& NC_salt_fab,
                        const Vector<FArrayBox>& NC_xvel_fab,
                        const Vector<FArrayBox>& NC_yvel_fab,
                        const Vector<FArrayBox>& NC_ubar_fab,
                        const Vector<FArrayBox>& NC_vbar_fab)
{
    int nboxes = NC_xvel_fab.size();
    for (int idx = 0; idx < nboxes; idx++)
    {
        //
        // FArrayBox to FArrayBox copy does "copy on intersection"
        // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
        //
        temp_fab.template copy<RunOn::Device>(NC_temp_fab[idx]);
        salt_fab.template copy<RunOn::Device>(NC_salt_fab[idx]);
        x_vel_fab.template copy<RunOn::Device>(NC_xvel_fab[idx]);
        y_vel_fab.template copy<RunOn::Device>(NC_yvel_fab[idx]);
        ubar_fab.template copy<RunOn::Device>(NC_ubar_fab[idx],0,0,1);
        vbar_fab.template copy<RunOn::Device>(NC_vbar_fab[idx],0,0,1);
    } // idx
}

/**
 * @param lev Integer specifying the current level
 */
void
REMORA::init_clim_nudg_coeff_from_netcdf (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_M2NC_fab     ; NC_M2NC_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_M3NC_fab     ; NC_M3NC_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_TempNC_fab   ; NC_TempNC_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_SaltNC_fab   ; NC_SaltNC_fab.resize(num_boxes_at_level[lev]);

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_clim_nudg_coeff_from_netcdf(lev,boxes_at_level[lev][idx], nc_clim_coeff_file,
                                         solverChoice.do_m2_clim_nudg,
                                         solverChoice.do_m3_clim_nudg,
                                         solverChoice.do_temp_clim_nudg,
                                         solverChoice.do_salt_clim_nudg,
                                         NC_M2NC_fab[idx],NC_M3NC_fab[idx],
                                         NC_TempNC_fab[idx],NC_SaltNC_fab[idx]);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
        // Don't tile this since we are operating on full FABs in this routine
        for ( MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi )
        {
            if (solverChoice.do_m2_clim_nudg) {
                FArrayBox &ubarNC_fab  = (*vec_nudg_coeff[BdyVars::ubar][lev])[mfi];
                ubarNC_fab.template    copy<RunOn::Device>(NC_M2NC_fab[idx]);
                FArrayBox &vbarNC_fab  = (*vec_nudg_coeff[BdyVars::vbar][lev])[mfi];
                vbarNC_fab.template    copy<RunOn::Device>(NC_M2NC_fab[idx]);
            }
            if (solverChoice.do_m3_clim_nudg) {
                FArrayBox &uNC_fab  = (*vec_nudg_coeff[BdyVars::u][lev])[mfi];
                uNC_fab.template    copy<RunOn::Device>(NC_M3NC_fab[idx]);
                FArrayBox &vNC_fab  = (*vec_nudg_coeff[BdyVars::v][lev])[mfi];
                vNC_fab.template    copy<RunOn::Device>(NC_M3NC_fab[idx]);
            }
            if (solverChoice.do_temp_clim_nudg) {
                FArrayBox &TempNC_fab  = (*vec_nudg_coeff[BdyVars::t][lev])[mfi];
                TempNC_fab.template    copy<RunOn::Device>(NC_TempNC_fab[idx]);
            }
            if (solverChoice.do_salt_clim_nudg) {
                FArrayBox &SaltNC_fab  = (*vec_nudg_coeff[BdyVars::s][lev])[mfi];
                SaltNC_fab.template    copy<RunOn::Device>(NC_SaltNC_fab[idx]);
            }

        } // mf
        } // omp
    } // idx

    if (solverChoice.do_m2_clim_nudg) {
        vec_nudg_coeff[BdyVars::ubar][lev]->FillBoundary(geom[lev].periodicity());
        vec_nudg_coeff[BdyVars::vbar][lev]->FillBoundary(geom[lev].periodicity());
        convert_inv_days_to_inv_s(vec_nudg_coeff[BdyVars::ubar][lev].get());
        convert_inv_days_to_inv_s(vec_nudg_coeff[BdyVars::vbar][lev].get());
    }
    if (solverChoice.do_m3_clim_nudg) {
        vec_nudg_coeff[BdyVars::u][lev]->FillBoundary(geom[lev].periodicity());
        vec_nudg_coeff[BdyVars::v][lev]->FillBoundary(geom[lev].periodicity());
        convert_inv_days_to_inv_s(vec_nudg_coeff[BdyVars::u][lev].get());
        convert_inv_days_to_inv_s(vec_nudg_coeff[BdyVars::v][lev].get());
    }
    if (solverChoice.do_temp_clim_nudg) {
        vec_nudg_coeff[BdyVars::t][lev]->FillBoundary(geom[lev].periodicity());
        convert_inv_days_to_inv_s(vec_nudg_coeff[BdyVars::t][lev].get());
    }
    if (solverChoice.do_salt_clim_nudg) {
        vec_nudg_coeff[BdyVars::s][lev]->FillBoundary(geom[lev].periodicity());
        convert_inv_days_to_inv_s(vec_nudg_coeff[BdyVars::s][lev].get());
    }
}

/**
 * @param[in   ] lev     level to read in river data
 */
void
REMORA::init_riv_pos_from_netcdf (int lev)
{
    amrex::Vector<int> river_pos_x;
    amrex::Vector<int> river_pos_y;
    amrex::Vector<int> river_direction_tmp;

    std::string river_x_name = "river_Xposition";
    std::string river_y_name = "river_Eposition";
    std::string river_dir_name = "river_direction";

    read_vec_from_netcdf(lev, nc_riv_file, river_x_name, river_pos_x);
    read_vec_from_netcdf(lev, nc_riv_file, river_y_name, river_pos_y);
    read_vec_from_netcdf(lev, nc_riv_file, river_dir_name, river_direction_tmp);

    int nriv = river_pos_x.size();
    amrex::Gpu::DeviceVector<int> xpos_d(nriv);
    amrex::Gpu::DeviceVector<int> ypos_d(nriv);
    river_direction.resize(nriv);
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy(xpos_d.data(), river_pos_x.data(), sizeof(int)*nriv);
    Gpu::htod_memcpy(ypos_d.data(), river_pos_y.data(), sizeof(int)*nriv);
    Gpu::htod_memcpy(river_direction.data(), river_direction_tmp.data(), sizeof(int)*nriv);
#else
    std::memcpy(xpos_d.data(), river_pos_x.data(), sizeof(int)*nriv);
    std::memcpy(ypos_d.data(), river_pos_y.data(), sizeof(int)*nriv);
    std::memcpy(river_direction.data(), river_direction_tmp.data(), sizeof(int)*nriv);
#endif
    const int* xpos_ptr = xpos_d.data();
    const int* ypos_ptr = ypos_d.data();

    for (amrex::MFIter mfi(*(vec_river_position[lev]).get(),true); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.growntilebox(amrex::IntVect(NGROW,NGROW,0));
        auto river_pos = vec_river_position[lev]->array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int ) {
            for (int iriv=0; iriv < nriv; iriv++) {
                int xriv = xpos_ptr[iriv]-1;
                int yriv = ypos_ptr[iriv]-1;
                if (i==xriv && j==yriv) {
                    river_pos(i,j,0) = iriv;
                }
            }
        });
    }
}

/**
 * @param[inout] mf    multifab of data to convert
 */
void
REMORA::convert_inv_days_to_inv_s (MultiFab* mf) {
    Real inv_days_to_inv_s = 1.0_rt / (3600._rt * 24._rt);

    for ( MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& arr = mf->array(mfi);
        Box bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            arr(i,j,k) *= inv_days_to_inv_s;
        });
    }

}

#endif // REMORA_USE_NETCDF
