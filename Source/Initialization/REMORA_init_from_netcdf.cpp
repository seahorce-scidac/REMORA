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
                       FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab);

/** \brief helper function for reading in full domain high-resolution initial state data from netcdf */
void
read_data_full_domain_from_netcdf (int /*lev*/, const Box& domain, const std::string& fname,
                       FArrayBox& NC_temp_fab, FArrayBox& NC_salt_fab,
                       FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                       IntVect ngrow);

/** \brief helper function for reading in initial biology data from netcdf */
void
read_biology_from_netcdf (int /*lev*/, const Box& domain, const std::string& fname,
                          const Vector<std::string>& biology_names,
                          Vector<FArrayBox>& NC_biology_fab);

/** \brief helper function for reading in full domain high-resolution initial biology data from netcdf */
void
read_biology_full_domain_from_netcdf (int /*lev*/, const Box& domain, const std::string& fname,
                                      const Vector<std::string>& biology_names,
                                      Vector<FArrayBox>& NC_biology_fab, IntVect ngrow);

/** \brief helper function for reading in land-sea masks from netcdf */
void
read_masks_from_netcdf (int /*lev*/, const Box& domain, const std::string& fname,
                       FArrayBox& NC_mskr_fab, FArrayBox& NC_msku_fab,
                       FArrayBox& NC_mskv_fab);

/** \brief helper function to initialize state from netcdf */
void
init_state_from_netcdf (int lev,
                        FArrayBox&  temp_fab, FArrayBox&  salt_fab,
                        FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                        const Vector<FArrayBox>& NC_temp_fab,
                        const Vector<FArrayBox>& NC_salt_fab,
                        const Vector<FArrayBox>& NC_xvel_fab,
                        const Vector<FArrayBox>& NC_yvel_fab);

/** \brief helper function to initialize biology tracer state from netcdf */
void
init_biology_state_from_netcdf (int lev, FArrayBox& biology_fab,
                                const Vector<Vector<FArrayBox>>& NC_biology_fab);

/** \brief helper function to read bathymetry from netcdf */
void
read_bathymetry_from_netcdf (int lev, const Box& domain, const std::string& fname,
                             FArrayBox& NC_h_fab);

/** \brief helper function to read grid variables from netcdf */
void
read_grid_vars_from_netcdf  (int lev, const Box& domain, const std::string& fname,
                             FArrayBox& NC_pm_fab, FArrayBox& NC_pn_fab,
                             FArrayBox& NC_xr_fab, FArrayBox& NC_yr_fab,
                             FArrayBox& NC_xu_fab, FArrayBox& NC_yu_fab,
                             FArrayBox& NC_xv_fab, FArrayBox& NC_yv_fab,
                             FArrayBox& NC_xp_fab, FArrayBox& NC_yp_fab);

/** \brief helper function to read optional spherical psi coordinates from netcdf */
bool
read_spherical_grid_vars_from_netcdf (int lev, const Box& domain, const std::string& fname,
                                      FArrayBox& NC_lonp_fab, FArrayBox& NC_latp_fab);

/** \brief helper function to read full-domain high resolution bathymetry from netcdf */
void
read_bathymetry_full_domain_from_netcdf (const Box& domain, const std::string& fname,
                                         FArrayBox& NC_h_fab, IntVect ngrow);

/** \brief helper function to read full-domain high resolution grid variables from netcdf */
void
read_grid_vars_full_domain_from_netcdf (const Box& domain, const std::string& fname,
                                         FArrayBox& NC_pm_fab, FArrayBox& NC_pn_fab,
                                         IntVect ngrow);

/** \brief helper function to read coriolis factor from netcdf */
void
read_coriolis_from_netcdf (int lev, const Box& domain, const std::string& fname, FArrayBox& NC_fcor_fab);

/** \brief helper function to read sea surface height from netcdf */
void
read_zeta_from_netcdf (int lev, const Box& domain, const std::string& fname,
                             FArrayBox& NC_zeta_fab);

/** \brief helper function to read high-resolution full-domain sea surface height from netcdf */
void
read_zeta_full_domain_from_netcdf (int lev, const Box& domain, const std::string& fname,
                                   FArrayBox& NC_zeta_fab, IntVect ngrow);

/** \brief helper function to read climatology nudging from netcdf */
void
read_clim_nudg_coeff_from_netcdf (int lev, const Box& domain, const std::string& fname,
                                  bool do_m2_clim_nudg,
                                  bool do_m3_clim_nudg,
                                  const amrex::Vector<int>& do_cons_clim_nudg,
                                  const amrex::Vector<std::string>& cons_names,
                                  FArrayBox& NC_M2NC_fab,
                                  FArrayBox& NC_M3NC_fab,
                                  amrex::Vector<FArrayBox>& NC_ConsNC_fab,
                                  amrex::Vector<int>& cons_coeff_in_file);

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

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_data_from_netcdf(lev, boxes_at_level[lev][idx], nc_init_file[lev][idx],
                              NC_temp_fab[idx], NC_salt_fab[idx],
                              NC_xvel_fab[idx], NC_yvel_fab[idx]);
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

        init_state_from_netcdf(lev, temp_fab, salt_fab,
                               xvel_fab, yvel_fab,
                               NC_temp_fab, NC_salt_fab,
                               NC_xvel_fab, NC_yvel_fab);
    } // mf
    } // omp

    // Zero every tracer past salinity. The passive scalars stay zero unless a problem
    // sets them; the biology block is overwritten by init_biology_ic immediately after,
    // so zeroing it here only guards against a partially filled biology IC.
    cons_new[lev]->setVal(zero, Tracer_comp, ncons - Tracer_comp, cons_new[lev]->nGrowVect());
}

void
REMORA::init_biology_from_netcdf (int lev)
{
    if (!REMORABiology::has_biology(biology_model)) {
        return;
    }

    Vector<std::string> biology_names;
    biology_names.reserve(nbio);
    for (int icomp = Bio_comp; icomp < ncons; ++icomp) {
        biology_names.push_back(cons_names[icomp]);
    }

    Vector<Vector<FArrayBox>> NC_biology_fab;
    NC_biology_fab.resize(num_boxes_at_level[lev]);
    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        NC_biology_fab[idx].resize(biology_names.size());
        read_biology_from_netcdf(lev, boxes_at_level[lev][idx], nc_init_file[lev][idx],
                                 biology_names, NC_biology_fab[idx]);
    }

    MultiFab mf_biology(*cons_new[lev], make_alias, Bio_comp, nbio);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
    // Dont tile this since we are operating on full FABs in this routine
    for (MFIter mfi(mf_biology, false); mfi.isValid(); ++mfi)
    {
        FArrayBox& biology_fab = mf_biology[mfi];
        init_biology_state_from_netcdf(lev, biology_fab, NC_biology_fab);
    }
    }
}

void
REMORA::init_data_full_domain_from_netcdf ()
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_temp_fab ; NC_temp_fab.resize(1);
    Vector<FArrayBox> NC_salt_fab ; NC_salt_fab.resize(1);
    Vector<FArrayBox> NC_xvel_fab ; NC_xvel_fab.resize(1);
    Vector<FArrayBox> NC_yvel_fab ; NC_yvel_fab.resize(1);

    read_data_full_domain_from_netcdf(hires_init_level, nc_hires_init_box, nc_init_file_hires,
                          NC_temp_fab[0], NC_salt_fab[0],
                          NC_xvel_fab[0], NC_yvel_fab[0],
                          cum_ref_ratios[hires_init_level]);

    MultiFab mf_temp(*vec_cons_full_domain[hires_init_level], make_alias, Temp_comp, 1);
    MultiFab mf_salt(*vec_cons_full_domain[hires_init_level], make_alias, Salt_comp, 1);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
    // Don't tile this since we are operating on full FABs in this routine
    for ( MFIter mfi(mf_temp, false); mfi.isValid(); ++mfi )
    {
        // Define fabs for holding the initial data
        FArrayBox &temp_fab = mf_temp[mfi];
        FArrayBox &salt_fab = mf_salt[mfi];
        FArrayBox &xvel_fab = (*vec_xvel_full_domain[hires_init_level])[mfi];
        FArrayBox &yvel_fab = (*vec_yvel_full_domain[hires_init_level])[mfi];

        init_state_from_netcdf(hires_init_level, temp_fab, salt_fab,
                               xvel_fab, yvel_fab,
                               NC_temp_fab, NC_salt_fab,
                               NC_xvel_fab, NC_yvel_fab);
    } // mf
    } // omp

    vec_cons_full_domain[hires_init_level]->setVal(zero, Tracer_comp, ncons - Tracer_comp, vec_cons_full_domain[hires_init_level]->nGrowVect());

    // Average down to fill levels below hires_grid_level. Use a special average_down so
    // grow cells get populated by averaged down fine data
    for (int lev=hires_init_level-1; lev >= 0; lev--) {
        average_down_with_grow_cells(lev, vec_cons_full_domain);
        average_down_with_grow_cells(lev, vec_xvel_full_domain);
        average_down_with_grow_cells(lev, vec_yvel_full_domain);
    }
}

/** \brief Biology initialization from full-domain NetCDF file */
void
REMORA::init_biology_full_domain_from_netcdf ()
{
    if (!REMORABiology::has_biology(biology_model)) {
        return;
    }

    Vector<std::string> biology_names;
    biology_names.reserve(nbio);
    for (int icomp = Bio_comp; icomp < ncons; ++icomp) {
        biology_names.push_back(cons_names[icomp]);
    }

    Vector<Vector<FArrayBox>> NC_biology_fab;
    NC_biology_fab.resize(1);
    NC_biology_fab[0].resize(biology_names.size());

    read_biology_full_domain_from_netcdf(hires_init_level, nc_hires_init_box, nc_init_file_hires,
                                         biology_names, NC_biology_fab[0],
                                         cum_ref_ratios[hires_init_level]);

    MultiFab mf_biology(*vec_cons_full_domain[hires_init_level], make_alias, Bio_comp, nbio);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
    // Dont tile this since we are operating on full FABs in this routine
    for (MFIter mfi(mf_biology, false); mfi.isValid(); ++mfi)
    {
        FArrayBox& biology_fab = mf_biology[mfi];
        init_biology_state_from_netcdf(hires_init_level, biology_fab, NC_biology_fab);
    }
    }

    for (int lev=hires_init_level-1; lev >= 0; lev--) {
        average_down_with_grow_cells(lev, vec_cons_full_domain);
    }
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
    (*physbcs[lev])(*vec_zeta[lev],*vec_mskr[lev].get(),0,1,vec_zeta[lev]->nGrowVect(),t_new[lev],zeta_bc(),0,*vec_zeta[lev],*vec_msku[lev],*vec_mskv[lev]);
//    (*physbcs[lev])(*vec_zeta[lev],*vec_mskr[lev].get(),1,1,vec_zeta[lev]->nGrowVect(),t_new[lev],BCVars::zeta_bc);
//    (*physbcs[lev])(*vec_zeta[lev],*vec_mskr[lev].get(),2,1,vec_zeta[lev]->nGrowVect(),t_new[lev],BCVars::zeta_bc);

    if (solverChoice.boundary_from_netcdf) {
        Real told = t_new[lev];
        fill_from_bdyfiles(lev, *vec_zeta[lev], *vec_mskr[lev], told, zeta_bc(),bdy_zeta(),0,0,
                           *vec_zeta[lev]);
    }
    if (lev>0) {
        FillPatch(lev, t_old[lev], *vec_zeta[lev], GetVecOfPtrs(vec_zeta), zeta_bc(), bdy_zeta(),
                  0, false,false,0,0,zero,*vec_zeta[lev]);
    }
//    fill_from_bdyfiles(lev, *vec_zeta[lev], *vec_mskr[lev], told, BCVars::zeta_bc,bdy_zeta(),1,1);
//    fill_from_bdyfiles(lev, *vec_zeta[lev], *vec_mskr[lev], told, BCVars::zeta_bc,bdy_zeta(),2,2);
}

void
REMORA::init_zeta_full_domain_from_netcdf ()
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_zeta_fab     ; NC_zeta_fab.resize(1);

    read_zeta_full_domain_from_netcdf(hires_init_level, nc_hires_init_box, nc_init_file_hires,
                          NC_zeta_fab[0], cum_ref_ratios[hires_init_level]);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
        // Don't tile this since we are operating on full FABs in this routine
        for ( MFIter mfi(*vec_zeta_full_domain[hires_init_level], false); mfi.isValid(); ++mfi )
        {
            FArrayBox &zeta_fab  = (*vec_zeta_full_domain[hires_init_level])[mfi];

            //
            // FArrayBox to FArrayBox copy does "copy on intersection"
            // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
            //

            zeta_fab.template    copy<RunOn::Device>(NC_zeta_fab[0],0,0,1);
        } // mf
        } // omp

    // Average down to fill levels below hires_grid_level. Use a special average_down so
    // grow cells get populated by averaged down fine data
    for (int lev=hires_init_level-1; lev >= 0; lev--) {
        average_down_with_grow_cells(lev, vec_zeta_full_domain);
    }
}


/**
 * @param lev Integer specifying the current level
 */
void
REMORA::init_grid_vars_from_netcdf (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
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

    // Optional spherical psi coordinates (see read_spherical_grid_vars_from_netcdf)
    Vector<FArrayBox> NC_lonp_fab  ; NC_lonp_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_latp_fab  ; NC_latp_fab.resize(num_boxes_at_level[lev]);
    bool have_spherical_psi = true;

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_grid_vars_from_netcdf(lev, boxes_at_level[lev][idx], nc_grid_file[lev][idx],
                                    NC_pm_fab[idx], NC_pn_fab[idx],
                                    NC_xr_fab[idx], NC_yr_fab[idx],
                                    NC_xu_fab[idx], NC_yu_fab[idx],
                                    NC_xv_fab[idx], NC_yv_fab[idx],
                                    NC_xp_fab[idx], NC_yp_fab[idx]);

        // All boxes at this level must agree: a partial set would leave holes in
        // the corner mesh that a conservative remap cannot detect.
        have_spherical_psi = read_spherical_grid_vars_from_netcdf(
                                 lev, boxes_at_level[lev][idx], nc_grid_file[lev][idx],
                                 NC_lonp_fab[idx], NC_latp_fab[idx]) && have_spherical_psi;

        // Mirror vec_xp exactly: same nodal psi BoxArray, DistributionMap and ngrow.
        if (have_spherical_psi && vec_lonp[lev] == nullptr) {
            vec_lonp[lev].reset(new MultiFab(vec_xp[lev]->boxArray(), vec_xp[lev]->DistributionMap(),
                                             1, vec_xp[lev]->nGrowVect()));
            vec_latp[lev].reset(new MultiFab(vec_xp[lev]->boxArray(), vec_xp[lev]->DistributionMap(),
                                             1, vec_xp[lev]->nGrowVect()));
        }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
        // Don't tile this since we are operating on full FABs in this routine
        for ( MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi )
        {
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

            if (have_spherical_psi) {
                (*vec_lonp[lev])[mfi].template copy<RunOn::Device>(NC_lonp_fab[idx]);
                (*vec_latp[lev])[mfi].template copy<RunOn::Device>(NC_latp_fab[idx]);
            }
        } // mf
        } // omp
    } // idx

    // Drop any partial allocation so the getter's null contract stays honest.
    if (!have_spherical_psi) {
        vec_lonp[lev].reset();
        vec_latp[lev].reset();
    }

    Real dummy_time = zero;
    if (lev > 0) {
        FillPatch(lev,dummy_time,*vec_pm[lev],GetVecOfPtrs(vec_pm),
                foextrap_periodic_bc(),
                BdyVars::null,0,false);
        FillPatch(lev,dummy_time,*vec_pn[lev],GetVecOfPtrs(vec_pn),
                foextrap_periodic_bc(),
                BdyVars::null,0,false);
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
    if (vec_lonp[lev] != nullptr) {
        vec_lonp[lev]->FillBoundary(geom[lev].periodicity());
        vec_latp[lev]->FillBoundary(geom[lev].periodicity());
    }

    extrapolate_metric_to_physical_boundaries(*vec_pm[lev], geom[lev]);
    extrapolate_metric_to_physical_boundaries(*vec_pn[lev], geom[lev]);
}

/**
 * @param lev Integer specifying the current level
 */
void
REMORA::init_bathymetry_from_netcdf (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_h_fab     ; NC_h_fab.resize(num_boxes_at_level[lev]);
    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_bathymetry_from_netcdf(lev, boxes_at_level[lev][idx], nc_grid_file[lev][idx],
                                    NC_h_fab[idx]);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
        // Don't tile this since we are operating on full FABs in this routine
        for ( MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi )
        {
            FArrayBox &h_fab     = (*vec_h[lev])[mfi];

            //
            // FArrayBox to FArrayBox copy does "copy on intersection"
            // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
            //

            // Copy into both components of h
            h_fab.template     copy<RunOn::Device>(NC_h_fab[idx],0,0,1);
            h_fab.template     copy<RunOn::Device>(NC_h_fab[idx],0,1,1);
        } // mf
        } // omp
    } // idx

    const double dummy_time = zero;
    // Unconditional foextrap will overwrite periodicity, but EnforcePeriodicity will
    // be called on h afterwards
    FillPatch(lev,dummy_time,*vec_h[lev],GetVecOfPtrs(vec_h),
            foextrap_periodic_bc(),
            BdyVars::null,0,false,false,1);
    FillPatch(lev,dummy_time,*vec_h[lev],GetVecOfPtrs(vec_h),
            foextrap_periodic_bc(),
            BdyVars::null,1,false,false,1);

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

    // One boundary series per BdyVars slot. Every cell-centered tracer gets one, named
    // for the variable itself ("temp", "salt", "NO3", ...), so the file is expected to
    // hold e.g. NO3_west following the same convention ROMS uses. A series whose
    // variable needs no boundary data reads nothing -- NCTimeSeriesBoundary only asks
    // for the sides flagged in phys_bc_need_data -- so tracers left on a local BC do
    // not require the file to contain anything for them.
    amrex::Vector<std::string> field_name (num_bdy_vars());
    amrex::Vector<IntVect    > index_types(num_bdy_vars());
    std::vector<bool         > is_2d      (num_bdy_vars());

    field_name[BdyVars::u] = "u"; index_types[BdyVars::u] = IntVect(1,0,0); is_2d[BdyVars::u] = false;
    field_name[BdyVars::v] = "v"; index_types[BdyVars::v] = IntVect(0,1,0); is_2d[BdyVars::v] = false;
    for (int icomp = 0; icomp < ncons; ++icomp) {
        const int ibdy = BdyVars::cons(icomp);
        field_name[ibdy]  = cons_names[icomp];
        index_types[ibdy] = IntVect(0,0,0);
        is_2d[ibdy]       = false;
    }
    field_name[bdy_ubar()] = "ubar"; index_types[bdy_ubar()] = IntVect(1,0,0); is_2d[bdy_ubar()] = true;
    field_name[bdy_vbar()] = "vbar"; index_types[bdy_vbar()] = IntVect(0,1,0); is_2d[bdy_vbar()] = true;
    field_name[bdy_zeta()] = "zeta"; index_types[bdy_zeta()] = IntVect(0,0,0); is_2d[bdy_zeta()] = true;

    amrex::Print() << "DOING INIT AT LEVEL " << lev << std::endl;
    int rx = 1; int ry = 1;
    if (lev > 0) {
        for (int k = lev-1; k >= 0; k--) {
            rx *= ref_ratio[k][0];
            ry *= ref_ratio[k][1];
        }
    }
    for (int ivar = 0; ivar < num_bdy_vars(); ivar++) {
        boundary_series[lev].push_back(std::unique_ptr<NCTimeSeriesBoundary>(new NCTimeSeriesBoundary(lev, geom, nc_bdry_file, field_name[ivar],
                                                                bdry_time_name_byvar[ivar],
                                                                index_types[ivar],
                                                                phys_bc_need_data[ivar], is_2d[ivar], rx, ry)));
        boundary_series[lev][ivar]->Initialize();
    }
}

/**
 * \brief Helper function to initialize state and velocity data in a Fab.
 *
 * @param lev Integer specifying current level
 * @param state_fab FArrayBox object holding the state data we initialize
 * @param temp_fab  FArrayBox object holding the temperature data we initialize
 * @param salt_fab  FArrayBox object holding the salt        data we initialize
 * @param x_vel_fab FArrayBox object holding the x-velocity data we initialize
 * @param y_vel_fab FArrayBox object holding the y-velocity data we initialize
 * @param NC_temp_fab Vector of FArrayBox objects with the REMORA dataset specifying temperature
 * @param NC_salt_fab Vector of FArrayBox objects with the REMORA dataset specifying salinity
 * @param NC_xvel_fab Vector of FArrayBox objects with the REMORA dataset specifying x-velocity
 * @param NC_yvel_fab Vector of FArrayBox objects with the REMORA dataset specifying y-velocity
 */
void
init_state_from_netcdf (int /*lev*/,
                        FArrayBox&  temp_fab, FArrayBox&  salt_fab,
                        FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                        const Vector<FArrayBox>& NC_temp_fab,
                        const Vector<FArrayBox>& NC_salt_fab,
                        const Vector<FArrayBox>& NC_xvel_fab,
                        const Vector<FArrayBox>& NC_yvel_fab)
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
    } // idx
}

void
init_biology_state_from_netcdf (int /*lev*/, FArrayBox& biology_fab,
                                const Vector<Vector<FArrayBox>>& NC_biology_fab)
{
    int nboxes = NC_biology_fab.size();
    for (int idx = 0; idx < nboxes; idx++)
    {
        AMREX_ALWAYS_ASSERT(static_cast<int>(NC_biology_fab[idx].size()) == biology_fab.nComp());
        for (int ibio = 0; ibio < biology_fab.nComp(); ++ibio) {
            biology_fab.template copy<RunOn::Device>(NC_biology_fab[idx][ibio], 0, ibio, 1);
        }
    }
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
    // One coefficient FAB per cons component per box
    Vector<Vector<FArrayBox>> NC_ConsNC_fab(num_boxes_at_level[lev]);
    // Per box, whether each tracer's coefficient was actually present in the file
    Vector<Vector<int>> cons_coeff_in_file(num_boxes_at_level[lev]);
    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++) {
        NC_ConsNC_fab[idx].resize(ncons);
        cons_coeff_in_file[idx].assign(ncons, 0);
    }
    // Aggregated over boxes: whether each tracer's coefficient came from the file
    Vector<int> cons_coeff_read(ncons, 0);

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_clim_nudg_coeff_from_netcdf(lev,boxes_at_level[lev][idx], nc_clim_coeff_file,
                                         solverChoice.do_m2_clim_nudg,
                                         solverChoice.do_m3_clim_nudg,
                                         solverChoice.do_cons_clim_nudg,
                                         cons_names,
                                         NC_M2NC_fab[idx],NC_M3NC_fab[idx],
                                         NC_ConsNC_fab[idx],cons_coeff_in_file[idx]);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
        // Don't tile this since we are operating on full FABs in this routine
        for ( MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi )
        {
            if (solverChoice.do_m2_clim_nudg) {
                FArrayBox &ubarNC_fab  = (*vec_nudg_coeff[bdy_ubar()][lev])[mfi];
                ubarNC_fab.template    copy<RunOn::Device>(NC_M2NC_fab[idx]);
                FArrayBox &vbarNC_fab  = (*vec_nudg_coeff[bdy_vbar()][lev])[mfi];
                vbarNC_fab.template    copy<RunOn::Device>(NC_M2NC_fab[idx]);
            }
            if (solverChoice.do_m3_clim_nudg) {
                FArrayBox &uNC_fab  = (*vec_nudg_coeff[BdyVars::u][lev])[mfi];
                uNC_fab.template    copy<RunOn::Device>(NC_M3NC_fab[idx]);
                FArrayBox &vNC_fab  = (*vec_nudg_coeff[BdyVars::v][lev])[mfi];
                vNC_fab.template    copy<RunOn::Device>(NC_M3NC_fab[idx]);
            }
            for (int icomp = 0; icomp < ncons; ++icomp) {
                if (!cons_coeff_in_file[idx][icomp]) { continue; }
                FArrayBox &ConsNC_fab  = (*vec_nudg_coeff[BdyVars::cons(icomp)][lev])[mfi];
                ConsNC_fab.template    copy<RunOn::Device>(NC_ConsNC_fab[idx][icomp]);
            }

        } // mf
        } // omp

        for (int icomp = 0; icomp < ncons; ++icomp) {
            cons_coeff_read[icomp] = cons_coeff_read[icomp] || cons_coeff_in_file[idx][icomp];
        }
    } // idx

    if (solverChoice.do_m2_clim_nudg) {
        vec_nudg_coeff[bdy_ubar()][lev]->FillBoundary(geom[lev].periodicity());
        vec_nudg_coeff[bdy_vbar()][lev]->FillBoundary(geom[lev].periodicity());
        convert_inv_days_to_inv_s(vec_nudg_coeff[bdy_ubar()][lev].get());
        convert_inv_days_to_inv_s(vec_nudg_coeff[bdy_vbar()][lev].get());
    }
    if (solverChoice.do_m3_clim_nudg) {
        vec_nudg_coeff[BdyVars::u][lev]->FillBoundary(geom[lev].periodicity());
        vec_nudg_coeff[BdyVars::v][lev]->FillBoundary(geom[lev].periodicity());
        convert_inv_days_to_inv_s(vec_nudg_coeff[BdyVars::u][lev].get());
        convert_inv_days_to_inv_s(vec_nudg_coeff[BdyVars::v][lev].get());
    }
    // Only convert the tracers whose coefficient actually came from the file: the rest
    // still hold the constant set by init_clim_nudg_coeff, which is already in 1/s.
    for (int icomp = 0; icomp < ncons; ++icomp) {
        if (!cons_coeff_read[icomp]) { continue; }
        vec_nudg_coeff[BdyVars::cons(icomp)][lev]->FillBoundary(geom[lev].periodicity());
        convert_inv_days_to_inv_s(vec_nudg_coeff[BdyVars::cons(icomp)][lev].get());
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

    if (river_pos_x.empty() ||
        river_pos_y.size() != river_pos_x.size() ||
        river_direction_tmp.size() != river_pos_x.size())
    {
        amrex::Abort("River metadata arrays must be nonempty and have matching lengths: " +
                     river_x_name + "=" + std::to_string(river_pos_x.size()) + ", " +
                     river_y_name + "=" + std::to_string(river_pos_y.size()) + ", " +
                     river_dir_name + "=" + std::to_string(river_direction_tmp.size()));
    }

    int nriv = river_pos_x.size();
    amrex::Gpu::DeviceVector<int> xpos_d(nriv);
    amrex::Gpu::DeviceVector<int> ypos_d(nriv);
    river_direction.resize(nriv);

    int rrx = (lev > 0) ? cum_ref_ratios[lev][0] : 1;
    int rry = (lev > 0) ? cum_ref_ratios[lev][1] : 1;

    // Map river source indices to the target AMR level while keeping one source
    // point per river so total prescribed transport is unchanged.
    amrex::Vector<int> river_pos_x_lev(nriv);
    amrex::Vector<int> river_pos_y_lev(nriv);
    for (int iriv = 0; iriv < nriv; ++iriv) {
        int x0 = river_pos_x[iriv] - 1;
        int y0 = river_pos_y[iriv] - 1;

        if (river_direction_tmp[iriv] == 0) {
            // u-face aligned in x, centered in y within the refined coarse cell.
            river_pos_x_lev[iriv] = x0 * rrx;
            river_pos_y_lev[iriv] = y0 * rry + (rry - 1) / 2;
        } else {
            // v-face aligned in y, centered in x within the refined coarse cell.
            river_pos_x_lev[iriv] = x0 * rrx + (rrx - 1) / 2;
            river_pos_y_lev[iriv] = y0 * rry;
        }
    }
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy(xpos_d.data(), river_pos_x_lev.data(), sizeof(int)*nriv);
    Gpu::htod_memcpy(ypos_d.data(), river_pos_y_lev.data(), sizeof(int)*nriv);
    Gpu::htod_memcpy(river_direction.data(), river_direction_tmp.data(), sizeof(int)*nriv);
#else
    std::memcpy(xpos_d.data(), river_pos_x_lev.data(), sizeof(int)*nriv);
    std::memcpy(ypos_d.data(), river_pos_y_lev.data(), sizeof(int)*nriv);
    std::memcpy(river_direction.data(), river_direction_tmp.data(), sizeof(int)*nriv);
#endif
    const int* xpos_ptr = xpos_d.data();
    const int* ypos_ptr = ypos_d.data();

    for (amrex::MFIter mfi(*(vec_river_position[lev]).get(),true); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.growntilebox(amrex::IntVect(NGROW,NGROW,0));
        auto river_pos = vec_river_position[lev]->array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int ) {
            for (int iriv=0; iriv < nriv; iriv++) {
                int xriv = xpos_ptr[iriv];
                int yriv = ypos_ptr[iriv];
                if (i==xriv && j==yriv) {
                    river_pos(i,j,0) = iriv;
                }
            }
        });
    }

    if (verbose) {
        amrex::Print() << "[river-debug] lev=" << lev
                       << " ref_ratio=(" << rrx << "," << rry << ")"
                       << " nriv=" << nriv << '\n';
        for (int iriv = 0; iriv < nriv; ++iriv) {
            amrex::Print() << "[river-debug]  river " << iriv
                           << " dir=" << river_direction_tmp[iriv]
                           << " nc=(" << river_pos_x[iriv] << "," << river_pos_y[iriv] << ")"
                           << " lev=(" << river_pos_x_lev[iriv] + 1 << "," << river_pos_y_lev[iriv] + 1 << ")"
                           << '\n';
        }
    }
}

/**
 * @param[inout] mf    multifab of data to convert
 */
void
REMORA::convert_inv_days_to_inv_s (MultiFab* mf) {
    Real inv_days_to_inv_s = one / (Real(3600.0) * Real(24.0));

    for ( MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& arr = mf->array(mfi);
        Box bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            arr(i,j,k) *= inv_days_to_inv_s;
        });
    }

}

void
REMORA::init_bathymetry_full_domain_from_netcdf ()
{
    if (nc_grid_file_hires.empty()) {
        Abort("Must specify high-resolution grid file when initializing from NetCDF and hires_grid_level > 0");
    }
    Vector<FArrayBox> NC_h_fab     ; NC_h_fab.resize(1);
    read_bathymetry_full_domain_from_netcdf(nc_hires_grid_box, nc_grid_file_hires, NC_h_fab[0],cum_ref_ratios[hires_grid_level]);

    // Don't tile this since we are operating on full FABs in this routine
    for ( MFIter mfi(*vec_h_full_domain[hires_grid_level], false); mfi.isValid(); ++mfi )
    {
        FArrayBox &h_fab     = (*vec_h_full_domain[hires_grid_level])[mfi];
        h_fab.template    copy<RunOn::Device>(NC_h_fab[0]);
    }

    // Average down to fill levels below hires_grid_level. Use a special average_down so
    // grow cells get populated by averaged down fine data
    for (int lev=hires_grid_level-1; lev >= 0; lev--) {
        average_down_with_grow_cells(lev, vec_h_full_domain);;
    }
}

void
REMORA::init_grid_vars_full_domain_from_netcdf ()
{
    if (nc_grid_file_hires.empty()) {
        Abort("Must specify high-resolution grid file when initializing from NetCDF and hires_grid_level > 0");
    }
    Vector<FArrayBox> NC_pm_fab    ; NC_pm_fab.resize(1);
    Vector<FArrayBox> NC_pn_fab    ; NC_pn_fab.resize(1);

    read_grid_vars_full_domain_from_netcdf(nc_hires_grid_box, nc_grid_file_hires,
                                    NC_pm_fab[0], NC_pn_fab[0],
                                    cum_ref_ratios[hires_grid_level]);

    // Don't tile this since we are operating on full FABs in this routine
    for ( MFIter mfi(*vec_h_full_domain[hires_grid_level], false); mfi.isValid(); ++mfi )
    {
        FArrayBox &pm_fab    = (*vec_pm_full_domain[hires_grid_level])[mfi];
        FArrayBox &pn_fab    = (*vec_pn_full_domain[hires_grid_level])[mfi];
        pm_fab.template    copy<RunOn::Device>(NC_pm_fab[0]);
        pn_fab.template    copy<RunOn::Device>(NC_pn_fab[0]);
    }

    // Average down to fill levels below hires_grid_level. Use a special average_down so
    // grow cells get populated by averaged down fine data
    for (int lev=hires_grid_level-1; lev >= 0; lev--) {
        average_down_with_grow_cells(lev, vec_pm_full_domain);
        average_down_with_grow_cells(lev, vec_pn_full_domain);

        int rrx = ref_ratio[lev][0];
        int rry = ref_ratio[lev][1];
        // pm and pn need to be rescaled by the refinement ratio
        for ( MFIter mfi(*vec_h_full_domain[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<Real> const& pm   = vec_pm_full_domain[lev]->array(mfi);
            Array4<Real> const& pn   = vec_pn_full_domain[lev]->array(mfi);
            Box ubx = mfi.growntilebox(cum_ref_ratios[lev] - IntVect(1,0,0));;
            Box vbx = mfi.growntilebox(cum_ref_ratios[lev] - IntVect(0,1,0));;
            ParallelFor(makeSlab(ubx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int ) {
                pm(i,j,0) = pm(i,j,0) / Real(rrx);
            });
            ParallelFor(makeSlab(vbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int ) {
                pn(i,j,0) = pn(i,j,0) / Real(rry);
            });
        }
    }

    for (int lev=0; lev<=hires_grid_level; lev++) {
        extrapolate_metric_to_physical_boundaries(*vec_pm_full_domain[lev], geom[lev]);
        extrapolate_metric_to_physical_boundaries(*vec_pn_full_domain[lev], geom[lev]);
    }
}

/**
 * @param[inout] mf    multifab of data to extrapolate on
 * @param[in   ] geom  geometry
 */
void
REMORA::extrapolate_metric_to_physical_boundaries (MultiFab& mf, const Geometry& geom)
{
    const IntVect ng = mf.nGrowVect();

    const auto& dom_lo = amrex::lbound(geom.Domain());
    const auto& dom_hi = amrex::ubound(geom.Domain());

    for ( MFIter mfi(mf); mfi.isValid(); ++mfi )
    {
        Box bx = mfi.tilebox();

        auto mf_arr = mf.array(mfi);

        Box gbx_lox = adjCellLo(bx,0,ng[0]); gbx_lox.grow(1,ng[1]); gbx_lox.setBig  (0,dom_lo.x-2);
        Box gbx_hix = adjCellHi(bx,0,ng[0]); gbx_hix.grow(1,ng[1]); gbx_hix.setSmall(0,dom_hi.x+2);
        Box gbx_loy = adjCellLo(bx,1,ng[1]); gbx_loy.grow(0,ng[0]); gbx_loy.setBig  (1,dom_lo.y-2);
        Box gbx_hiy = adjCellHi(bx,1,ng[1]); gbx_hiy.grow(0,ng[0]); gbx_hiy.setSmall(1,dom_hi.y+2);

        if (gbx_lox.ok()) {
            ParallelFor(gbx_lox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                mf_arr(i,j,k,0) = mf_arr(dom_lo.x-1,j,k,0);
            });
        }
        if (gbx_hix.ok()) {
            ParallelFor(gbx_hix, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                mf_arr(i,j,k,0) = mf_arr(dom_hi.x+1,j,k,0);
            });
        }
        if (gbx_loy.ok()) {
            ParallelFor(gbx_loy, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                mf_arr(i,j,k,0) = mf_arr(i,dom_lo.y-1,k,0);
            });
        }
        if (gbx_hiy.ok()) {
            ParallelFor(gbx_hiy, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                mf_arr(i,j,k,0) = mf_arr(i,dom_hi.y+1,k,0);
            });
        }
    } // mfi
}

#endif // REMORA_USE_NETCDF
