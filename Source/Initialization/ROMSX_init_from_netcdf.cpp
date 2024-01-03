/**
 * \file ROMSX_init_from_netcdf.cpp
 */

#include <ROMSX.H>
#include <EOS.H>
#include <ROMSX_Constants.H>
#include <prob_common.H>
#include <DataStruct.H>

using namespace amrex;

#ifdef ROMSX_USE_NETCDF

void
read_data_from_netcdf (int /*lev*/, const Box& domain, const std::string& fname,
                       FArrayBox& NC_temp_fab, FArrayBox& NC_salt_fab,
                       FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                       FArrayBox& NC_ubar_fab, FArrayBox& NC_vbar_fab,
                       FArrayBox& NC_zeta_fab);

Real
read_bdry_from_netcdf (const Box& domain, const std::string& fname,
                       Vector<Vector<FArrayBox>>& bdy_data_xlo,
                       Vector<Vector<FArrayBox>>& bdy_data_xhi,
                       Vector<Vector<FArrayBox>>& bdy_data_ylo,
                       Vector<Vector<FArrayBox>>& bdy_data_yhi,
                       int& width, amrex::Real& start_bdy_time);

void
init_state_from_netcdf (int lev,
                        FArrayBox&  temp_fab, FArrayBox&  salt_fab,
                        FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                        FArrayBox&  ubar_fab, FArrayBox&  vbar_fab,
                        FArrayBox& zeta_fab,
                        const Vector<FArrayBox>& NC_temp_fab,
                        const Vector<FArrayBox>& NC_salt_fab,
                        const Vector<FArrayBox>& NC_xvel_fab,
                        const Vector<FArrayBox>& NC_yvel_fab,
                        const Vector<FArrayBox>& NC_ubar_fab,
                        const Vector<FArrayBox>& NC_vbar_fab,
                        const Vector<FArrayBox>& NC_zeta_fab);

void
read_bathymetry_from_netcdf (int lev, const Box& domain, const std::string& fname,
                             FArrayBox& NC_h_fab,
                             FArrayBox& NC_pm_fab, FArrayBox& NC_pn_fab);

void
init_bathymetry_from_netcdf (int lev);

/**
 * ROMSX function that initializes solution data from a netcdf file
 *
 * @param lev Integer specifying the current level
 */
void
ROMSX::init_data_from_netcdf (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_temp_fab ; NC_temp_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_salt_fab ; NC_salt_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_xvel_fab ; NC_xvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_yvel_fab ; NC_yvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_ubar_fab ; NC_ubar_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_vbar_fab ; NC_vbar_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_zeta_fab ; NC_zeta_fab.resize(num_boxes_at_level[lev]);

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_data_from_netcdf(lev, boxes_at_level[lev][idx], nc_init_file[lev][idx],
                              NC_temp_fab[idx], NC_salt_fab[idx],
                              NC_xvel_fab[idx], NC_yvel_fab[idx],
                              NC_ubar_fab[idx], NC_vbar_fab[idx],
                              NC_zeta_fab[idx]);
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
        FArrayBox &zeta_fab = (*vec_zeta[lev])[mfi];

        init_state_from_netcdf(lev, temp_fab, salt_fab,
                               xvel_fab, yvel_fab,
                               ubar_fab, vbar_fab, zeta_fab,
                               NC_temp_fab, NC_salt_fab,
                               NC_xvel_fab, NC_yvel_fab,
                               NC_ubar_fab, NC_vbar_fab,
                               NC_zeta_fab);
    } // mf
    } // omp
}

/**
 * ROMSX function that initializes bathymetry from a netcdf file
 *
 * @param lev Integer specifying the current level
 */
void
ROMSX::init_bathymetry_from_netcdf (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_h_fab     ; NC_h_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_pm_fab    ; NC_pm_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_pn_fab    ; NC_pn_fab.resize(num_boxes_at_level[lev]);

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_bathymetry_from_netcdf(lev, boxes_at_level[lev][idx], nc_grid_file[lev][idx],
                                    NC_h_fab[idx],
                                    NC_pm_fab[idx], NC_pn_fab[idx]);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
        // Don't tile this since we are operating on full FABs in this routine
        for ( MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi )
        {
            FArrayBox &h_fab     = (*vec_hOfTheConfusingName[lev])[mfi];
            FArrayBox &pm_fab    = (*vec_pm[lev])[mfi];
            FArrayBox &pn_fab    = (*vec_pn[lev])[mfi];

            //
            // FArrayBox to FArrayBox copy does "copy on intersection"
            // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
            //

            // Copy into both components of h
            h_fab.template     copy<RunOn::Device>(NC_h_fab[idx],0,0,1);
            h_fab.template     copy<RunOn::Device>(NC_h_fab[idx],0,1,1);

            pm_fab.template    copy<RunOn::Device>(NC_pm_fab[idx]);
            pn_fab.template    copy<RunOn::Device>(NC_pn_fab[idx]);
        } // mf
        } // omp
    } // idx
    vec_hOfTheConfusingName[lev]->FillBoundary(geom[lev].periodicity());
    vec_pm[lev]->FillBoundary(geom[lev].periodicity());
    vec_pn[lev]->FillBoundary(geom[lev].periodicity());
}

/**
 * ROMSX function that initializes time series of boundary data from a netcdf file
 *
 * @param lev Integer specifying the current level
 */
void
ROMSX::init_bdry_from_netcdf ()
{
    if (nc_bdry_file.empty()) {
        amrex::Error("NetCDF boundary file name must be provided via input");
    }

    bdy_time_interval = read_bdry_from_netcdf(geom[0].Domain(), nc_bdry_file,
                                              bdy_data_xlo,bdy_data_xhi,bdy_data_ylo,bdy_data_yhi,
                                              bdy_width, start_bdy_time);

    if (bdy_width-1 <= bdy_set_width) bdy_set_width = bdy_width;
    amrex::Print() << "Read in boundary data with width "  << bdy_width << std::endl;
    amrex::Print() << "Running with specification width: " << bdy_set_width
                   << " and relaxation width: " << bdy_width - bdy_set_width << std::endl;

    // NOTE: Last bdy cell is a ghost cell for Laplacian relaxation.
    //       Without relaxation zones, we must augment this value by 1.
    if (bdy_width == bdy_set_width) bdy_width += 1;
}

/**
 * Helper function to initialize state and velocity data in a Fab from a ROMS-x dataset.
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
 * @param NC_temp_fab Vector of FArrayBox objects with the ROMS-x dataset specifying temperature
 * @param NC_salt_fab Vector of FArrayBox objects with the ROMS-x dataset specifying salinity
 * @param NC_xvel_fab Vector of FArrayBox objects with the ROMS-x dataset specifying x-velocity
 * @param NC_yvel_fab Vector of FArrayBox objects with the ROMS-x dataset specifying y-velocity
 * @param NC_ubar_fab Vector of FArrayBox objects with the ROMS-x dataset specifying ubar
 * @param NC_vbar_fab Vector of FArrayBox objects with the ROMS-x dataset specifying vbar
 * @param NC_zeta_fab Vector of FArrayBox objects with the ROMS-x dataset specifying zeta
 */
void
init_state_from_netcdf (int /*lev*/,
                        FArrayBox&  temp_fab, FArrayBox&  salt_fab,
                        FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                        FArrayBox&  ubar_fab, FArrayBox&  vbar_fab,
                        FArrayBox& zeta_fab,
                        const Vector<FArrayBox>& NC_temp_fab,
                        const Vector<FArrayBox>& NC_salt_fab,
                        const Vector<FArrayBox>& NC_xvel_fab,
                        const Vector<FArrayBox>& NC_yvel_fab,
                        const Vector<FArrayBox>& NC_ubar_fab,
                        const Vector<FArrayBox>& NC_vbar_fab,
                        const Vector<FArrayBox>& NC_zeta_fab)
{
    int nboxes = NC_xvel_fab.size();
    for (int idx = 0; idx < nboxes; idx++)
    {
        //
        // FArrayBox to FArrayBox copy does "copy on intersection"
        // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
        //

        // This copies temperature
        temp_fab.template copy<RunOn::Device>(NC_temp_fab[idx]);

        // This copies salinity
        salt_fab.template copy<RunOn::Device>(NC_salt_fab[idx]);

        // This copies x-vel
        x_vel_fab.template copy<RunOn::Device>(NC_xvel_fab[idx]);

        // This copies y-vel
        y_vel_fab.template copy<RunOn::Device>(NC_yvel_fab[idx]);

        // This copies ubar -- from a single value to fill all three values
        ubar_fab.template copy<RunOn::Device>(NC_ubar_fab[idx],0,0,1);
        // ubar_fab.template copy<RunOn::Device>(NC_ubar_fab[idx],0,1,1);
        // ubar_fab.template copy<RunOn::Device>(NC_ubar_fab[idx],0,2,1);

        // This copies vbar
        vbar_fab.template copy<RunOn::Device>(NC_vbar_fab[idx],0,0,1);
        // vbar_fab.template copy<RunOn::Device>(NC_vbar_fab[idx],0,1,1);
        // vbar_fab.template copy<RunOn::Device>(NC_vbar_fab[idx],0,2,1);

        // This copies zeta
        zeta_fab.template copy<RunOn::Device>(NC_zeta_fab[idx],0,0,1);
        // zeta_fab.template copy<RunOn::Device>(NC_zeta_fab[idx],0,1,1);
        // zeta_fab.template copy<RunOn::Device>(NC_zeta_fab[idx],0,2,1);
    } // idx
}

#endif // ROMSX_USE_NETCDF
