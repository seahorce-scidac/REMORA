/**
 * \file ROMSX_init.cpp
 */

#include <ROMSX.H>
#include <EOS.H>
#include <ROMSX_Constants.H>
#include <prob_common.H>

using namespace amrex;

#ifdef ROMSX_USE_NETCDF

Real
read_from_wrfbdy(std::string nc_bdy_file, const Box& domain,
                 Vector<Vector<FArrayBox>>& bdy_data_xlo,
                 Vector<Vector<FArrayBox>>& bdy_data_xhi,
                 Vector<Vector<FArrayBox>>& bdy_data_ylo,
                 Vector<Vector<FArrayBox>>& bdy_data_yhi);

void
convert_wrfbdy_data (int which, const Box& domain,
                     Vector<Vector<FArrayBox>>& bdy_data,
                     const FArrayBox& NC_MUB_fab,
                     const FArrayBox& NC_MSFU_fab,
                     const FArrayBox& NC_MSFV_fab,
                     const FArrayBox& NC_MSFM_fab,
                     const FArrayBox& NC_PH_fab,
                     const FArrayBox& NC_PHB_fab,
                     const FArrayBox& NC_C1H_fab,
                     const FArrayBox& NC_C2H_fab,
                     const FArrayBox& NC_RDNW_fab,
                     const FArrayBox& NC_xvel_fab,
                     const FArrayBox& NC_yvel_fab,
                     const FArrayBox& NC_rho_fab,
                     const FArrayBox& NC_rhoth_fab);

void
ROMSX::init_from_wrfinput (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_xvel_fab ; NC_xvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_yvel_fab ; NC_yvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_zvel_fab ; NC_zvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_rho_fab  ; NC_rho_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_rhop_fab ; NC_rhop_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_rhoth_fab; NC_rhoth_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MUB_fab  ; NC_MUB_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFU_fab ; NC_MSFU_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFV_fab ; NC_MSFV_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFM_fab ; NC_MSFM_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_SST_fab  ; NC_SST_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_C1H_fab  ; NC_C1H_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_C2H_fab  ; NC_C2H_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_RDNW_fab ; NC_RDNW_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_PH_fab   ; NC_PH_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_PHB_fab  ; NC_PHB_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_ALB_fab  ; NC_ALB_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_PB_fab   ; NC_PB_fab.resize(num_boxes_at_level[lev]);

    if (nc_init_file.size() == 0)
        amrex::Error("NetCDF initialization file name must be provided via input");

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_from_wrfinput(lev,idx,NC_xvel_fab,NC_yvel_fab,NC_zvel_fab,NC_rho_fab,
                           NC_rhop_fab,NC_rhoth_fab,NC_MUB_fab,
                           NC_MSFU_fab,NC_MSFV_fab,NC_MSFM_fab,
                           NC_SST_fab,
                           NC_C1H_fab,NC_C2H_fab,NC_RDNW_fab,
                           NC_PH_fab,NC_PHB_fab,NC_ALB_fab,NC_PB_fab);
    }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    // INITIAL DATA common for "ideal" as well as "real" simulation
    for ( MFIter mfi(cons_new, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Define fabs for holding the initial data
        FArrayBox &cons_fab = cons_new[lev][mfi];
        FArrayBox &xvel_fab = xvel_new[lev][mfi];
        FArrayBox &yvel_fab = yvel_new[lev][mfi];
        FArrayBox &zvel_fab = zvel_new[lev][mfi];

        init_state_from_wrfinput(lev, cons_fab, xvel_fab, yvel_fab, zvel_fab,
                                 NC_xvel_fab, NC_yvel_fab, NC_zvel_fab,
                                 NC_rho_fab, NC_rhoth_fab);
    } // mf

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    // Map scale factors common for "ideal" as well as "real" simulation
    for ( MFIter mfi(*mapfac_u[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Define fabs for holding the initial data
        FArrayBox &msfu_fab = (*mapfac_u[lev])[mfi];
        FArrayBox &msfv_fab = (*mapfac_v[lev])[mfi];
        FArrayBox &msfm_fab = (*mapfac_m[lev])[mfi];

        init_msfs_from_wrfinput(lev, msfu_fab, msfv_fab, msfm_fab,
                                NC_MSFU_fab, NC_MSFV_fab, NC_MSFM_fab);
    } // mf

    if (init_type == "real" && (lev == 0)) {
        if (nc_bdy_file.empty())
            amrex::Error("NetCDF boundary file name must be provided via input");
        bdy_time_interval = read_from_wrfbdy(nc_bdy_file,geom[0].Domain(),bdy_data_xlo,bdy_data_xhi,bdy_data_ylo,bdy_data_yhi);

        const Box& domain = geom[lev].Domain();

        convert_wrfbdy_data(0,domain,bdy_data_xlo,
                            NC_MUB_fab[0], NC_MSFU_fab[0], NC_MSFV_fab[0], NC_MSFM_fab[0],
                            NC_PH_fab[0] , NC_PHB_fab[0],
                            NC_C1H_fab[0], NC_C2H_fab[0], NC_RDNW_fab[0],
                            NC_xvel_fab[0],NC_yvel_fab[0],NC_rho_fab[0],NC_rhoth_fab[0]);
        convert_wrfbdy_data(1,domain,bdy_data_xhi,
                            NC_MUB_fab[0], NC_MSFU_fab[0], NC_MSFV_fab[0], NC_MSFM_fab[0],
                            NC_PH_fab[0] , NC_PHB_fab[0],
                            NC_C1H_fab[0], NC_C2H_fab[0], NC_RDNW_fab[0],
                            NC_xvel_fab[0],NC_yvel_fab[0],NC_rho_fab[0],NC_rhoth_fab[0]);
        convert_wrfbdy_data(2,domain,bdy_data_ylo,
                            NC_MUB_fab[0], NC_MSFU_fab[0], NC_MSFV_fab[0], NC_MSFM_fab[0],
                            NC_PH_fab[0] , NC_PHB_fab[0],
                            NC_C1H_fab[0], NC_C2H_fab[0], NC_RDNW_fab[0],
                            NC_xvel_fab[0],NC_yvel_fab[0],NC_rho_fab[0],NC_rhoth_fab[0]);
        convert_wrfbdy_data(3,domain,bdy_data_yhi,
                            NC_MUB_fab[0], NC_MSFU_fab[0], NC_MSFV_fab[0], NC_MSFM_fab[0],
                            NC_PH_fab[0] , NC_PHB_fab[0],
                            NC_C1H_fab[0], NC_C2H_fab[0], NC_RDNW_fab[0],
                            NC_xvel_fab[0],NC_yvel_fab[0],NC_rho_fab[0],NC_rhoth_fab[0]);
    }
}

void
ROMSX::init_state_from_wrfinput (int lev, FArrayBox& state_fab,
                                 FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                                 FArrayBox& z_vel_fab,
                                 const Vector<FArrayBox>& NC_xvel_fab,
                                 const Vector<FArrayBox>& NC_yvel_fab,
                                 const Vector<FArrayBox>& NC_zvel_fab,
                                 const Vector<FArrayBox>& NC_rho_fab,
                                 const Vector<FArrayBox>& NC_temp_fab)
{
    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        //
        // FArrayBox to FArrayBox copy does "copy on intersection"
        // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
        //
        // This copies x-vel
        x_vel_fab.template copy<RunOn::Device>(NC_xvel_fab[idx]);

        // This copies y-vel
        y_vel_fab.template copy<RunOn::Device>(NC_yvel_fab[idx]);

        // This copies z-vel
        z_vel_fab.template copy<RunOn::Device>(NC_zvel_fab[idx]);

        // We first initialize all state_fab variables to zero
        state_fab.template setVal<RunOn::Device>(0.);

        // This copies temperature
        state_fab.template copy<RunOn::Device>(NC_temp_fab[idx], 0, Temp_comp, 1);

#ifdef ROMSX_USE_SALINITY
        // This copies salt
        state_fab.template copy<RunOn::Device>(NC_temp_fab[idx], 0, Salt_comp, 1);
#endif
    } // idx
}

void
ROMSX::init_msfs_from_wrfinput (int lev, FArrayBox& msfu_fab,
                                FArrayBox& msfv_fab, FArrayBox& msfm_fab,
                                const Vector<FArrayBox>& NC_MSFU_fab,
                                const Vector<FArrayBox>& NC_MSFV_fab,
                                const Vector<FArrayBox>& NC_MSFM_fab)
{
    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        //
        // FArrayBox to FArrayBox copy does "copy on intersection"
        // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
        //
        // This copies mapfac_u
        msfu_fab.template copy<RunOn::Device>(NC_MSFU_fab[idx]);

        // This copies mapfac_v
        msfv_fab.template copy<RunOn::Device>(NC_MSFV_fab[idx]);

        // This copies mapfac_m
        msfm_fab.template copy<RunOn::Device>(NC_MSFM_fab[idx]);
    } // idx
}
#endif // ROMSX_USE_NETCDF

void
ROMSX::init_custom(int lev)
{
    std::unique_ptr<MultiFab>& mf_z_w = vec_z_w[lev];
    std::unique_ptr<MultiFab>& mf_z_r = vec_z_r[lev];
    std::unique_ptr<MultiFab>& mf_Hz  = vec_Hz[lev];
    std::unique_ptr<MultiFab>& mf_h  = vec_hOfTheConfusingName[lev];
    std::unique_ptr<MultiFab>& mf_Zt_avg1  = vec_Zt_avg1[lev];

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
        const auto &cons_arr = cons_new[lev]->array(mfi);
        const auto &xvel_arr = xvel_new[lev]->array(mfi);
        const auto &yvel_arr = yvel_new[lev]->array(mfi);
        const auto &zvel_arr = zvel_new[lev]->array(mfi);

        Array4<const Real> const& z_w_arr = (mf_z_w)->array(mfi);
        Array4<const Real> const& z_r_arr = (mf_z_r)->array(mfi);
        Array4<const Real> const& Hz_arr  = (mf_Hz)->array(mfi);
        Array4<const Real> const& h_arr  = (mf_h)->array(mfi);
        Array4<const Real> const& Zt_avg1_arr  = (mf_Zt_avg1)->array(mfi);

        init_custom_prob(bx, cons_arr, xvel_arr, yvel_arr, zvel_arr,
                         z_w_arr, z_r_arr, Hz_arr, h_arr, Zt_avg1_arr, geom[lev].data(),
                         solverChoice);

    } //mfi

    set_2darrays(lev);

}

void
ROMSX::set_2darrays (int lev)
{
    std::unique_ptr<MultiFab>& mf_x_r = vec_x_r[lev];
    std::unique_ptr<MultiFab>& mf_y_r = vec_y_r[lev];
    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    for ( MFIter mfi(*(mf_x_r), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

      Array4<Real> const& x_r = (mf_x_r)->array(mfi);
      Array4<Real> const& y_r = (mf_y_r)->array(mfi);
      const Box& bx = mfi.growntilebox();
      const auto & geomdata = Geom(lev).data();
      Gpu::synchronize();
      amrex::ParallelFor(amrex::makeSlab(bx,2,0),
      [=] AMREX_GPU_DEVICE (int i, int j, int  )
      {
        const auto prob_lo         = geomdata.ProbLo();
        const auto dx              = geomdata.CellSize();

        x_r(i,j,0) = prob_lo[0] + (i + 0.5) * dx[0];
        y_r(i,j,0) = prob_lo[1] + (j + 0.5) * dx[1];
        //        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

      });
    }

    MultiFab* U_old = xvel_new[lev];
    MultiFab* V_old = yvel_new[lev];
    std::unique_ptr<MultiFab>& mf_ubar = vec_ubar[lev];
    std::unique_ptr<MultiFab>& mf_vbar = vec_vbar[lev];
    std::unique_ptr<MultiFab>& mf_Hz  = vec_Hz[lev];
    int nstp = 0;
    int kstp = 0;
    int knew = 0;

    for ( MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& ubar = (mf_ubar)->array(mfi);
        Array4<Real> const& vbar = (mf_vbar)->array(mfi);

        Array4<const Real> const& Hz       = mf_Hz->const_array(mfi);
        Array4<const Real> const& u        = U_old->const_array(mfi);
        Array4<const Real> const& v        = V_old->const_array(mfi);

        Box ubx2 = mfi.nodaltilebox(0); ubx2.grow(IntVect(NGROW  ,NGROW  ,0)); // x-face-centered, grown by 2
        Box vbx2 = mfi.nodaltilebox(1); vbx2.grow(IntVect(NGROW  ,NGROW  ,0)); // y-face-centered, grown by 2

        amrex::ParallelFor(makeSlab(ubx2,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real CF = 0.;
                Real sum_of_hz = 0.;

                for (int k=0; k<=N; k++) {
                    Real avg_hz = 0.5*(Hz(i,j,k)+Hz(i-1,j,k));
                    sum_of_hz += avg_hz;
                    CF += avg_hz*u(i,j,k,nstp);
                }
                ubar(i,j,0,kstp) = CF / sum_of_hz;
                ubar(i,j,0,knew) = CF / sum_of_hz;
            });

        amrex::ParallelFor(makeSlab(vbx2,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real CF = 0.;
                Real sum_of_hz = 0.;

                for(int k=0; k<=N; k++) {
                    Real avg_hz = 0.5*(Hz(i,j,k)+Hz(i,j-1,k));
                    sum_of_hz += avg_hz;
                    CF += avg_hz*v(i,j,k,nstp);
                }
                vbar(i,j,0,kstp) = CF / sum_of_hz;
                vbar(i,j,0,knew) = CF / sum_of_hz;
            });
    }

    vec_ubar[lev]->FillBoundary(geom[lev].periodicity());
    vec_vbar[lev]->FillBoundary(geom[lev].periodicity());
}
