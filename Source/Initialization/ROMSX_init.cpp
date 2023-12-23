/**
 * \file ROMSX_init.cpp
 */

#include <ROMSX.H>
#include <EOS.H>
#include <ROMSX_Constants.H>
#include <prob_common.H>

using namespace amrex;

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
        Array4<const Real> const& Zt_avg1_arr  = mf_Zt_avg1->const_array(mfi);

        init_custom_prob(bx, cons_arr, xvel_arr, yvel_arr, zvel_arr,
                         z_w_arr, z_r_arr, Hz_arr, h_arr, Zt_avg1_arr, geom[lev].data(),
                         solverChoice);

    } //mfi

    // Initialize the "pm" and "pn" arrays
    const auto dxi = Geom(lev).InvCellSize();
    vec_pm[lev]->setVal(dxi[0]);
    vec_pn[lev]->setVal(dxi[1]);

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

    // DEBUGGING NOTE -- DoublyPeriodic fails if these are commented out
    const Real time = 0.0;
    FillPatch(lev,time, *vec_ubar[lev], GetVecOfPtrs(vec_ubar));
    FillPatch(lev,time, *vec_vbar[lev], GetVecOfPtrs(vec_vbar));
}
