/**
 * \file REMORA_init.cpp
 */

#include <REMORA.H>
#include <EOS.H>
#include <REMORA_Constants.H>
#include <prob_common.H>

using namespace amrex;

void
REMORA::init_custom(int lev)
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
    set_pm_pn(lev);

}

void
REMORA::init_beta_plane_coriolis (int lev)
{
    std::unique_ptr<MultiFab>& mf_fcor = vec_fcor[lev];
    auto geomdata  = Geom(lev).data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto fcor_arr = (mf_fcor)->array(mfi);
        Real coriolis_f0 = solverChoice.coriolis_f0;
        Real coriolis_beta = solverChoice.coriolis_beta;
        Real Esize = geomdata.ProbHi()[1] - geomdata.ProbLo()[1];
        Real prob_lo = geomdata.ProbLo()[1];
        Real dx = geomdata.CellSize()[1];

        ParallelFor(Box(fcor_arr), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real y = prob_lo + (j + 0.5_rt) * dx;
            fcor_arr(i,j,0) = coriolis_f0 + coriolis_beta * (y - 0.5_rt * Esize);
        });
    } //mfi

    vec_fcor[lev]->FillBoundary(geom[lev].periodicity());
}

void
REMORA::set_zeta_average (int lev)
{
    std::unique_ptr<MultiFab>& mf_zeta = vec_zeta[lev];
    std::unique_ptr<MultiFab>& mf_Zt_avg1  = vec_Zt_avg1[lev];
    for ( MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        int nstp = 0;
        Array4<Real> const& Zt_avg1 = (mf_Zt_avg1)->array(mfi);
        Array4<const Real> const& zeta     = mf_zeta->const_array(mfi);

        Box  bx3 = mfi.tilebox()      ;  bx3.grow(IntVect(NGROW+1,NGROW+1,0)); //   cell-centered, grown by 3

        ParallelFor(makeSlab(bx3,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Zt_avg1(i,j,0) = zeta(i,j,0,nstp);
        });

    }
}

void
REMORA::set_2darrays (int lev)
{
    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    vec_ubar[lev]->setVal(0.0_rt);
    vec_vbar[lev]->setVal(0.0_rt);

    MultiFab* U_old = xvel_new[lev];
    MultiFab* V_old = yvel_new[lev];
    std::unique_ptr<MultiFab>& mf_ubar = vec_ubar[lev];
    std::unique_ptr<MultiFab>& mf_vbar = vec_vbar[lev];
    std::unique_ptr<MultiFab>& mf_Hz  = vec_Hz[lev];
    int nstp = 0;

    for ( MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& ubar = (mf_ubar)->array(mfi);
        Array4<Real> const& vbar = (mf_vbar)->array(mfi);

        Array4<const Real> const& Hz       = mf_Hz->const_array(mfi);
        Array4<const Real> const& u        = U_old->const_array(mfi);
        Array4<const Real> const& v        = V_old->const_array(mfi);

        Box  bx2 = mfi.tilebox()      ;  bx2.grow(IntVect(NGROW  ,NGROW  ,0)); //   cell-centered, grown by 2
        Box ubx2 = mfi.nodaltilebox(0); ubx2.grow(IntVect(NGROW  ,NGROW  ,0)); // x-face-centered, grown by 2
        Box vbx2 = mfi.nodaltilebox(1); vbx2.grow(IntVect(NGROW  ,NGROW  ,0)); // y-face-centered, grown by 2

        ParallelFor(makeSlab(ubx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real CF = 0.;
            Real sum_of_hz = 0.;

            for (int k=0; k<=N; k++) {
                Real avg_hz = 0.5_rt*(Hz(i,j,k)+Hz(i-1,j,k));
                sum_of_hz += avg_hz;
                CF += avg_hz*u(i,j,k,nstp);
            }
            ubar(i,j,0,0) = CF / sum_of_hz;
        });

        ParallelFor(makeSlab(vbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real CF = 0.;
            Real sum_of_hz = 0.;

            for(int k=0; k<=N; k++) {
                Real avg_hz = 0.5_rt*(Hz(i,j,k)+Hz(i,j-1,k));
                sum_of_hz += avg_hz;
                CF += avg_hz*v(i,j,k,nstp);
            }
            vbar(i,j,0,0) = CF / sum_of_hz;
        });
    }

    FillPatch(lev, t_new[lev], *vec_ubar[lev], GetVecOfPtrs(vec_ubar), BCVars::ubar_bc, BdyVars::ubar,0,false,false);
    FillPatch(lev, t_new[lev], *vec_vbar[lev], GetVecOfPtrs(vec_vbar), BCVars::vbar_bc, BdyVars::vbar,0,false,false);
}

void
REMORA::init_gls_vmix (int lev, SolverChoice solver_choice)
{
    vec_tke[lev]->setVal(solver_choice.gls_Kmin);
    vec_gls[lev]->setVal(solver_choice.gls_Pmin);
    vec_Lscale[lev]->setVal(0.0_rt);
    vec_Akk[lev]->setVal(solver_choice.Akk_bak);
    vec_Akp[lev]->setVal(solver_choice.Akp_bak);
    vec_Akv[lev]->setVal(solver_choice.Akv_bak);
    vec_Akt[lev]->setVal(solver_choice.Akt_bak);

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    for (MFIter mfi(*vec_Akk[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real> const& Akk = vec_Akk[lev]->array(mfi);
        Array4<Real> const& Akp = vec_Akp[lev]->array(mfi);
        Array4<Real> const& Akt = vec_Akt[lev]->array(mfi);
        Array4<Real> const& Akv = vec_Akv[lev]->array(mfi);

        ParallelFor(makeSlab(Box(Akk),2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Akk(i,j, 0) = 0.0_rt;
            Akk(i,j, N+1) = 0.0_rt;

            Akp(i,j, 0) = 0.0_rt;
            Akp(i,j, N+1) = 0.0_rt;

            Akv(i,j, 0) = 0.0_rt;
            Akv(i,j, N+1) = 0.0_rt;

            Akt(i,j, 0) = 0.0_rt;
            Akt(i,j, N+1) = 0.0_rt;
        });
    }
}
