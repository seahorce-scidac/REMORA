#include <REMORA.H>

using namespace amrex;
/** Nonlinear shallow-water rpimitive equations predictor (Leap-frog) and
 * corrector (Adams-Moulton) time-stepping engine. Corresponds to Nonlinear/step2d_LF_AM3.h
 * in ROMS.
 *
 * @param[in   ] lev        level of refinement (coarsest level is 0)
 * @param[in   ] mf_rhoS    density perturbation
 * @param[in   ] mf_rhoA    vertically-averaged density
 * @param[inout] mf_ru2d    RHS contributions to 2D u-momentum
 * @param[inout] mf_rv2d    RHS contribtuions to 2D v-momentum
 * @param[inout] mf_rufrc   before first predictor, vertical integral of 3D RHS for uvel, converted to forcing terms
 * @param[inout] mf_rvfrc   before first predictor, vertical integral of 3D RHS for vvel, converted to forcing term
 * @param[inout] mf_Zt_avg1 average of sea surface height over all fast steps
 * @param[inout] mf_DU_avg1 time-averaged u-flux for 2D equations
 * @param[inout] mf_DU_avg2 time-averaged u-flux for 3D equation coupling
 * @param[inout] mf_DV_avg1 time-averaged v-flux for 2D equations
 * @param[inout] mf_DV_avg2 time-averaged v-flux for 3D equation coupling
 * @param[inout] mf_rubar   RHS of vertically integrated u-momentum
 * @param[inout] mf_rvbar   RHS of vertically integrated v-momentum
 * @param[inout] mf_rzeta   RHS of sea surface height
 * @param[inout] mf_ubar    vertically integrated u-momentum
 * @param[inout] mf_vbar    vertically integrated v-momentum
 * @param[inout] mf_zeta    Sea-surface height
 * @param[in   ] mf_h       Bathymetry
 * @param[in   ] mf_pm      1 / dx
 * @param[in   ] mf_pn      1 / dy
 * @param[in   ] mf_fcor    Coriolis factor
 * @param[inout] mf_visc2_p Harmonic viscosity at psi points
 * @param[inout] mf_visc2_r Harmoic viscosity at rho points
 * @param[in   ] mf_mskr    Land-sea mask at rho-points
 * @param[in   ] mf_msku    Land-sea mask at u-points
 * @param[in   ] mf_mskv    Land-sea mask at v-points
 * @param[in   ] mf_mskp    Land-sea mask at psi-points
 * @param[in   ] dtfast_lev Length of current barotropic step
 * @param[in   ] predictor_2d_step Is this a predictor step?
 * @param[in   ] first_2d_step Is this the first barotropic step?
 * @param[in   ] my_iif     Which barotropic predictor-corrector pair?
 * @param[inout] next_indx1 Cached index for
 */

void
REMORA::advance_2d (int lev,
                   MultiFab const* mf_rhoS,
                   MultiFab const* mf_rhoA,
                   MultiFab      * mf_ru2d,
                   MultiFab      * mf_rv2d,
                   MultiFab      * mf_rufrc,
                   MultiFab      * mf_rvfrc,
                   MultiFab      * mf_Zt_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg2,
                   std::unique_ptr<MultiFab>& mf_DV_avg1,
                   std::unique_ptr<MultiFab>& mf_DV_avg2,
                   std::unique_ptr<MultiFab>& mf_rubar,
                   std::unique_ptr<MultiFab>& mf_rvbar,
                   std::unique_ptr<MultiFab>& mf_rzeta,
                   std::unique_ptr<MultiFab>& mf_ubar,
                   std::unique_ptr<MultiFab>& mf_vbar,
                   MultiFab      * mf_zeta,
                   MultiFab const* mf_h,
                   MultiFab const* mf_pm,
                   MultiFab const* mf_pn,
                   MultiFab const* mf_fcor,
                   MultiFab const* mf_visc2_p,
                   MultiFab const* mf_visc2_r,
                   MultiFab const* mf_mskr,
                   MultiFab const* mf_msku,
                   MultiFab const* mf_mskv,
                   MultiFab const* mf_mskp,
                   Real dtfast_lev,
                   bool predictor_2d_step,
                   bool first_2d_step, int my_iif,
                   int & next_indx1)
{
    int iic = istep[lev];
    const int nnew  = 0;
    const int nstp  = 0;
    int ntfirst = 0;

    int knew = 3;
    int krhs = (my_iif + iic) % 2 + 1;
    int kstp = my_iif <=1 ? iic % 2 + 1 : (iic % 2 + my_iif % 2 + 1) % 2 + 1;
    int indx1 = krhs;
    if (predictor_2d_step) {
        next_indx1 = 3 - indx1;
    } else {
        knew = next_indx1;
        kstp = 3 - knew;
        krhs = 3;
        //If it's not the auxiliary time step, set indx1 to next_indx1
        // NOTE: should this ever not execute?
        // Include indx1 updates for diagnostic purposes?
        //        if (my_iif<nfast+1)
        //            indx1=next_indx1;
    }
    int ptsk = 3-kstp;
    knew-=1;
    krhs-=1;
    kstp-=1;
    // Include indx1 updates for diagnostic purposes?
    //indx1-=1;
    ptsk-=1;
    auto ba = mf_h->boxArray();
    auto dm = mf_h->DistributionMap();

    MultiFab mf_DUon(convert(ba,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_DVom(convert(ba,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0));

    const auto dlo = amrex::lbound(Geom(lev).Domain());
    const auto dhi = amrex::ubound(Geom(lev).Domain());

    for ( MFIter mfi(*mf_rhoS, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real      > const& ubar = mf_ubar->array(mfi);
        Array4<Real      > const& vbar = mf_vbar->array(mfi);
        Array4<Real      > const& zeta = mf_zeta->array(mfi);
        Array4<Real const> const& h    = mf_h->const_array(mfi);

        Array4<Real const> const& pm   = mf_pm->const_array(mfi);
        Array4<Real const> const& pn   = mf_pn->const_array(mfi);

        Box bx = mfi.tilebox();
        Box gbx = mfi.growntilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box xgbx2 = mfi.grownnodaltilebox(0, IntVect(NGROW,NGROW,0));
        Box ygbx2 = mfi.grownnodaltilebox(1, IntVect(NGROW,NGROW,0));

        Box tbxp1 = bx;
        Box tbxp11 = bx;
        Box tbxp2 = bx;
        Box tbxp3 = bx;
        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp2.grow(IntVect(NGROW,NGROW,0));
        tbxp11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));
        tbxp3.grow(IntVect(NGROW+1,NGROW+1,0));

        Box bxD   = bx  ;   bxD.makeSlab(2,0);
        Box gbxD  = gbx ;  gbxD.makeSlab(2,0);
        Box gbx1D = gbx1; gbx1D.makeSlab(2,0);
        Box gbx2D = gbx2; gbx2D.makeSlab(2,0);

        Box tbxp2D = tbxp2;
        tbxp2D.makeSlab(2,0);

        // step2d work arrays
        FArrayBox fab_Drhs(makeSlab(tbxp3,2,0),1,The_Async_Arena());
        auto Drhs=fab_Drhs.array();

        auto DUon = mf_DUon.array(mfi);
        auto DVom = mf_DVom.array(mfi);

        ParallelFor(makeSlab(tbxp3,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Drhs(i,j,0)=zeta(i,j,0,krhs)+h(i,j,0);
        });

        ParallelFor(makeSlab(xgbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real on_u = 2.0_rt / (pn(i,j,0)+pn(i-1,j,0));
            Real cff1= 0.5_rt * on_u *(Drhs(i,j,0)+Drhs(i-1,j,0));
            DUon(i,j,0)=ubar(i,j,0,krhs)*cff1;
        });

        ParallelFor(makeSlab(ygbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real om_v = 2.0_rt / (pm(i,j,0)+pm(i,j-1,0));
            Real cff1= 0.5_rt * om_v * (Drhs(i,j,0)+Drhs(i,j-1,0));
            DVom(i,j,0)=vbar(i,j,0,krhs)*cff1;
         });
    }

    // These are needed to pass the tests with bathymetry but I don't quite see why
    mf_DUon.FillBoundary(geom[lev].periodicity());
    mf_DVom.FillBoundary(geom[lev].periodicity());

#ifdef REMORA_USE_NETCDF
    if (solverChoice.do_m2_clim_nudg) {
        ubar_clim_data_from_file->update_interpolated_to_time(t_new[lev]);
        vbar_clim_data_from_file->update_interpolated_to_time(t_new[lev]);
    }
#endif

    for ( MFIter mfi(*mf_rhoS, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real const> const& rhoS = mf_rhoS->const_array(mfi);
        Array4<Real const> const& rhoA = mf_rhoA->const_array(mfi);
        Array4<Real const> const& h    = mf_h->const_array(mfi);

        Array4<Real      > const& rufrc   = mf_rufrc->array(mfi);
        Array4<Real      > const& rvfrc   = mf_rvfrc->array(mfi);
        Array4<Real      > const& Zt_avg1 = mf_Zt_avg1->array(mfi);
        Array4<Real      > const& ubar    = mf_ubar->array(mfi);
        Array4<Real      > const& vbar    = mf_vbar->array(mfi);
        Array4<Real      > const& ubar_krhs    = mf_ubar->array(mfi, krhs);
        Array4<Real      > const& vbar_krhs    = mf_vbar->array(mfi, krhs);
        Array4<Real      > const& zeta = mf_zeta->array(mfi);
        Array4<Real      > const& DU_avg1 = (mf_DU_avg1)->array(mfi);
        Array4<Real      > const& DU_avg2 = (mf_DU_avg2)->array(mfi);
        Array4<Real      > const& DV_avg1 = (mf_DV_avg1)->array(mfi);
        Array4<Real      > const& DV_avg2 = (mf_DV_avg2)->array(mfi);
        Array4<Real      > const& ru2d = (mf_ru2d)->array(mfi);
        Array4<Real      > const& rv2d = (mf_rv2d)->array(mfi);
        Array4<Real      > const& rubar = (mf_rubar)->array(mfi);
        Array4<Real      > const& rvbar = (mf_rvbar)->array(mfi);
        Array4<Real      > const& rzeta = (mf_rzeta)->array(mfi);
        Array4<Real const> const& visc2_p = mf_visc2_p->const_array(mfi);
        Array4<Real const> const& visc2_r = mf_visc2_r->const_array(mfi);

        Array4<Real const> const& pm   = mf_pm->const_array(mfi);
        Array4<Real const> const& pn   = mf_pn->const_array(mfi);
        Array4<Real const> const& fcor = mf_fcor->const_array(mfi);

        Array4<Real const> const& mskr = mf_mskr->const_array(mfi);
        Array4<Real const> const& msku = mf_msku->const_array(mfi);
        Array4<Real const> const& mskv = mf_mskv->const_array(mfi);
        Array4<Real const> const& mskp = mf_mskp->const_array(mfi);

        Box bx = mfi.tilebox();
        Box gbx = mfi.growntilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box gbx3 = mfi.growntilebox(IntVect(NGROW+1,NGROW+1,0));
        Box xgbx2 = mfi.grownnodaltilebox(0, IntVect(NGROW,NGROW,0));
        Box ygbx2 = mfi.grownnodaltilebox(1, IntVect(NGROW,NGROW,0));

        Box xbxD = mfi.nodaltilebox(0);
        xbxD.makeSlab(2,0);

        Box ybxD = mfi.nodaltilebox(1);
        ybxD.makeSlab(2,0);

        Box xbxD_adj = mfi.nodaltilebox(0);
        xbxD_adj.makeSlab(2,0);
        Box ybxD_adj = mfi.nodaltilebox(1);
        ybxD_adj.makeSlab(2,0);

        auto xbxD_lo = lbound(xbxD_adj);
        auto xbxD_hi = ubound(xbxD_adj);
        auto ybxD_lo = lbound(ybxD_adj);
        auto ybxD_hi = ubound(ybxD_adj);

        if (xbxD_lo.x == dlo.x) {
            xbxD_adj.growLo(0,-1);
        } else if (xbxD_hi.x == dhi.x) {
            xbxD_adj.growHi(0,-1);
        }

        if (ybxD_lo.y == dlo.y) {
            ybxD_adj.growLo(1,-1);
        } else if (ybxD_hi.y == dhi.y) {
            ybxD_adj.growHi(1,-1);
        }

        Box tbxp1  = bx;  tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        Box tbxp11 = bx; tbxp11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));
        Box tbxp2  = bx;  tbxp2.grow(IntVect(NGROW,NGROW,0));
        Box tbxp3  = bx;  tbxp3.grow(IntVect(NGROW+1,NGROW+1,0));

        Box bxD   = bx;   bxD.makeSlab(2,0);
        Box gbxD  = gbx;  gbxD.makeSlab(2,0);
        Box gbx1D = gbx1; gbx1D.makeSlab(2,0);
        Box gbx2D = gbx2; gbx2D.makeSlab(2,0);

        Box tbxp2D = tbxp2;
        tbxp2D.makeSlab(2,0);

        FArrayBox fab_fomn(tbxp2,1,The_Async_Arena());

        //step2d work arrays
        FArrayBox fab_Drhs(tbxp3,1,The_Async_Arena());
        FArrayBox fab_Dnew(tbxp2,1,The_Async_Arena());
        FArrayBox fab_zwrk(tbxp1,1,The_Async_Arena());
        FArrayBox fab_gzeta(tbxp1,1,The_Async_Arena());
        FArrayBox fab_gzeta2(tbxp1,1,The_Async_Arena());
        FArrayBox fab_gzetaSA(tbxp1,1,The_Async_Arena());
        FArrayBox fab_Dstp(tbxp3,1,The_Async_Arena());
        FArrayBox & fab_DUon=mf_DUon[mfi];
        FArrayBox & fab_DVom=mf_DVom[mfi];
        FArrayBox fab_rhs_ubar(xbxD,1,The_Async_Arena());
        FArrayBox fab_rhs_vbar(ybxD,1,The_Async_Arena());
        FArrayBox fab_rhs_zeta(tbxp1,1,The_Async_Arena());
        FArrayBox fab_zeta_new(tbxp1,1,The_Async_Arena());

        auto fomn=fab_fomn.array();

        auto Drhs = fab_Drhs.array();
        auto Drhs_const = fab_Drhs.const_array();
        auto Dnew=fab_Dnew.array();
        auto zwrk=fab_zwrk.array();
        auto gzeta=fab_gzeta.array();
        auto gzeta2=fab_gzeta2.array();
        auto gzetaSA=fab_gzetaSA.array();
        auto Dstp=fab_Dstp.array();
        auto DUon=fab_DUon.array();
        auto DVom=fab_DVom.array();
        auto rhs_ubar=fab_rhs_ubar.array();
        auto rhs_vbar=fab_rhs_vbar.array();
        auto rhs_zeta=fab_rhs_zeta.array();
        auto zeta_new=fab_zeta_new.array();

        auto weight1 = vec_weight1.dataPtr();
        auto weight2 = vec_weight2.dataPtr();

        //From ana_grid.h and metrics.F
        ParallelFor(xbxD, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            rhs_ubar(i,j,0)=0.0_rt;
        });

        ParallelFor(ybxD, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            rhs_vbar(i,j,0)=0.0_rt;
        });

        if (solverChoice.use_coriolis) {
            ParallelFor(tbxp2D, [=] AMREX_GPU_DEVICE (int i, int j, int  )
            {
                fomn(i,j,0) = fcor(i,j,0)*(1.0_rt/(pm(i,j,0)*pn(i,j,0)));
            });
        }

        ParallelFor(makeSlab(tbxp3,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Drhs(i,j,0)=zeta(i,j,0,krhs)+h(i,j,0);
        });

        if(predictor_2d_step)
        {
            if(first_2d_step) {
                Real cff2=(-1.0_rt/12.0_rt)*weight2[my_iif+1];
                ParallelFor(makeSlab(gbx3,2,0),
                [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    Zt_avg1(i,j,0)=0.0_rt;
                });
                ParallelFor(makeSlab(xgbx2,2,0),
                [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    DU_avg1(i,j,0)=0.0_rt;
                    DU_avg2(i,j,0)=cff2*DUon(i,j,0);
                });
                ParallelFor(makeSlab(ygbx2,2,0),
                [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    DV_avg1(i,j,0)=0.0_rt;
                    DV_avg2(i,j,0)=cff2*DVom(i,j,0);
                });
            }
            else {
                Real cff1_wt1 = weight1[my_iif-1];
                Real cff2_wt1 = (8.0_rt/12.0_rt)*weight2[my_iif]-
                                (1.0_rt/12.0_rt)*weight2[my_iif+1];

                ParallelFor(makeSlab(gbx3,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    Zt_avg1(i,j,0) += cff1_wt1*zeta(i,j,0,krhs);
                });

                ParallelFor(makeSlab(xgbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    DU_avg1(i,j,0) += cff1_wt1*DUon(i,j,0);
                    DU_avg2(i,j,0) += cff2_wt1*DUon(i,j,0);
                });

                ParallelFor(makeSlab(ygbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    DV_avg1(i,j,0) += cff1_wt1*DVom(i,j,0);
                    DV_avg2(i,j,0) += cff2_wt1*DVom(i,j,0);
                });
            }
        }
        else {
            Real cff2_wt2;

            if (first_2d_step) {
                cff2_wt2=weight2[my_iif];
            } else {
                cff2_wt2=5.0_rt/12.0_rt*weight2[my_iif];
            }

            ParallelFor(makeSlab(xgbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                DU_avg2(i,j,0)=DU_avg2(i,j,0)+cff2_wt2*DUon(i,j,0);
            });

            ParallelFor(makeSlab(ygbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                DV_avg2(i,j,0)=DV_avg2(i,j,0)+cff2_wt2*DVom(i,j,0);
            });
        }
        //
        //  Do not perform the actual time stepping during the auxiliary
        //  (nfast(ng)+1) time step. Jump to next box
        //

        if (my_iif>=nfast) {
            continue; }
        //Load new free-surface values into shared array at both predictor
        //and corrector steps
        //
        //=======================================================================
        //  Time step free-surface equation.
        //=======================================================================
        //
        //  During the first time-step, the predictor step is Forward-Euler
        //  and the corrector step is Backward-Euler. Otherwise, the predictor
        //  step is Leap-frog and the corrector step is Adams-Moulton.
        //

        // todo: gzeta

        // todo: HACKHACKHACK Should use rho0 from prob.H
        Real fac=1000.0_rt/1025.0_rt;

        if (my_iif==0) {
            Real cff1=dtfast_lev;

            ParallelFor(makeSlab(tbxp1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rhs_zeta(i,j,0) = (DUon(i,j,0)-DUon(i+1,j,0))+
                                  (DVom(i,j,0)-DVom(i,j+1,0));
                zeta_new(i,j,0) = (zeta(i,j,0,kstp)+ pm(i,j,0)*pn(i,j,0)*cff1*rhs_zeta(i,j,0)) * mskr(i,j,0);
                Dnew(i,j,0) = zeta_new(i,j,0)+h(i,j,0);

                //Pressure gradient terms:
                zwrk(i,j,0)=0.5_rt*(zeta(i,j,0,kstp)+zeta_new(i,j,0));
                gzeta(i,j,0)=(fac+rhoS(i,j,0))*zwrk(i,j,0);
                gzeta2(i,j,0)=gzeta(i,j,0)*zwrk(i,j,0);
                gzetaSA(i,j,0)=zwrk(i,j,0)*(rhoS(i,j,0)-rhoA(i,j,0));
            });

        } else if (predictor_2d_step) {

            Real cff1=2.0_rt * dtfast_lev;
            Real cff4=4.0_rt / 25.0_rt;
            Real cff5=1.0_rt - 2.0_rt*cff4;

            ParallelFor(makeSlab(tbxp1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rhs_zeta(i,j,0)=(DUon(i,j,0)-DUon(i+1,j,0))+
                                (DVom(i,j,0)-DVom(i,j+1,0));
                zeta_new(i,j,0)=(zeta(i,j,0,kstp)+
                                pm(i,j,0)*pn(i,j,0)*cff1*rhs_zeta(i,j,0)) * mskr(i,j,0);
                Dnew(i,j,0)=zeta_new(i,j,0)+h(i,j,0);
                //Pressure gradient terms
                zwrk(i,j,0)=cff5*zeta(i,j,0,krhs)+
                    cff4*(zeta(i,j,0,kstp)+zeta_new(i,j,0));
                gzeta(i,j,0)=(fac+rhoS(i,j,0))*zwrk(i,j,0);
                gzeta2(i,j,0)=gzeta(i,j,0)*zwrk(i,j,0);
                gzetaSA(i,j,0)=zwrk(i,j,0)*(rhoS(i,j,0)-rhoA(i,j,0));
            });

        } else if (!predictor_2d_step) { //AKA if(corrector_2d_step)

            Real cff1=dtfast_lev * 5.0_rt/12.0_rt;
            Real cff2=dtfast_lev * 8.0_rt/12.0_rt;
            Real cff3=dtfast_lev * 1.0_rt/12.0_rt;
            Real cff4=2.0_rt/5.0_rt;
            Real cff5=1.0_rt-cff4;

            ParallelFor(makeSlab(tbxp1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=cff1*((DUon(i,j,0)-DUon(i+1,j,0))+
                               (DVom(i,j,0)-DVom(i,j+1,0)));
                zeta_new(i,j,0)=zeta(i,j,0,kstp)+
                    pm(i,j,0)*pn(i,j,0)*(cff+
                                         cff2*rzeta(i,j,0,kstp)-
                                         cff3*rzeta(i,j,0,ptsk));
                zeta_new(i,j,0) *= mskr(i,j,0);
                Dnew(i,j,0)=zeta_new(i,j,0)+h(i,j,0);
                //Pressure gradient terms
                zwrk(i,j,0)=cff5*zeta_new(i,j,0)+cff4*zeta(i,j,0,krhs);
                gzeta(i,j,0)=(fac+rhoS(i,j,0))*zwrk(i,j,0);
                gzeta2(i,j,0)=gzeta(i,j,0)*zwrk(i,j,0);
                gzetaSA(i,j,0)=zwrk(i,j,0)*(rhoS(i,j,0)-rhoA(i,j,0));
            });
        }

        //
        //  Load new free-surface values into shared array at both predictor
        //  and corrector steps.
        //
        //// zeta(knew) only valid at zeta_new, i.e. tbxp1
        ParallelFor(makeSlab(gbx1,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            zeta(i,j,0,knew) = zeta_new(i,j,0);
        });

        //
        //  If predictor step, load right-side-term into shared array.
        //
        if (predictor_2d_step) {
            ParallelFor(makeSlab(gbx1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rzeta(i,j,0,krhs)=rhs_zeta(i,j,0);
            });
        }

        //
        //=======================================================================
        //  Compute right-hand-side for the 2D momentum equations.
        //=======================================================================
        //
/*
!
!-----------------------------------------------------------------------
!  Compute pressure gradient terms.
!-----------------------------------------------------------------------
!
*/

        Real cff1 = 0.5_rt * g;
        Real cff2 = 1.0_rt / 3.0_rt;
        ParallelFor(xbxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real on_u = 2.0_rt / (pn(i,j,0)+pn(i-1,j,0));
            rhs_ubar(i,j,0)=cff1 * on_u *
                          ((    h(i-1,j,0) +     h(i,j,0))*
                           (gzeta(i-1,j,0) - gzeta(i,j,0))+
                           (    h(i-1,j,0) -     h(i,j,0))*
                           (      gzetaSA(i-1,j,0) + gzetaSA(i,j,0)+
                            cff2*(   rhoA(i-1,j,0) -    rhoA(i,j,0))*
                                 (   zwrk(i-1,j,0) -    zwrk(i,j,0)))+
                           (gzeta2(i-1,j,0)- gzeta2(i  ,j,0)));
        });

        ParallelFor(ybxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real om_v = 2.0_rt / (pm(i,j,0)+pm(i,j-1,0));
            rhs_vbar(i,j,0) = cff1*om_v *
                          ((    h(i,j-1,0) +     h(i,j,0))*
                           (gzeta(i,j-1,0) - gzeta(i,j,0))+
                           (    h(i,j-1,0) -     h(i,j,0))*
                           (gzetaSA(i,j-1,0)+ gzetaSA(i,j  ,0)+
                            cff2*(rhoA(i,j-1,0)- rhoA(i,j  ,0))*
                                 (zwrk(i,j-1,0)- zwrk(i,j  ,0)))+
                           (gzeta2(i,j-1,0)- gzeta2(i,j  ,0)));
        });

        // Advection terms for 2d ubar, vbar added to rhs_ubar and rhs_vbar
        //
        //-----------------------------------------------------------------------
        // rhs_uv_2d
        //-----------------------------------------------------------------------
        //
        Array4<Real const> const& ubar_const = mf_ubar->const_array(mfi);
        Array4<Real const> const& vbar_const = mf_vbar->const_array(mfi);

        rhs_uv_2d(xbxD, ybxD, ubar_const, vbar_const, rhs_ubar, rhs_vbar, DUon, DVom, krhs);

        //-----------------------------------------------------------------------
        // Add Coriolis forcing
        //-----------------------------------------------------------------------
        if (solverChoice.use_coriolis) {
            // Coriolis terms for 2d ubar, vbar added to rhs_ubar and rhs_vbar
            //
            //-----------------------------------------------------------------------
            // coriolis
            //-----------------------------------------------------------------------
            //
            coriolis(xbxD, ybxD, ubar_const, vbar_const, rhs_ubar, rhs_vbar, Drhs, fomn, krhs, 0);
        }
        //-----------------------------------------------------------------------
        //Add in horizontal harmonic viscosity.
        // Consider generalizing or copying uv3dmix, where Drhs is used instead of Hz and u=>ubar v=>vbar, drop dt terms
        //-----------------------------------------------------------------------
        uv3dmix(xbxD, ybxD, ubar, vbar, ubar, vbar, rhs_ubar, rhs_vbar,
                visc2_p, visc2_r, Drhs_const,
                pm, pn, mskp, krhs, nnew, 0.0_rt);

#ifdef REMORA_USE_NETCDF
        if (solverChoice.do_m2_clim_nudg) {
            Array4<const Real> const& ubar_clim = ubar_clim_data_from_file->mf_interpolated->const_array(mfi);
            Array4<const Real> const& vbar_clim = vbar_clim_data_from_file->mf_interpolated->const_array(mfi);
            Array4<const Real> const& ubar_nudg_coeff = vec_nudg_coeff[BdyVars::ubar][lev]->const_array(mfi);
            Array4<const Real> const& vbar_nudg_coeff = vec_nudg_coeff[BdyVars::vbar][lev]->const_array(mfi);
            // Boxes are like this to match ROMS
            apply_clim_nudg(xbxD_adj, 1, 0, rhs_ubar, ubar_krhs, ubar_clim, ubar_nudg_coeff, Drhs_const, pm, pn);
            apply_clim_nudg(ybxD_adj, 0, 1, rhs_vbar, vbar_krhs, vbar_clim, vbar_nudg_coeff, Drhs_const, pm, pn);
        }
#endif

        //-----------------------------------------------------------------------
        // Coupling from 3d to 2d
        //-----------------------------------------------------------------------
        if (first_2d_step&&predictor_2d_step)
        {
            if (iic==ntfirst) {
                ParallelFor(xbxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rufrc(i,j,0)    -= rhs_ubar(i,j,0);
                    rhs_ubar(i,j,0) += rufrc(i,j,0);
                    ru2d(i,j,0,nstp)  = rufrc(i,j,0);
                });

                ParallelFor(ybxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rvfrc(i,j,0)    -= rhs_vbar(i,j,0);
                    rhs_vbar(i,j,0) += rvfrc(i,j,0);
                    rv2d(i,j,0,nstp)  = rvfrc(i,j,0);
                });

            } else if (iic==(ntfirst+1)) {

                ParallelFor(xbxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rufrc(i,j,0)=rufrc(i,j,0)-rhs_ubar(i,j,0);
                    rhs_ubar(i,j,0)=rhs_ubar(i,j,0)+1.5_rt*rufrc(i,j,0)-0.5_rt*ru2d(i,j,0,0);
                    ru2d(i,j,0,1)=rufrc(i,j,0);
                    Real r_swap= ru2d(i,j,0,1);
                    ru2d(i,j,0,1) = ru2d(i,j,0,0);
                    ru2d(i,j,0,0) = r_swap;
                });

                ParallelFor(ybxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rvfrc(i,j,0)=rvfrc(i,j,0)-rhs_vbar(i,j,0);
                    rhs_vbar(i,j,0)=rhs_vbar(i,j,0)+1.5_rt*rvfrc(i,j,0)-0.5_rt*rv2d(i,j,0,0);
                    rv2d(i,j,0,1)=rvfrc(i,j,0);
                    Real r_swap= rv2d(i,j,0,1);
                    rv2d(i,j,0,1) = rv2d(i,j,0,0);
                    rv2d(i,j,0,0) = r_swap;
                });

            } else {
                cff1=23.0_rt/12.0_rt;
                cff2=16.0_rt/12.0_rt;
                Real cff3= 5.0_rt/12.0_rt;

                ParallelFor(xbxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rufrc(i,j,0)=rufrc(i,j,0)-rhs_ubar(i,j,0);
                    rhs_ubar(i,j,0)=rhs_ubar(i,j,0)+
                        cff1*rufrc(i,j,0)-
                        cff2*ru2d(i,j,0,0)+
                        cff3*ru2d(i,j,0,1);
                    ru2d(i,j,0,1)=rufrc(i,j,0);
                    Real r_swap= ru2d(i,j,0,1);
                    ru2d(i,j,0,1) = ru2d(i,j,0,0);
                    ru2d(i,j,0,0) = r_swap;
                });

                ParallelFor(ybxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rvfrc(i,j,0)=rvfrc(i,j,0)-rhs_vbar(i,j,0);
                    rhs_vbar(i,j,0)=rhs_vbar(i,j,0)+
                          cff1*rvfrc(i,j,0)-
                          cff2*rv2d(i,j,0,0)+
                          cff3*rv2d(i,j,0,1);
                    rv2d(i,j,0,1)=rvfrc(i,j,0);

                    Real r_swap= rv2d(i,j,0,1);
                    rv2d(i,j,0,1) = rv2d(i,j,0,0);
                    rv2d(i,j,0,0) = r_swap;
                });
            }
        } else {
            ParallelFor(xbxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rhs_ubar(i,j,0) += rufrc(i,j,0);
            });

            ParallelFor(ybxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rhs_vbar(i,j,0) += rvfrc(i,j,0);
            });
        }

        //
        //=======================================================================
        //  Time step 2D momentum equations.
        //=======================================================================
        //
        //  Compute total water column depth.
        //
        ParallelFor(makeSlab(tbxp3,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
              Dstp(i,j,0)=zeta(i,j,0,kstp)+h(i,j,0);
        });

        //
        //  During the first time-step, the predictor step is Forward-Euler
        //  and the corrector step is Backward-Euler. Otherwise, the predictor
        //  step is Leap-frog and the corrector step is Adams-Moulton.
        //
        if (my_iif==0) {
            cff1=0.5_rt*dtfast_lev;
            ParallelFor(xbxD,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i-1,j,0))*(pn(i,j,0)+pn(i-1,j,0));
                Real Dnew_avg =1.0_rt/(Dnew(i,j,0)+Dnew(i-1,j,0));
                ubar(i,j,0,knew)=(ubar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i-1,j,0))+
                                  cff*cff1*rhs_ubar(i,j,0))*Dnew_avg * msku(i,j,0);
            });
            ParallelFor(ybxD,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
                Real Dnew_avg=1.0_rt/(Dnew(i,j,0)+Dnew(i,j-1,0));
                vbar(i,j,0,knew)=(vbar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i,j-1,0))+
                                  cff*cff1*rhs_vbar(i,j,0))*Dnew_avg * mskv(i,j,0);
            });

        } else if (predictor_2d_step) {

            cff1=dtfast_lev;
            ParallelFor(xbxD,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i-1,j,0))*(pn(i,j,0)+pn(i-1,j,0));
                Real Dnew_avg=1.0_rt/(Dnew(i,j,0)+Dnew(i-1,j,0));
                ubar(i,j,0,knew)=(ubar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i-1,j,0))+
                                  cff*cff1*rhs_ubar(i,j,0))*Dnew_avg * msku(i,j,0);
            });
            ParallelFor(ybxD,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
                Real Dnew_avg=1.0_rt/(Dnew(i,j,0)+Dnew(i,j-1,0));
                vbar(i,j,0,knew)=(vbar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i,j-1,0))+
                                  cff*cff1*rhs_vbar(i,j,0))*Dnew_avg * mskv(i,j,0);
            });

        } else if ((!predictor_2d_step)) {

            cff1=0.5_rt*dtfast_lev*5.0_rt/12.0_rt;
            cff2=0.5_rt*dtfast_lev*8.0_rt/12.0_rt;
            Real cff3=0.5_rt*dtfast_lev*1.0_rt/12.0_rt;
            ParallelFor(xbxD,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i-1,j,0))*(pn(i,j,0)+pn(i-1,j,0));
                Real Dnew_avg=1.0_rt/(Dnew(i,j,0)+Dnew(i-1,j,0));
                ubar(i,j,0,knew)=(ubar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i-1,j,0))+
                                 cff*(cff1*rhs_ubar(i,j,0)+
                                      cff2*rubar(i,j,0,kstp)-
                                      cff3*rubar(i,j,0,ptsk)))*Dnew_avg * msku(i,j,0);
            });
            ParallelFor(ybxD,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
                Real Dnew_avg=1.0_rt/(Dnew(i,j,0)+Dnew(i,j-1,0));
                vbar(i,j,0,knew)=(vbar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i,j-1,0))+
                                 cff*(cff1*rhs_vbar(i,j,0)+
                                      cff2*rvbar(i,j,0,kstp)-
                                      cff3*rvbar(i,j,0,ptsk)))*Dnew_avg * mskv(i,j,0);
            });
        }

        //store rhs_ubar and rhs_vbar to save later
        //
        //  If predictor step, load right-side-term into shared arrays for
        //  future use during the subsequent corrector step.
        //

        if (predictor_2d_step) {
            ParallelFor(xbxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rubar(i,j,0,krhs)=rhs_ubar(i,j,0);
            });
            ParallelFor(ybxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rvbar(i,j,0,krhs)=rhs_vbar(i,j,0);
            });
        }
    }

    // Don't do the FillPatch at the last truncated predictor step.
    // We may need to move the zeta FillPatch further up
    if (my_iif<nfast) {
        int know;
        Real dt2d;
        if (my_iif==0) {
            know = krhs;
            dt2d = dtfast_lev;
        } else if (predictor_2d_step) {
            know = krhs;
            dt2d = 2.0_rt * dtfast_lev;
        } else {
            know = kstp;
            dt2d = dtfast_lev;
        }

        FillPatch(lev, t_old[lev], *vec_ubar[lev], GetVecOfPtrs(vec_ubar), BCVars::ubar_bc, BdyVars::ubar,
                  knew, false,true, 0,know, dt2d);
        FillPatch(lev, t_old[lev], *vec_vbar[lev], GetVecOfPtrs(vec_vbar), BCVars::vbar_bc, BdyVars::vbar,
                  knew, false,true, 0,know, dt2d);
        FillPatch(lev, t_old[lev], *vec_zeta[lev], GetVecOfPtrs(vec_zeta), BCVars::zeta_bc, BdyVars::zeta,
                  knew, false,false, 0,know, dt2d);
    }
}
