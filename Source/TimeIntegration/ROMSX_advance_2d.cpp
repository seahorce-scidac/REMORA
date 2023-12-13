#include <ROMSX.H>

using namespace amrex;
//
// Start 2d step
//
/**
 * Function that coordinates the evolution across levels -- this calls Advance to do the
 * actual advance at this level,  then recursively calls itself at finer levels
 *
 * @param[in] lev level of refinement (coarsest level is 0)
 * @param[in] mf_u
 * @param[in] mf_v
 * @param[in] mf_rhoS
 * @param[in] mf_rhoA
 * @param[in] mf_ru
 * @param[in] mf_rv
 * @param[inout] mf_rufrc
 * @param[inout] mf_rvfrc
 * @param[inout] mf_Zt_avg1
 * @param[inout] mf_DU_avg1
 * @param[inout] mf_DU_avg2
 * @param[inout] mf_DV_avg1
 * @param[inout] mf_DV_avg2
 * @param[inout] mf_rubar
 * @param[inout] mf_rvbar
 * @param[inout] mf_rbar
 * @param[inout] mf_rbar
 * @param[inout] mf_zeta
 * @param[inout] mf_h
 * @param[inout] mf_visc2_p
 * @param[inout] mf_visc2_f
 * @param[in   ] ncomp
 * @param[in   ] dtfast_lev
 * @param[in   ] predictor_2d_step
 * @param[in   ] first_2d_step
 * @param[in   ] my_iif
 * @param[in   ] next_indx1
 */

void
ROMSX::advance_2d (int lev,
                   MultiFab& mf_u, MultiFab& mf_v,
                   std::unique_ptr<MultiFab>& mf_rhoS,
                   std::unique_ptr<MultiFab>& mf_rhoA,
                   std::unique_ptr<MultiFab>& mf_ru,
                   std::unique_ptr<MultiFab>& mf_rv,
                   std::unique_ptr<MultiFab>& mf_rufrc,
                   std::unique_ptr<MultiFab>& mf_rvfrc,
                   std::unique_ptr<MultiFab>& mf_Zt_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg2,
                   std::unique_ptr<MultiFab>& mf_DV_avg1,
                   std::unique_ptr<MultiFab>& mf_DV_avg2,
                   std::unique_ptr<MultiFab>& mf_rubar,
                   std::unique_ptr<MultiFab>& mf_rvbar,
                   std::unique_ptr<MultiFab>& mf_rzeta,
                   std::unique_ptr<MultiFab>& mf_ubar,
                   std::unique_ptr<MultiFab>& mf_vbar,
                   std::unique_ptr<MultiFab>& mf_zeta,
                   std::unique_ptr<MultiFab>& mf_h,
                   std::unique_ptr<MultiFab>& mf_visc2_p,
                   std::unique_ptr<MultiFab>& mf_visc2_r,
                   const int ncomp, Real dtfast_lev,
                   bool predictor_2d_step,
                   bool first_2d_step, int my_iif,
                   int & next_indx1)
{
    auto geomdata  = Geom(lev).data();
    const auto dxi = Geom(lev).InvCellSizeArray();

    const int Mm = Geom(lev).Domain().size()[1];

    int iic = istep[lev];
    const int nnew  = ncomp-1;
    const int nstp  = ncomp-1;
    int ntfirst = 0;
    //bool predictor_2d_step = true;
    //    int my_iif = 1; //substep index
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
        if (my_iif<nfast+1)
            indx1=next_indx1;
    }
    int ptsk = 3-kstp;
    knew-=1;
    krhs-=1;
    kstp-=1;
    indx1-=1;
    ptsk-=1;
    auto ba = mf_h->boxArray();
    auto dm = mf_h->DistributionMap();

    MultiFab mf_DUon(convert(ba,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_DVom(convert(ba,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0));

    for ( MFIter mfi(*mf_rhoS, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);
        Array4<Real> const& ubar = (mf_ubar)->array(mfi);
        Array4<Real> const& vbar = (mf_vbar)->array(mfi);
        Array4<Real> const& zeta = (mf_zeta)->array(mfi);
        Array4<Real> const& h = (mf_h)->array(mfi);
        Array4<Real> const& ru = (mf_ru)->array(mfi);
        Array4<Real> const& rv = (mf_rv)->array(mfi);

        Box bx = mfi.tilebox();
        Box gbx = mfi.growntilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box xbxD = mfi.nodaltilebox(0).makeSlab(2,0);
        Box ybxD = mfi.nodaltilebox(1).makeSlab(2,0);
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

        Box bxD = bx;
        bxD.makeSlab(2,0);
        Box gbxD = gbx;
        gbxD.makeSlab(2,0);
        Box gbx1D = gbx1;
        gbx1D.makeSlab(2,0);
        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);

        Box tbxp2D = tbxp2;
        tbxp2D.makeSlab(2,0);

        FArrayBox fab_pn(tbxp2,1,The_Async_Arena());
        FArrayBox fab_pm(tbxp2,1,The_Async_Arena());
        FArrayBox fab_on_u(tbxp3,1,The_Async_Arena());
        FArrayBox fab_om_v(tbxp3,1,The_Async_Arena());
        FArrayBox fab_om_u(tbxp3,1,amrex::The_Async_Arena());
        FArrayBox fab_on_v(tbxp3,1,amrex::The_Async_Arena());
        FArrayBox fab_om_r(tbxp3,1,amrex::The_Async_Arena());
        FArrayBox fab_on_r(tbxp3,1,amrex::The_Async_Arena());
        FArrayBox fab_om_p(tbxp3,1,amrex::The_Async_Arena());
        FArrayBox fab_on_p(tbxp3,1,amrex::The_Async_Arena());

        //step2d work arrays
        FArrayBox fab_Drhs(tbxp3,1,The_Async_Arena());
        FArrayBox fab_Dnew(tbxp1,1,The_Async_Arena());
        FArrayBox fab_gzeta(tbxp1,1,The_Async_Arena());
        FArrayBox fab_gzeta2(tbxp1,1,The_Async_Arena());
        FArrayBox fab_gzetaSA(tbxp1,1,The_Async_Arena());
        FArrayBox fab_Dstp(tbxp2,1,The_Async_Arena());
        FArrayBox & fab_DUon=mf_DUon[mfi];
        FArrayBox & fab_DVom=mf_DVom[mfi];
        FArrayBox fab_rhs_ubar(xbxD,1,The_Async_Arena());
        FArrayBox fab_rhs_vbar(ybxD,1,The_Async_Arena());
        FArrayBox fab_rhs_zeta(tbxp1,1,The_Async_Arena());
        FArrayBox fab_zeta_new(tbxp1,1,The_Async_Arena());

        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto om_u=fab_om_u.array();
        auto on_v=fab_on_v.array();
        auto om_r=fab_om_r.array();
        auto on_r=fab_on_r.array();
        auto om_p=fab_om_p.array();
        auto on_p=fab_on_p.array();

        auto Drhs=fab_Drhs.array();
        auto DUon=fab_DUon.array();
        auto DVom=fab_DVom.array();

        if ((verbose > 2) && predictor_2d_step && my_iif == 0) {
            amrex::PrintToFile("ru_startadvance2").SetPrecision(18)<<FArrayBox(ru)<<std::endl;
            amrex::PrintToFile("rv_startadvance2").SetPrecision(18)<<FArrayBox(rv)<<std::endl;
            amrex::PrintToFile("u_startadvance2").SetPrecision(18)<<FArrayBox(u)<<std::endl;
            amrex::PrintToFile("v_startadvance2").SetPrecision(18)<<FArrayBox(v)<<std::endl;
            amrex::PrintToFile("ubar_startadvance2").SetPrecision(18)<<FArrayBox(ubar)<<std::endl;
            amrex::PrintToFile("vbar_startadvance2").SetPrecision(18)<<FArrayBox(vbar)<<std::endl;
        }
        if (verbose > 2) {
            ParallelFor(gbx2D,
            [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                printf("%d %d  %15.15g %15.15g %15.15g all zeta start adv2 %d %d\n", i,j,zeta(i,j,0,0),zeta(i,j,0,1), zeta(i,j,0,2),my_iif,predictor_2d_step);
            });
        }

        ParallelFor(makeSlab(xgbx2,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
              DUon(i,j,0)=0.0;
        });
        ParallelFor(makeSlab(ygbx2,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
              DVom(i,j,0)=0.0;
        });

        // NOTE: these will eventually be on smaller boxes, since these
        // are all centered on different points
        ParallelFor(makeSlab(tbxp3,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
              //Note: are the comment definitons right? Don't seem to match metrics.f90
              om_v(i,j,0)=1.0/dxi[0]; // 2/(pm(i,j-1)+pm(i,j))
              on_u(i,j,0)=1.0/dxi[1]; // 2/(pm(i,j-1)+pm(i,j))
              om_r(i,j,0)=1.0/dxi[0]; // 1/pm(i,j)
              on_r(i,j,0)=1.0/dxi[1]; // 1/pn(i,j)
              //todo: om_p on_p
              om_p(i,j,0)=1.0/dxi[0]; // 4/(pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))
              on_p(i,j,0)=1.0/dxi[1]; // 4/(pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))
              on_v(i,j,0)=1.0/dxi[1]; // 2/(pn(i-1,j)+pn(i,j))
              om_u(i,j,0)=1.0/dxi[0]; // 2/(pm(i-1,j)+pm(i,j))
        });


        if (verbose > 1) {
            Print() << "gbx2D advance_2d " << gbx2D << std::endl;
            Print() << "tbxp2D advance_2d " << tbxp2D << std::endl;
        }

        if (verbose > 2) {
            ParallelFor(tbxp2D,
            [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                printf("%d %d  %15.15g zeta start adv2  %d %d \n", i,j, zeta(i,j,0,krhs), my_iif, predictor_2d_step);
            });
        }

        ParallelFor(makeSlab(tbxp3,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Drhs(i,j,0)=zeta(i,j,0,krhs)+h(i,j,0);
        });
        ParallelFor(makeSlab(xgbx2,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real cff1= 0.5 * on_u(i,j,0) *(Drhs(i,j,0)+Drhs(i-1,j,0));
            DUon(i,j,0)=ubar(i,j,0,krhs)*cff1;
        });
        ParallelFor(makeSlab(ygbx2,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real cff1= 0.5*om_v(i,j,0)*(Drhs(i,j,0)+Drhs(i,j-1,0));
            DVom(i,j,0)=vbar(i,j,0,krhs)*cff1;
        });
    }
    mf_DUon.FillBoundary(geom[lev].periodicity());
    mf_DVom.FillBoundary(geom[lev].periodicity());


    for ( MFIter mfi(*mf_rhoS, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);
        Array4<Real> const& rhoS = (mf_rhoS)->array(mfi);
        Array4<Real> const& rhoA = (mf_rhoA)->array(mfi);
        Array4<Real> const& ubar = (mf_ubar)->array(mfi);
        Array4<Real> const& vbar = (mf_vbar)->array(mfi);
        Array4<Real> const& zeta = (mf_zeta)->array(mfi);
        Array4<Real> const& h = (mf_h)->array(mfi);
        Array4<Real> const& Zt_avg1 = (mf_Zt_avg1)->array(mfi);
        Array4<Real> const& DU_avg1 = (mf_DU_avg1)->array(mfi);
        Array4<Real> const& DU_avg2 = (mf_DU_avg2)->array(mfi);
        Array4<Real> const& DV_avg1 = (mf_DV_avg1)->array(mfi);
        Array4<Real> const& DV_avg2 = (mf_DV_avg2)->array(mfi);
        Array4<Real> const& ru = (mf_ru)->array(mfi);
        Array4<Real> const& rv = (mf_rv)->array(mfi);
        Array4<Real> const& rufrc = (mf_rufrc)->array(mfi);
        Array4<Real> const& rvfrc = (mf_rvfrc)->array(mfi);
        Array4<Real> const& rubar = (mf_rubar)->array(mfi);
        Array4<Real> const& rvbar = (mf_rvbar)->array(mfi);
        Array4<Real> const& rzeta = (mf_rzeta)->array(mfi);
        Array4<Real> const& visc2_p = (mf_visc2_p)->array(mfi);
        Array4<Real> const& visc2_r = (mf_visc2_r)->array(mfi);

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

        Box tbxp1 = bx;
        Box tbxp11 = bx;
        Box tbxp2 = bx;
        Box tbxp3 = bx;
        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp2.grow(IntVect(NGROW,NGROW,0));
        tbxp11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));
        tbxp3.grow(IntVect(NGROW+1,NGROW+1,0));

        Box bxD = bx;
        bxD.makeSlab(2,0);
        Box gbxD = gbx;
        gbxD.makeSlab(2,0);
        Box gbx1D = gbx1;
        gbx1D.makeSlab(2,0);
        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);

        Box tbxp2D = tbxp2;
        tbxp2D.makeSlab(2,0);

        FArrayBox fab_pn(tbxp3,1,The_Async_Arena());
        FArrayBox fab_pm(tbxp3,1,The_Async_Arena());
        FArrayBox fab_on_u(tbxp3,1,The_Async_Arena());
        FArrayBox fab_om_v(tbxp3,1,The_Async_Arena());
        FArrayBox fab_fomn(tbxp2,1,The_Async_Arena());
        FArrayBox fab_om_u(tbxp3,1,amrex::The_Async_Arena());
        FArrayBox fab_on_v(tbxp3,1,amrex::The_Async_Arena());
        FArrayBox fab_om_r(tbxp3,1,amrex::The_Async_Arena());
        FArrayBox fab_on_r(tbxp3,1,amrex::The_Async_Arena());
        FArrayBox fab_om_p(tbxp3,1,amrex::The_Async_Arena());
        FArrayBox fab_on_p(tbxp3,1,amrex::The_Async_Arena());

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

        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto fomn=fab_fomn.array();
        auto om_u=fab_om_u.array();
        auto on_v=fab_on_v.array();
        auto om_r=fab_om_r.array();
        auto on_r=fab_on_r.array();
        auto om_p=fab_om_p.array();
        auto on_p=fab_on_p.array();

        auto Drhs=fab_Drhs.array();
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
       if ((verbose > 2) && predictor_2d_step && my_iif == 0) {
           amrex::PrintToFile("ru_startadvance").SetPrecision(18)<<FArrayBox(ru)<<std::endl;
           amrex::PrintToFile("rv_startadvance").SetPrecision(18)<<FArrayBox(rv)<<std::endl;
           amrex::PrintToFile("u_startadvance2").SetPrecision(18)<<FArrayBox(u)<<std::endl;
           amrex::PrintToFile("v_startadvance2").SetPrecision(18)<<FArrayBox(v)<<std::endl;
       }

       if (verbose > 2) {
           PrintToFile("DUon_mid") << "step  " << my_iif << " " << predictor_2d_step << std::endl;
           PrintToFile("DUon_mid").SetPrecision(18) << FArrayBox(DUon) << std::endl;
           PrintToFile("DVom_mid") << "step  " << my_iif << " " << predictor_2d_step << std::endl;
           PrintToFile("DVom_mid").SetPrecision(18) << FArrayBox(DVom) << std::endl;
       }

        //From ana_grid.h and metrics.F
        ParallelFor(makeSlab(tbxp3,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
              pm(i,j,0)=dxi[0];
              pn(i,j,0)=dxi[1];
        });
        ParallelFor(xbxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
              rhs_ubar(i,j,0)=0.0;
        });
        ParallelFor(ybxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
              rhs_vbar(i,j,0)=0.0;
        });

        ParallelFor(makeSlab(tbxp3,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
              //Note: are the comment definitons right? Don't seem to match metrics.f90
              om_v(i,j,0)=1.0/dxi[0]; // 2/(pm(i,j-1)+pm(i,j))
              on_u(i,j,0)=1.0/dxi[1]; // 2/(pm(i,j-1)+pm(i,j))
              om_r(i,j,0)=1.0/dxi[0]; // 1/pm(i,j)
              on_r(i,j,0)=1.0/dxi[1]; // 1/pn(i,j)
              //todo: om_p on_p
              om_p(i,j,0)=1.0/dxi[0]; // 4/(pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))
              on_p(i,j,0)=1.0/dxi[1]; // 4/(pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))
              on_v(i,j,0)=1.0/dxi[1]; // 2/(pn(i-1,j)+pn(i,j))
              om_u(i,j,0)=1.0/dxi[0]; // 2/(pm(i-1,j)+pm(i,j))
        });
        ParallelFor(tbxp2D,
        [=] AMREX_GPU_DEVICE (int i, int j, int  )
        {

              const auto prob_lo         = geomdata.ProbLo();
              const auto dx              = geomdata.CellSize();

              //defined UPWELLING
              Real Esize=geomdata.ProbHi()[1] - geomdata.ProbLo()[1];
              Real y = prob_lo[1] + (j + 0.5) * dx[1];
              Real f=solverChoice.coriolis_f0 + solverChoice.coriolis_beta*(y-.5*Esize);
              fomn(i,j,0)=f*(1.0/(pm(i,j,0)*pn(i,j,0)));
        });

        ParallelFor(makeSlab(tbxp3,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Drhs(i,j,0)=zeta(i,j,0,krhs)+h(i,j,0);
        });

        if(predictor_2d_step)
        {
            if(first_2d_step) {
                Real cff2=(Real(-1.0)/Real(12.0))*weight2[my_iif+1];
                ParallelFor(makeSlab(gbx3,2,0),
                [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    Zt_avg1(i,j,0)=0.0;
                });
                ParallelFor(makeSlab(xgbx2,2,0),
                [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    DU_avg1(i,j,0)=0.0;
                    DU_avg2(i,j,0)=cff2*DUon(i,j,0);
                });
                ParallelFor(makeSlab(ygbx2,2,0),
                [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    DV_avg1(i,j,0)=0.0;
                    DV_avg2(i,j,0)=cff2*DVom(i,j,0);
                });
            }
            else {
                Real cff1_wt1 = weight1[my_iif-1];
                Real cff2_wt1 = (Real(8.0)/Real(12.0))*weight2[my_iif]-
                                (Real(1.0)/Real(12.0))*weight2[my_iif+1];
                ParallelFor(makeSlab(gbx3,2,0),
                [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    Zt_avg1(i,j,0)=Zt_avg1(i,j,0)+cff1_wt1*zeta(i,j,0,krhs);
                });
                ParallelFor(makeSlab(xgbx2,2,0),
                [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    DU_avg1(i,j,0)=DU_avg1(i,j,0)+cff1_wt1*DUon(i,j,0);
                    DU_avg2(i,j,0)=DU_avg2(i,j,0)+cff2_wt1*DUon(i,j,0);
                });
                ParallelFor(makeSlab(ygbx2,2,0),
                [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    DV_avg1(i,j,0)=DV_avg1(i,j,0)+cff1_wt1*DVom(i,j,0);
                    DV_avg2(i,j,0)=DV_avg2(i,j,0)+cff2_wt1*DVom(i,j,0);
                });
            }
        }
        else {
            Real cff2_wt2;
            if(first_2d_step)
                cff2_wt2=weight2[my_iif];
            else
                cff2_wt2=Real(5.0)/Real(12.0)*weight2[my_iif];
            ParallelFor(makeSlab(xgbx2,2,0),
            [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                DU_avg2(i,j,0)=DU_avg2(i,j,0)+cff2_wt2*DUon(i,j,0);
            });
            ParallelFor(makeSlab(ygbx2,2,0),
            [=] AMREX_GPU_DEVICE (int i, int j, int)
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
        Real fac=1000.0/1025.0;

        if(my_iif==0) {
            Real cff1=dtfast_lev;
            ParallelFor(tbxp1,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rhs_zeta(i,j,0)=(DUon(i,j,0)-DUon(i+1,j,0))+
                                (DVom(i,j,0)-DVom(i,j+1,0));
                zeta_new(i,j,0)=zeta(i,j,0,kstp)+
                                pm(i,j,0)*pn(i,j,0)*cff1*rhs_zeta(i,j,0);
                Dnew(i,j,0)=zeta_new(i,j,0)+h(i,j,0);

                //Pressure gradient terms:
                zwrk(i,j,0)=0.5_rt*(zeta(i,j,0,kstp)+zeta_new(i,j,0));
                gzeta(i,j,0)=(fac+rhoS(i,j,0))*zwrk(i,j,0);
                gzeta2(i,j,0)=gzeta(i,j,0)*zwrk(i,j,0);
                gzetaSA(i,j,0)=zwrk(i,j,0)*(rhoS(i,j,0)-rhoA(i,j,0));
            });
        } else if (predictor_2d_step) {
            Real cff1=2.0_rt*dtfast_lev;
            Real cff4=4.0/25.0;
            Real cff5=1.0-2.0*cff4;
            ParallelFor(tbxp1,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                rhs_zeta(i,j,0)=(DUon(i,j,0)-DUon(i+1,j,0))+
                                (DVom(i,j,0)-DVom(i,j+1,0));
                zeta_new(i,j,0)=zeta(i,j,0,kstp)+
                                pm(i,j,0)*pn(i,j,0)*cff1*rhs_zeta(i,j,0);
                Dnew(i,j,0)=zeta_new(i,j,0)+h(i,j,0);
                //Pressure gradient terms
                zwrk(i,j,0)=cff5*zeta(i,j,0,krhs)+
                    cff4*(zeta(i,j,0,kstp)+zeta_new(i,j,0));
                gzeta(i,j,0)=(fac+rhoS(i,j,0))*zwrk(i,j,0);
                gzeta2(i,j,0)=gzeta(i,j,0)*zwrk(i,j,0);
                gzetaSA(i,j,0)=zwrk(i,j,0)*(rhoS(i,j,0)-rhoA(i,j,0));
            });
        } else if (!predictor_2d_step) { //AKA if(corrector_2d_step)
            Real cff1=dtfast_lev*5.0_rt/12.0_rt;
            Real cff2=dtfast_lev*8.0_rt/12.0_rt;
            Real cff3=dtfast_lev*1.0_rt/12.0_rt;
            Real cff4=2.0_rt/5.0_rt;
            Real cff5=1.0_rt-cff4;
            ParallelFor(tbxp1,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=cff1*((DUon(i,j,0)-DUon(i+1,j,0))+
                               (DVom(i,j,0)-DVom(i,j+1,0)));
                zeta_new(i,j,0)=zeta(i,j,0,kstp)+
                    pm(i,j,0)*pn(i,j,0)*(cff+
                                         cff2*rzeta(i,j,0,kstp)-
                                         cff3*rzeta(i,j,0,ptsk));
                Dnew(i,j,0)=zeta_new(i,j,0)+h(i,j,0);
                //Pressure gradient terms
                zwrk(i,j,0)=cff5*zeta_new(i,j,0)+cff4*zeta(i,j,0,krhs);
                gzeta(i,j,0)=(fac+rhoS(i,j,0))*zwrk(i,j,0);
                gzeta2(i,j,0)=gzeta(i,j,0)*zwrk(i,j,0);
                gzetaSA(i,j,0)=zwrk(i,j,0)*(rhoS(i,j,0)-rhoA(i,j,0));
            });
        }

        if (verbose > 2) {
            PrintToFile("zeta_new").SetPrecision(18) << FArrayBox(zeta_new) << std::endl;
        }
        //
        //  Load new free-surface values into shared array at both predictor
        //  and corrector steps.
        //
        //// zeta(knew) only valid at zeta_new, i.e. tbxp1
        ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            zeta(i,j,0,knew)=zeta_new(i,j,0);
        });

        //
        //  If predictor step, load right-side-term into shared array.
        //
        if (predictor_2d_step) {
            ParallelFor(gbx1,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
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

        Real cff1 = 0.5 * solverChoice.g; // Should be the variable gravitational field strength
        Real cff2 = 1.0 / 3.0;
        ParallelFor(xbxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
          rhs_ubar(i,j,0)=cff1*on_u(i,j,0)*
                        ((h(i-1,j,0)+
                          h(i ,j,0))*
                         (gzeta(i-1,j,0)-
                          gzeta(i  ,j,0))+
                         (h(i-1,j,0)-
                          h(i  ,j,0))*
                         (gzetaSA(i-1,j,0)+
                          gzetaSA(i  ,j,0)+
                          cff2*(rhoA(i-1,j,0)-
                                rhoA(i  ,j,0))*
                               (zwrk(i-1,j,0)-
                                zwrk(i  ,j,0)))+
                         (gzeta2(i-1,j,0)-
                          gzeta2(i  ,j,0)));
        });

        ParallelFor(ybxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            rhs_vbar(i,j,0)=cff1*om_v(i,j,0)*
                          ((h(i,j-1,0)+
                            h(i,j  ,0))*
                           (gzeta(i,j-1,0)-
                            gzeta(i,j  ,0))+
                           (h(i,j-1,0)-
                            h(i,j  ,0))*
                           (gzetaSA(i,j-1,0)+
                            gzetaSA(i,j  ,0)+
                            cff2*(rhoA(i,j-1,0)-
                                  rhoA(i,j  ,0))*
                                 (zwrk(i,j-1,0)-
                                  zwrk(i,j  ,0)))+
                           (gzeta2(i,j-1,0)-
                            gzeta2(i,j  ,0)));
        });

       // Advection terms for 2d ubar, vbar added to rhs_ubar and rhs_vbar
       //
       //-----------------------------------------------------------------------
       // rhs_uv_2d
       //-----------------------------------------------------------------------
       //
       rhs_uv_2d(xbxD, ybxD, ubar, vbar, rhs_ubar, rhs_vbar, DUon, DVom, krhs);

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
            coriolis(xbxD, ybxD, ubar, vbar, rhs_ubar, rhs_vbar, Drhs, fomn, krhs, 0);
       }

       //-----------------------------------------------------------------------
       //Add in horizontal harmonic viscosity.
       // Consider generalizing or copying uv3dmix, where Drhs is used instead of Hz and u=>ubar v=>vbar, drop dt terms
       //-----------------------------------------------------------------------

       uv3dmix(xbxD, ybxD, ubar, vbar, ubar, vbar, rhs_ubar, rhs_vbar, visc2_p, visc2_r, Drhs, on_r, om_r, on_p, om_p, pn, pm, krhs, nnew, 0.0);

       //-----------------------------------------------------------------------
       // Coupling from 3d to 2d
       //-----------------------------------------------------------------------
       if (first_2d_step&&predictor_2d_step)
       {
            if (iic==ntfirst) {
                ParallelFor(xbxD,
                [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rufrc(i,j,0)=rufrc(i,j,0)-rhs_ubar(i,j,0);
                    rhs_ubar(i,j,0)=rhs_ubar(i,j,0)+rufrc(i,j,0);
                    ru(i,j,-1,nstp)=rufrc(i,j,0);
                });
                ParallelFor(ybxD,
                [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rvfrc(i,j,0)=rvfrc(i,j,0)-rhs_vbar(i,j,0);
                    rhs_vbar(i,j,0)=rhs_vbar(i,j,0)+rvfrc(i,j,0);
                    rv(i,j,-1,nstp)=rvfrc(i,j,0);
                });

            } else if (iic==(ntfirst+1)) {

                ParallelFor(xbxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rufrc(i,j,0)=rufrc(i,j,0)-rhs_ubar(i,j,0);
                    rhs_ubar(i,j,0)=rhs_ubar(i,j,0)+1.5_rt*rufrc(i,j,0)-0.5_rt*ru(i,j,-1,0);
                    ru(i,j,-1,1)=rufrc(i,j,0);
                    Real r_swap= ru(i,j,-1,1);
                    ru(i,j,-1,1) = ru(i,j,-1,0);
                    ru(i,j,-1,0) = r_swap;
                });

                ParallelFor(ybxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rvfrc(i,j,0)=rvfrc(i,j,0)-rhs_vbar(i,j,0);
                    rhs_vbar(i,j,0)=rhs_vbar(i,j,0)+1.5_rt*rvfrc(i,j,0)-0.5_rt*rv(i,j,-1,0);
                    rv(i,j,-1,1)=rvfrc(i,j,0);
                    Real r_swap= rv(i,j,-1,1);
                    rv(i,j,-1,1) = rv(i,j,-1,0);
                    rv(i,j,-1,0) = r_swap;
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
                        cff2*ru(i,j,-1,0)+
                        cff3*ru(i,j,-1,1);
                    ru(i,j,-1,1)=rufrc(i,j,0);
                    Real r_swap= ru(i,j,-1,1);
                    ru(i,j,-1,1) = ru(i,j,-1,0);
                    ru(i,j,-1,0) = r_swap;
                });

                ParallelFor(ybxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    rvfrc(i,j,0)=rvfrc(i,j,0)-rhs_vbar(i,j,0);
                    rhs_vbar(i,j,0)=rhs_vbar(i,j,0)+
                          cff1*rvfrc(i,j,0)-
                          cff2*rv(i,j,-1,0)+
                          cff3*rv(i,j,-1,1);
                    rv(i,j,-1,1)=rvfrc(i,j,0);

                    Real r_swap= rv(i,j,-1,1);
                    rv(i,j,-1,1) = rv(i,j,-1,0);
                    rv(i,j,-1,0) = r_swap;
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
                                  cff*cff1*rhs_ubar(i,j,0))*Dnew_avg;
            });
            ParallelFor(ybxD,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
                Real Dnew_avg=1.0_rt/(Dnew(i,j,0)+Dnew(i,j-1,0));
                vbar(i,j,0,knew)=(vbar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i,j-1,0))+
                                  cff*cff1*rhs_vbar(i,j,0))*Dnew_avg;
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
                                  cff*cff1*rhs_ubar(i,j,0))*Dnew_avg;
            });
            ParallelFor(ybxD,
            [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real cff=(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
                Real Dnew_avg=1.0_rt/(Dnew(i,j,0)+Dnew(i,j-1,0));
                vbar(i,j,0,knew)=(vbar(i,j,0,kstp)*
                                 (Dstp(i,j,0)+Dstp(i,j-1,0))+
                                  cff*cff1*rhs_vbar(i,j,0))*Dnew_avg;
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
                                      cff3*rubar(i,j,0,ptsk)))*Dnew_avg;
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
                                      cff3*rvbar(i,j,0,ptsk)))*Dnew_avg;
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
    mf_DU_avg1->FillBoundary(geom[lev].periodicity());
    mf_DU_avg2->FillBoundary(geom[lev].periodicity());
    mf_DV_avg1->FillBoundary(geom[lev].periodicity());
    mf_DV_avg2->FillBoundary(geom[lev].periodicity());
    mf_ru->FillBoundary(geom[lev].periodicity());
    mf_rv->FillBoundary(geom[lev].periodicity());
    mf_rubar->FillBoundary(geom[lev].periodicity());
    mf_rvbar->FillBoundary(geom[lev].periodicity());
    mf_rzeta->FillBoundary(geom[lev].periodicity());
    mf_ubar->FillBoundary(geom[lev].periodicity());
    mf_vbar->FillBoundary(geom[lev].periodicity());
    mf_zeta->FillBoundary(geom[lev].periodicity());

}
