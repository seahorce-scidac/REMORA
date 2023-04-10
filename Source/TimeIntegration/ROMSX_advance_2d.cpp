#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;
//
// Start 2d step
//
void
ROMSX::advance_2d (int lev,
                   MultiFab& mf_u, MultiFab& mf_v,
                   std::unique_ptr<MultiFab>& /*mf_ru*/,
                   std::unique_ptr<MultiFab>& /*mf_rv*/,
                   std::unique_ptr<MultiFab>& mf_Zt_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg2,
                   std::unique_ptr<MultiFab>& mf_DV_avg1,
                   std::unique_ptr<MultiFab>& mf_DV_avg2,
                   std::unique_ptr<MultiFab>& mf_rubar,
                   std::unique_ptr<MultiFab>& mf_rvbar,
                   std::unique_ptr<MultiFab>& /*mf_rzeta*/,
                   std::unique_ptr<MultiFab>& mf_ubar,
                   std::unique_ptr<MultiFab>& mf_vbar,
                   std::unique_ptr<MultiFab>& mf_zeta,
                   std::unique_ptr<MultiFab>& mf_h,
                   const int ncomp, Real dt_lev)
{
    auto geomdata  = Geom(lev).data();
    const auto dxi = Geom(lev).InvCellSizeArray();

    const int Lm = Geom(lev).Domain().size()[0];
    const int Mm = Geom(lev).Domain().size()[1];

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    int iic = istep[lev];
    bool predictor_2d_step = true;
    for(int my_iif = 0; my_iif <=1; my_iif++) {
        //    int my_iif = 1; //substep index
    int knew = 3;
    int krhs = (my_iif + iic) % 2 + 1;
    int kstp = my_iif <=1 ? iic % 2 + 1 : (iic % 2 + my_iif % 2 + 1) % 2 + 1;
    int indx1 = krhs;
    //    Print()<<knew<<"\t"<<krhs<<"\t"<<kstp<<"\t"<<indx1<<std::endl;
    knew-=1;
    krhs-=1;
    kstp-=1;
    indx1-=1;
    //Hardcode for 1 fast timestep (predictor+corrector only)
    Real weighta = 0.0;
    Real weightb = 1.0;
    Real weightc = 0.0;
    Real weightd = 0.0;
    Real Fgamma = 0.28400;
    Real  gamma = 0.00000;

    for ( MFIter mfi(mf_u, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);
        Array4<Real> const& ubar = (mf_ubar)->array(mfi);
        Array4<Real> const& vbar = (mf_vbar)->array(mfi);
        Array4<Real> const& zeta = (mf_zeta)->array(mfi);
        Array4<Real> const& h = (mf_h)->array(mfi);
        Array4<Real> const& Zt_avg1 = (mf_Zt_avg1)->array(mfi);
        Array4<Real> const& DU_avg1 = (mf_DU_avg1)->array(mfi);
        Array4<Real> const& DU_avg2 = (mf_DU_avg2)->array(mfi);
        Array4<Real> const& DV_avg1 = (mf_DV_avg1)->array(mfi);
        Array4<Real> const& DV_avg2 = (mf_DV_avg2)->array(mfi);
        Array4<Real> const& rubar = (mf_rubar)->array(mfi);
        Array4<Real> const& rvbar = (mf_rvbar)->array(mfi);

        Box bx = mfi.tilebox();
        //copy the tilebox
        Box gbx1 = bx;
        Box gbx11 = bx;
        Box gbx2 = bx;
        //make only gbx be grown to match multifabs
        gbx2.grow(IntVect(2,2,0));
        gbx1.grow(IntVect(1,1,0));
        gbx11.grow(IntVect(1,1,1));
        Box bxD = bx;
        Box ubxD = surroundingNodes(bx,0);
        Box vbxD = surroundingNodes(bx,1);
        bxD.makeSlab(2,0);
        ubxD.makeSlab(2,0);
        vbxD.makeSlab(2,0);
        //AKA
        //ubxD.setRange(2,0);
        //vbxD.setRange(2,0);

        FArrayBox fab_pn(gbx2,1,The_Async_Arena());
        FArrayBox fab_pm(gbx2,1,The_Async_Arena());
        FArrayBox fab_on_u(gbx2,1,The_Async_Arena());
        FArrayBox fab_om_v(gbx2,1,The_Async_Arena());
        FArrayBox fab_fomn(gbx2,1,The_Async_Arena());
        FArrayBox fab_Huon(gbx2,1,The_Async_Arena()); fab_Huon.setVal(0.0);
        FArrayBox fab_Hvom(gbx2,1,The_Async_Arena()); fab_Hvom.setVal(0.0);
        FArrayBox fab_oHz(gbx11,1,The_Async_Arena()); fab_oHz.setVal(0.0);

        //step2d work arrays
        FArrayBox fab_Drhs(gbx2,1,The_Async_Arena());
        FArrayBox fab_DUon(gbx2,1,The_Async_Arena());
        FArrayBox fab_DVom(gbx2,1,The_Async_Arena());

        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto fomn=fab_fomn.array();

        auto Drhs=fab_Drhs.array();
        auto DUon=fab_DUon.array();
        auto DVom=fab_DVom.array();

        //From ana_grid.h and metrics.F
        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
              pm(i,j,0)=dxi[0];
              pn(i,j,0)=dxi[1];
        });

        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
              om_v(i,j,0)=1.0/dxi[0];
              on_u(i,j,0)=1.0/dxi[1];
        });
        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int  )
            {

              const auto prob_lo         = geomdata.ProbLo();
              const auto dx              = geomdata.CellSize();

              pm(i,j,0)=dxi[0];
              pn(i,j,0)=dxi[1];
              //defined UPWELLING
              Real f0=-8.26e-5;
              Real beta=0.0;
              Real Esize=1000*(Mm);
              Real y = prob_lo[1] + (j + 0.5) * dx[1];
              Real f=fomn(i,j,0)=f0+beta*(y-.5*Esize);
              fomn(i,j,0)=f*(1.0/(pm(i,j,0)*pn(i,j,0)));
            });

        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Drhs(i,j,0)=zeta(i,j,0,krhs)+h(i,j,0);
        });
        amrex::ParallelFor(ubxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real cff=.5*on_u(i,j,0);
            Real cff1=cff*(Drhs(i,j,0)+Drhs(i-1,j,0));
            DUon(i,j,0)=ubar(i,j,0,krhs)*cff1;
        });
        amrex::ParallelFor(vbxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real cff=.5*om_v(i,j,0);
            Real cff1=cff*(Drhs(i,j,0)+Drhs(i,j-1,0));
            DVom(i,j,0)=vbar(i,j,0,krhs)*cff1;
        });
        if(predictor_2d_step)
        {
        if(my_iif==0) {
        Real cff2=(Real(-1.0)/Real(12.0))*weighta;
        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Zt_avg1(i,j,0)=0.0;
        });
        amrex::ParallelFor(ubxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            DU_avg1(i,j,0)=0.0;
            DU_avg2(i,j,0)=cff2*DUon(i,j,0);
        });
        amrex::ParallelFor(vbxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            DV_avg1(i,j,0)=0.0;
            DV_avg2(i,j,0)=cff2*DVom(i,j,0);
        });
        }
        else {
        Real cff1=weightb;
        Real cff2=(Real(8.0)/Real(12.0))*weightc-
                  (Real(1.0)/Real(12.0))*weightd;
        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if(k==0)
            Zt_avg1(i,j,0)=Zt_avg1(i,j,0)+cff1*zeta(i,j,0,krhs);
        });
        amrex::ParallelFor(ubxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            DU_avg1(i,j,0)=DU_avg1(i,j,0)+cff1*DUon(i,j,0);
            DU_avg2(i,j,0)=DU_avg2(i,j,0)+cff2*DUon(i,j,0);
        });
        amrex::ParallelFor(vbxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            DV_avg1(i,j,0)=DV_avg1(i,j,0)+cff1*DVom(i,j,0);
            DV_avg2(i,j,0)=DV_avg2(i,j,0)+cff2*DVom(i,j,0);
        });
        }
        } else {
        Real cff2;
        if(my_iif==0)
            cff2=weightc;
        else
            cff2=Real(5.0)/Real(12.0)*weightc;
        amrex::ParallelFor(ubxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            DU_avg2(i,j,0)=DU_avg2(i,j,0)+cff2*DUon(i,j,0);
        });
        amrex::ParallelFor(vbxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            DV_avg2(i,j,0)=DV_avg2(i,j,0)+cff2*DVom(i,j,0);
        });
        }
#ifdef UV_COR
        //
        //-----------------------------------------------------------------------
        // coriolis
        //-----------------------------------------------------------------------
        //
        coriolis(bxD, ubar, vbar, rubar, rvbar, Drhs, fomn, krhs);
#endif
    }
    }
}
