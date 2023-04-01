#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;
//
// Start 2d step
//
void
ROMSX::advance_2d (int lev,
                   MultiFab& mf_u, MultiFab& mf_v,
                   std::unique_ptr<MultiFab>& mf_ru,
                   std::unique_ptr<MultiFab>& mf_rv,
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
                   const int ncomp, Real dt_lev)
{
    auto geomdata  = Geom(lev).data();
    const auto dxi = Geom(lev).InvCellSizeArray();

    const int Lm = Geom(lev).Domain().size()[0];
    const int Mm = Geom(lev).Domain().size()[1];

    const int nrhs = ncomp-1;
    const int nnew = ncomp-1;
    const int nstp = ncomp-1;

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    int iic = istep[lev];
    int ntfirst = 0;
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

        Box bx = mfi.tilebox();
        //copy the tilebox
        Box gbx1 = bx;
        Box gbx11 = bx;
        Box gbx2 = bx;
        //make only gbx be grown to match multifabs
        gbx2.grow(IntVect(2,2,0));
        gbx1.grow(IntVect(1,1,0));
        gbx11.grow(IntVect(1,1,1));
        Box gbx=gbx2;

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_CF(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_pn(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_pm(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_on_u(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_om_v(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_fomn(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_Huon(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_Hvom(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena);
        //rhs3d work arrays
        FArrayBox fab_Huxx(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_Huee(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_Hvxx(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_Hvee(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_uxx(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_uee(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_vxx(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_vee(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_UFx(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_UFe(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_VFx(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_VFe(gbx2,1,amrex::The_Async_Arena);
        //step2d work arrays
        FArrayBox fab_Drhs(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_DUon(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_DVom(gbx2,1,amrex::The_Async_Arena);

        auto FC=fab_FC.array();
        auto BC=fab_BC.array();
        auto CF=fab_CF.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();
        auto fomn=fab_fomn.array();
        auto Huon=fab_Huon.array();
        auto Hvom=fab_Hvom.array();
        auto oHz_arr=fab_oHz.array();
        auto Huxx=fab_Huxx.array();
        auto Huee=fab_Huee.array();
        auto Hvxx=fab_Hvxx.array();
        auto Hvee=fab_Hvee.array();
        auto uxx=fab_uxx.array();
        auto uee=fab_uee.array();
        auto vxx=fab_vxx.array();
        auto vee=fab_vee.array();
        auto UFx=fab_UFx.array();
        auto UFe=fab_UFe.array();
        auto VFx=fab_VFx.array();
        auto VFe=fab_VFe.array();

        auto Drhs=fab_Drhs.array();
        auto DUon=fab_DUon.array();
        auto DVom=fab_DVom.array();

        //From ana_grid.h and metrics.F
        amrex::ParallelFor(gbx2, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
              const auto prob_lo         = geomdata.ProbLo();
              const auto dx              = geomdata.CellSize();

              pm(i,j,0)=dxi[0];
              pn(i,j,0)=dxi[1];
        });

        amrex::ParallelFor(gbx2, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
              om_v(i,j,0)=1.0/dxi[0];
              on_u(i,j,0)=1.0/dxi[1];
        });

        fab_Huon.setVal(0.0);
        fab_Hvom.setVal(0.0);
        fab_Huxx.setVal(0.0);
        fab_Huee.setVal(0.0);
        fab_Hvxx.setVal(0.0);
        fab_Hvee.setVal(0.0);
        fab_uxx.setVal(0.0);
        fab_uee.setVal(0.0);
        fab_UFx.setVal(0.0);
        fab_UFe.setVal(0.0);
        fab_vxx.setVal(0.0);
        fab_vee.setVal(0.0);
        fab_VFx.setVal(0.0);
        fab_VFe.setVal(0.0);

        amrex::ParallelFor(gbx2, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            Drhs(i,j,0)=zeta(i,j,0,krhs)+h(i,j,0);
        });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            Real cff=.5*on_u(i,j,0);
            Real cff1=cff*(Drhs(i,j,0)+Drhs(i-1,j,0));
            DUon(i,j,0)=ubar(i,j,0,krhs)*cff1;
        });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            Real cff=.5*om_v(i,j,0);
            Real cff1=cff*(Drhs(i,j,0)+Drhs(i,j-1,0));
            DVom(i,j,0)=vbar(i,j,0,krhs)*cff1;
        });
        if(iic!=0||my_iif>0)
        {
        if(my_iif==0) {
        Real cff2=(Real(-1.0)/Real(12.0))*weighta;
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            Zt_avg1(i,j,0)=0.0;
        });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            DU_avg1(i,j,0)=0.0;
            DU_avg2(i,j,0)=cff2*DUon(i,j,0);
        });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            DV_avg1(i,j,0)=0.0;
            DV_avg2(i,j,0)=cff2*DVom(i,j,0);
        });
        }
        else {
        Real cff1=weightb;
        Real cff2=(Real(8.0)/Real(12.0))*weightc-
                  (Real(1.0)/Real(12.0))*weightd;
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            Zt_avg1(i,j,0)=Zt_avg1(i,j,0)+cff1*zeta(i,j,0,krhs);
        });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            DU_avg1(i,j,0)=DU_avg1(i,j,0)+cff1*DUon(i,j,0);
            DU_avg2(i,j,0)=DU_avg2(i,j,0)+cff2*DUon(i,j,0);
        });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
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
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            DU_avg2(i,j,0)=DU_avg2(i,j,0)+cff2*DUon(i,j,0);
        });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            DV_avg2(i,j,0)=DV_avg2(i,j,0)+cff2*DVom(i,j,0);
        });
        }
    }
    }
}
