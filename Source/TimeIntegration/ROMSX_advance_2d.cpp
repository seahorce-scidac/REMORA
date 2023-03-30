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
                   Real dt_lev)
{
    auto geomdata  = Geom(lev).data();
    const auto dxi = Geom(lev).InvCellSizeArray();

    const int Lm = Geom(lev).Domain().size()[0];
    const int Mm = Geom(lev).Domain().size()[1];

    const int ncomp = 1;
    const int nrhs = ncomp-1;
    const int nnew = ncomp-1;
    const int nstp = ncomp-1;

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    int iic = istep[lev];
    int ntfirst = 0;

    for ( MFIter mfi(mf_u, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);

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

        //From ana_grid.h and metrics.F
        amrex::ParallelFor(gbx2, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
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
    }
}
