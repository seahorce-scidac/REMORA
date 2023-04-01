#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// Start 3d step
//
void
ROMSX::advance_3d (int lev,
                   MultiFab& mf_u , MultiFab& mf_v ,
                   std::unique_ptr<MultiFab>& mf_ru,
                   std::unique_ptr<MultiFab>& mf_rv,
                   std::unique_ptr<MultiFab>& mf_DU_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg2,
                   std::unique_ptr<MultiFab>& mf_DV_avg1,
                   std::unique_ptr<MultiFab>& mf_DV_avg2,
                   std::unique_ptr<MultiFab>& mf_ubar,
                   std::unique_ptr<MultiFab>& mf_vbar,
                   MultiFab& mf_AK, MultiFab& mf_DC,
                   MultiFab& mf_Hzk,
                   std::unique_ptr<MultiFab>& mf_Akv,
                   std::unique_ptr<MultiFab>& mf_Hz,
                   const int ncomp, const int N, Real dt_lev)
{
    // Need to include uv3dmix

    auto geomdata  = Geom(lev).data();
    const auto dxi = Geom(lev).InvCellSizeArray();

    const int Mm = Geom(lev).Domain().size()[1];

    const int nrhs  = ncomp-1;
    const int nnew  = ncomp-1;

    int iic = istep[lev];
    int ntfirst = 0;

    for ( MFIter mfi(mf_u, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& u = mf_u.array(mfi);
        Array4<Real> const& v = mf_v.array(mfi);

        Array4<Real> const& ru_arr = mf_ru->array(mfi);
        Array4<Real> const& rv_arr = mf_rv->array(mfi);

        Array4<Real> const& AK_arr = mf_AK.array(mfi);
        Array4<Real> const& DC_arr = mf_DC.array(mfi);

        Array4<Real> const& Hzk_arr = mf_Hzk.array(mfi);
        Array4<Real> const& Akv_arr = mf_Akv->array(mfi);
        Array4<Real> const& Hz_arr  = mf_Hz->array(mfi);

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

        Box ubx = surroundingNodes(bx,0);
        Box vbx = surroundingNodes(bx,1);
        amrex::Print() << " BX " <<  bx << std::endl;
        amrex::Print() << "UBX " << ubx << std::endl;
        amrex::Print() << "VBX " << vbx << std::endl;

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_CF(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena);
        FArrayBox fab_pn(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_pm(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_fomn(gbx2,1,amrex::The_Async_Arena);
#if 0
        FArrayBox fab_Huon(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_Hvom(gbx2,1,amrex::The_Async_Arena);
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
#endif

        auto FC_arr = fab_FC.array();
        auto BC_arr = fab_BC.array();
        auto CF_arr = fab_CF.array();
        auto oHz_arr= fab_oHz.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto fomn=fab_fomn.array();
#if 0
        auto Huon=fab_Huon.array();
        auto Hvom=fab_Hvom.array();
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
#endif
        //From ana_grid.h and metrics.F

        //
        // Update to u
        //
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

        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real cff;
                if(iic==ntfirst) {
                  cff=0.25*dt_lev;
                } else if(iic==ntfirst+1) {
                  cff=0.25*dt_lev*3.0/2.0;
                } else {
                  cff=0.25*dt_lev*23.0/12.0;
                }

                DC_arr(i,j,k)=cff*(pm(i,j,0)+pm(i-1,j,0))*(pn(i,j,0)+pn(i-1,j,0));

                u(i,j,k) += DC_arr(i,j,k)*ru_arr(i,j,k,nrhs);
                v(i,j,k) += DC_arr(i,j,k)*rv_arr(i,j,k,nrhs);

                //ifdef SPLINES_VVISC is true
                u(i,j,k) *= 2.0 / (Hz_arr(i-1,j,k) + Hz_arr(i,j,k));

                //if(j>0&&j<Mm-1)
                v(i,j,k) *= 2.0 / (Hz_arr(i,j-1,k) + Hz_arr(i,j,k));
            });
        // End previous

       vert_visc_3d(bx,1,0,u,Hz_arr,Hzk_arr,oHz_arr,AK_arr,Akv_arr,BC_arr,DC_arr,FC_arr,CF_arr,nnew,N,dt_lev);
       vert_visc_3d(bx,0,1,v,Hz_arr,Hzk_arr,oHz_arr,AK_arr,Akv_arr,BC_arr,DC_arr,FC_arr,CF_arr,nnew,N,dt_lev);

    } // MFiter
}