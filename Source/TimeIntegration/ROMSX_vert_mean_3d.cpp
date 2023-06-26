#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// vert_mean_3d
//

void
ROMSX::vert_mean_3d (const Box& phi_bx, const int ioff, const int joff,
                     Array4<Real> phi,
                     Array4<Real> Hz, Array4<Real> Hzk,
                     Array4<Real> Dphi_avg1,
                     Array4<Real> /*AK*/, Array4<Real> /*Akv*/,
                     Array4<Real> /*BC*/, Array4<Real> DC,
                     Array4<Real> /*FC*/, Array4<Real> CF,
                     Array4<Real> dxlen,
                     const int nnew, const int N, const Real dt_lev)
{
    // Interior points vertical mean correction
    //
    // Put Hzk on the x- or y-face as appropriate, or leave on cell center for tracers
    //
    amrex::Print() << "updating on box in vert_mean_3d: " << phi_bx << std::endl;
    auto phi_bxD=phi_bx;
    phi_bxD.makeSlab(2,0);

    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Hzk(i,j,k)=0.5*(Hz(i-ioff,j-joff,k)+Hz(i,j,k));
    });
    Gpu::streamSynchronize();
    auto the_arena = amrex::The_Arena();
    FArrayBox fab_CF_tmp(Box(CF),1,the_arena); //fab_grad.setVal(0.0);
    FArrayBox fab_DC_tmp(Box(DC),1,the_arena); //fab_grad.setVal(0.0);
    FArrayBox fab_Dphi_avg1_tmp(Box(Dphi_avg1),1,the_arena); //fab_grad.setVal(0.0);
    FArrayBox fab_dxlen_tmp(Box(dxlen),1,the_arena); //fab_grad.setVal(0.0);
    Array4<Real> const & CF_tmp = fab_CF_tmp.array();
    Array4<Real> const & DC_tmp = fab_DC_tmp.array();
    Array4<Real> const & Dphi_avg1_tmp = fab_Dphi_avg1_tmp.array();
    Array4<Real> const & dxlen_tmp = fab_dxlen_tmp.array();

    amrex::ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
      CF(i,j,-1)=Hzk(i,j,0);
      DC(i,j,-1)=phi(i,j,0,nnew)*Hzk(i,j,0);
      for(int k=1; k<=N; k++) {
        CF(i,j,-1)+=Hzk(i,j,k);
        DC(i,j,-1)+=phi(i,j,k,nnew)*Hzk(i,j,k);
        CF_tmp(i,j,-1)=CF(i,j,-1);
        DC_tmp(i,j,-1)=DC(i,j,-1);	
      }
      Dphi_avg1_tmp(i,j,0)=Dphi_avg1_tmp(i,j,0);
      dxlen_tmp(i,j,0)=dxlen(i,j,0);
    });
    Gpu::streamSynchronize();
    for(int i=(unsigned long) phi_bxD.smallEnd(0); i<=phi_bxD.bigEnd(0);i++)
            for(int j=(unsigned long) phi_bxD.smallEnd(1); j<=phi_bxD.bigEnd(1);j++)
              for(int k=(unsigned long) phi_bxD.smallEnd(2); k<=phi_bxD.bigEnd(2);k++)
          {
	      /*amrex::ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {*/
        Real cff1=1.0/(CF_tmp(i,j,-1)*(1.0/dxlen_tmp(i,j,0)));
        DC_tmp(i,j,-1)=(DC_tmp(i,j,-1)*(1.0/dxlen_tmp(i,j,0))-Dphi_avg1_tmp(i,j,0))*cff1; // recursive
	//    });
	  }
#ifdef AMREX_USE_GPU
    Gpu::synchronize();
    amrex::PrintToFile("DC4_cuda")<<FArrayBox(DC)<<std::endl;
    amrex::PrintToFile("DCtmp4_cuda")<<FArrayBox(DC_tmp)<<std::endl;
#else
    amrex::PrintToFile("DC4_nocuda")<<FArrayBox(DC)<<std::endl;
    amrex::PrintToFile("DCtmp4_nocuda")<<FArrayBox(DC_tmp)<<std::endl;
#endif
    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        phi(i,j,k) -= DC(i,j,-1);
    });
}
