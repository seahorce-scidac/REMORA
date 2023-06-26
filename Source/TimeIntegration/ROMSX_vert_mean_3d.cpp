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
    Gpu::synchronize();
    Real CF_tmp=0.0;
    Real DC_tmp=0.0;
    amrex::LoopOnCpu(phi_bx,
    [=] (int i, int j, int k)
    {
      if(k==0) {
        CF(i,j,-1)=Hzk(i,j,k);
        DC(i,j,-1)=phi(i,j,k,nnew)*Hzk(i,j,k);
      } else {
        CF(i,j,-1)+=Hzk(i,j,k);
        DC(i,j,-1)+=phi(i,j,k,nnew)*Hzk(i,j,k);
      }
    });
    Gpu::synchronize();
    amrex::LoopOnCpu(phi_bxD,
    [=] (int i, int j, int k)
    {
        Real cff1=1.0/(CF(i,j,-1)*(1.0/dxlen(i,j,0)));
        DC(i,j,-1)=(DC(i,j,-1)*(1.0/dxlen(i,j,0))-Dphi_avg1(i,j,0))*cff1; // recursive
    });

    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        phi(i,j,k) -= DC(i,j,-1);
    });
}
