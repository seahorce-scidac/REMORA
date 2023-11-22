#include <ROMSX.H>

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
                     const int nnew, const int N, const Real /*dt_lev*/)
{
    // Interior points vertical mean correction
    //
    // Put Hzk on the x- or y-face as appropriate, or leave on cell center for tracers
    //
    if (verbose > 1)
        amrex::Print() << "updating on box in vert_mean_3d: " << phi_bx << std::endl;
    auto phi_bxD=phi_bx;
    phi_bxD.makeSlab(2,0);

    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Hzk(i,j,k)=0.5*(Hz(i-ioff,j-joff,k)+Hz(i,j,k));
    });
    if ((verbose > 2) && (ioff==1)) {
        PrintToFile("Hzk_vertmean").SetPrecision(18) << FArrayBox(Hzk) << std::endl;
    }

    Gpu::synchronize();

    ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
      for(int k=0; k<=N; k++) {
      if(k==0) {
        CF(i,j,-1)=Hzk(i,j,k);
        DC(i,j,-1)=phi(i,j,k,nnew)*Hzk(i,j,k);
      } else {
        CF(i,j,-1)+=Hzk(i,j,k);
        DC(i,j,-1)+=phi(i,j,k,nnew)*Hzk(i,j,k);
      }
      }
    });
    Gpu::synchronize();
    ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        Real cff1=1.0/(CF(i,j,-1)*(1.0/dxlen(i,j,0)));
        DC(i,j,-1)=(DC(i,j,-1)*(1.0/dxlen(i,j,0))-Dphi_avg1(i,j,0))*cff1; // recursive
    });

    if ((verbose > 2) && (ioff==1)) {
        PrintToFile("u_vertmean").SetPrecision(18) << FArrayBox(phi) << std::endl;
        PrintToFile("DC_vertmean").SetPrecision(18) << FArrayBox(DC) << std::endl;
        PrintToFile("CF_vertmean").SetPrecision(18) << FArrayBox(CF) << std::endl;
    }
    ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        phi(i,j,k) -= DC(i,j,-1);
    });
}
