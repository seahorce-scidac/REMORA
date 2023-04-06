#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// vert_mean_3d
//
void
ROMSX::vert_mean_3d (const Box& phi_bx, const int ioff, const int joff,
                     Array4<Real> phi_arr,
                     Array4<Real> Hz_arr, Array4<Real> Hzk_arr,
                     Array4<Real> Dphi_avg1_arr,
                     Array4<Real> /*AK_arr*/, Array4<Real> /*Akv_arr*/,
                     Array4<Real> /*BC_arr*/, Array4<Real> DC_arr,
                     Array4<Real> /*FC_arr*/, Array4<Real> CF_arr,
                     Array4<Real> dxlen_arr,
                     const int nnew, const int N, const Real dt_lev)
{
    //
    // Put Hzk on the x- or y-face as appropriate, or leave on cell center for tracers
    //
    amrex::Print() << "updating on box in vert_mean_3d: " << phi_bx << std::endl;

    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Hzk_arr(i,j,k)=0.5*(Hz_arr(i-ioff,j-joff,k)+Hz_arr(i,j,k));
    });

    Real CF_tmp=0.0;
    Real DC_tmp=0.0;
    amrex::LoopOnCpu(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
      if(k==0) {
        CF_arr(i,j,-1)=Hzk_arr(i,j,k);
        DC_arr(i,j,-1)=phi_arr(i,j,k,nnew)*Hzk_arr(i,j,k);
      } else {
        CF_arr(i,j,-1)+=Hzk_arr(i,j,k);
        DC_arr(i,j,-1)+=phi_arr(i,j,k,nnew)*Hzk_arr(i,j,k);
      }
    });

    amrex::LoopOnCpu(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff1=1.0/(CF_arr(i,j,-1)*(1.0/dxlen_arr(i,j,0)));
        DC_arr(i,j,-1)=(DC_arr(i,j,-1)*(1.0/dxlen_arr(i,j,0))-Dphi_avg1_arr(i,j,0))*cff1; // recursive
    });

    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        phi_arr(i,j,k) -= DC_arr(i,j,-1);
    });
}
