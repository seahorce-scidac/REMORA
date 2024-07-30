#include <REMORA.H>

using namespace amrex;

//
// vert_mean_3d
//

void
REMORA::vert_mean_3d (const Box& phi_bx, const int ioff, const int joff,
                     const Array4<Real      >& phi,
                     const Array4<Real const>& Hz,
                     const Array4<Real const>& Dphi_avg1,
                     const Array4<Real      >& DC,
                     const Array4<Real      >& CF,
                     const Array4<Real const>& dxlen,
                     const Array4<Real const>& msk,
                     const int nnew, const int N)
{

    ParallelFor(makeSlab(phi_bx,2,0),
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        Real Hzk_on_face = 0.5_rt*(Hz(i-ioff,j-joff,0)+Hz(i,j,0));
        CF(i,j,-1) =                 Hzk_on_face;
        DC(i,j,-1) = phi(i,j,0,nnew)*Hzk_on_face;

        for (int k=1; k<=N; k++) {
            Hzk_on_face = 0.5_rt*(Hz(i-ioff,j-joff,k)+Hz(i,j,k));
            CF(i,j,-1) +=                 Hzk_on_face;
            DC(i,j,-1) += phi(i,j,k,nnew)*Hzk_on_face;
        }
    });

    ParallelFor(makeSlab(phi_bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        Real on_u_or_om_v = 2.0_rt / (dxlen(i-ioff,j-joff,0) + dxlen(i,j,0));
        Real cff1=1.0_rt/(CF(i,j,-1)*(on_u_or_om_v));
        DC(i,j,-1) = (DC(i,j,-1)*(on_u_or_om_v) - Dphi_avg1(i,j,0))*cff1; // recursive
    });

    ParallelFor(phi_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        phi(i,j,k) -= DC(i,j,-1);
        phi(i,j,k) *= msk(i,j,0);
    });
}
