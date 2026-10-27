#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] phi_bx     box to update on
 * @param[in   ] ioff       x-direction offset
 * @param[in   ] joff       y-direction offset
 * @param[inout] phi        velocity (u or v)
 * @param[in   ] Dphi_avg1  time average of barotropic velocity
 * @param        DC         temporary
 * @param        CF         temporary
 * @param[in   ] on_u_or_om_v  on_u for u, om_v for v
 * @param[in   ] msk        land-sea mask
 * @param[in   ] nnew       index of time step to update
 * @param[in   ] N          number of vertical levels
 */
void
REMORA::vert_mean_3d (const Box& phi_bx, const int ioff, const int joff,
                     const Array4<Real      >& phi,
                     const Array4<Real const>& Hz,
                     const Array4<Real const>& Dphi_avg1,
                     const Array4<Real      >& DC,
                     const Array4<Real      >& CF,
                     const Array4<Real const>& on_u_or_om_v,
                     const Array4<Real const>& msk,
                     const int nnew, const int N)
{
    BL_PROFILE("REMORA::vert_mean_3d()");

    ParallelFor(makeSlab(phi_bx,2,0),
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        Real Hzk_on_face = Real(0.5)*(Hz(i-ioff,j-joff,0)+Hz(i,j,0));
        CF(i,j,-1) =                 Hzk_on_face;
        DC(i,j,-1) = phi(i,j,0,nnew)*Hzk_on_face;

        for (int k=1; k<=N; k++) {
            Hzk_on_face = Real(0.5)*(Hz(i-ioff,j-joff,k)+Hz(i,j,k));
            CF(i,j,-1) +=                 Hzk_on_face;
            DC(i,j,-1) += phi(i,j,k,nnew)*Hzk_on_face;
        }
    });

    ParallelFor(makeSlab(phi_bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        Real cff1=one/(CF(i,j,-1)*(on_u_or_om_v(i,j,0)));
        DC(i,j,-1) = (DC(i,j,-1)*(on_u_or_om_v(i,j,0)) - Dphi_avg1(i,j,0))*cff1; // recursive
    });

    ParallelFor(phi_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        phi(i,j,k) -= DC(i,j,-1);
        phi(i,j,k) *= msk(i,j,0);
    });
}
