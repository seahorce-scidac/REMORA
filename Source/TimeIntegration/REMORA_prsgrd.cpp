#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] phi_bx   tilebox
 * @param[in   ] phi_gbx  grownbox
 * @param[in   ] utbx     u-nodal tilebox
 * @param[in   ] vtbx     v-nodal tilebox
 * @param[  out] ru       u-velocity RHS
 * @param[  out] rv       v-velocity RHS
 * @param[in   ] pm       1/dx
 * @param[in   ] pn       1/dy
 * @param[in   ] rho      density
 * @param        FC       temporary
 * @param[in   ] Hz       vertical cell height
 * @param[in   ] z_r      z coordinates at rho points
 * @param[in   ] z_w      z coordinates at w points
 * @param[in   ] msku     land-sea mask on u-points
 * @param[in   ] mskv     land-sea mask on v-points
 * @param[in   ] nrhs     index of RHS component
 * @param[in   ] N        number of vertical levels
 */
void
REMORA::prsgrd (const Box& phi_bx, const Box& phi_gbx,
               const Box& utbx, const Box& vtbx,
               const Array4<Real      >& ru,
               const Array4<Real      >& rv,
               const Array4<Real const>& pn,
               const Array4<Real const>& pm,
               const Array4<Real const>& rho,
               const Array4<Real      >& FC,
               const Array4<Real const>& Hz,
               const Array4<Real const>& z_r,
               const Array4<Real const>& z_w,
               const Array4<Real const>& msku,
               const Array4<Real const>& mskv,
               const int nrhs, const int N)
{
    BL_PROFILE("REMORA::prsgrd()");
    auto phi_bxD=phi_bx;
    phi_bxD.makeSlab(2,0);
    auto phi_gbxD=phi_gbx & phi_bx;
    phi_gbxD.makeSlab(2,0);
    Box phi_ubx = surroundingNodes(phi_bx,0);
    Box phi_vbx = surroundingNodes(phi_bx,1);
    auto utbxD = utbx;
    auto vtbxD = vtbx;
    utbxD.makeSlab(2,0);
    vtbxD.makeSlab(2,0);

    const Real OneFifth = Real(0.2);
    const Real OneTwelfth = one/Real(12.0);
    const Real eps = Real(1.0e-10);
    Real GRho     = g/solverChoice.rho0;
    Real GRho0    = Real(1000.0) * GRho;
    Real HalfGRho = half    * GRho;

    int ncomp = 0;
    int P_comp = ncomp++;
    int aux_comp = ncomp++;
    int dR_comp = ncomp++;
    int dZ_comp = ncomp++;
    int dRx_comp = ncomp++;
    int dZx_comp = ncomp++;

    FArrayBox fab(grow(phi_bx,IntVect(1,1,0)),ncomp,The_Async_Arena());

    auto P=  fab.array(P_comp);
    auto aux=fab.array(aux_comp);
    auto dR= fab.array(dR_comp);
    auto dZ= fab.array(dZ_comp);
    auto dRx=fab.array(dRx_comp);
    auto dZx=fab.array(dZx_comp);

    // Derivatives in the z direction
    ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if(k>=0&&k<N) {
            dR(i,j,k)=rho(i,j,k+1)-rho(i,j,k);
            dZ(i,j,k)=z_r(i,j,k+1)-z_r(i,j,k);
        } else {
            dR(i,j,N)=rho(i,j,N)-rho(i,j,N-1);
            dZ(i,j,N)=z_r(i,j,N)-z_r(i,j,N-1);
        }
    });

    ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        for(int k=N;k>=0;k--) {
            Real cff= k>0 ? two*dR(i,j,k)*dR(i,j,k-1) : two*dR(i,j,k)*dR(i,j,k);
            if (cff>eps) {
                dR(i,j,k)= k>0 ? cff/(dR(i,j,k)+dR(i,j,k-1)) : cff/(dR(i,j,k)+dR(i,j,k));
            } else {
                dR(i,j,k)=zero;
            }
            dZ(i,j,k)= k>0 ? two*dZ(i,j,k)*dZ(i,j,k-1)/(dZ(i,j,k)+dZ(i,j,k-1)) :
                             two*dZ(i,j,k)*dZ(i,j,k)/(dZ(i,j,k)+dZ(i,j,k));
        }
    });

    ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        Real cff1=one/(z_r(i,j,N)-z_r(i,j,N-1));
        Real cff2=half*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N+1)-z_r(i,j,N))*cff1;

        P(i,j,N)=GRho0*z_w(i,j,N+1)+GRho*(rho(i,j,N)+cff2)*(z_w(i,j,N+1)-z_r(i,j,N));

        for (int k=N-1;k>=0;k--)
        {
            Real rho_diff = rho(i,j,k+1)-rho(i,j,k) - OneTwelfth* (dR(i,j,k+1)+dR(i,j,k));
            Real   z_diff = z_r(i,j,k+1)-z_r(i,j,k) - OneTwelfth* (dZ(i,j,k+1)+dZ(i,j,k));
            Real   rz_avg = (rho(i,j,k+1)+rho(i,j,k)) * (z_r(i,j,k+1)-z_r(i,j,k));

            P(i,j,k) = P(i,j,k+1) + HalfGRho * ( rz_avg -
                          OneFifth* ( (dR(i,j,k+1)-dR(i,j,k)) *  z_diff -
                                      (dZ(i,j,k+1)-dZ(i,j,k)) * rho_diff ) );
         }
    });

    //This should be nodal
    // Derivatives in the x direction
    ParallelFor(phi_ubx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        FC(i,j,k)=(rho(i,j,k)-rho(i-1,j,k)) * msku(i,j,0);
        aux(i,j,k)=(z_r(i,j,k)-z_r(i-1,j,k)) * msku(i,j,0);
    });

    //This should be nodal aux and FC need wider boxes above
    //dZx and dRx may have index mismatch issues at k=2 and k=N
    ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        for(int k=N;k>=0;k--) {
            Real cff= two*aux(i,j,k)*aux(i+1,j,k);
            if (cff>eps) {
                Real cff1= one/(aux(i+1,j,k)+aux(i,j,k));
                dZx(i,j,k)=cff*cff1;
            } else {
                dZx(i,j,k)=zero;
            }
            Real cff1= two*FC(i,j,k)*FC(i+1,j,k);
            if (cff1>eps) {
                Real cff2= one/(FC(i,j,k)+FC(i+1,j,k));
                dRx(i,j,k)=cff1*cff2;
            } else {
                dRx(i,j,k)=zero;
            }
        }
    });

    //This should be nodal aux and FC need wider boxes above
    ParallelFor(utbxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        for(int k=N;k>=0;k--)
        {
            Real rho_diff   = rho(i,j,k)-rho(i-1,j,k)- OneTwelfth* (dRx(i,j,k)+dRx(i-1,j,k));
            Real z_r_diff   = z_r(i,j,k)-z_r(i-1,j,k)- OneTwelfth* (dZx(i,j,k)+dZx(i-1,j,k));
            Real   Hz_avg   = half * (Hz(i,j,k)+Hz(i-1,j,k));

            Real on_u = two / (pn(i-1,j,0)+pn(i,j,0));
            ru(i,j,k,nrhs) = on_u * Hz_avg * (
                            P(i-1,j,k) - P(i,j,k) - HalfGRho *
                            ( (rho(i,j,k)+rho(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k))-
                              OneFifth * ( (dRx(i,j,k)-dRx(i-1,j,k)) * z_r_diff -
                                           (dZx(i,j,k)-dZx(i-1,j,k)) * rho_diff ) ) );
        }
    });

    //This should be nodal
    ParallelFor(phi_vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        FC(i,j,k)= (rho(i,j,k)-rho(i,j-1,k)) * mskv(i,j,0);
        aux(i,j,k)= (z_r(i,j,k)-z_r(i,j-1,k)) * mskv(i,j,0);
    });

    //This should be nodal aux and FC need wider boxes above
    //dZx and dRx may have index mismatch issues at k=2 and k=N
    ParallelFor(phi_bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        for(int k=N;k>=0;k--) {
            Real cff= two*aux(i,j,k)*aux(i,j+1,k);
            if (cff>eps) {
                Real cff1= one/(aux(i,j+1,k)+aux(i,j,k));
                dZx(i,j,k)=cff*cff1;
            } else {
                dZx(i,j,k)=zero;
            }
            Real cff1= two*FC(i,j,k)*FC(i,j+1,k);
            if (cff1>eps) {
                Real cff2= one/(FC(i,j,k)+FC(i,j+1,k));
                dRx(i,j,k)=cff1*cff2;
            } else {
                dRx(i,j,k)=zero;
            }
        }
    });

    //This should be nodal aux and FC need wider boxes above
    ParallelFor(vtbxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        for (int k=N;k>=0;k--)
        {
            Real rho_diff   = rho(i,j,k)-rho(i,j-1,k)- OneTwelfth* (dRx(i,j,k)+dRx(i,j-1,k));
            Real z_r_diff   = z_r(i,j,k)-z_r(i,j-1,k)- OneTwelfth* (dZx(i,j,k)+dZx(i,j-1,k));
            Real   Hz_avg   = half * (Hz(i,j,k)+Hz(i,j-1,k));

            Real om_v = two / (pm(i,j-1,0)+pm(i,j,0));
            rv(i,j,k,nrhs) = om_v * Hz_avg * (
                            P(i,j-1,k) - P(i,j,k) - HalfGRho *
                            ( (rho(i,j,k)+rho(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))-
                              OneFifth * ( (dRx(i,j,k)-dRx(i,j-1,k)) * z_r_diff -
                                           (dZx(i,j,k)-dZx(i,j-1,k)) * rho_diff ) ) );
        } // k
    });
}
