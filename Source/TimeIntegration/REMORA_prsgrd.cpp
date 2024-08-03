#include <REMORA.H>

using namespace amrex;

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

    const Real OneFifth = 0.2_rt;
    const Real OneTwelfth = 1.0_rt/12.0_rt;
    const Real eps = 1.0E-10_rt;
    Real GRho     = solverChoice.g/solverChoice.rho0;
    Real GRho0    = 1000.0_rt * GRho;
    Real HalfGRho = 0.5_rt    * GRho;

    FArrayBox fab_P(phi_bx,1,The_Async_Arena());
    FArrayBox fab_aux(Box(z_r),1,The_Async_Arena());
    FArrayBox fab_dR(phi_bx,1,The_Async_Arena());
    FArrayBox fab_dZ(phi_bx,1,The_Async_Arena());
    FArrayBox fab_dRx(phi_bx,1,The_Async_Arena());
    FArrayBox fab_dZx(phi_bx,1,The_Async_Arena());

    auto P=fab_P.array();
    auto aux=fab_aux.array();
    auto dR=fab_dR.array();
    auto dZ=fab_dZ.array();
    auto dRx=fab_dRx.array();
    auto dZx=fab_dZx.array();

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
            Real cff= k>0 ? 2.0_rt*dR(i,j,k)*dR(i,j,k-1) : 2.0_rt*dR(i,j,k)*dR(i,j,k);
            if (cff>eps) {
                dR(i,j,k)= k>0 ? cff/(dR(i,j,k)+dR(i,j,k-1)) : cff/(dR(i,j,k)+dR(i,j,k));
            } else {
                dR(i,j,k)=0.0_rt;
            }
            dZ(i,j,k)= k>0 ? 2.0_rt*dZ(i,j,k)*dZ(i,j,k-1)/(dZ(i,j,k)+dZ(i,j,k-1)) :
                             2.0_rt*dZ(i,j,k)*dZ(i,j,k)/(dZ(i,j,k)+dZ(i,j,k));
        }
    });

    ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        Real cff1=1.0_rt/(z_r(i,j,N)-z_r(i,j,N-1));
        Real cff2=0.5_rt*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N+1)-z_r(i,j,N))*cff1;

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
            Real cff= 2.0_rt*aux(i,j,k)*aux(i+1,j,k);
            if (cff>eps) {
                Real cff1= 1.0_rt/(aux(i+1,j,k)+aux(i,j,k));
                dZx(i,j,k)=cff*cff1;
            } else {
                dZx(i,j,k)=0.0_rt;
            }
            Real cff1= 2.0_rt*FC(i,j,k)*FC(i+1,j,k);
            if (cff1>eps) {
                Real cff2= 1.0_rt/(FC(i,j,k)+FC(i+1,j,k));
                dRx(i,j,k)=cff1*cff2;
            } else {
                dRx(i,j,k)=0.0_rt;
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
            Real   Hz_avg   = 0.5_rt * (Hz(i,j,k)+Hz(i-1,j,k));

            Real on_u = 2.0_rt / (pn(i-1,j,0)+pn(i,j,0));
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
            Real cff= 2.0_rt*aux(i,j,k)*aux(i,j+1,k);
            if (cff>eps) {
                Real cff1= 1.0_rt/(aux(i,j+1,k)+aux(i,j,k));
                dZx(i,j,k)=cff*cff1;
            } else {
                dZx(i,j,k)=0.0_rt;
            }
            Real cff1= 2.0_rt*FC(i,j,k)*FC(i,j+1,k);
            if (cff1>eps) {
                Real cff2= 1.0_rt/(FC(i,j,k)+FC(i,j+1,k));
                dRx(i,j,k)=cff1*cff2;
            } else {
                dRx(i,j,k)=0.0_rt;
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
            Real   Hz_avg   = 0.5_rt * (Hz(i,j,k)+Hz(i,j-1,k));

            Real om_v = 2.0_rt / (pm(i,j-1,0)+pm(i,j,0));
            rv(i,j,k,nrhs) = om_v * Hz_avg * (
                            P(i,j-1,k) - P(i,j,k) - HalfGRho *
                            ( (rho(i,j,k)+rho(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))-
                              OneFifth * ( (dRx(i,j,k)-dRx(i,j-1,k)) * z_r_diff -
                                           (dZx(i,j,k)-dZx(i,j-1,k)) * rho_diff ) ) );
        } // k
    });
}
