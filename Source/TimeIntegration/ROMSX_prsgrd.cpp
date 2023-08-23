#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

void
ROMSX::prsgrd (const Box& phi_bx,
               Array4<Real> ru , Array4<Real> rv,
               Array4<Real> rho,
               Array4<Real> FC,
               Array4<Real> Hz,
               Array4<Real> z_r,
               Array4<Real> z_w,
               const int nrhs,
               const int N)
{
    const int Mn = Geom(0).Domain().size()[0];
    const int Mm = Geom(0).Domain().size()[1];
    auto phi_bxD=phi_bx;
    phi_bxD.makeSlab(2,0);
    //hardcode these for now instead of reading them from inputs
    Real T0=14.0;
    Real S0=35.0;
    Real R0=1027;
    Real Tcoef=1.7e-4;
    Real Scoef=0.0;
    Real rho0=1025.0;

    const Real OneFifth = 0.2_rt;
    const Real OneTwelfth = 1.0_rt/12.0_rt;
    const Real eps = 1.0E-10_rt;
    Real g=9.81;
    Real GRho=g/rho0;
    Real GRho0=1000.0_rt*GRho;
    Real HalfGRho=0.5_rt*GRho;

    FArrayBox fab_zwrk(phi_bx,1,The_Async_Arena());
    FArrayBox fab_P(phi_bx,1,The_Async_Arena());
    FArrayBox fab_aux(phi_bx,1,The_Async_Arena());
    FArrayBox fab_dR(phi_bx,1,The_Async_Arena());
    FArrayBox fab_dZ(phi_bx,1,The_Async_Arena());
    FArrayBox fab_dRx(phi_bx,1,The_Async_Arena());
    FArrayBox fab_dZx(phi_bx,1,The_Async_Arena());

    auto zwrk=fab_zwrk.array();
    auto P=fab_P.array();
    auto aux=fab_aux.array();
    auto dR=fab_dR.array();
    auto dZ=fab_dZ.array();
    auto dRx=fab_dRx.array();
    auto dZx=fab_dZx.array();

    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if(k>=0||k<N) {
            dR(i,j,k)=rho(i,j,k+1)-rho(i,j,k);
            dZ(i,j,k)=z_r(i,j,k+1)-z_r(i,j,k);
        } else {
            dR(i,j,N)=dR(i,j,N-1);
            dZ(i,j,N)=dZ(i,j,N-1);
            //This is really k=-1
            //dR(i,j,0)=dR(i,j,1);
            //dZ(i,j,0)=dZ(i,j,1);
        }
    });

    amrex::ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        for(int k=N;k>=0;k--) {
            Real cff= k>0 ? 2.0*dR(i,j,k)*dR(i,j,k-1) : 2.0*dR(i,j,k)*dR(i,j,k);
            if (cff>eps) {
                dR(i,j,k)= k>0 ? cff/(dR(i,j,k)+dR(i,j,k-1)) : cff/(dR(i,j,k)+dR(i,j,k));
            } else {
                dR(i,j,k)=0.0;
            }
            dZ(i,j,k)= k>0 ? 2.0_rt*dZ(i,j,k)*dZ(i,j,k-1)/(dZ(i,j,k)+dZ(i,j,k-1)) :
                             2.0_rt*dZ(i,j,k)*dZ(i,j,k)/(dZ(i,j,k)+dZ(i,j,k));
        }
        Real cff1=1.0_rt/(z_r(i,j,N)-z_r(i,j,N-1));
        Real cff2=0.5_rt*(rho(i,j,N)-rho(i,j,N-1))*
            (z_w(i,j,N)-z_r(i,j,N))*cff1;
        P(i,j,N)=GRho0*z_w(i,j,N)+GRho*(rho(i,j,N)+cff2)*
                                       (z_w(i,j,N)-z_r(i,j,N));
        for(int k=N-1;k>=0;k--) {
            P(i,j,k)=P(i,j,k+1)+
                     HalfGRho*((rho(i,j,k+1)+rho(i,j,k))*
                               (z_r(i,j,k+1)-z_r(i,j,k))-
                               OneFifth*
                               ((dR(i,j,k+1)-dR(i,j,k))*
                                (z_r(i,j,k+1)-z_r(i,j,k)-
                                 OneTwelfth*
                                 (dZ(i,j,k+1)+dZ(i,j,k)))-
                                (dZ(i,j,k+1)-dZ(i,j,k))*
                                (rho(i,j,k+1)-rho(i,j,k)-
                                 OneTwelfth*
                                 (dR(i,j,k+1)+dR(i,j,k)))));
         }
    });
}
