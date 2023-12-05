#include <ROMSX.H>

using namespace amrex;

/**
 * rho_eos
 *
 * @param[in ] bx    box for calculation
 * @param[in ] temp  temperature
 * @param[out] rho   density
 * @param[out] rhoA  vertically-averaged density
 * @param[out] rhoS  density perturbation
 * @param[in ] Hz
 * @param[in ] z_w
 * @param[in ] h
 * @param[in ] nrhs
 * @param[in ] N
 */

void
ROMSX::rho_eos (const Box& bx,
                Array4<Real> temp , Array4<Real> salt,
                Array4<Real> rho,
                Array4<Real> rhoA,
                Array4<Real> rhoS,
                Array4<Real> Hz,
                Array4<Real> z_w,
                Array4<Real> h,
                const int nrhs,
                const int N)
{
    auto bxD=bx;
    bxD.makeSlab(2,0);
    //hardcode these for now instead of reading them from inputs
    Real T0=14.0;
    Real S0=35.0;
    Real R0=1027;
    Real Tcoef=1.7e-4;
    Real Scoef=0.0;
    Real rho0=1025.0;

//
//=======================================================================
//  Linear equation of state.
//=======================================================================
//
//-----------------------------------------------------------------------
//  Compute "in situ" density anomaly (kg/m3 - 1000) using the linear
//  equation of state.
//-----------------------------------------------------------------------
//
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        rho(i,j,k)=R0-
            R0*Tcoef*(temp(i,j,k,nrhs)-T0);
        rho(i,j,k)=rho(i,j,k)+
            R0*Scoef*(salt(i,j,k,nrhs)-S0);
        rho(i,j,k)=rho(i,j,k)-1000.0_rt;
//        pden(i,j,k)=rho(i,j,k);
    });

//
//-----------------------------------------------------------------------
//  Compute vertical averaged density (rhoA) and density perturbation
//  used (rhoS) in barotropic pressure gradient.
//-----------------------------------------------------------------------
//
    amrex::ParallelFor(bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        //printf("%d %d  %15.15g Hzstart\n", i,j, Hz(i,j,N));
        Real cff1=rho(i,j,N)*Hz(i,j,N);
        rhoS(i,j,0)=0.5_rt*cff1*Hz(i,j,N);
        rhoA(i,j,0)=cff1;
        //printf("%d %d %d   %15.15g cff1\n", i,j,N, cff1);
        //printf("%d %d %d   %15.15g rhoN\n", i,j,N, rho(i,j,N));
        //printf("%d %d %d   %15.15g rhoA\n", i,j,N, rhoA(i,j,0));
        //printf("%d %d %d   %15.15g rhoS\n", i,j,N, rhoS(i,j,0));
        //printf("%d %d %d   %15.15g Hzend  \n", i,j,N, Hz(i,j,N));
        //printf("%d %d  %15.15g Hzend\n", i,j, Hz(i,j,N));
        //printf("%d %d  %15.15g %15.15g %15.15g %15.15g %15.15g cff1 rhoN rhoA rhoS Hz rhoeos\n",
        //        i,j, cff1, rho(i,j,N), rhoA(i,j,0), rhoS(i,j,0), Hz(i,j,N));
    });
    AMREX_ASSERT(bx.smallEnd(2) == 0 &&
                 bx.bigEnd(2) == N);
    amrex::ParallelFor(bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
      for (int k = 1; k <= N; ++k) {
        Real cff1=rho(i,j,N-k)*Hz(i,j,N-k);
        rhoS(i,j,0)=rhoS(i,j,0)+Hz(i,j,N-k)*(rhoA(i,j,0)+0.5_rt*cff1);
        rhoA(i,j,0)=rhoA(i,j,0)+cff1;
        //printf("%d %d %d  %15.15g %15.15g %15.15g %15.15g   cff1 rhoA rhoS Hz rhoeos\n", i,j,k, cff1, rhoA(i,j,0), rhoS(i,j,0), Hz(i,j,N-k));
      }
    });
    amrex::ParallelFor(bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        Real cff2=1.0_rt/rho0;
        Real cff1=1.0_rt/(z_w(i,j,N)+h(i,j,0,0));
        rhoA(i,j,0)=cff2*cff1*rhoA(i,j,0);
        rhoS(i,j,0)=2.0_rt*cff1*cff1*cff2*rhoS(i,j,0);
        //printf("%d %d  %15.15g %15.15g %15.15g %15.15g z_wN z_w0\n", i,j, z_w(i,j,N), z_w(i,j,0), h(i,j,0,0),h(i,j,0,1));
        //printf("%d %d  %15.15g %15.15g %15.15g %15.15g  cff2 cff1 rhoA rhoS rhoeos\n", i,j,cff2, cff1, rhoA(i,j,0), rhoS(i,j,0));
    });

}
