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
                Array4<Real const> temp,
                Array4<Real const> salt,
                Array4<Real      > rho,
                Array4<Real      > rhoA,
                Array4<Real      > rhoS,
                Array4<Real const> Hz,
                Array4<Real const> z_w,
                Array4<Real const> h,
                const int nrhs,
                const int N)
{
    // Hardcode these for now instead of reading them from inputs
    Real T0=14.0;
    Real S0=35.0;
    Real R0=1027;
    Real Tcoef=1.7e-4;
    Real Scoef=0.0;
    Real rho0=1025.0;

    AMREX_ASSERT(bx.smallEnd(2) == 0 && bx.bigEnd(2) == N);
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
        rho(i,j,k)  = R0 - R0*Tcoef*(temp(i,j,k,nrhs)-T0);
        rho(i,j,k) += R0*Scoef*(salt(i,j,k,nrhs)-S0) - 1000.0_rt;
    });

//
//-----------------------------------------------------------------------
//  Compute vertical averaged density (rhoA) and density perturbation
//  used (rhoS) in barotropic pressure gradient.
//-----------------------------------------------------------------------
//
    Real cff2 =1.0_rt/rho0;

    amrex::ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        Real cff0 = rho(i,j,N)*Hz(i,j,N);
        rhoS(i,j,0) = 0.5_rt*cff0*Hz(i,j,N);
        rhoA(i,j,0) = cff0;

        for (int k = 1; k <= N; ++k) {
            Real cff1=rho(i,j,N-k)*Hz(i,j,N-k);
            rhoS(i,j,0) += Hz(i,j,N-k)*(rhoA(i,j,0)+0.5_rt*cff1);
            rhoA(i,j,0) += cff1;
        }

        Real cff11 =1.0_rt/(z_w(i,j,N)+h(i,j,0,0));

        rhoA(i,j,0) *= cff2*cff11;

        rhoS(i,j,0) *= 2.0_rt*cff11*cff11*cff2;
    });
}
