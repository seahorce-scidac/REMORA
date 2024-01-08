#include <REMORA.H>

using namespace amrex;

/**
 * rho_eos
 *
 * @param[in ] bx    box for calculation
 * @param[in ] state state holds temp, salt
 * @param[out] rho   density
 * @param[out] rhoA  vertically-averaged density
 * @param[out] rhoS  density perturbation
 * @param[in ] Hz
 * @param[in ] z_w
 * @param[in ] h
 * @param[in ] N
 */

void
REMORA::rho_eos (const Box& bx,
                const Array4<Real const>& state,
                const Array4<Real      >& rho,
                const Array4<Real      >& rhoA,
                const Array4<Real      >& rhoS,
                const Array4<Real const>& Hz,
                const Array4<Real const>& z_w,
                const Array4<Real const>& h,
                const int N)
{
//
    AMREX_ASSERT(bx.smallEnd(2) == 0 && bx.bigEnd(2) == N);
//
//=======================================================================
//  Linear equation of state.
//=======================================================================
//
//
//-----------------------------------------------------------------------
//  Compute "in situ" density anomaly (kg/m3 - 1000) using the linear
//  equation of state.
//-----------------------------------------------------------------------
//
    Real R0 = solverChoice.R0;
    Real S0 = solverChoice.S0;
    Real T0 = solverChoice.T0;
    Real Tcoef = solverChoice.Tcoef;
    Real Scoef = solverChoice.Scoef;

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        rho(i,j,k)  = R0 - R0*Tcoef*(state(i,j,k,Temp_comp)-T0)
                         + R0*Scoef*(state(i,j,k,Salt_comp)-S0)
                         - 1000.0_rt;
    });

//
//-----------------------------------------------------------------------
//  Compute vertical averaged density (rhoA) and density perturbation
//  used (rhoS) in barotropic pressure gradient.
//-----------------------------------------------------------------------
//
    Real cff2 =1.0_rt/solverChoice.rho0;

    ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
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
