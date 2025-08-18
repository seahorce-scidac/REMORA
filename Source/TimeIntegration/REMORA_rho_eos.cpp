#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] bx       box for calculation
 * @param[in   ] state    state holds temp, salt
 * @param[  out] rho      density
 * @param[  out] rhoA     vertically-averaged density
 * @param[  out] rhoS     density perturbation
 * @param[  out] bvf      Brunt-Vaisala frequency
 * @param[  out] alpha    thermal expansion coefficient
 * @param[  out] beta     saline contraction coefficient
 * @param[in   ] Hz       vertical cell height
 * @param[in   ] z_r      z coordinates at rho points
 * @param[in   ] z_w      z coordinates at w points
 * @param[in   ] h        bathymetry
 * @param[in   ] mskr     land-sea mask on rho points
 * @param[in   ] N        number of vertical levels
 */
void
REMORA::rho_eos (const Box& bx,
                const Array4<Real const>& state,
                const Array4<Real      >& rho,
                const Array4<Real      >& rhoA,
                const Array4<Real      >& rhoS,
                const Array4<Real      >& bvf,
                const Array4<Real      >& alpha,
                const Array4<Real      >& beta,
                const Array4<Real const>& Hz,
                const Array4<Real const>& z_w,
                const Array4<Real const>& z_r,
                const Array4<Real const>& h,
                const Array4<Real const>& mskr,
                const int N)
{
    if (solverChoice.eos_type == EOSType::linear) {
        lin_eos(bx, state, rho, rhoA, rhoS, bvf, Hz, z_w, z_r, h, mskr, N);
    } else if (solverChoice.eos_type == EOSType::nonlinear) {
        nonlin_eos(bx, state, rho, rhoA, rhoS, bvf, alpha, beta, Hz, z_w, z_r, h, mskr, N);
    } else {
        Abort("Unknown EOS type in rho_eos");
    }
}

/**
 * @param[in   ] bx       box for calculation
 * @param[in   ] state    state holds temp, salt
 * @param[  out] rho      density
 * @param[  out] rhoA     vertically-averaged density
 * @param[  out] rhoS     density perturbation
 * @param[  out] bvf      Brunt-Vaisala frequency
 * @param[in   ] Hz       vertical cell height
 * @param[in   ] z_r      z coordinates at rho points
 * @param[in   ] z_w      z coordinates at w points
 * @param[in   ] h        bathymetry
 * @param[in   ] mskr     land-sea mask on rho points
 * @param[in   ] N        number of vertical levels
 */
void
REMORA::lin_eos (const Box& bx,
                const Array4<Real const>& state,
                const Array4<Real      >& rho,
                const Array4<Real      >& rhoA,
                const Array4<Real      >& rhoS,
                const Array4<Real      >& bvf,
                const Array4<Real const>& Hz,
                const Array4<Real const>& z_w,
                const Array4<Real const>& z_r,
                const Array4<Real const>& h,
                const Array4<Real const>& mskr,
                const int N)
{
    BL_PROFILE("REMORA::lin_eos()");
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
        rho(i,j,k)  = (R0 - R0*Tcoef*(state(i,j,k,Temp_comp)-T0)
                         + R0*Scoef*(state(i,j,k,Salt_comp)-S0)
                         - 1000.0_rt) * mskr(i,j,0);
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

        Real cff11 =1.0_rt/(z_w(i,j,N+1)+h(i,j,0,0));

        rhoA(i,j,0) *= cff2*cff11;

        rhoS(i,j,0) *= 2.0_rt*cff11*cff11*cff2;
    });

    // Compute Brunt-Vaisala frequency (1/s2)
    Real gorho0 = g / solverChoice.rho0;
    // Really want enclosed nodes or something similar
    Box box_w = bx;
    box_w.surroundingNodes(2);
    box_w.grow(IntVect(0,0,-1));

    ParallelFor(box_w, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        bvf(i,j,k) = -gorho0 * (rho(i,j,k)-rho(i,j,k-1)) / (z_r(i,j,k)-z_r(i,j,k-1));
    });
}

/**
 * @param[in   ] bx       box for calculation
 * @param[in   ] state    state holds temp, salt
 * @param[  out] rho      density
 * @param[  out] rhoA     vertically-averaged density
 * @param[  out] rhoS     density perturbation
 * @param[  out] bvf      Brunt-Vaisala frequency
 * @param[  out] alpha    thermal expansion coefficient
 * @param[  out] beta     saline contraction coefficient
 * @param[in   ] Hz       vertical cell height
 * @param[in   ] z_r      z coordinates at rho points
 * @param[in   ] z_w      z coordinates at w points
 * @param[in   ] h        bathymetry
 * @param[in   ] mskr     land-sea mask on rho points
 * @param[in   ] N        number of vertical levels
 */
void
REMORA::nonlin_eos (const Box& bx,
                const Array4<Real const>& state,
                const Array4<Real      >& rho,
                const Array4<Real      >& rhoA,
                const Array4<Real      >& rhoS,
                const Array4<Real      >& bvf,
                const Array4<Real      >& alpha,
                const Array4<Real      >& beta,
                const Array4<Real const>& Hz,
                const Array4<Real const>& z_w,
                const Array4<Real const>& z_r,
                const Array4<Real const>& h,
                const Array4<Real const>& mskr,
                const int N)
{
    BL_PROFILE("REMORA::nonlin_eos()");
//
    AMREX_ASSERT(bx.smallEnd(2) == 0 && bx.bigEnd(2) == N);

    FArrayBox fab_den1 (bx,1,amrex::The_Async_Arena()); auto den1  = fab_den1.array();
    FArrayBox fab_den  (bx,1,amrex::The_Async_Arena()); auto den   = fab_den.array();
    FArrayBox fab_bulk (bx,1,amrex::The_Async_Arena()); auto bulk  = fab_bulk.array();
    FArrayBox fab_bulk0(bx,1,amrex::The_Async_Arena()); auto bulk0 = fab_bulk0.array();
    FArrayBox fab_bulk1(bx,1,amrex::The_Async_Arena()); auto bulk1 = fab_bulk1.array();
    FArrayBox fab_bulk2(bx,1,amrex::The_Async_Arena()); auto bulk2 = fab_bulk2.array();
    FArrayBox fab_DbulkDT(bx,1,amrex::The_Async_Arena()); auto DbulkDT = fab_DbulkDT.array();
    FArrayBox fab_Dden1DT(bx,1,amrex::The_Async_Arena()); auto Dden1DT = fab_Dden1DT.array();
    FArrayBox fab_DbulkDS(bx,1,amrex::The_Async_Arena()); auto DbulkDS = fab_DbulkDS.array();
    FArrayBox fab_Dden1DS(bx,1,amrex::The_Async_Arena()); auto Dden1DS = fab_Dden1DS.array();

    bool bulk_fluxes = solverChoice.bulk_fluxes;
//
//=======================================================================
//  Non-linear equation of state.
//=======================================================================
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real Tt = std::max(-2.5_rt, state(i,j,k,Temp_comp));
        Tt = std::min(40.0_rt, Tt);
        Real Ts = std::max(0.0_rt, state(i,j,k,Salt_comp));
        Ts = std::min(100.0_rt, Ts);
        Real sqrtTs = std::sqrt(Ts);

        Real Tp = z_r(i,j,k);
        Real Tpr10 = 0.1_rt*Tp;

        Real C0 = Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))));
        Real C1 = U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)));
        Real C2 = V00+Tt*(V01+Tt*V02);

        if (bulk_fluxes) {
            Real dCdT0=Q01+Tt*(2.0_rt*Q02+Tt*(3.0_rt*Q03+Tt*(4.0_rt*Q04+
                         Tt*5.0_rt*Q05)));
            Real dCdT1=U01+Tt*(2.0_rt*U02+Tt*(3.0_rt*U03+Tt*4.0_rt*U04));
            Real dCdT2=V01+Tt*2.0_rt*V02;
            //  Compute d(den1)/d(S) and d(den1)/d(T) derivatives used in the
            //  computation of thermal expansion and saline contraction
            //  coefficients.

            Dden1DS(i,j,k)=C1+1.5_rt*C2*sqrtTs+2.0_rt*W00*Ts;
            Dden1DT(i,j,k)=dCdT0+Ts*(dCdT1+sqrtTs*dCdT2);
        }

        den1(i,j,k) = C0 + Ts*(C1+sqrtTs*C2+Ts*W00);

        //-----------------------------------------------------------------------
        //  Compute secant bulk modulus.
        //-----------------------------------------------------------------------
        Real C3=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)));
        Real C4=B00+Tt*(B01+Tt*(B02+Tt*B03));
        Real C5=D00+Tt*(D01+Tt*D02);
        Real C6=E00+Tt*(E01+Tt*(E02+Tt*E03));
        Real C7=F00+Tt*(F01+Tt*F02);
        Real C8=G01+Tt*(G02+Tt*G03);
        Real C9=H00+Tt*(H01+Tt*H02);

        bulk0(i,j,k)=C3+Ts*(C4+sqrtTs*C5);
        bulk1(i,j,k)=C6+Ts*(C7+sqrtTs*G00);
        bulk2(i,j,k)=C8+Ts*C9;
        bulk (i,j,k)=bulk0(i,j,k)-Tp*(bulk1(i,j,k)-Tp*bulk2(i,j,k));

        Real cff = 1.0_rt / (bulk(i,j,k) + Tpr10);
        den(i,j,k) = (den1(i,j,k)*bulk(i,j,k)*cff - 1000.0_rt) * mskr(i,j,0);
        // This line may need to move once bluk fluxes are added
        rho(i,j,k) = den(i,j,k);

        if (bulk_fluxes) {
            Real dCdT3=A01+Tt*(2.0_rt*A02+Tt*(3.0_rt*A03+Tt*4.0_rt*A04));
            Real dCdT4=B01+Tt*(2.0_rt*B02+Tt*3.0_rt*B03);
            Real dCdT5=D01+Tt*2.0_rt*D02;
            Real dCdT6=E01+Tt*(2.0_rt*E02+Tt*3.0_rt*E03);
            Real dCdT7=F01+Tt*2.0_rt*F02;
            Real dCdT8=G02+Tt*2.0_rt*G03;
            Real dCdT9=H01+Tt*2.0_rt*H02;

            DbulkDS(i,j,k)=C4+sqrtTs*1.5_rt*C5-
                        Tp*(C7+sqrtTs*1.5_rt*G00-Tp*C9);
            DbulkDT(i,j,k)=dCdT3+Ts*(dCdT4+sqrtTs*dCdT5)-
                        Tp*(dCdT6+Ts*dCdT7-
                            Tp*(dCdT8+Ts*dCdT9));
        }
    });

//
//-----------------------------------------------------------------------
//  Compute vertical averaged density (rhoA) and density perturbation
//  used (rhoS) in barotropic pressure gradient.
//-----------------------------------------------------------------------
//
    Real cff_rho =1.0_rt/solverChoice.rho0;

    ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        Real cff0 = den(i,j,N)*Hz(i,j,N);
        rhoS(i,j,0) = 0.5_rt*cff0*Hz(i,j,N);
        rhoA(i,j,0) = cff0;

        for (int k = 1; k <= N; ++k) {
            Real cff1=den(i,j,N-k)*Hz(i,j,N-k);
            rhoS(i,j,0) += Hz(i,j,N-k)*(rhoA(i,j,0)+0.5_rt*cff1);
            rhoA(i,j,0) += cff1;
        }
        Real cff11 =1.0_rt/(z_w(i,j,N+1)+h(i,j,0,0));
        rhoA(i,j,0) *= cff_rho*cff11;
        rhoS(i,j,0) *= 2.0_rt*cff11*cff11*cff_rho;
    });

    // Compute Brunt-Vaisala frequency (1/s2)
    Box bxD = bx; bxD.makeSlab(2,0);

    ParallelFor(bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        bvf(i,j,0) = 0.0_rt;
        bvf(i,j,N+1) = 0.0_rt;
        for (int k=0; k<=N-1; k++) {
            Real bulk_up = bulk0(i,j,k+1) - (z_w(i,j,k+1) * (bulk1(i,j,k+1) - bulk2(i,j,k+1)*z_w(i,j,k+1)));
            Real bulk_dn = bulk0(i,j,k  ) - (z_w(i,j,k+1) * (bulk1(i,j,k  ) - bulk2(i,j,k  )*z_w(i,j,k+1)));
            Real cff1 = 1.0_rt / (bulk_up + 0.1_rt * z_w(i,j,k+1));
            Real cff2 = 1.0_rt / (bulk_dn + 0.1_rt * z_w(i,j,k+1));
            Real den_up = cff1 * (den1(i,j,k+1) * bulk_up);
            Real den_dn = cff2 * (den1(i,j,k  ) * bulk_dn);
            bvf(i,j,k+1) = -g * (den_up - den_dn) / (0.5_rt * (den_up+den_dn) * (z_r(i,j,k+1) - z_r(i,j,k)));
        }
    });

    if (solverChoice.bulk_fluxes) {
        ParallelFor(bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real Tp = z_r(i,j,N);
            Real Tpr10 = 0.1_rt*Tp;

            // Compute thermal expansion and saline contraction coefficients
            Real cff = bulk(i,j,N) + Tpr10;
            Real cff1 = Tpr10 * den1(i,j,N);
            Real cff2 = bulk(i,j,N) * cff;
            Real wrk = (den(i,j,N) + 1000.0_rt) * cff * cff;
            Real Tcof = -(DbulkDT(i,j,N) * cff1 + Dden1DT(i,j,N) * cff2);
            Real Scof = (DbulkDS(i,j,N) * cff1 + Dden1DS(i,j,N) * cff2);
            alpha(i,j,0) = Tcof / wrk;
            beta(i,j,0)  = Scof / wrk;
        });
    }
}
