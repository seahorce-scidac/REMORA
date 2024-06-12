#include <REMORA.H>

using namespace amrex;

//
// vert_visc_3d
//

void
REMORA::vert_visc_3d (const Box& phi_bx, const int ioff, const int joff,
                     const Array4<Real      >& phi,
                     const Array4<Real const>& Hz,
                     const Array4<Real      >& Hzk, /*temp var */
                     const Array4<Real      >& AK,
                     const Array4<Real      >& Akv,
                     const Array4<Real      >& BC,
                     const Array4<Real      >& DC,
                     const Array4<Real      >& FC,
                     const Array4<Real      >& CF,
                     const int nnew, const int N, const Real dt_lev)
{
    //
    // Put Hzk on the x- or y-face as appropriate, or leave on cell center for tracers
    //
    if (verbose > 1)
        amrex::Print() << "updating on box in vert_visc_3d: " << phi_bx << std::endl;

    ParallelFor(phi_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Hzk(i,j,k)=0.5_rt*(Hz(i-ioff,j-joff,k)+Hz(i,j,k));
    });

    //
    // Put Akv on the x- or y-face as appropriate, or leave on cell center for tracers
    //
    ParallelFor(surroundingNodes(phi_bx,2), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        AK(i,j,k) = 0.5_rt * (Akv(i-ioff,j-joff,k)+Akv(i,j,k));
    });

    Gpu::streamSynchronize();
#ifdef AMREX_USE_GPU
    Gpu::synchronize();
#else
#endif
    /////////////////// This and the following loop is the first non-matching thing that affects plotfile comparison for cuda
    // NOTE: vertical viscosity term for tracers is identical except AK=Akt
    const Real sixth = 1.0_rt / 6.0_rt;
    const Real third = 1.0_rt / 3.0_rt;

    ParallelFor(makeSlab(phi_bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        DC(i,j,0) = 0.0_rt;
        FC(i,j,0) = 0.0_rt;
        CF(i,j,0) = 0.0_rt;
        for (int k=1; k<=N; k++)
        {
            //
            //  Use conservative, parabolic spline reconstruction of vertical
            //  viscosity derivatives.  Then, time step vertical viscosity term
            //  implicitly by solving a tridiagonal system.
            //
            const Real oHzkm1 = 1.0_rt/ Hzk(i,j,k-1);
            const Real oHz = 1.0_rt/ Hzk(i,j,k);
            //const Real oHzkp1 = 1.0_rt/ Hzk(i,j,k+1);

            FC(i,j,k) = sixth * Hzk(i,j,k-1) - dt_lev * AK(i,j,k-1) / Hzk(i,j,k-1);
            CF(i,j,k) = sixth * Hzk(i,j,k  ) - dt_lev * AK(i,j,k+1) / Hzk(i,j,k  );

            BC(i,j,k) = third * (Hzk(i,j,k-1) + Hzk(i,j,k  )) + dt_lev * AK(i,j,k) * (oHzkm1 + oHz);
            Real cff = 1.0_rt / (BC(i,j,k) - FC(i,j,k) * CF(i,j,k-1));
            CF(i,j,k) *= cff;
            DC(i,j,k) = cff * (phi(i,j,k  ,nnew) - phi(i,j,k-1,nnew) - FC(i,j,k)*DC(i,j,k-1));
        } // k
    });
#ifdef AMREX_USE_GPU
    Gpu::synchronize();
#else
#endif
    Gpu::streamSynchronize();
    ParallelFor(makeSlab(phi_bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        DC(i,j,N+1) = 0.0_rt;
        for(int k=N; k>=1; k--) {
            //
            //  Backward substitution.
            //
            DC(i,j,k) -= CF(i,j,k)*DC(i,j,k+1);
        }
    });

#ifdef AMREX_USE_GPU
    Gpu::synchronize();
#endif

    ParallelFor(surroundingNodes(phi_bx,2), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        DC(i,j,k) *= AK(i,j,k);
    });

    ParallelFor(phi_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const Real oHz = 1.0_rt/ Hzk(i,j,k);
        Real cff = dt_lev*oHz*(DC(i,j,k+1)-DC(i,j,k));
        phi(i,j,k) += cff;
    });
}
