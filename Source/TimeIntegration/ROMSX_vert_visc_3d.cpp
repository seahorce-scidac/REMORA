#include <ROMSX.H>

using namespace amrex;

//
// vert_visc_3d
//

void
ROMSX::vert_visc_3d (const Box& phi_bx, const int ioff, const int joff,
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
        Hzk(i,j,k)=0.5*(Hz(i-ioff,j-joff,k)+Hz(i,j,k));
    });

    //
    // Put Akv on the x- or y-face as appropriate, or leave on cell center for tracers
    //
    ParallelFor(phi_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        AK(i,j,k) = 0.5 * (Akv(i-ioff,j-joff,k)+Akv(i,j,k));
    });

    Gpu::streamSynchronize();
#ifdef AMREX_USE_GPU
    Gpu::synchronize();
#else
#endif
    /////////////////// This and the following loop is the first non-matching thing that affects plotfile comparison for cuda
    // NOTE: vertical viscosity term for tracers is identical except AK=Akt
    const Real sixth = Real(1.0) / Real(6.0);
    const Real third = Real(1.0) / Real(3.0);

    ParallelFor(makeSlab(phi_bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        for (int k=0; k<=N; k++)
        {
            //
            //  Use conservative, parabolic spline reconstruction of vertical
            //  viscosity derivatives.  Then, time step vertical viscosity term
            //  implicitly by solving a tridiagonal system.
            //
            Real cff;

            const Real oHz = 1.0/ Hzk(i,j,k);

            FC(i,j,k) = (k >= 1) ? sixth*Hzk(i,j,k  )-dt_lev*AK(i,j,k-1)*oHz :
                                   sixth*Hzk(i,j,k);

            if (k < N)
            {
                const Real oHzkp1 = 1.0/ Hzk(i,j,k+1);

                CF(i,j,k)=sixth*Hzk(i,j,k+1)-dt_lev*AK(i,j,k+1)*oHzkp1;

                if (k==0)
                {
                    BC(i,j,k)=third*(Hzk(i,j,k)+Hzk(i,j,k+1)) + dt_lev*AK(i,j,k)*(oHz+oHzkp1);

                    cff=1.0/(BC(i,j,k)-FC(i,j,k)*0.0);
                    CF(i,j,k) *= cff;
                    DC(i,j,k) = cff*(phi(i,j,k+1,nnew)-phi(i,j,k,nnew)-FC(i,j,k)*0.0);

                } else {

                    BC(i,j,k)=third*(Hzk(i,j,k)+Hzk(i,j,k+1)) + dt_lev*AK(i,j,k)*(oHz+oHzkp1);
                    cff=1.0/(BC(i,j,k)-FC(i,j,k)*CF(i,j,k-1));
                    CF(i,j,k) *= cff;
                    DC(i,j,k) = cff*(phi(i,j,k+1,nnew)-phi(i,j,k,nnew)-FC(i,j,k)*DC(i,j,k-1));
                }
            } // k < N
        } // k
    });
#ifdef AMREX_USE_GPU
    Gpu::synchronize();
#else
#endif
    Gpu::streamSynchronize();
    ParallelFor(makeSlab(phi_bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
       for(int k=0; k<=N; k++) {
           //
           //  Backward substitution.
           //
           DC(i,j,N)=0.0;

           if(N-k+1<=N && N>=k) //-N,1,-1 => kidx =N-k+1
           {
               if(N+1<k || N+2<k) amrex::Abort("-1 here");
               DC(i,j,N-k) -= CF(i,j,N-k)*DC(i,j,N-k+1);
           }
       }
    });

#ifdef AMREX_USE_GPU
    Gpu::synchronize();
#endif

    ParallelFor(phi_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        DC(i,j,k) *= AK(i,j,k);
    });

    ParallelFor(phi_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff;
        if(k-1>=0) {
            const Real oHz = 1.0/ Hzk(i,j,k);
            cff = dt_lev*oHz*(DC(i,j,k)-DC(i,j,k-1));
        } else {
            const Real oHz = 1.0/ Hzk(i,j,k);
            cff = dt_lev*oHz*(DC(i,j,k));
        }
        phi(i,j,k) += cff;
     });
}
