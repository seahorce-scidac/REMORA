#include <ROMSX.H>

using namespace amrex;

//
// prestep_uv_3d
//

void
ROMSX::vert_visc_3d (const Box& phi_bx, const int ioff, const int joff,
                     Array4<Real> phi,
                     Array4<Real> Hz, Array4<Real> Hzk,
                     Array4<Real> oHz,
                     Array4<Real> AK, Array4<Real> Akv,
                     Array4<Real> BC, Array4<Real> DC,
                     Array4<Real> FC, Array4<Real> CF,
                     const int nnew, const int N, const Real dt_lev)
{
    //
    // Put Hzk on the x- or y-face as appropriate, or leave on cell center for tracers
    //
    if (verbose > 1)
        amrex::Print() << "updating on box in vert_visc_3d: " << phi_bx << std::endl;
    auto phi_bxD = phi_bx;
    phi_bxD.makeSlab(2,0);
    ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Hzk(i,j,k)=0.5*(Hz(i-ioff,j-joff,k)+Hz(i,j,k));
    });
    if (verbose > 2) {
    ParallelFor(phi_bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            printf("%d %d %d  %15.15g %15.15g %15.15g  Hzk Hz2\n", i,j,k, Hzk(i,j,k), Hz(i-ioff, j-joff,k),Hz(i,j,k));
        });
    }

    //
    // Define oHz = (1/Hz)
    //
    ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        oHz(i,j,k) = 1.0/ Hzk(i,j,k);
    });

    //
    // Put Akv on the x- or y-face as appropriate, or leave on cell center for tracers
    //
    ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        AK(i,j,k) = 0.5 * (Akv(i-ioff,j-joff,k)+Akv(i,j,k));
    });
    if (verbose > 2) {
        ParallelFor(phi_bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            printf("%d %d %d  %15.15g %15.15g %15.15g  AKs\n", i,j,k, AK(i,j,k), Akv(i-ioff,j-joff,k), Akv(i,j,k));
        });
    }
    Gpu::streamSynchronize();
#ifdef AMREX_USE_GPU
    Gpu::synchronize();
#else
#endif
    /////////////////// This and the following loop is the first non-matching thing that affects plotfile comparison for cuda
    // Begin vertical viscosity term
    // NOTE: vertical viscosity term for tracers is identical except AK=Akt
    ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        for(int k=0; k<=N; k++) {
        //
        //  Use conservative, parabolic spline reconstruction of vertical
        //  viscosity derivatives.  Then, time step vertical viscosity term
        //  implicitly by solving a tridiagonal system.
        //
        Real cff;
        Real cff1=1.0/6.0;

        if(k-1>=0)
            FC(i,j,k)=cff1*Hzk(i,j,k  )-dt_lev*AK(i,j,k-1)*oHz(i,j,k  );
        else
            FC(i,j,k)=cff1*Hzk(i,j,k  );
        //printf("%d %d %d %d %15.15g FC\n",i,j,k,0,FC(i,j,k));
        //      amrex::Abort("first index");
        if(k<=N-1)
         {
             CF(i,j,k)=cff1*Hzk(i,j,k+1)-dt_lev*AK(i,j,k+1)*oHz(i,j,k+1);
         }
         //printf("%d %d %d %d %15.15g CF\n",i,j,k,0,CF(i,j,k));
         //
         //  LU decomposition and forward substitution.
         //
         cff1=1.0/3.0;
         if(k==0)
         {
             BC(i,j,k)=cff1*(Hzk(i,j,k)+Hzk(i,j,k+1))+
                     dt_lev*AK(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1));
             cff=1.0/(BC(i,j,k)-FC(i,j,k)*0.0);
             CF(i,j,k) *= cff;
             DC(i,j,k) = cff*(phi(i,j,k+1,nnew)-phi(i,j,k,nnew)-FC(i,j,k)*0.0);
         }
         if(k+1<=N&&k>=1)
         {
                 BC(i,j,k)=cff1*(Hzk(i,j,k)+Hzk(i,j,k+1))+
                     dt_lev*AK(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1));
                 cff=1.0/(BC(i,j,k)-FC(i,j,k)*CF(i,j,k-1));
             CF(i,j,k) *= cff;
             DC(i,j,k) = cff*(phi(i,j,k+1,nnew)-phi(i,j,k,nnew)-FC(i,j,k)*DC(i,j,k-1));
         }
         //if (verbose > 2) {
         //   printf("%d %d %d  %15.15g %15.15g %15.15g %15.15g %15.15g vvisc Hzk2 AK oHz2\n", i,j,k,Hzk(i,j,k+1), Hzk(i,j,k), AK(i,j,k), oHz(i,j,k+1), oHz(i,j,k));
         //   printf("%d %d %d %d %15.15g %15.15g %15.15g %15.15g %15.15g vvisc BC cff CF DC\n",i,j,k,0,BC(i,j,k),cff,CF(i,j,k),FC(i,j,k),DC(i,j,k));
         //   printf("%d %d %d %d %15.15g %15.15g %15.15g %15.15g %15.15g vvisc cff u(k+1) u(k)FC DC(k-1)\n",i,j,k,0,cff,phi(i,j,k+1,nnew),phi(i,j,k,nnew),FC(i,j,k),DC(i,j,k-1));
        //}
        }
    });
#ifdef AMREX_USE_GPU
    Gpu::synchronize();
#else
#endif
    Gpu::streamSynchronize();
    ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
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
            cff = dt_lev*oHz(i,j,k)*(DC(i,j,k)-DC(i,j,k-1));
        } else {
            cff = dt_lev*oHz(i,j,k)*(DC(i,j,k));
        }
        //if (verbose > 2) {
        //    printf("%d %d %d  %15.15g %15.15g %15.15g %15.15g tracer phi cff1\n",
        //        i,j,k,phi(i,j,k),cff, DC(i,j,k), DC(i,j,k-1));
        //}
        phi(i,j,k) += cff;
     });
}
