#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// prestep_uv_3d
//
void
ROMSX::vert_visc_3d (const Box& phi_bx, const int ioff, const int joff,
                     Array4<Real> phi_arr,
                     Array4<Real> Hz_arr, Array4<Real> Hzk_arr,
                     Array4<Real> oHz_arr,
                     Array4<Real> AK_arr, Array4<Real> Akv_arr,
                     Array4<Real> BC_arr, Array4<Real> DC_arr,
                     Array4<Real> FC_arr, Array4<Real> CF_arr,
                     const int nnew, const int N, const Real dt_lev)
{
    //
    // Put Hzk on the x- or y-face as appropriate, or leave on cell center for tracers
    //
    amrex::Print() << "updating on box in vert_visc_3d: " << phi_bx << std::endl;
    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Hzk_arr(i,j,k)=0.5*(Hz_arr(i-ioff,j-joff,k)+Hz_arr(i,j,k));
    });

    //
    // Define oHz = (1/Hz)
    //
    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        oHz_arr(i,j,k) = 1.0/ Hzk_arr(i,j,k);
    });

    //
    // Put Akv on the x- or y-face as appropriate, or leave on cell center for tracers
    //
    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        AK_arr(i,j,k) = 0.5 * (Akv_arr(i-ioff,j-joff,k)+Akv_arr(i,j,k));
    });

    // Begin vertical viscosity term
    // NOTE: vertical viscosity term for tracers is identical except AK=Akt
    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //
        //  Use conservative, parabolic spline reconstruction of vertical
        //  viscosity derivatives.  Then, time step vertical viscosity term
        //  implicitly by solving a tridiagonal system.
        //
        Real cff;
        Real cff1=1.0/6.0;

        if(k>=1)
            FC_arr(i,j,k)=cff1*Hzk_arr(i,j,k  )-dt_lev*AK_arr(i,j,k-1)*oHz_arr(i,j,k  );
        else
            FC_arr(i,j,k)=cff1*Hzk_arr(i,j,k  );
        if(k<=N-1)
         {
             CF_arr(i,j,k)=cff1*Hzk_arr(i,j,k+1)-dt_lev*AK_arr(i,j,k+1)*oHz_arr(i,j,k+1);
         }

         //
         //  LU decomposition and forward substitution.
         //
         cff1=1.0/3.0;
         if(k==0)
         {
             BC_arr(i,j,k)=cff1*(Hzk_arr(i,j,k)+Hzk_arr(i,j,k+1))+
                     dt_lev*AK_arr(i,j,k)*(oHz_arr(i,j,k)+oHz_arr(i,j,k+1));
             cff=1.0/(BC_arr(i,j,k)-FC_arr(i,j,k)*0.0);
             CF_arr(i,j,k) *= cff;
             DC_arr(i,j,k) = cff*(phi_arr(i,j,k+1,nnew)-phi_arr(i,j,k,nnew)-FC_arr(i,j,k)*0.0);
         }
         if(k+1<=N&&k>=1)
         {
                 BC_arr(i,j,k)=cff1*(Hzk_arr(i,j,k)+Hzk_arr(i,j,k+1))+
                     dt_lev*AK_arr(i,j,k)*(oHz_arr(i,j,k)+oHz_arr(i,j,k+1));
                 cff=1.0/(BC_arr(i,j,k)-FC_arr(i,j,k)*CF_arr(i,j,k-1));
             CF_arr(i,j,k) *= cff;
             DC_arr(i,j,k) = cff*(phi_arr(i,j,k+1,nnew)-phi_arr(i,j,k,nnew)-FC_arr(i,j,k)*DC_arr(i,j,k-1));
         }

    });

    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
       //
       //  Backward substitution.
       //
       DC_arr(i,j,N)=0.0;

       if(N-k+1<=N && N>=k) //-N,1,-1 => kidx =N-k+1
       {
           if(N+1<k || N+2<k) amrex::Abort("-1 here");
           DC_arr(i,j,N-k) -= CF_arr(i,j,N-k)*DC_arr(i,j,N-k+1);
       }
    });

    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        DC_arr(i,j,k) *= AK_arr(i,j,k);
    });

    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff;
        if(k-1>=0) {
            cff = dt_lev*oHz_arr(i,j,k)*(DC_arr(i,j,k)-DC_arr(i,j,k-1));
        } else {
            cff = dt_lev*oHz_arr(i,j,k)*(DC_arr(i,j,k));
        }
        phi_arr(i,j,k) += cff;
     });
}
