#include <ROMSX.H>

using namespace amrex;

//
// Start 3d step
//

void
ROMSX::coriolis (const Box& xbx, const Box& ybx,
                 Array4<Real> uold  , Array4<Real> vold,
                 Array4<Real> ru, Array4<Real> rv,
                 Array4<Real> Hz, Array4<Real> fomn,
                 int nrhs, int nr)
{
    //
    //-----------------------------------------------------------------------
    //  Add in Coriolis terms.
    //-----------------------------------------------------------------------
    //

    ParallelFor(xbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real UFx_i   = 0.5 * Hz(i  ,j,k) * fomn(i  ,j,0) * (vold(i  ,j,k,nrhs)+vold(i  ,j+1,k,nrhs));
        Real UFx_im1 = 0.5 * Hz(i-1,j,k) * fomn(i-1,j,0) * (vold(i-1,j,k,nrhs)+vold(i-1,j+1,k,nrhs));
        ru(i,j,k,nr) += 0.5*(UFx_i + UFx_im1);
    });

    ParallelFor(ybx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real VFe_j   = 0.5 * Hz(i,j  ,k) * fomn(i,j  ,0) * (uold(i,j  ,k,nrhs)+uold(i+1,j  ,k,nrhs));
        Real VFe_jm1 = 0.5 * Hz(i,j-1,k) * fomn(i,j-1,0) * (uold(i,j-1,k,nrhs)+uold(i+1,j-1,k,nrhs));
        rv(i,j,k,nr) -= 0.5*(VFe_j + VFe_jm1);
    });
}
