#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

#ifdef UV_COR
//
// Start 3d step
//

void
ROMSX::coriolis (const Box& bx,
                 Array4<Real> uold  , Array4<Real> vold,
                 Array4<Real> ru, Array4<Real> rv,
                 Array4<Real> Hz, Array4<Real> fomn,
                 int nrhs)
{
    // Need to include uv3dmix
    //
    //-----------------------------------------------------------------------
    //  Add in Coriolis terms.
    //-----------------------------------------------------------------------
    //

    Box ubx = surroundingNodes(bx,0);
    Box vbx = surroundingNodes(bx,1);

    amrex::ParallelFor(ubx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real UFx_i   = 0.5 * Hz(i  ,j,k) * fomn(i  ,j,0) * (vold(i  ,j,k,nrhs)+vold(i  ,j+1,k,nrhs));
        Real UFx_im1 = 0.5 * Hz(i-1,j,k) * fomn(i-1,j,0) * (vold(i-1,j,k,nrhs)+vold(i-1,j+1,k,nrhs));
        ru(i,j,k,nrhs) += 0.5*(UFx_i + UFx_im1);
    });

    amrex::ParallelFor(vbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real VFe_j   = 0.5 * Hz(i,j  ,k) * fomn(i,j  ,0) * (uold(i,j  ,k,nrhs)+uold(i+1,j  ,k,nrhs));
        Real VFe_jm1 = 0.5 * Hz(i,j-1,k) * fomn(i,j-1,0) * (uold(i,j-1,k,nrhs)+uold(i+1,j-1,k,nrhs));
        rv(i,j,k,nrhs) -= 0.5*(VFe_j + VFe_jm1);
    });
}
#endif
