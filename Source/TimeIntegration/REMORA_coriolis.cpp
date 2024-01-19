#include <REMORA.H>

using namespace amrex;

//
// Start 3d step
//

void
REMORA::coriolis (const Box& xbx, const Box& ybx,
                 const Array4<Real const>& uold,
                 const Array4<Real const>& vold,
                 const Array4<Real      >& ru,
                 const Array4<Real      >& rv,
                 const Array4<Real const>& Hz,
                 const Array4<Real const>& fomn,
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
        Real UFx_i   = 0.5_rt * Hz(i  ,j,k) * fomn(i  ,j,0) * (vold(i  ,j,k,nrhs)+vold(i  ,j+1,k,nrhs));
        Real UFx_im1 = 0.5_rt * Hz(i-1,j,k) * fomn(i-1,j,0) * (vold(i-1,j,k,nrhs)+vold(i-1,j+1,k,nrhs));
        ru(i,j,k,nr) += 0.5_rt*(UFx_i + UFx_im1);
    });

    ParallelFor(ybx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real VFe_j   = 0.5_rt * Hz(i,j  ,k) * fomn(i,j  ,0) * (uold(i,j  ,k,nrhs)+uold(i+1,j  ,k,nrhs));
        Real VFe_jm1 = 0.5_rt * Hz(i,j-1,k) * fomn(i,j-1,0) * (uold(i,j-1,k,nrhs)+uold(i+1,j-1,k,nrhs));
        rv(i,j,k,nr) -= 0.5_rt*(VFe_j + VFe_jm1);
    });
}
