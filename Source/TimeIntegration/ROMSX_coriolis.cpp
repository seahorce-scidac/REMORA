#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// Start 3d step
//

void
ROMSX::coriolis (const Box& bx, const Box& gbx,
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

    //Box ubx = surroundingNodes(bx,0);
    //Box vbx = surroundingNodes(bx,1);
    Box ubxShift = bx;
    ((((ubxShift.growLo(0,NGROW-1)).growLo(1,NGROW)).growHi(0,NGROW)).growHi(1,NGROW-1));
    Box vbxShift = bx;
    ((((vbxShift.growLo(0,NGROW)).growLo(1,NGROW-1)).growHi(0,NGROW-1)).growHi(1,NGROW));

    BoxArray ba_ubxShift = intersect(BoxArray(ubxShift), gbx);
    AMREX_ASSERT((ba_ubxShift.size() == 1));
    ubxShift = ba_ubxShift[0];

    BoxArray ba_vbxShift = intersect(BoxArray(vbxShift), gbx);
    AMREX_ASSERT((ba_vbxShift.size() == 1));
    vbxShift = ba_vbxShift[0];

    amrex::ParallelFor(ubxShift,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real UFx_i   = 0.5 * Hz(i  ,j,k) * fomn(i  ,j,0) * (vold(i  ,j,k,nrhs)+vold(i  ,j+1,k,nrhs));
        Real UFx_im1 = 0.5 * Hz(i-1,j,k) * fomn(i-1,j,0) * (vold(i-1,j,k,nrhs)+vold(i-1,j+1,k,nrhs));
        if ((verbose>=2) && (printinloop > 0))
            printf("cor %d %d %d %25.25g %25.25g %25.25g\n", i,j,k,UFx_i, UFx_im1, ru(i,j,k,nr));
        ru(i,j,k,nr) += 0.5*(UFx_i + UFx_im1);
    });

    amrex::ParallelFor(vbxShift,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real VFe_j   = 0.5 * Hz(i,j  ,k) * fomn(i,j  ,0) * (uold(i,j  ,k,nrhs)+uold(i+1,j  ,k,nrhs));
        Real VFe_jm1 = 0.5 * Hz(i,j-1,k) * fomn(i,j-1,0) * (uold(i,j-1,k,nrhs)+uold(i+1,j-1,k,nrhs));
        rv(i,j,k,nr) -= 0.5*(VFe_j + VFe_jm1);
    });
}
