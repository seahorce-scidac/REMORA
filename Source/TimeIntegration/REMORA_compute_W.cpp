#include <REMORA.H>

using namespace amrex;

/**
 * Vertical volume flux across the w-levels (ROMS omega times the cell area), from the
 * vertical integral of the horizontal mass flux divergence. W is zero at the bed and at
 * the free surface.
 *
 * @param[in   ] bx2d   k = 0 slab of the columns to compute
 * @param[inout] W      vertical volume flux; W(i,j,0) must already be zero
 * @param[in   ] Huon   u-volume flux
 * @param[in   ] Hvom   v-volume flux
 * @param[in   ] z_w    z coordinates at w points
 * @param[in   ] h      bathymetry
 * @param[in   ] N      number of vertical levels minus one
 */
void
REMORA::compute_W (const Box& bx2d,
                   const Array4<Real      >& W,
                   const Array4<Real const>& Huon,
                   const Array4<Real const>& Hvom,
                   const Array4<Real const>& z_w,
                   const Array4<Real const>& h,
                   const int N)
{
    ParallelFor(bx2d, N+1,
    [=] AMREX_GPU_DEVICE (int i, int j, int , int kk)
    {
        //  Starting with zero vertical velocity at the bottom, integrate
        //  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
        //  contains the vertical velocity at the free-surface, d(zeta)/d(t).
        //  Notice that barotropic mass flux divergence is not used directly.
        //
        int k = kk + 1;
        W(i,j,k) = W(i,j,k-1) - (Huon(i+1,j,k-1)-Huon(i,j,k-1)) - (Hvom(i,j+1,k-1)-Hvom(i,j,k-1));
    });
    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        W(i,j,N+1)=W(i,j,N+1)/(z_w(i,j,N+1)+h(i,j,0,0)); // wrk_i
    });
    ParallelFor(bx2d, N,
    [=] AMREX_GPU_DEVICE (int i, int j, int , int kk)
    {
        int k = kk + 1;
        W(i,j,k) = W(i,j,k)- W(i,j,N+1)*(z_w(i,j,k)+h(i,j,0,0));
    });
    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        W(i,j,N+1) = zero;
    });
}
