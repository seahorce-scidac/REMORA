#ifndef REMORA_MASKED_AVERAGE_DOWN_H_
#define REMORA_MASKED_AVERAGE_DOWN_H_

#include <AMReX_Algorithm.H>
#include <AMReX_Array4.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>

/**
 * Two-way nesting's fine-to-coarse update, weighted by the land/sea mask.
 *
 * Ported from ROMS fine2coarse2d/fine2coarse3d in Nonlinear/nesting.F, the AreaAvg = .FALSE.
 * branch that every call site there uses:
 *
 *     num  = sum of (fine value * fine mask) over the fine points under the coarse one
 *     den  = sum of min(1, fine mask) over the same
 *     crse = num/den * coarse mask
 *
 * Dividing by the count of wet fine points rather than by how many there are means a coarse
 * point only partly covered by water takes the mean of that water, not a value diluted toward
 * zero by land. That is not conservative; ROMS makes the same trade deliberately. The coarse
 * mask has the final say on wet or dry.
 *
 * Where every fine point is wet the mask-weighted kernels have to reproduce amrex_avgdown and
 * amrex_avgdown_faces bit for bit, or every existing multilevel answer moves in its last bit.
 * avgdown_weighted with a volume weight (cf_avgdown_volume_weighted) is the exception: it
 * changes the answer by design.
 */
namespace REMORAMaskedAvgDown {

/** Weight of a fine cell: its 2D mask, capped at 1 and read at k=0 */
struct MaskWeight {
    amrex::Array4<amrex::Real const> msk;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int i, int j, int) const noexcept {
        return amrex::min(amrex::Real(1.0), msk(i,j,0));
    }
};

/** Weight of a fine cell read from a 3D array, e.g. mask times volume */
struct FieldWeight {
    amrex::Array4<amrex::Real const> w;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int i, int j, int k) const noexcept { return w(i,j,k); }
};

/** Cell-centered fields: the weighted mean of every fine cell under the coarse one. */
template <typename Weight>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void avgdown_weighted (int i, int j, int k, int n,
                       amrex::Array4<amrex::Real>       const& crse,
                       amrex::Array4<amrex::Real const> const& fine,
                       Weight const& weight,
                       amrex::Array4<amrex::Real const> const& cmsk,
                       int ccomp, int fcomp, amrex::IntVect const& ratio) noexcept
{
    using amrex::Real;
    const int facx = ratio[0];
    const int facy = ratio[1];
    const int facz = ratio[2];
    const int ii = i*facx;
    const int jj = j*facy;
    const int kk = k*facz;

    Real num = Real(0.0);
    Real den = Real(0.0);
    // AMReX's loop order, and the mask multiplying the value before it is accumulated rather
    // than after (x * 1.0 == x), are what make the all-wet case bit-identical to amrex_avgdown.
    for (int kref = 0; kref < facz; ++kref) {
    for (int jref = 0; jref < facy; ++jref) {
    for (int iref = 0; iref < facx; ++iref) {
        const Real m = weight(ii+iref,jj+jref,kk+kref);
        num += fine(ii+iref,jj+jref,kk+kref,n+fcomp) * m;
        den += m;
    }}}

    // Multiplication by the reciprocal, not num/den: this is AMReX's volfrac, and dividing
    // would differ in the last bit. den == 0 under a wet coarse point cannot happen once
    // check_mask_consistency passes, so write zero rather than diagnose it from a kernel.
    crse(i,j,k,n+ccomp) = (den > Real(0.0))
                        ? num * (Real(1.0)/den) * cmsk(i,j,0)
                        : Real(0.0);
}

/** Cell-centered fields weighted by the fine mask alone, as ROMS's fine2coarse3d does. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void avgdown_masked (int i, int j, int k, int n,
                     amrex::Array4<amrex::Real>       const& crse,
                     amrex::Array4<amrex::Real const> const& fine,
                     amrex::Array4<amrex::Real const> const& fmsk,
                     amrex::Array4<amrex::Real const> const& cmsk,
                     int ccomp, int fcomp, amrex::IntVect const& ratio) noexcept
{
    avgdown_weighted(i,j,k,n,crse,fine,MaskWeight{fmsk},cmsk,ccomp,fcomp,ratio);
}

/** Face-centered fields: only the fine faces in the coarse face's own plane, the stencil
 *  amrex_avgdown_faces uses, so "wet" means the flow really does cross here. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void avgdown_faces_masked (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real>       const& crse,
                           amrex::Array4<amrex::Real const> const& fine,
                           amrex::Array4<amrex::Real const> const& fmsk,
                           amrex::Array4<amrex::Real const> const& cmsk,
                           int ccomp, int fcomp, amrex::IntVect const& ratio, int idir) noexcept
{
    using amrex::Real;
    const int facx = ratio[0];
    const int facy = ratio[1];
    const int facz = ratio[2];
    const int ii = i*facx;
    const int jj = j*facy;
    const int kk = k*facz;

    Real num = Real(0.0);
    Real den = Real(0.0);

    if (idir == 0) {
        for (int kref = 0; kref < facz; ++kref) {
        for (int jref = 0; jref < facy; ++jref) {
            const Real m = amrex::min(Real(1.0), fmsk(ii,jj+jref,0));
            num += fine(ii,jj+jref,kk+kref,n+fcomp) * m;
            den += m;
        }}
    } else if (idir == 1) {
        for (int kref = 0; kref < facz; ++kref) {
        for (int iref = 0; iref < facx; ++iref) {
            const Real m = amrex::min(Real(1.0), fmsk(ii+iref,jj,0));
            num += fine(ii+iref,jj,kk+kref,n+fcomp) * m;
            den += m;
        }}
    } else {
        for (int jref = 0; jref < facy; ++jref) {
        for (int iref = 0; iref < facx; ++iref) {
            const Real m = amrex::min(Real(1.0), fmsk(ii+iref,jj+jref,0));
            num += fine(ii+iref,jj+jref,kk,n+fcomp) * m;
            den += m;
        }}
    }

    // As in avgdown_masked: reciprocal rather than divide, and zero when nothing is wet.
    crse(i,j,k,n+ccomp) = (den > Real(0.0))
                        ? num * (Real(1.0)/den) * cmsk(i,j,0)
                        : Real(0.0);
}

} // namespace REMORAMaskedAvgDown

#endif
