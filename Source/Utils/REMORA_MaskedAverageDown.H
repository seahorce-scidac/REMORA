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
 * Where every fine point is wet these have to reproduce amrex_avgdown and amrex_avgdown_faces
 * bit for bit, or every existing multilevel answer moves in its last bit.
 */
namespace REMORAMaskedAvgDown {

/** Cell-centered fields: every fine cell under the coarse one. The mask is 2D, read at k=0. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void avgdown_masked (int i, int j, int k, int n,
                     amrex::Array4<amrex::Real>       const& crse,
                     amrex::Array4<amrex::Real const> const& fine,
                     amrex::Array4<amrex::Real const> const& fmsk,
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
        const Real m = amrex::min(Real(1.0), fmsk(ii+iref,jj+jref,0));
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

/** Cell-centered tracers weighted by each fine cell's volume, Hz/(pm*pn), as well as its
 *  mask, so the coarse cell takes the content of the fine cells under it rather than their
 *  mean concentration. Not what ROMS does; see cf_avgdown_volume_weighted. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void avgdown_masked_volume (int i, int j, int k, int n,
                            amrex::Array4<amrex::Real>       const& crse,
                            amrex::Array4<amrex::Real const> const& fine,
                            amrex::Array4<amrex::Real const> const& fmsk,
                            amrex::Array4<amrex::Real const> const& cmsk,
                            amrex::Array4<amrex::Real const> const& Hz,
                            amrex::Array4<amrex::Real const> const& pm,
                            amrex::Array4<amrex::Real const> const& pn,
                            amrex::IntVect const& ratio) noexcept
{
    using amrex::Real;
    const int ii = i*ratio[0];
    const int jj = j*ratio[1];
    const int kk = k*ratio[2];

    Real num = Real(0.0);
    Real den = Real(0.0);
    for (int kref = 0; kref < ratio[2]; ++kref) {
    for (int jref = 0; jref < ratio[1]; ++jref) {
    for (int iref = 0; iref < ratio[0]; ++iref) {
        const int fi = ii+iref;
        const int fj = jj+jref;
        const int fk = kk+kref;
        const Real w = amrex::min(Real(1.0), fmsk(fi,fj,0)) * Hz(fi,fj,fk)
                     / (pm(fi,fj,0) * pn(fi,fj,0));
        num += fine(fi,fj,fk,n) * w;
        den += w;
    }}}

    crse(i,j,k,n) = (den > Real(0.0)) ? num / den * cmsk(i,j,0) : Real(0.0);
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
