#ifndef REMORA_BOUNDARYINTERPMATH_H_
#define REMORA_BOUNDARYINTERPMATH_H_

/**
 * \file REMORA_BoundaryInterpMath.H
 *
 * Index arithmetic for interpolating boundary data along a boundary strip onto a refined level
 * (NCTimeSeriesBoundary::interp_fab), for data that is cell-centered or nodal along the strip.
 * Nothing here depends on REMORA_USE_NETCDF, so Tests/Unit can check it directly.
 */

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace boundary_interp {

/**
 * For fine cell f along a cell-centered strip refined by r, find the coarse cell c below it
 * and the weight w of coarse cell c+1, so that the fine value is
 * crse[c] + w*(crse[c+1] - crse[c]). Fine cell f is centered at (f+0.5)/r coarse cell widths
 * and coarse cell c at c+0.5. When w is zero, crse[c+1] need not exist.
 *
 * @param[in   ] f  fine cell index
 * @param[in   ] r  refinement ratio along the strip
 * @param[  out] c  lower coarse cell index
 * @param[  out] w  weight of coarse cell c+1, in [0, 1)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
cell_centered_weight (int f, int r, int& c, amrex::Real& w)
{
    // 2r times the fine center's position, measured in coarse cells from the center of
    // coarse cell 0. Integer arithmetic gives the coarse index exactly, and the weight with a
    // single rounding.
    const int num = 2*f + 1 - r;
    const int den = 2*r;
    c = (num >= 0) ? num / den : -((-num + den - 1) / den);
    w = static_cast<amrex::Real>(num - c*den) / static_cast<amrex::Real>(den);
}

/**
 * For fine node f along a nodal strip refined by r, find the coarse node c at or below it and
 * the weight w of coarse node c+1, as for cell_centered_weight. Fine node f lies at f/r coarse
 * cell widths and coarse node c at c, so w is zero where the nodes coincide.
 *
 * @param[in   ] f  fine node index
 * @param[in   ] r  refinement ratio along the strip
 * @param[  out] c  lower coarse node index
 * @param[  out] w  weight of coarse node c+1, in [0, 1)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
nodal_weight (int f, int r, int& c, amrex::Real& w)
{
    c = (f >= 0) ? f / r : -((-f + r - 1) / r);
    w = static_cast<amrex::Real>(f - c*r) * (amrex::Real(1.0) / static_cast<amrex::Real>(r));
}

} // namespace boundary_interp

#endif
