#ifndef REMORA_BOUNDARYINTERPMATH_H_
#define REMORA_BOUNDARYINTERPMATH_H_

/**
 * \file REMORA_BoundaryInterpMath.H
 *
 * Index arithmetic for interpolating cell-centered boundary data along a boundary strip onto
 * a refined level (NCTimeSeriesBoundary::interp_fab). Nothing here depends on
 * REMORA_USE_NETCDF, so Tests/Unit can check it directly.
 */

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace boundary_interp {

/**
 * For fine cell f along a cell-centered strip refined by r, find the coarse cell c below it
 * and the weight w of coarse cell c+1, so that the fine value is
 * (1-w)*crse[c] + w*crse[c+1]. Fine cell f is centered at (f+0.5)/r coarse cell widths and
 * coarse cell c at c+0.5. When w is zero, crse[c+1] need not exist.
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
    // coarse cell 0. Integers keep the weights exact.
    const int num = 2*f + 1 - r;
    const int den = 2*r;
    c = (num >= 0) ? num / den : -((-num + den - 1) / den);
    w = static_cast<amrex::Real>(num - c*den) / static_cast<amrex::Real>(den);}

} // namespace boundary_interp

#endif
