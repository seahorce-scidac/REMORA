#include <REMORA.H>

using namespace amrex;

/**
 * Convert the RHS momentum terms ru, rv, ru2d or rv2d between the form the kernels keep
 * them in and the per-area form, which is that times cff = (pm+pm')*(pn+pn'), the
 * inverse cell area of the face up to a factor of four. Coarse and fine faces differ in
 * area, so data moved between levels is interpolated in the per-area form. The whole
 * array is converted, ghost cells included.
 *
 * mf must share its level's boxes and distribution, so that the metrics of the same box
 * index cover it, and carry NGROW ghost cells laterally.
 *
 * @param[in   ] lev          level whose pm and pn apply
 * @param[inout] mf           x-face (dir 0) or y-face (dir 1) data, 3D or 2D
 * @param[in   ] dir          face direction
 * @param[in   ] to_per_area  true to multiply by cff, false to divide by it
 */
void
REMORA::scale_rhs_by_area (int lev, MultiFab& mf, int dir, bool to_per_area)
{
    const int ioff = (dir == 0) ? 1 : 0;
    const int joff = (dir == 1) ? 1 : 0;
    const int ncomp = mf.nComp();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real const> const& pm = vec_pm[lev]->const_array(mfi);
        Array4<Real const> const& pn = vec_pn[lev]->const_array(mfi);
        Array4<Real      > const& r  = mf.array(mfi);

        Box bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));

        if (to_per_area) {
            ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff = (pm(i,j,0)+pm(i-ioff,j-joff,0)) * (pn(i,j,0)+pn(i-ioff,j-joff,0));
                r(i,j,k,n) = r(i,j,k,n) * cff;
            });
        } else {
            ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff = (pm(i,j,0)+pm(i-ioff,j-joff,0)) * (pn(i,j,0)+pn(i-ioff,j-joff,0));
                r(i,j,k,n) = r(i,j,k,n) / cff;
            });
        }
    }
}

/**
 * A per-area copy of one level's ru, rv, ru2d or rv2d, ghost cells included, for a fill
 * of another level to interpolate from.
 *
 * @param[in] lev  level of src
 * @param[in] src  RHS data in the kernels' form
 * @param[in] dir  face direction of src
 */
MultiFab
REMORA::rhs_per_area_copy (int lev, const MultiFab& src, int dir)
{
    MultiFab dst(src.boxArray(), src.DistributionMap(), src.nComp(), src.nGrowVect());
    MultiFab::Copy(dst, src, 0, 0, src.nComp(), src.nGrowVect());
    scale_rhs_by_area(lev, dst, dir, true);
    return dst;
}
