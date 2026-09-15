#include <REMORA.H>

using namespace amrex;

/**
 * Scale one level's RHS momentum terms by 1/cell area.
 *
 * @param[in] lev            level of refinement
 */
void
REMORA::scale_rhs_vars (int lev)
{
    // Per level, not over all levels: /cff then *cff is not an identity in floating point,
    // so scaling every level on every level's step perturbs untouched ones.
    MultiFab& mf_cons = *cons_new[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real const> const& pm   = vec_pm[lev]->array(mfi);
        Array4<Real const> const& pn   = vec_pn[lev]->array(mfi);
        Array4<Real      > const& ru   = vec_ru[lev]->array(mfi);
        Array4<Real      > const& rv   = vec_rv[lev]->array(mfi);
        Array4<Real      > const& ru2d = vec_ru2d[lev]->array(mfi);
        Array4<Real      > const& rv2d = vec_rv2d[lev]->array(mfi);

        Box ubx = mfi.grownnodaltilebox(0,IntVect(NGROW,NGROW,0));
        Box vbx = mfi.grownnodaltilebox(1,IntVect(NGROW,NGROW,0));
        Box ubx2d = ubx; ubx2d.makeSlab(2,0);
        Box vbx2d = vbx; vbx2d.makeSlab(2,0);

        ParallelFor(ubx, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            Real cff = (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0));
            ru(i,j,k,n) = ru(i,j,k,n) / cff;
        });

        ParallelFor(vbx, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            Real cff = (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0));
            rv(i,j,k,n) = rv(i,j,k,n) / cff;
        });

        ParallelFor(ubx2d, 2, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
        {
            Real cff = (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0));
            ru2d(i,j,0,n) = ru2d(i,j,0,n) / cff;
        });

        ParallelFor(vbx2d, 2, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
        {
            Real cff = (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0));
            rv2d(i,j,0,n) = rv2d(i,j,0,n) / cff;
        });
    }
}

/**
 * Undo scale_rhs_vars on one level.
 *
 * @param[in] lev            level of refinement
 */
void
REMORA::scale_rhs_vars_inv (int lev)
{
    MultiFab& mf_cons = *cons_new[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real const> const& pm   = vec_pm[lev]->array(mfi);
        Array4<Real const> const& pn   = vec_pn[lev]->array(mfi);
        Array4<Real      > const& ru   = vec_ru[lev]->array(mfi);
        Array4<Real      > const& rv   = vec_rv[lev]->array(mfi);
        Array4<Real      > const& ru2d = vec_ru2d[lev]->array(mfi);
        Array4<Real      > const& rv2d = vec_rv2d[lev]->array(mfi);

        Box ubx = mfi.grownnodaltilebox(0,IntVect(NGROW,NGROW,0));
        Box vbx = mfi.grownnodaltilebox(1,IntVect(NGROW,NGROW,0));
        Box ubx2d = ubx; ubx2d.makeSlab(2,0);
        Box vbx2d = vbx; vbx2d.makeSlab(2,0);

        ParallelFor(ubx, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            Real cff = (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0));
            ru(i,j,k,n) = ru(i,j,k,n) * cff;
        });

        ParallelFor(vbx, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            Real cff = (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0));
            rv(i,j,k,n) = rv(i,j,k,n) * cff;
        });

        ParallelFor(ubx2d, 2, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
        {
            Real cff = (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0));
            ru2d(i,j,0,n) = ru2d(i,j,0,n) * cff;
        });

        ParallelFor(vbx2d, 2, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
        {
            Real cff = (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0));
            rv2d(i,j,0,n) = rv2d(i,j,0,n) * cff;
        });
    }
}
