#include <REMORA.H>

using namespace amrex;

void
REMORA::scale_rhs_vars ()
{
    for (int lev=0; lev<=finest_level;lev++) {
        MultiFab& mf_cons = *cons_new[lev];
        for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<Real const> const& pm   = vec_pm[lev]->array(mfi);
            Array4<Real const> const& pn   = vec_pn[lev]->array(mfi);
            Array4<Real      > const& ru   = vec_ru[lev]->array(mfi);
            Array4<Real      > const& rv   = vec_rv[lev]->array(mfi);
            Array4<Real      > const& ru2d = vec_ru2d[lev]->array(mfi);
            Array4<Real      > const& rv2d = vec_rv2d[lev]->array(mfi);

            ParallelFor(Box(ru), 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff = (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0));
                ru(i,j,k,n) = ru(i,j,k,n) / cff;
            });

            ParallelFor(Box(rv), 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff = (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0));
                rv(i,j,k,n) = rv(i,j,k,n) / cff;
            });

            ParallelFor(Box(ru2d), 2, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
            {
                Real cff = (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0));
                ru2d(i,j,0,n) = ru2d(i,j,0,n) / cff;
            });

            ParallelFor(Box(rv2d), 2, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
            {
                Real cff = (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0));
                rv2d(i,j,0,n) = rv2d(i,j,0,n) / cff;
            });
        }
    }
}

void
REMORA::scale_rhs_vars_inv ()
{
    for (int lev=0; lev<=finest_level;lev++) {
        MultiFab& mf_cons = *cons_new[lev];
        for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<Real const> const& pm   = vec_pm[lev]->array(mfi);
            Array4<Real const> const& pn   = vec_pn[lev]->array(mfi);
            Array4<Real      > const& ru   = vec_ru[lev]->array(mfi);
            Array4<Real      > const& rv   = vec_rv[lev]->array(mfi);
            Array4<Real      > const& ru2d = vec_ru2d[lev]->array(mfi);
            Array4<Real      > const& rv2d = vec_rv2d[lev]->array(mfi);

            ParallelFor(Box(ru), 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff = (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0));
                ru(i,j,k,n) = ru(i,j,k,n) * cff;
            });

            ParallelFor(Box(rv), 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff = (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0));
                rv(i,j,k,n) = rv(i,j,k,n) * cff;
            });

            ParallelFor(Box(ru2d), 2, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
            {
                Real cff = (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0));
                ru2d(i,j,0,n) = ru2d(i,j,0,n) * cff;
            });

            ParallelFor(Box(rv2d), 2, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
            {
                Real cff = (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0));
                rv2d(i,j,0,n) = rv2d(i,j,0,n) * cff;
            });
        }
    }

}
