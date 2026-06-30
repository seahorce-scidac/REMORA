    mf_zeta.setVal(0.0_rt);

    for (MFIter mfi(mf_zeta, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.growntilebox(IntVect(1,1,0));

        Array4<      Real> const& zeta = mf_zeta.array(mfi);
        Array4<const Real> const& x_r  = remora.vec_xr[lev]->const_array(mfi);

        ParallelFor(bx, 3, [=] AMREX_GPU_DEVICE(int i, int j, int , int n) noexcept
        {
            if (x_r(i,j,0) < 1100.0_rt) {
                zeta(i,j,0,n) = -0.00125_rt * (x_r(i,j,0) - 100.0_rt) + 1.25_rt;
            }
        });
    }
