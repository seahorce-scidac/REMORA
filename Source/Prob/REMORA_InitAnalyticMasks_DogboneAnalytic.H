    for (MFIter mfi(mf_mskr, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));

        Array4<      Real> const& mskr = mf_mskr.array(mfi);
        Array4<const Real> const& x_r  = remora.vec_xr[lev]->const_array(mfi);
        Array4<const Real> const& y_r  = remora.vec_yr[lev]->const_array(mfi);

        // mskr is 1 everywhere already
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
        {
            if ((x_r(i,j,0) > 2800.0_rt && x_r(i,j,0) < 5600.0_rt) &&
                (y_r(i,j,0) < 250.0_rt || y_r(i,j,0) > 500.0_rt)) {
                mskr(i,j,0) = 0.0_rt;
            }
        });
    }
