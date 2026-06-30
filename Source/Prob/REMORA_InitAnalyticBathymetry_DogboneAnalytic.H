    ParmParse pp_prob("remora.prob");

    bool traditional = true;
    pp_prob.query("traditional", traditional);

    if (traditional) {
        mf_h.setVal(10.0_rt);
    } else {
        mf_h.setVal(0.0_rt);
        for (MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));

            Array4<      Real> const& h = mf_h.array(mfi);
            Array4<const Real> const& x_r  = remora.vec_xr[lev]->const_array(mfi);
            Array4<const Real> const& y_r  = remora.vec_yr[lev]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
            {
                h(i,j,0) = 10.0_rt + (x_r(i,j,0) - 4200._rt) * (x_r(i,j,0) - 4200._rt) * 5.0_rt / (4200.0_rt * 4200.0_rt) -
                            (y_r(i,j,0) - 375._rt) * (y_r(i,j,0) - 375._rt) * 5.0_rt / (375._rt * 375._rt);
            });
        }
    }
