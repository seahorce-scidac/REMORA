    ParmParse pp_prob("remora.prob");

    bool traditional = true;
    pp_prob.query("traditional", traditional);

    if (traditional) {
        mf_h.setVal(Real(10.0));
    } else {
        mf_h.setVal(zero);
        for (MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));

            Array4<      Real> const& h = mf_h.array(mfi);
            Array4<const Real> const& x_r  = remora.vec_xr[lev]->const_array(mfi);
            Array4<const Real> const& y_r  = remora.vec_yr[lev]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
            {
                h(i,j,0) = Real(10.0) + (x_r(i,j,0) - Real(4200.0)) * (x_r(i,j,0) - Real(4200.0)) * Real(5.0) / (Real(4200.0) * Real(4200.0)) -
                            (y_r(i,j,0) - Real(375.0)) * (y_r(i,j,0) - Real(375.0)) * Real(5.0) / (Real(375.0) * Real(375.0));
            });
        }
    }
