    mf_zeta.setVal(zero);

    // Off leaves the free surface flat, which a test needs if its exact solution is rest.
    ParmParse pp_zeta("remora.prob");
    bool zeta_bump = true;
    pp_zeta.query("zeta_bump", zeta_bump);

    if (zeta_bump) {
        for (MFIter mfi(mf_zeta, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.growntilebox(IntVect(1,1,0));

            Array4<      Real> const& zeta = mf_zeta.array(mfi);
            Array4<const Real> const& x_r  = coords.x_r->const_array(mfi);

            ParallelFor(bx, mf_zeta.nComp(), [=] AMREX_GPU_DEVICE(int i, int j, int , int n) noexcept
            {
                if (x_r(i,j,0) < Real(1100.0)) {
                    zeta(i,j,0,n) = Real(-0.00125) * (x_r(i,j,0) - Real(100.0)) + Real(1.25);
                }
            });
        }
    }
