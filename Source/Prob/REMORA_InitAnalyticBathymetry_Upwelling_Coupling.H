
    // Must not be periodic in E-W, and must allow bathymetry
    AMREX_ALWAYS_ASSERT(!EWPeriodic);
    AMREX_ALWAYS_ASSERT(!m_solverChoice.flat_bathymetry);

    mf_h.setVal(geomdata.ProbHi(2));

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& h  = (mf_h).array(mfi);

        Box bx = mfi.tilebox();
        Box gbx = bx;
        gbx.grow(mf_h.nGrowVect());

        Box gbxD = gbx;
        gbxD.makeSlab(2,0);

        ParallelFor(gbxD, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const auto prob_lo = geomdata.ProbLo();
            const auto dx = geomdata.CellSize();

            const Real x0 = prob_lo[0];
            const Real x  = x0 + (i + 0.5_rt) * dx[0];
            const Real x1 = geomdata.ProbHi(0);
            const Real h0 = 300.0_rt;
            const Real h1 = 10.0_rt;

            Real depth = h0 + (h1 - h0) * (x - x0) / (x1 - x0);
            depth = amrex::min(amrex::max(depth, h1), h0);
            h(i,j,0) = depth;
        });
    } // mfi
