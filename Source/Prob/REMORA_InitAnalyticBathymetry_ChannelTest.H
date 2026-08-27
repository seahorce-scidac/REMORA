
    // Must not be doubly periodic, and must have terrain
    AMREX_ALWAYS_ASSERT( !NSPeriodic || !EWPeriodic);

    // mf_h.setVal(geomdata.ProbHi(2));
    mf_h.setVal(Real(50.0));

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& h  = (mf_h).array(mfi);

        Box bx = mfi.tilebox();

        Box bxD = bx;
        bxD.makeSlab(2,0);

        Gpu::streamSynchronize();

        ParallelForRNG(bxD, [=] AMREX_GPU_DEVICE (int i, int j, int , auto /*engine*/)
        {
            h(i,j,0) = Real(50.0);
        });

    } // mfi
