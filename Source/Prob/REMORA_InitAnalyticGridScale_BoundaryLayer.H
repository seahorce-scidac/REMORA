
    auto geomdata = geom.data();
    const int Lm = geomdata.Domain().size()[0];
    const auto dxi = geom.InvCellSize();
    for ( MFIter mfi(mf_pm, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& pm = (mf_pm).array(mfi);

        Box bx = mfi.tilebox();
        Box gbx2 = bx;
        gbx2.grow(IntVect(NGROW+1,NGROW+1,0));

        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);

        Gpu::streamSynchronize();
        ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real dx = half * (Real(4000.0) / (Lm + 1)) * (i+1) + Real(675.0);
            pm(i,j,0) = one / dx;
        });
    }

    mf_pn.setVal(dxi[1]);
