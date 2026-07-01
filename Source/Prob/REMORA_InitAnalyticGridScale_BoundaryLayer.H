
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
            Real dx = 0.5_rt * (4000.0_rt / (Lm + 1)) * (i+1) + 675.0_rt;
            pm(i,j,0) = 1.0_rt / dx;
        });
    }

    mf_pn.setVal(dxi[1]);
