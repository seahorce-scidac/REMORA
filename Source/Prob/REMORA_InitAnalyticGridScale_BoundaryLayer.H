    // Cell width grows linearly across the domain in x; pn is uniform. Taken as the
    // difference of a face-position mapping X(xi) = c1*xi + c2*xi^2, xi = i / n_cell_x, so
    // the widths telescope and a level's cells subdivide its parent's exactly. Setting the
    // width from the cell index instead would give a refined level cells the physical size
    // of its parent's, and no conservative flux at the interface (see check_cf_metrics).
    //
    // The constants reproduce dx(i) = first_width + width_growth*i at 39 cells, exactly.
    constexpr Real first_width  = Real(725.0);
    constexpr Real width_growth = Real(50.0);
    constexpr Real n_base       = Real(39.0);

    constexpr Real c1 = n_base * (first_width - Real(0.5) * width_growth);
    constexpr Real c2 = Real(0.5) * width_growth * n_base * n_base;

    const int n_cell_x = geom.Domain().size()[0];
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
            const Real xi_lo = Real(i  ) / Real(n_cell_x);
            const Real xi_hi = Real(i+1) / Real(n_cell_x);
            const Real dx = (c1 * xi_hi + c2 * xi_hi * xi_hi)
                          - (c1 * xi_lo + c2 * xi_lo * xi_lo);
            pm(i,j,0) = one / dx;
        });
    }

    mf_pn.setVal(dxi[1]);
