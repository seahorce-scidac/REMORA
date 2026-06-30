
    // Must not be doubly periodic, and must have terrain
    AMREX_ALWAYS_ASSERT( !NSPeriodic || !EWPeriodic);
    AMREX_ALWAYS_ASSERT( !m_solverChoice.flat_bathymetry);

    mf_h.setVal(geomdata.ProbHi(2));

    [[maybe_unused]] const int Lm = geomdata.Domain().size()[0];
    [[maybe_unused]] const int Mm = geomdata.Domain().size()[1];

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& h  = (mf_h).array(mfi);
      Array4<const Real> const& xr  = remora.vec_xr[lev]->array(mfi);
      [[maybe_unused]] Array4<const Real> const& yr  = remora.vec_yr[lev]->array(mfi);

      Box bx = mfi.tilebox();
      Box gbx2 = bx;
      gbx2.grow(IntVect(NGROW,NGROW,0));

      Box gbx2D = gbx2;
      gbx2D.makeSlab(2,0);

      Gpu::streamSynchronize();

      ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int )
      {
          Real val1 = (xr(i,j,0) + 500.0_rt) / 15000.0_rt;
          h(i,j,0) = 14.0_rt + 25.0_rt * (1.0_rt - std::exp(-PI * xr(i,j,0) * 1.0e-5_rt)) - 8.0_rt * std::exp(-val1 * val1);
      });

    } // mfi
