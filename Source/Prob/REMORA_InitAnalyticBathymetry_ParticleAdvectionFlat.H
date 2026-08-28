    // This problem is run over a flat bottom at the base of the domain. The tracer
    // particles are seeded through the full water column and there is no bed
    // reflection in REMORAPC::AdvectWithFlow, so a shoaling bottom would leave
    // particles below the sea floor, outside the terrain-following grid.
    mf_h.setVal(geomdata.ProbHi(2));

    // vec_h has two components, but vec_h_full_domain (used with hires_grid_level) has one
    const int h_ncomp = mf_h.nComp();

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& h  = (mf_h).array(mfi);

      Gpu::streamSynchronize();
      ParallelFor(Box(h), [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
          h(i,j,0,0) = -geomdata.ProbLo(2);
          if (k==0 && h_ncomp > 1) {
              h(i,j,0,1) = h(i,j,0,0);
          }
      });
    } // mfi
