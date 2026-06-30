
    mf_h.setVal(geomdata.ProbHi(2));

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& h  = (mf_h).array(mfi);

      Gpu::streamSynchronize();
      ParallelFor(Box(h), [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
          h(i,j,0,0) = -geomdata.ProbLo(2);
          if (k==0) {
              h(i,j,0,1) = h(i,j,0,0);
          }
      });
    } // mfi
