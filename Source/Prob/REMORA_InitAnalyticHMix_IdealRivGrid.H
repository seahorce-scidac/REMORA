    for ( MFIter mfi((mf_visc2_p), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& visc2_p_arr = mf_visc2_p.array(mfi);
      Array4<Real> const& visc2_r_arr = mf_visc2_r.array(mfi);
      Array4<Real> const& diff2_arr   = mf_diff2.array(mfi);

      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));

      int ncomp = mf_diff2.nComp(); // temperature and salt and tracer
      Gpu::streamSynchronize();

      ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
      {
          visc2_p_arr(i,j,0) = Real(5.0);
          visc2_r_arr(i,j,0) = Real(5.0);

          for (int n = 0; n < ncomp; n++) {
              diff2_arr(i,j,0,n) = Real(5.0e3);
          }
      });
    } // mfi
