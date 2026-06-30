    for ( MFIter mfi((mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& Akv = (mf_Akv).array(mfi);
      Array4<Real> const& Akt = (mf_Akt).array(mfi);
      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));
      Gpu::streamSynchronize();
      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
        Akv(i,j,k) = 1.0_rt;

        Akt(i,j,k,Temp_comp) = 1.0_rt;
        Akt(i,j,k,Salt_comp) = 1.0_rt;
        Akt(i,j,k,Tracer_comp) = 0.0_rt;
      });
    }
