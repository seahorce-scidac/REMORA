    for ( MFIter mfi((mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& Akv = (mf_Akv).array(mfi);
      Array4<Real> const& Akt = (mf_Akt).array(mfi);
      Array4<const Real> const& z_w = remora.vec_z_w[lev]->const_array(mfi);
      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));
      Gpu::streamSynchronize();
      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
        Akv(i,j,k) = 0.0; //2.0e-03_rt+8.0e-03_rt*std::exp(z_w(i,j,k)/150.0_rt);

        Akt(i,j,k,Temp_comp) = 0.0; // 1.0e-6_rt;
        Akt(i,j,k,Salt_comp) = 0.0; // 1.0e-6_rt;
        Akt(i,j,k,Tracer_comp) = 0.0_rt;
      });
    }
