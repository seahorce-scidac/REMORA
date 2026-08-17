    for ( MFIter mfi((mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& Akv = (mf_Akv).array(mfi);
        Array4<Real> const& Akt = (mf_Akt).array(mfi);
        Array4<const Real> const& z_w = remora.vec_z_w[lev]->array(mfi);
        Box bx = mfi.tilebox();
        bx.grow(IntVect(NGROW,NGROW,0));
        Gpu::streamSynchronize();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Akv(i,j,k) = Real(2.0e-03)+Real(8.0e-03)*std::exp(z_w(i,j,k)/Real(150.0));

            Akt(i,j,k,Temp_comp) = Real(1.0e-6);
            Akt(i,j,k,Salt_comp) = Real(1.0e-6);
        });
    }
