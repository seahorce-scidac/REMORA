    // One Akt for every passive (dye) scalar, so they all diffuse alike
    const int l_nscalar = remora.num_scalar();
    for ( MFIter mfi((mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& Akv = (mf_Akv).array(mfi);
        Array4<Real> const& Akt = (mf_Akt).array(mfi);

        Box bx = mfi.tilebox();
        bx.grow(IntVect(NGROW,NGROW,0));

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Akv(i,j,k) = Real(1.0e-5);

            Akt(i,j,k,Temp_comp) = Real(1.0e-6);
            Akt(i,j,k,Salt_comp) = Real(1.0e-6);
            for (int n = Tracer_comp; n < Tracer_comp + l_nscalar; ++n) { Akt(i,j,k,n) = zero; }
        });
    }
