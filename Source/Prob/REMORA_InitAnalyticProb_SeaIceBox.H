    // Unstratified water just above its freezing point (-0.43 C at S = 8), at rest
    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& state = mf_cons.array(mfi);
        Array4<Real> const& x_vel = mf_xvel.array(mfi);
        Array4<Real> const& y_vel = mf_yvel.array(mfi);

        ParallelFor(grow(bx,IntVect(1,1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            state(i,j,k,Temp_comp) = Real(-0.35);
            state(i,j,k,Salt_comp) = Real(8.0);
        });
        ParallelFor(surroundingNodes(bx,0), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            x_vel(i,j,k) = zero;
        });
        ParallelFor(surroundingNodes(bx,1), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            y_vel(i,j,k) = zero;
        });
    }
