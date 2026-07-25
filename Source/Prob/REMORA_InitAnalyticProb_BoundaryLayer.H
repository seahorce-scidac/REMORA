
    auto T0 = m_solverChoice.T0;
    auto S0 = m_solverChoice.S0;

    auto geomdata = geom.data();

    const int khi = geomdata.Domain().bigEnd()[2];

    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
        AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

        Array4<      Real> const& state = mf_cons.array(mfi);
        Array4<      Real> const& x_vel = mf_xvel.array(mfi);
        Array4<      Real> const& y_vel = mf_yvel.array(mfi);

        Array4<const Real> const& z_r   = remora.vec_z_r[lev]->const_array(mfi);

        ParallelFor(grow(bx,IntVect(1,1,0)), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real z = z_r(i,j,k);
            Real val1 = std::tanh(Real(1.1) * z + Real(11.0));
            state(i,j,k,Temp_comp)=T0+Real(6.25)*val1;
            state(i,j,k,Salt_comp)=S0-Real(0.75)*val1;

            state(i, j, k, Tracer_comp) = zero;
        });

        const Box& xbx = surroundingNodes(bx,0);
        const Box& ybx = surroundingNodes(bx,1);
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            x_vel(i, j, k) = zero;
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            y_vel(i, j, k) = zero;
        });
    }
    Gpu::streamSynchronize();
