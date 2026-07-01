
    bool l_use_salt = m_solverChoice.use_salt;

    auto geomdata = geom.data();
    const int khi = geomdata.Domain().bigEnd()[2];

    auto T0 = m_solverChoice.T0;
    Real val1 = (44.69_rt / 39.382_rt) * (44.69_rt / 39.382_rt);
    Real val2 = val1 * (m_solverChoice.rho0 * 100.0_rt/g) * (5.0e-5_rt/((42.689_rt/44.69_rt) * (42.689_rt/44.69_rt)));

    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
        AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

        Array4<      Real> const& state = mf_cons.array(mfi);
        Array4<      Real> const& x_vel = mf_xvel.array(mfi);
        Array4<      Real> const& y_vel = mf_yvel.array(mfi);

        Array4<const Real> const& z_r = remora.vec_z_r[lev]->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const auto prob_lo         = geomdata.ProbLo();
            const auto prob_hi         = geomdata.ProbHi();
            const auto dx              = geomdata.CellSize();

            const Real y  = prob_lo[1] + (j + 0.5_rt) * dx[1];
            const Real z = z_r(i,j,k);
            const Real yextent = prob_hi[1] - prob_lo[1];

            const Real val3 = T0 + val2 * std::exp(z/100.0_rt) * (10.0_rt - 0.4_rt * std::tanh(z / 100.0_rt));
            const Real val4 = y / yextent;

            state(i,j,k,Temp_comp)=val3 - 3.0_rt * val4;
            if (l_use_salt) {
                state(i,j,k,Salt_comp)=34.5_rt - 0.001_rt * z - val4;
            }

            // Set scalar = 0 everywhere
            state(i, j, k, Tracer_comp) = 0.0_rt;
        });

        const Box& xbx = surroundingNodes(bx,0);
        const Box& ybx = surroundingNodes(bx,1);
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            x_vel(i,j,k) = 0.0_rt;
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            y_vel(i, j, k) = 0.0_rt;
        });
    }
    Gpu::streamSynchronize();
