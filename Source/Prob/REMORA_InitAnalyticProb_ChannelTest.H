
    bool l_use_salt = m_solverChoice.use_salt;

    auto geomdata = geom.data();
    const int khi = geomdata.Domain().bigEnd()[2];

    [[maybe_unused]] bool EWPeriodic = geomdata.isPeriodic(1);
    [[maybe_unused]] bool NSPeriodic = geomdata.isPeriodic(0);

    auto T0 = m_solverChoice.T0;
    auto S0 = m_solverChoice.S0;
    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
        AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

        Array4<      Real> const& state = mf_cons.array(mfi);
        Array4<      Real> const& x_vel = mf_xvel.array(mfi);
        Array4<      Real> const& y_vel = mf_yvel.array(mfi);

        Array4<const Real> const& z_r = remora.vec_z_r[lev]->const_array(mfi);

        ParallelFor(grow(grow(bx,0,1),1,1), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const auto prob_lo         = geomdata.ProbLo();
            const auto prob_hi         = geomdata.ProbHi();
            const auto dx              = geomdata.CellSize();

            const Real z = z_r(i,j,k);

            // Create bounding box for x and y to make spatially-dependent T and S
            const Real xcent = Real(0.5)*(prob_lo[0] + prob_hi[0]);
            [[maybe_unused]] const Real ycent = Real(0.5)*(prob_lo[1] + prob_hi[1]);

            [[maybe_unused]] const Real x  = prob_lo[0] + (i + Real(0.5)) * dx[0] - xcent;
            const Real y  = prob_lo[1] + (j + Real(0.5)) * dx[1];

            state(i,j,k,Temp_comp)=T0 + z / (Real(9.8) * Real(1.7));
            if (l_use_salt) {
                state(i,j,k,Salt_comp)=S0 + y / (Real(9.8) * Real(760.0));
            }

            // Set tracer = 0 everywhere
            state(i, j, k, Tracer_comp) = zero;
        });

        const Box& xbx = surroundingNodes(bx,0);
        const Box& ybx = surroundingNodes(bx,1);
        ParallelFor(grow(grow(xbx,1,1),0,1), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Create bounding box for x and y to make spatially-dependent T and S
            const auto prob_lo         = geomdata.ProbLo();
            const auto prob_hi         = geomdata.ProbHi();
            const auto dx              = geomdata.CellSize();

            [[maybe_unused]] const Real xcent = Real(0.5)*(prob_lo[0] + prob_hi[0]);
            [[maybe_unused]] const Real ycent = Real(0.5)*(prob_lo[1] + prob_hi[1]);

            const Real z = z_r(i,j,k);
            // h(i,j,0) = -geomdata.ProbLo(2);
            const Real depth = Real(50.0);

            const Real x  = prob_lo[0] + (i + Real(0.5)) * dx[0];// - xcent;
            const Real y  = prob_lo[1] + (j + Real(0.5)) * dx[1];// - ycent;

            const Real val1  = x * Real(6.28318530718) * Real(1.0e-5);
            const Real val2  = y * Real(6.28318530718) * Real(1.0e-5);
            x_vel(i, j, k) = (depth + z) * Real(0.01) + Real(0.05) * std::sin(val1) * std::sin(val2);
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            y_vel(i, j, k) = zero;
        });
    }
    Gpu::streamSynchronize();
