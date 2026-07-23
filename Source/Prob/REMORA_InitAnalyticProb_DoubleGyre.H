
    bool l_use_salt = m_solverChoice.use_salt;

    auto geomdata = geom.data();
    const int khi = geomdata.Domain().bigEnd()[2];

    auto T0 = m_solverChoice.T0;
    Real val1 = (Real(44.69) / Real(39.382)) * (Real(44.69) / Real(39.382));
    Real val2 = val1 * (m_solverChoice.rho0 * Real(100.0)/g) * (Real(5.0e-5)/((Real(42.689)/Real(44.69)) * (Real(42.689)/Real(44.69))));

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

            const Real y  = prob_lo[1] + (j + Real(half)) * dx[1];
            const Real z = z_r(i,j,k);
            const Real yextent = prob_hi[1] - prob_lo[1];

            const Real val3 = T0 + val2 * std::exp(z/Real(100.0)) * (Real(10.0) - Real(0.4) * std::tanh(z / Real(100.0)));
            const Real val4 = y / yextent;

            state(i,j,k,Temp_comp)=val3 - Real(3.0) * val4;
            if (l_use_salt) {
                state(i,j,k,Salt_comp)=Real(34.5) - Real(0.001) * z - val4;
            }

            // Set scalar = 0 everywhere
            state(i, j, k, Tracer_comp) = zero;
        });

        const Box& xbx = surroundingNodes(bx,0);
        const Box& ybx = surroundingNodes(bx,1);
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            x_vel(i,j,k) = zero;
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            y_vel(i, j, k) = zero;
        });
    }
    Gpu::streamSynchronize();
