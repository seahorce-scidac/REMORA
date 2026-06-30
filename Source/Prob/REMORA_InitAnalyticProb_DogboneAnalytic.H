
    ParmParse pp_prob("remora.prob");
    bool traditional = true;
    pp_prob.query("traditional", traditional);

    bool l_use_salt = m_solverChoice.use_salt;

    auto geomdata = geom.data();
    const int khi = geomdata.Domain().bigEnd()[2];
    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
        AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

        Array4<      Real> const& state = mf_cons.array(mfi);
        Array4<      Real> const& x_vel = mf_xvel.array(mfi);
        Array4<      Real> const& y_vel = mf_yvel.array(mfi);

        Array4<const Real> const& x_r = remora.vec_xr[lev]->const_array(mfi);
        if (traditional) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (l_use_salt) {
                    state(i, j, k, Salt_comp) = 35.0_rt;
                }

                state(i,j,k,Temp_comp)= (x_r(i,j,0) < 1000._rt) ? 5.0_rt : 10.0_rt;

                // Set tracer = 0 everywhere
                state(i, j, k, Tracer_comp) = 0.0_rt;
            });
        } else {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                [[maybe_unused]] const Real x = x_r(i,j,0);

                if (l_use_salt) {
                    state(i, j, k, Salt_comp) = 10.0_rt - 10.0_rt * (x_r(i,j,0) - 8400.0_rt) / 8400.0_rt;
                }

                state(i,j,k,Temp_comp)= 10.0_rt;

                // Set tracer = 0 everywhere
                state(i, j, k, Tracer_comp) = 0.0_rt;
            });
        }

        const Box& xbx = surroundingNodes(bx,0);
        const Box& ybx = surroundingNodes(bx,1);
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            x_vel(i, j, k) = 0.0_rt;
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            y_vel(i, j, k) = 0.0_rt;
        });
    }
    Gpu::streamSynchronize();
