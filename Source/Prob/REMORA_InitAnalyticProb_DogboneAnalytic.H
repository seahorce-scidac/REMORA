
    // Whether this run carries a passive dye tracer. nComp() cannot answer this: with
    // biology active and no dye, Tracer_comp is the first biology component, so testing
    // nComp() would write the dye initial value into NO3.
    const bool l_do_dye = (remora.num_passive_scalars() > 0);
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
                    state(i, j, k, Salt_comp) = Real(35.0);
                }

                state(i,j,k,Temp_comp)= (x_r(i,j,0) < Real(1000.)) ? Real(5.0) : Real(10.0);

                // Set tracer = 0 everywhere
                if (l_do_dye) { state(i, j, k,Tracer_comp) = zero; }
            });
        } else {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                [[maybe_unused]] const Real x = x_r(i,j,0);

                if (l_use_salt) {
                    state(i, j, k, Salt_comp) = Real(10.0) - Real(10.0) * (x_r(i,j,0) - Real(8400.0)) / Real(8400.0);
                }

                state(i,j,k,Temp_comp)= Real(10.0);

                // Set tracer = 0 everywhere
                if (l_do_dye) { state(i, j, k,Tracer_comp) = zero; }
            });
        }

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
