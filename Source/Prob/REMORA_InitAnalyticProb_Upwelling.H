
    // Whether this run carries a passive dye tracer. nComp() cannot answer this: with
    // biology active and no dye, Tracer_comp is the first biology component, so testing
    // nComp() would write the dye initial value into NO3.
    const bool l_do_dye = (remora.num_passive_scalars() > 0);
    bool l_use_salt = m_solverChoice.use_salt;

    auto geomdata = geom.data();
    const int khi = geomdata.Domain().bigEnd()[2];

    [[maybe_unused]] bool EWPeriodic = geomdata.isPeriodic(0);
    [[maybe_unused]] bool NSPeriodic = geomdata.isPeriodic(1);

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
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real z = z_r(i,j,k);

            state(i, j, k, Temp_comp) = one;

            state(i,j,k,Temp_comp)=T0+Real(8.0)*std::exp(z/Real(50.0));
            if (l_use_salt) {
                state(i,j,k,Salt_comp)=S0;
            }

            // Set tracer = 0 everywhere
            if (l_do_dye) { state(i, j, k,Tracer_comp) = zero; }
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
