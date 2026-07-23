
    const bool l_use_salt = m_solverChoice.use_salt;

    auto geomdata = geom.data();
    const int khi = geomdata.Domain().bigEnd()[2];

    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
    amrex::Print() << "Box is " << bx << "comparing to khi+1 " << khi+1 << std::endl;
        AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

        Array4<Real> const& state = mf_cons.array(mfi);
        Array4<Real> const& x_vel = mf_xvel.array(mfi);
        Array4<Real> const& y_vel = mf_yvel.array(mfi);

        Array4<const Real> const& z_r = remora.vec_z_r[lev]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Convert model z coordinate to positive depth below the surface.
            const Real depth = amrex::max(zero, -z_r(i, j, k));

            Real temp = Real(6.5);
            if (depth <= Real(20.0)) {
                temp = Real(14.0);
            } else if (depth <= Real(50.0)) {
                const Real alpha = (depth - Real(20.0)) / (Real(50.0) - Real(20.0));
                temp = Real(14.0) + alpha * (Real(10.0) - Real(14.0));
            } else if (depth <= Real(300.0)) {
                const Real alpha = (depth - Real(50.0)) / (Real(300.0) - Real(50.0));
                temp = Real(10.0) + alpha * (Real(6.5) - Real(10.0));
            }

            Real salt = Real(33.5);
            if (depth <= Real(50.0)) {
                salt = Real(33.5);
            } else if (depth <= Real(150.0)) {
                const Real alpha = (depth - Real(50.0)) / (Real(150.0) - Real(50.0));
                salt = Real(33.5) + alpha * (Real(34.3) - Real(33.5));
            } else {
                const Real alpha = amrex::min((depth - Real(150.0)) / (Real(300.0) - Real(150.0)), one);
                salt = Real(34.3) + alpha * (Real(34.6) - Real(34.3));
            }

            state(i, j, k, Temp_comp) = temp;
            if (l_use_salt) {
                state(i, j, k, Salt_comp) = salt;
            }

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
