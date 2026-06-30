
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
            const Real depth = amrex::max(0.0_rt, -z_r(i, j, k));

            Real temp = 6.5_rt;
            if (depth <= 20.0_rt) {
                temp = 14.0_rt;
            } else if (depth <= 50.0_rt) {
                const Real alpha = (depth - 20.0_rt) / (50.0_rt - 20.0_rt);
                temp = 14.0_rt + alpha * (10.0_rt - 14.0_rt);
            } else if (depth <= 300.0_rt) {
                const Real alpha = (depth - 50.0_rt) / (300.0_rt - 50.0_rt);
                temp = 10.0_rt + alpha * (6.5_rt - 10.0_rt);
            }

            Real salt = 33.5_rt;
            if (depth <= 50.0_rt) {
                salt = 33.5_rt;
            } else if (depth <= 150.0_rt) {
                const Real alpha = (depth - 50.0_rt) / (150.0_rt - 50.0_rt);
                salt = 33.5_rt + alpha * (34.3_rt - 33.5_rt);
            } else {
                const Real alpha = amrex::min((depth - 150.0_rt) / (300.0_rt - 150.0_rt), 1.0_rt);
                salt = 34.3_rt + alpha * (34.6_rt - 34.3_rt);
            }

            state(i, j, k, Temp_comp) = temp;
            if (l_use_salt) {
                state(i, j, k, Salt_comp) = salt;
            }

            state(i, j, k, Tracer_comp) = 0.0_rt;
        });

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
