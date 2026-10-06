    // Coordinates come from the Geometry rather than remora.vec_xr/vec_yr so that this hook
    // also works on the full-domain MultiFab used for hires_grid_level, which is
    // not built on grids[lev] and so has no per-level coordinate arrays to index. This
    // matches set_grid_coords_from_grid_scale, which puts the rho point of cell i at
    // (i + 1/2) * dx measured from index 0.
    const Real dx = geom.CellSize(0);
    const Real dy = geom.CellSize(1);

    // Defaults put the coast on coarse cell faces at every refinement ratio the tests use,
    // so every block is entirely land or entirely water. Moving y_lo/y_hi off a coarse face
    // is what produces a partially-masked coarse cell.
    ParmParse pp_mask("remora.prob");
    Real mask_x_lo = Real(2800.0);
    Real mask_x_hi = Real(5600.0);
    Real mask_y_lo = Real(250.0);
    Real mask_y_hi = Real(500.0);
    pp_mask.query("mask_x_lo", mask_x_lo);
    pp_mask.query("mask_x_hi", mask_x_hi);
    pp_mask.query("mask_y_lo", mask_y_lo);
    pp_mask.query("mask_y_hi", mask_y_hi);

    // A ROMS closed basin has land on its boundary ring of rho points, next to interior
    // water. mask_ring makes the grow cells outside the non-periodic domain edges land,
    // so the wall boundary conditions meet land beside water.
    bool mask_ring = false;
    pp_mask.query("mask_ring", mask_ring);
    const Box& dom = geom.Domain();
    const auto dlo = amrex::lbound(dom);
    const auto dhi = amrex::ubound(dom);
    const bool ring_x = mask_ring && !geom.isPeriodic(0);
    const bool ring_y = mask_ring && !geom.isPeriodic(1);

    for (MFIter mfi(mf_mskr, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Every grow cell, not just NGROW of them: the full-domain array carries as many as
        // the cumulative refinement ratio needs, and any it left unwritten would coarsen
        // down as water.
        const Box &bx = mfi.growntilebox();

        Array4<      Real> const& mskr = mf_mskr.array(mfi);

        // mskr is 1 everywhere already
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
        {
            const Real x_r = (i + Real(0.5)) * dx;
            const Real y_r = (j + Real(0.5)) * dy;
            if ((x_r > mask_x_lo && x_r < mask_x_hi) &&
                (y_r < mask_y_lo || y_r > mask_y_hi)) {
                mskr(i,j,0) = zero;
            }
            if ((ring_x && (i < dlo.x || i > dhi.x)) ||
                (ring_y && (j < dlo.y || j > dhi.y))) {
                mskr(i,j,0) = zero;
            }
        });
    }
