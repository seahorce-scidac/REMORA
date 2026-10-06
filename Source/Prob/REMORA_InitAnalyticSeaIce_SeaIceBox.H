    // Ice over the western 60% of the basin: concentration 0.9, 0.5 m thick
    // (0.45 m as mass per area), with 3 cm of snow, at rest. Open water
    // elsewhere, with ROMS's open-water ice temperature of -2 C.
    using namespace REMORASeaIce;
    const Real x_edge = Real(0.6) * (geom.ProbHi(0) - geom.ProbLo(0)) + geom.ProbLo(0);

    for (MFIter mfi(mf_ice_r, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& si = mf_ice_r.array(mfi);
        Array4<const Real> const& xr = remora.vec_xr[lev]->const_array(mfi);
        const Box gbx = mfi.growntilebox(IntVect(1,1,0));

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const bool ice = xr(i,j,0) < x_edge;
            si(i,j,0,rcomp(IceR::Aice,0)) = ice ? Real(0.9)  : zero;
            si(i,j,0,rcomp(IceR::Hice,0)) = ice ? Real(0.45) : zero;
            si(i,j,0,rcomp(IceR::Hsno,0)) = ice ? Real(0.03) : zero;
            si(i,j,0,rcomp(IceR::Tice,0)) = ice ? Real(-4.0) : Real(-2.0);
            si(i,j,0,rcomp(IceR::Iage,0)) = ice ? Real(10.0) : zero;
        });
    }
    mf_ice_u.setVal(zero);
    mf_ice_v.setVal(zero);
