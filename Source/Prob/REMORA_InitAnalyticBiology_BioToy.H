
    Real cff1=Real(20.0)/Real(3.0);
    Real cff2=Real( 2.0)/Real(3.0);
    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<const Real> const& temp = mf_cons.array(mfi,Temp_comp);
        Array4<      Real> const& cons = mf_cons.array(mfi);

        Box bx = mfi.growntilebox(mf_cons.nGrowVect());

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real SiO4;
            if (temp(i,j,k) < Real(8.0)) {
                SiO4 = Real(30.0);
            } else if ((temp(i,j,k) >= Real(8.0)) && (temp(i,j,k) <= Real(11.0))) {
                SiO4 = Real(30.0) - ((temp(i,j,k) - Real(8.0))*cff1);
            } else if ((temp(i,j,k) > Real(11.0)) && (temp(i,j,k) <= Real(13.0))) {
                SiO4 = Real(10.0) - ((temp(i,j,k) - Real(11.0))*Real(4.0));
            } else if ((temp(i,j,k) > Real(13.0)) && (temp(i,j,k) <= Real(16.0))) {
                SiO4 = Real(2.0) - ((temp(i,j,k) - Real(13.0))*cff2);
            } else { //else if (temp(i,j,k) > Real(16.0))
                SiO4 = Real(0.0);
            }
            cons(i,j,k,bio_comp.no3)=Real(1.67)+Real(0.5873)*SiO4+
                                     Real(0.0144)*SiO4*SiO4+
                                     Real(0.0003099)*SiO4*SiO4*SiO4;
            cons(i,j,k,bio_comp.phyt)=Real(0.08);
            cons(i,j,k,bio_comp.zoop)=Real(0.06);
            cons(i,j,k,bio_comp.nh4)=Real(0.1);
            cons(i,j,k,bio_comp.lden)=Real(0.02);
            cons(i,j,k,bio_comp.sden)=Real(0.04);
            cons(i,j,k,bio_comp.chlo)=Real(0.02);
            if (m_fennel_params.oxygen) {
                cons(i,j,k,bio_comp.oxyg)=Real(10.0)/Real(0.02241);
            }
            if (m_fennel_params.carbon) {
                cons(i,j,k,bio_comp.tic)=Real(2100.0);
                cons(i,j,k,bio_comp.talk)=Real(2350.0);
                cons(i,j,k,bio_comp.ldec)=Real(0.002);
                cons(i,j,k,bio_comp.sdec)=Real(0.06);
            }
            if (m_fennel_params.river_don) {
                // ROMS ana_biology.h has no river DON profile, so pick values
                // that make the pool visible: 0.05 mmol_N/m3 is the same order
                // as SdetritusN above, so its remineralization to NH4 is a
                // real contribution rather than a rounding artifact, and a
                // zero start would leave the whole option a no-op no matter
                // whether the terms are right.
                cons(i,j,k,bio_comp.rden)=Real(0.05);
                if (m_fennel_params.carbon) {
                    cons(i,j,k,bio_comp.rdec)=Real(0.3);
                }
            }
            if (m_fennel_params.po4) {
                // 0.5 mmol_P/m3: a plausible coastal subsurface value, far
                // above the 1e-6 MinVal floor that would otherwise mask the
                // uptake terms. With the default K_PO4=32 this gives
                // L_PO4 ~ 0.94, deliberately the same order as LTOT so that
                // MIN(LTOT, L_PO4) is a real comparison rather than a
                // foregone conclusion.
                //
                // Note the PO4 branch of biology_tile is structurally
                // different from the default one (mu-based cff4/cff5, and
                // LMIN rather than LTOT in the chlorophyll denominator), and
                // it is taken for any non-zero PO4 regardless of which side
                // the MIN selects. Which side it selects does depend on the
                // NO3/NH4 field, so a run that changes those should re-check
                // that both sides still get exercised.
                cons(i,j,k,bio_comp.po4)=Real(0.5);
            }
            if (m_fennel_params.odu) {
                // ODU is inert on both paths: ROMS fennel.h never references
                // iODU_ (zero occurrences), and the native kernel only clamps
                // it on entry and writes it back on exit. This value exists so
                // the tracer-count, index-mapping and pack/unpack comparison
                // is not vacuous -- initialising it to zero would compare
                // 0 against 0 and prove nothing about the plumbing.
                cons(i,j,k,bio_comp.odu)=Real(1.0);
            }
        });
    }
