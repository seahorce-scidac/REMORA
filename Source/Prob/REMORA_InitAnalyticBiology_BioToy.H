
    Real cff1=Real(20.0)/Real(3.0);
    Real cff2=Real( 2.0)/Real(3.0);
    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<const Real> const& temp = mf_cons.array(mfi,Temp_comp);
        Array4<      Real> const& cons = mf_cons.array(mfi);

        bx = mfi.growntilebox(mf_cons.nGrowVect());

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real SiO4;
            if (temp(i,j,k) < Real(8.0)) {
                SiO4 = Real(30.0);
            } else if ((temp(i,j,k) >= Real(8.0)) && (temp(i,j,k) <= Real(11.0))) {
                SiO4 = Real(30.0) - ((temp(i,j,k) - Real(8.0))*cff1);
            } else if ((temp(i,j,k) > Real(11.0)) && (temp(i,j,k) <= Real(13.0))) {
                SiO4 = Real(10.0) - ((temp(i,j,k) - Real(11.0))*Real(4.0));
            } else if ((temp(i,j,k) > Real(13.0)) && (temp(i,j,k) <= Real(16.0))) {
                SiO4 = Real(2.0) - ((temp(i,j,k) - Real(13.0))*cff2);
            } else if (temp(i,j,k) > Real(16.0)) {
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
                t(i,j,k,bio_comp.tic)=Real(2100.0);
                t(i,j,k,bio_comp.talk)=Real(2350.0);
                t(i,j,k,bio_comp.ldec)=Real(0.002);
                t(i,j,k,bio_comp.sdec)=Real(0.06);
            }
        });
