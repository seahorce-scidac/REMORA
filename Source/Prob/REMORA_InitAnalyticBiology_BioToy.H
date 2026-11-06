// Analytic biology profile for BioToy and Upwelling, from ROMS ana_biology.h.
//
// A code fragment included from Problem::init_analytic_biology in
// Exec/REMORA_Prob.cpp (and from the biology unit test). It needs only
// `biology` (REMORABiology::Biology const&) and `context`
// (REMORABiology::BiologyInitContext&) in scope, and locates its fields by
// tracer name, so it never sees FennelParameters. Optional Fennel pools that
// the run does not carry have index -1 and are skipped; every write is
// guarded by its own index.
//
// The values are Fennel's, in Fennel's units, so the profile refuses any
// other model even if it reuses the same tracer names.
//
// The profile is a function of temperature alone, so it carries over unchanged
// to any problem that stratifies temperature.
{
    if (biology.name() != "fennel") {
        amrex::Error("The BioToy analytic biology profile is Fennel's;"
                     " remora.biology_model = " + biology.name() +
                     " needs its own profile or"
                     " remora.biology_ic_type = constant or netcdf.");
    }
    auto& mf_cons = context.state;
    auto const& bio_declaration = biology.declaration();
    auto comp = [&] (const std::string& name) {
        return bio_declaration.component_index(name, context.component_start);
    };
    const int no3 = comp("NO3");
    const int nh4 = comp("NH4");
    const int chlo = comp("chlorophyll");
    const int phyt = comp("phytoplankton");
    const int zoop = comp("zooplankton");
    const int lden = comp("LdetritusN");
    const int sden = comp("SdetritusN");
    const int rden = comp("RdetritusN");
    const int po4 = comp("PO4");
    const int ldec = comp("LdetritusC");
    const int sdec = comp("SdetritusC");
    const int tic = comp("TIC");
    const int talk = comp("alkalinity");
    const int rdec = comp("RdetritusC");
    const int oxyg = comp("oxygen");
    const int odu = comp("ODU");
    // Fennel always declares the seven nitrogen tracers.
    AMREX_ALWAYS_ASSERT(no3 >= 0 && nh4 >= 0 && chlo >= 0 && phyt >= 0 &&
                        zoop >= 0 && lden >= 0 && sden >= 0);

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
            cons(i,j,k,no3)=Real(1.67)+Real(0.5873)*SiO4+
                                     Real(0.0144)*SiO4*SiO4+
                                     Real(0.0003099)*SiO4*SiO4*SiO4;
            cons(i,j,k,phyt)=Real(0.08);
            cons(i,j,k,zoop)=Real(0.06);
            cons(i,j,k,nh4)=Real(0.1);
            cons(i,j,k,lden)=Real(0.02);
            cons(i,j,k,sden)=Real(0.04);
            cons(i,j,k,chlo)=Real(0.02);
            if (oxyg >= 0) {
                cons(i,j,k,oxyg)=Real(10.0)/Real(0.02241);
            }
            if (tic >= 0) {
                cons(i,j,k,tic)=Real(2100.0);
            }
            if (talk >= 0) {
                cons(i,j,k,talk)=Real(2350.0);
            }
            if (ldec >= 0) {
                cons(i,j,k,ldec)=Real(0.002);
            }
            if (sdec >= 0) {
                cons(i,j,k,sdec)=Real(0.06);
            }
            if (rden >= 0) {
                // ROMS ana_biology.h has no river DON profile, so pick values
                // that make the pool visible: 0.05 mmol_N/m3 is the same order
                // as SdetritusN above, so its remineralization to NH4 is a
                // real contribution rather than a rounding artifact, and a
                // zero start would leave the whole option a no-op no matter
                // whether the terms are right.
                cons(i,j,k,rden)=Real(0.05);
            }
            if (rdec >= 0) {
                cons(i,j,k,rdec)=Real(0.3);
            }
            if (po4 >= 0) {
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
                cons(i,j,k,po4)=Real(0.5);
            }
            if (odu >= 0) {
                // ODU is inert on both paths: ROMS fennel.h never references
                // iODU_ (zero occurrences), and the native kernel only clamps
                // it on entry and writes it back on exit. This value exists so
                // the tracer-count, index-mapping and pack/unpack comparison
                // is not vacuous -- initialising it to zero would compare
                // 0 against 0 and prove nothing about the plumbing.
                cons(i,j,k,odu)=Real(1.0);
            }
        });
    }
}
