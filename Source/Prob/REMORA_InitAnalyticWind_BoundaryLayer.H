    Real time = remora.get_t_old(lev);
    Real tdays = time/Real(24.0*60.0*60.0);
    Real dstart = 0.0_rt;

    mf_Uwind.setVal(0.0_rt);
    if (tdays - dstart <= 6.0_rt) {
        mf_Vwind.setVal(10.0_rt);
    } else {
        mf_Vwind.setVal(0.0_rt);
    }
