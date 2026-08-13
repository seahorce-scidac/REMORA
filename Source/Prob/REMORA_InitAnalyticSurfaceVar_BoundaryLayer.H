    Real time = remora.get_t_old(lev);
    Real tdays = time/(Real(24.0)*Real(60.0)*Real(60.0));
    Real dstart = zero;

    mf_Uwind.setVal(zero);
    if (tdays - dstart <= Real(6.0)) {
        mf_Vwind.setVal(Real(10.0));
    } else {
        mf_Vwind.setVal(zero);
    }
