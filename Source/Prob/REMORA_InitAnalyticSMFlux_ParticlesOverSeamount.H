    //If we had wind stress and bottom stress we would need to set these:
    Real pi = Real(3.14159265359);
    Real tdays=remora.get_t_old(lev)/(Real(24.0)*Real(60.0)*Real(60.0));
    //this is a hack because time is off by dt. this needs to be fixed for non-fixed dt
    Real rho0=m_solverChoice.rho0;
    Real windamp;
    Real dstart = zero;
    //It's possible these should be set to be nonzero only at the boundaries they affect
    if(NSPeriodic) {
        mf_sustr.setVal(zero);
    }
    else if(EWPeriodic) {
        if ((tdays-dstart)<=two)
            windamp=Real(-0.1)*sin(pi*(tdays-dstart)/Real(4.0))/rho0;
        else
            windamp=Real(-0.1)/rho0;
        mf_sustr.setVal(windamp);
    }
    if(NSPeriodic) {
        if ((tdays-dstart)<=two)
            windamp=Real(-0.1)*sin(pi*(tdays-dstart)/Real(4.0))/rho0;
        else
            windamp=Real(-0.1)/rho0;
        mf_svstr.setVal(windamp);
    }
    else if(EWPeriodic) {
        mf_svstr.setVal(zero);
    }

