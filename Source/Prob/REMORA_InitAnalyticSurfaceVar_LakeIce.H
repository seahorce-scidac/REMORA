    // ROMS LAKE_ICE analytic forcing (roms_test lake_ice/Forward ana_*.h), in
    // the order set_data.F calls them: cloud, Tair, humidity, then shortwave,
    // which reads those three, then winds, rain and pressure.
    const double tdays = remora.model_time(remora.get_t_old(lev)) / 86400.0;
    const Real pi = Real(M_PI);
    const Real deg2rad = pi / Real(180.0);

    const Real cloud_val = Real(0.65);
    const Real tair_val = Real(5.0) + Real(15.0) *
        std::cos(Real(2.0) * pi * (Real(tdays) - Real(220.0)) / Real(365.25));
    const Real hair_val = Real(0.8);
    const Real pair_val = Real(1013.48);
    const Real rain_val = Real(2.45e-5) - Real(1.47e-5) * std::cos(Real(2.0) * pi * Real(tdays) / Real(365.25));
    const Real vwind_val = (Real(5.0) + Real(2.0) * std::cos(Real(tdays) * Real(2.0) * pi / Real(365.25))) *
                           std::cos(Real(tdays) * Real(0.5));

    // Shortwave (ALBEDO branch of ana_srflux.h) at 42N, 81W, reduced by the
    // water albedo. ROMS stores srflx in Celsius m/s; REMORA keeps W/m2.
    int year = 0;
    double yday = 0.0;
    remora_caldate(m_solverChoice.time_ref, tdays, year, yday);
    const double hour = (yday - std::trunc(yday)) * 24.0;
    const double Dangle = 23.44 * std::cos((172.0 - yday) * 2.0 * M_PI / 365.2425) * (M_PI / 180.0);
    const Real Hangle = (Real(12.0) - Real(hour)) * pi / Real(12.0);
    const Real Csolar = Real(1353.0);
    const Real alb_w = Real(0.06);
    const Real LatRad = Real(42.0) * deg2rad;
    const Real LonRad = Real(-81.0) * deg2rad;
    const Real cff1 = std::sin(LatRad) * std::sin(Real(Dangle));
    const Real cff2 = std::cos(LatRad) * std::cos(Real(Dangle));
    const Real zenith = cff1 + cff2 * std::cos(Hangle - LonRad);
    Real srflx_val = zero;
    if (zenith > zero) {
        const Real cff = (Real(0.7859) + Real(0.03477) * tair_val) / (one + Real(0.00412) * tair_val);
        const Real e_sat = std::pow(Real(10.0), cff);
        const Real vap_p = e_sat * hair_val;
        srflx_val = Csolar * zenith * zenith * (one - Real(0.6) * cloud_val * cloud_val * cloud_val) /
                    ((zenith + Real(2.7)) * vap_p * Real(1.0e-3) + Real(1.085) * zenith + Real(0.1));
    }
    srflx_val = (one - alb_w) * srflx_val;

    mf_cloud.setVal(cloud_val);
    mf_Tair.setVal(tair_val);
    mf_qair.setVal(hair_val);
    mf_srflx.setVal(srflx_val);
    mf_Uwind.setVal(zero);
    mf_Vwind.setVal(vwind_val);
    mf_rain.setVal(rain_val);
    mf_Pair.setVal(pair_val);
