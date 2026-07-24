#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] lev            level to operate on
 * @param[in   ] mf_cons        scalar data: temperature, salinity, passsive scalar, etc
 * @param[in   ] mf_uwind       u-direction wind dvelocity
 * @param[in   ] mf_vwind       v-direction wind dvelocity
 * @param[in   ] mf_Tair        air temperature [°C]
 * @param[in   ] mf_qair        specific humidity [kg/kg]
 * @param[in   ] mf_Pair        air pressure [mb]
 * @param[in   ] mf_srflx       shortwave radiation flux [W/m²]
 * @param[in   ] mf_longwave_down external longwave radiation flux [W/m²]
 * @param[inout] mf_evap        evaporation rate
 * @param[  out] mf_sustr       u-direction surface momentum stress
 * @param[  out] mf_svstr       v-direction surface momentum stress
 * @param[  out] mf_stflux      surface scalar flux (temperature, salinity)
 * @param[  out] mf_lrflx       longwave radiation flux
 * @param[inout] mf_lhflx       latent heat flux
 * @param[inout] mf_shflx       sensible heat flux
 * @param[in   ] N              number of vertical levels
 */
void
REMORA::bulk_fluxes (int lev, MultiFab* mf_cons, MultiFab* mf_uwind, MultiFab* mf_vwind,
                     MultiFab* mf_Tair, MultiFab* mf_qair, MultiFab* mf_Pair,
                     MultiFab* mf_srflx,
                     MultiFab* mf_longwave_down,
                     MultiFab* mf_evap, MultiFab* mf_sustr, MultiFab* mf_svstr,
                     MultiFab* mf_stflux, MultiFab* mf_lrflx, MultiFab* mf_lhflx,
                     MultiFab* mf_shflx,
                     const int N)
{
    BL_PROFILE("REMORA::bulk_fluxes()");
    const int IterMax = 3;
    const BoxArray& ba = mf_cons->boxArray();
    const DistributionMapping& dm = mf_cons->DistributionMap();
    MultiFab mf_Taux(ba, dm, 1, IntVect(NGROW,NGROW,0));
    MultiFab mf_Tauy(ba, dm, 1, IntVect(NGROW,NGROW,0));

    // temps: Taux, Tauy,
    for ( MFIter mfi(*mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real const> const& uwind = mf_uwind->const_array(mfi);
        Array4<Real const> const& vwind = mf_vwind->const_array(mfi);
        Array4<Real const> const& Tair_arr = mf_Tair->const_array(mfi);
        Array4<Real const> const& qair_arr = mf_qair->const_array(mfi);
        Array4<Real const> const& Pair_arr = mf_Pair->const_array(mfi);
        Array4<Real const> const& srflx_arr = mf_srflx->const_array(mfi);
        Array4<Real const> longwave_down_arr;
        if (mf_longwave_down != nullptr) {
            longwave_down_arr = mf_longwave_down->const_array(mfi);
        }
        Array4<Real const> const& cons = mf_cons->const_array(mfi);
        Array4<Real> const& sustr = mf_sustr->array(mfi);
        Array4<Real> const& svstr = mf_svstr->array(mfi);
        Array4<Real> const& stflux = mf_stflux->array(mfi);
        Array4<Real> const& lrflx = mf_lrflx->array(mfi);
        Array4<Real> const& lhflx = mf_lhflx->array(mfi);
        Array4<Real> const& shflx = mf_shflx->array(mfi);
        Array4<Real> const& evap  = mf_evap->array(mfi);
        Array4<Real> const& Taux = mf_Taux.array(mfi);
        Array4<Real> const& Tauy = mf_Tauy.array(mfi);

        Array4<const Real> const& mskr = vec_mskr[lev]->const_array(mfi);
        Array4<const Real> const& msku  = vec_msku[lev]->const_array(mfi);
        Array4<const Real> const& mskv  = vec_mskv[lev]->const_array(mfi);
        Array4<const Real> const& rain  = vec_rain[lev]->const_array(mfi);
        Array4<const Real> const& EminusP = vec_EminusP[lev]->const_array(mfi);
        Array4<const Real> const& cloud_arr = vec_cloud[lev]->const_array(mfi);

        Real Hscale = solverChoice.rho0 * Cp;
        Real Hscale2 = one / (solverChoice.rho0 * Cp);
        Real blk_ZQ = solverChoice.blk_ZQ;
        Real blk_ZT = solverChoice.blk_ZT;
        Real blk_ZW = solverChoice.blk_ZW;

        bool use_longwave_down = solverChoice.longwave_down;
        bool longwave_is_net = solverChoice.longwave_is_net;
        bool have_external_longwave = (mf_longwave_down != nullptr);
        bool use_EminusP_from_input = solverChoice.eminusp &&
            solverChoice.bulk_flux_type[BulkFlux::EminusP] != BulkForcingType::computed;

        Real eps = Real(1e-20);

        Box bx = mfi.tilebox();
        Box ubx = mfi.grownnodaltilebox(0,IntVect(NGROW-1,NGROW-1,0));
        Box vbx = mfi.grownnodaltilebox(1,IntVect(NGROW-1,NGROW-1,0));
        Box gbx1 = bx; gbx1.grow(IntVect(NGROW,NGROW,0));

        ParallelFor(makeSlab(gbx1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int ) {
            // Get spatially-varying atmospheric forcing from input arrays
            Real PairM = Pair_arr(i,j,0);  // Air pressure [mb]
            Real TairC = Tair_arr(i,j,0);  // Air temperature [°C]
            Real TairK = TairC + Real(273.16); // Air temperature [K]
            Real Hair = qair_arr(i,j,0);   // Specific humidity [kg/kg] or RH [fraction]
            Real RH = Hair;
            Real srflux = srflx_arr(i,j,0); // Shortwave radiation flux [W/m²]
            Real cloud = cloud_arr(i,j,0);  // Cloud cover fraction [0-1]

            // Input bulk parametrization fields
            Real wind_mag = std::sqrt(uwind(i,j,0)*uwind(i,j,0) + vwind(i,j,0) * vwind(i,j,0)) + eps;
            Real TseaK = cons(i,j,N,Temp_comp) + Real(273.16);

            // Initialize
            Real delTc = zero;
            Real delQc = zero;
            Real cff = zero;

            Real LHeat = lhflx(i,j,0) * Hscale;
            Real SHeat = shflx(i,j,0) * Hscale;
            Real Taur = zero;
            Taux(i,j,0) = zero;
            Tauy(i,j,0) = zero;
            Real LRad;

            /*-----------------------------------------------------------------------
               Compute outward or net longwave radiation (W/m2), LRad.
             -----------------------------------------------------------------------
               If external longwave radiation is net, use it directly. If it is
               downward, compute net longwave radiation as Ldown - Lemit, where
               Lemit is computed from the model SST and an emissivity. Otherwise,
               use Berliand (1952) formula to calculate net longwave radiation.
               The equation for saturation vapor pressure is from Gill (Atmosphere-
               Ocean Dynamics, pp 606). Here the coefficient in the cloud term
               is assumed constant, but it is a function of latitude varying from
               1.0 at poles to 0.5 at the Equator).

            */
            if (have_external_longwave && longwave_is_net) {
                // External forcing provides net longwave directly (W/m2).
                LRad = longwave_down_arr(i,j,0);
            } else if (have_external_longwave && use_longwave_down) {
                Real Ldown = longwave_down_arr(i,j,0);
                Real Lemit = emmiss * StefBo * std::pow(TseaK,4);
                LRad = Ldown - Lemit;
            } else {
                // Original Berliand parameterization
                cff=(Real(0.7859)+Real(0.03477)*TairC)/(one+Real(0.00412)*TairC);
                Real e_sat=std::pow(Real(10.0),cff);
                Real vap_p=e_sat*RH;
                Real cff2=TairK*TairK*TairK;
                Real cff1=cff2*TairK;

                LRad=-emmiss*StefBo*
                    (cff1*(Real(0.39)-Real(0.05)*std::sqrt(vap_p))*
                        (one-Real(0.6823)*cloud*cloud)+
                        cff2*Real(4.0)*(TseaK-TairK));
            }
           /*
            -----------------------------------------------------------------------
              Compute specific humidities (kg/kg).

                note that Qair is the saturation specific humidity at Tair
                             Q is the actual specific humidity
                          Qsea is the saturation specific humidity at Tsea

                      Saturation vapor pressure in mb is first computed and then
                      converted to specific humidity in kg/kg

                      The saturation vapor pressure is computed from Teten formula
                      using the approach of Buck (1981):

                      Esat(mb) = (Real(1.0007)+3.46E-Real(6)*PairM(mb))*Real(6.1121)*
                              EXP(Real(17.502)*TairC(C)/(Real(240.97)+TairC(C)))

                      The ambient vapor is found from the definition of the
                      Relative humidity:

                      RH = W/Ws*100 ~ E/Esat*100   E = RH/100*Esat if RH is in %
                                                   E = RH*Esat     if RH fractional

                      The specific humidity is then found using the relationship:

                      Q = 0.622 E/(P + (0.622-1)e)

                      Q(kg/kg) = Real(0.62197)*(E(mb)/(PairM(mb)Real(-0.378)*E(mb)))

            -----------------------------------------------------------------------
            */

            //  Compute air saturation vapor pressure (mb), using Teten formula.

            Real cff_saturation_air=(Real(1.0007)+Real(3.46e-6)*PairM)*Real(6.1121)*
                                std::exp(Real(17.502)*TairC/(Real(240.97)+TairC));

            //  Compute specific humidity at Saturation, Qair (kg/kg).

            Real Qair = Real(0.62197)*(cff_saturation_air/(PairM-Real(0.378)*cff_saturation_air+eps));

            //  Compute specific humidity, Q (kg/kg).
            Real Q;
            if (RH < 2.0) {
                Real cff_Q = cff_saturation_air*RH;                  //Vapor pressure (mb)
                Q=Real(0.62197)*(cff_Q/(PairM-Real(0.378)*cff_Q+eps)); //Spec hum (kg/kg)
            } else { // RH input was actually specific humidity in g/kg
                Q=RH/Real(1000.0);                          //!Spec Hum (kg/kg)
            }

            //  Compute water saturation vapor pressure (mb), using Teten formula.

            Real cff_saturation_water=(Real(1.0007)+Real(3.46e-6)*PairM)*Real(6.1121)*
                                std::exp(Real(17.502)*cons(i,j,N,Temp_comp)/(Real(240.97)+cons(i,j,N,Temp_comp)));

            //  Compute water saturation vapor pressure (mb), using Teten formula.
            //   Vapor Pressure reduced for salinity (Kraus and Businger, 1994, pp42).
            Real cff_vp=cff_saturation_water*Real(0.98);

            //   Compute Qsea (kg/kg) from vapor pressure.
            //   NOTE: ROMS does not have the small-value guard here, but does for
            //   Q and Qair

            Real Qsea=Real(0.62197)*(cff_vp/(PairM-Real(0.378)*cff_vp+eps));
            //
            // -----------------------------------------------------------------------
            //   Compute Monin-Obukhov similarity parameters for wind (Wstar),
            //   heat (Tstar), and moisture (Qstar), Liu et al. (1979).
            // -----------------------------------------------------------------------
            //
            //   Moist air density (kg/m3).

            Real rhoAir=PairM*Real(100.0)/(blk_Rgas*TairK*(one+Real(0.61)*Q));

            //  Kinematic viscosity of dry air (m2/s), Andreas (1989).

            Real VisAir=Real(1.326e-5)*(one+TairC*(Real(6.542e-3)+TairC*
                                     (Real(8.301e-6)-Real(4.84e-9)*TairC)));


            //  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.

            Real Hlv = (Real(2.501)-Real(0.00237)*cons(i,j,N,Temp_comp))*Real(1.0e6);

            //  Assume that wind is measured relative to sea surface and include
            //  gustiness.

            Real Wgus=half;
            Real delW=std::sqrt(wind_mag*wind_mag+Wgus*Wgus);
            Real delQ=Qsea-Q;
            Real delT=cons(i,j,N,Temp_comp)-TairC;

            // Neutral coefficients.
            Real ZoW=Real(0.0001);
            Real u10=delW*std::log(Real(10.0)/ZoW)/std::log(blk_ZW/ZoW);
            Real Wstar=Real(0.035) * u10;
            Real Zo10=Real(0.011)*Wstar*Wstar/g+Real(0.11)*VisAir/Wstar;
            Real Cd10 =(vonKar/std::log(Real(10.0)/Zo10));
            Cd10 = Cd10 * Cd10;
            Real Ch10 =Real(0.00115);
            Real Ct10 = Ch10/std::sqrt(Cd10);
            Real ZoT10=Real(10.0)/std::exp(vonKar/Ct10);
            Real Cd=(vonKar/std::log(blk_ZW/Zo10));
            Cd = Cd * Cd;

            //  Compute Richardson number.
            Real Ct=vonKar/std::log(blk_ZT/ZoT10);  // T transfer coefficient
            Real CC=vonKar*Ct/Cd;

            Real Ribcu = -blk_ZW/(blk_Zabl*Real(0.004)*blk_beta*blk_beta*blk_beta);
            Real Ri = -g*blk_ZW*((delT-delTc)+Real(0.61)*TairK*delQ)/
                                 (TairK*delW*delW+eps);
            Real Zetu;
            if (Ri < zero) {
                Zetu=CC*Ri/(one+Ri/Ribcu);       // Unstable
            } else {
                Zetu=CC*Ri/(one+Real(3.0)*Ri/CC);   // Stable
            }
            Real L10 = blk_ZW/Zetu;

            // First guesses for Monon-Obukhov similarity scales.
            Wstar=delW*vonKar/(std::log(blk_ZW/Zo10)-
                                bulk_psiu(blk_ZW/L10));
            Real Tstar=-(delT-delTc)*vonKar/(std::log(blk_ZT/ZoT10)-
                                         bulk_psit(blk_ZT/L10));
            Real Qstar=-(delQ-delQc)*vonKar/(std::log(blk_ZQ/ZoT10)-
                                         bulk_psit(blk_ZQ/L10));

            //  Modify Charnock for high wind speeds. The 0.125 factor below is for
            //  1.0/(18.0-10.0).

            Real charn;
            if (delW > Real(18.0)) {
                charn=Real(0.018);
            } else if ((Real(10.0) < delW) and (delW <= Real(18.0))) {
                charn=Real(0.011)+Real(0.125)*(Real(0.018)-Real(0.011))*(delW-Real(10.0));
            } else {
                charn=Real(0.011);
            }

            //  Iterate until convergence. It usually converges within 3 iterations.
            for (int it=0; it<IterMax; it++) {
                ZoW=charn*Wstar*Wstar/g+Real(0.11)*VisAir/(Wstar+eps);
                Real Rr=ZoW*Wstar/VisAir;
                //  Compute Monin-Obukhov stability parameter, Z/L.
                Real ZoQ=std::min(Real(1.15e-4),Real(5.5e-5)/std::pow(Rr,Real(0.6)));
                Real ZoT=ZoQ;
                Real ZoL=vonKar*g*blk_ZW*(Tstar*(one+Real(0.61)*Q)+
                             Real(0.61)*TairK*Qstar)/
                            (TairK*Wstar*Wstar*(one+Real(0.61)*Q)+eps);
                Real L=blk_ZW/(ZoL+eps);

                //  Evaluate stability functions at Z/L.
                Real Wpsi=bulk_psiu(ZoL);
                Real Tpsi=bulk_psit(blk_ZT/L);
                Real Qpsi=bulk_psit(blk_ZQ/L);

                //  Compute wind scaling parameters, Wstar.
                Wstar=std::max(eps,delW*vonKar/(std::log(blk_ZW/ZoW)-Wpsi));
                Tstar=-(delT-delTc)*vonKar/(std::log(blk_ZT/ZoT)-Tpsi);
                Qstar=-(delQ-delQc)*vonKar/(std::log(blk_ZQ/ZoQ)-Qpsi);

                //  Compute gustiness in wind speed.
                Real Bf=-g/TairK*Wstar*(Tstar+Real(0.61)*TairK*Qstar);
                if (Bf>zero) {
                    Wgus=blk_beta*std::pow(Bf*blk_Zabl,one/Real(3.0));
                } else {
                    Wgus=Real(0.2);
                }
                delW=std::sqrt(wind_mag*wind_mag+Wgus*Wgus);
            }

            // Compute transfer coefficients for momentum (Cd).
            Real Wspeed=std::sqrt(wind_mag*wind_mag+Wgus*Wgus);
            Cd=Wstar*Wstar/(Wspeed*Wspeed+eps);

            // Compute turbulent sensible heat flux (W/m2), Hs.
            Real Hs=-blk_Cpa*rhoAir*Wstar*Tstar;

            //  Compute sensible heat flux (W/m2) due to rainfall (kg/m2/s), Hsr.
            Real diffw=Real(2.11e-5)*std::pow(TairK/Real(273.16),Real(1.94));
            Real diffh=Real(0.02411)*(one+TairC*
                               (Real(3.309e-3)-Real(1.44e-6)*TairC))/
                               (rhoAir*blk_Cpa+eps);
            cff=Qair*Hlv/(blk_Rgas*TairK*TairK);
            Real wet_bulb=one/(one+Real(0.622)*(cff*Hlv*diffw)/
                                                  (blk_Cpa*diffh));
            Real Hsr=rain(i,j,0)*wet_bulb*blk_Cpw*
                              ((cons(i,j,N,Temp_comp)-TairC)+(Qsea-Q)*Hlv/blk_Cpa);
            SHeat=(Hs+Hsr) * mskr(i,j,0);

            // Compute turbulent latent heat flux (W/m2), Hl.

            Real Hl=-Hlv*rhoAir*Wstar*Qstar;

            // Compute Webb correction (Webb effect) to latent heat flux, Hlw.
            Real upvel=Real(-1.61)*Wstar*Qstar-
                        (one+Real(1.61)*Q)*Wstar*Tstar/TairK;
            Real Hlw=rhoAir*Hlv*upvel*Q;
            LHeat=(Hl+Hlw) * mskr(i,j,0);

            // Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
            Taur=Real(0.85)*rain(i,j,0)*wind_mag;

            // Compute wind stress components (N/m2), Tau.
            cff=rhoAir*Cd*Wspeed;
            // amrex::Print() << "rhoAir: " << rhoAir << " Cd: " << Cd << " Wspeed: " << Wspeed << " cff: " << cff << "\n";
            Real sign_u = (uwind(i,j,0) >= zero) ? 1 : -1;
            Real sign_v = (vwind(i,j,0) >= zero) ? 1 : -1;
            Taux(i,j,0)=(cff*uwind(i,j,0)+Taur*sign_u) * mskr(i,j,0);
            Tauy(i,j,0)=(cff*vwind(i,j,0)+Taur*sign_v) * mskr(i,j,0);
            // amrex::Print() << "Taux: " << Taux(i,j,0) << " Tauy: " << Tauy(i,j,0) << "\n";

            //=======================================================================
            //  Compute surface net heat flux and surface wind stress.
            //=======================================================================
            //
            //  Compute kinematic, surface, net heat flux (degC m/s).  Notice that
            //  the signs of latent and sensible fluxes are reversed because fluxes
            //  calculated from the bulk formulations above are positive out of the
            //  ocean.
            //
            //  For EMINUSP option,  EVAP = LHeat (W/m2) / Hlv (J/kg) = kg/m2/s
            //                       PREC = rain = kg/m2/s
            //
            //  To convert these rates to m/s divide by freshwater density, rhow.
            //
            //  Note that when the air is undersaturated in water vapor (Q < Qsea)
            //  the model will evaporate and LHeat > 0:
            //
            //                   LHeat positive out of the ocean
            //                    evap positive out of the ocean
            //
            //  Note that if evaporating, the salt flux is positive
            //        and if     raining, the salt flux is negative
            //
            //  Note that stflux(:,:,isalt) is the E-P flux. The fresh water flux
            //  is positive out of the ocean and the salt flux is positive into the
            //  ocean. It is  multiplied by surface salinity when computing state
            //  variable stflx(:,:,isalt) in "set_vbc.F".

//            Real one_over_rhow=one/rhow;
            lrflx(i,j,0) = LRad*Hscale2;
            lhflx(i,j,0) = -LHeat*Hscale2;
            shflx(i,j,0) = -SHeat*Hscale2;
            // Note: srflx from NetCDF is in W/m², convert to degC m/s by multiplying by Hscale2
            stflux(i,j,0,Temp_comp)=(srflux*Hscale2 + lrflx(i,j,0) + lhflx(i,j,0) + shflx(i,j,0)) * mskr(i,j,0);
            evap(i,j,0) = (LHeat / Hlv+eps) * mskr(i,j,0);
            if (use_EminusP_from_input) {
                // Use prescribed E-P directly
                stflux(i,j,0,Salt_comp) = mskr(i,j,0) * EminusP(i,j,0);
            } else {
                stflux(i,j,0,Salt_comp) = mskr(i,j,0) * (evap(i,j,0)-rain(i,j,0)) / rhow;
            }
        });

        Real cff_rho = half / solverChoice.rho0;
        ParallelFor(makeSlab(ubx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int ) {
            sustr(i,j,0) = cff_rho*(Taux(i-1,j,0) + Taux(i,j,0)) * msku(i,j,0);
        });
        ParallelFor(makeSlab(vbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int ) {
            svstr(i,j,0) = cff_rho*(Tauy(i,j-1,0) + Tauy(i,j,0)) * mskv(i,j,0);
        });
    }

}
