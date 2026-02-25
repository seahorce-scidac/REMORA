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
 * @param[in   ] mf_longwave_down longwave radiation flux [W/m²]
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
        Array4<const Real> const& cloud_arr = vec_cloud[lev]->const_array(mfi);

        Real Hscale = solverChoice.rho0 * Cp;
        Real Hscale2 = 1.0_rt / (solverChoice.rho0 * Cp);
        Real blk_ZQ = solverChoice.blk_ZQ;
        Real blk_ZT = solverChoice.blk_ZT;
        Real blk_ZW = solverChoice.blk_ZW;

        bool use_longwave_down = solverChoice.longwave_down;

        Real eps = 1e-20_rt;

        Box bx = mfi.tilebox();
        Box ubx = mfi.grownnodaltilebox(0,IntVect(NGROW-1,NGROW-1,0));
        Box vbx = mfi.grownnodaltilebox(1,IntVect(NGROW-1,NGROW-1,0));
        Box gbx1 = bx; gbx1.grow(IntVect(NGROW,NGROW,0));

        ParallelFor(makeSlab(gbx1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int ) {
            // Get spatially-varying atmospheric forcing from input arrays
            Real PairM = Pair_arr(i,j,0);  // Air pressure [mb]
            Real TairC = Tair_arr(i,j,0);  // Air temperature [°C]
            Real TairK = TairC + 273.16_rt; // Air temperature [K]
            Real Hair = qair_arr(i,j,0);   // Specific humidity [kg/kg] or RH [fraction]
            Real RH = Hair;
            Real srflux = srflx_arr(i,j,0); // Shortwave radiation flux [W/m²]
            Real cloud = cloud_arr(i,j,0);  // Cloud cover fraction [0-1]

            // Input bulk parametrization fields
            Real wind_mag = std::sqrt(uwind(i,j,0)*uwind(i,j,0) + vwind(i,j,0) * vwind(i,j,0)) + eps;
            Real TseaK = cons(i,j,N,Temp_comp) + 273.16_rt;

            // Initialize
            Real cff = 0.0_rt;
            Real delTc = 0.0_rt;
            Real delQc = 0.0_rt;

            Real LHeat = lhflx(i,j,0) * Hscale;
            Real SHeat = shflx(i,j,0) * Hscale;
            Real Taur = 0.0_rt;
            Taux(i,j,0) = 0.0_rt;
            Tauy(i,j,0) = 0.0_rt;


            /*-----------------------------------------------------------------------
               Compute outward or net longwave radiation (W/m2), LRad.
             -----------------------------------------------------------------------
               If given downward longwave radiation, compute net longwave radiation as 
               Ldown - Lemit, where Lemit is computed from the model SST and an emissivity.
               Or use Berliand (1952) formula to calculate net longwave radiation.
               The equation for saturation vapor pressure is from Gill (Atmosphere-
               Ocean Dynamics, pp 606). Here the coefficient in the cloud term
               is assumed constant, but it is a function of latitude varying from
               1.0 at poles to 0.5 at the Equator).

            */
            Real LRad;
            // If given downward longwave radiation, compute net longwave radiation as 
            // Ldown - Lemit, where Lemit is computed from the model SST and an assumed emissivity.

            if (use_longwave_down) {

                Real Ldown = longwave_down_arr(i,j,0);
                Real Lemit = emmiss * StefBo * std::pow(TseaK,4);

                LRad = Ldown - Lemit;

            } else {

                // Original Berliand parameterization
                Real cff=(0.7859_rt+0.03477_rt*TairC)/(1.0_rt+0.00412_rt*TairC);
                Real e_sat=std::pow(10.0_rt,cff);
                Real vap_p=e_sat*RH;
                Real cff2=TairK*TairK*TairK;
                Real cff1=cff2*TairK;

                LRad=-emmiss*StefBo*
                    (cff1*(0.39_rt-0.05_rt*std::sqrt(vap_p))*
                        (1.0_rt-0.6823_rt*cloud*cloud)+
                        cff2*4.0_rt*(TseaK-TairK));
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

                      Esat(mb) = (1.0007_rt+3.46E-6_rt*PairM(mb))*6.1121_rt*
                              EXP(17.502_rt*TairC(C)/(240.97_rt+TairC(C)))

                      The ambient vapor is found from the definition of the
                      Relative humidity:

                      RH = W/Ws*100 ~ E/Esat*100   E = RH/100*Esat if RH is in %
                                                   E = RH*Esat     if RH fractional

                      The specific humidity is then found using the relationship:

                      Q = 0.622 E/(P + (0.622-1)e)

                      Q(kg/kg) = 0.62197_rt*(E(mb)/(PairM(mb)-0.378_rt*E(mb)))

            -----------------------------------------------------------------------
            */

            //  Compute air saturation vapor pressure (mb), using Teten formula.

            Real cff_saturation_air=(1.0007_rt+3.46E-6_rt*PairM)*6.1121_rt*
                                std::exp(17.502_rt*TairC/(240.97_rt+TairC));

            //  Compute specific humidity at Saturation, Qair (kg/kg).

            Real Qair = 0.62197_rt*(cff_saturation_air/(PairM-0.378_rt*cff_saturation_air+eps));

            //  Compute specific humidity, Q (kg/kg).
            Real Q;
            if (RH < 2.0) {
                Real cff_Q = cff_saturation_air*RH;                  //Vapor pressure (mb)
                Q=0.62197_rt*(cff_Q/(PairM-0.378_rt*cff_Q+eps)); //Spec hum (kg/kg)
            } else { // RH input was actually specific humidity in g/kg
                Q=RH/1000.0_rt;                          //!Spec Hum (kg/kg)
            }

            //  Compute water saturation vapor pressure (mb), using Teten formula.

            Real cff_saturation_water=(1.0007_rt+3.46E-6_rt*PairM)*6.1121_rt*
                                std::exp(17.502_rt*cons(i,j,N,Temp_comp)/(240.97_rt+cons(i,j,N,Temp_comp)));

            //  Compute water saturation vapor pressure (mb), using Teten formula.
            //   Vapor Pressure reduced for salinity (Kraus and Businger, 1994, pp42).
            Real cff_vp=cff_saturation_water*0.98_rt;

            //   Compute Qsea (kg/kg) from vapor pressure.

            Real Qsea=0.62197_rt*(cff_vp/(PairM-0.378_rt*cff_vp));
            //
            // -----------------------------------------------------------------------
            //   Compute Monin-Obukhov similarity parameters for wind (Wstar),
            //   heat (Tstar), and moisture (Qstar), Liu et al. (1979).
            // -----------------------------------------------------------------------
            //
            //   Moist air density (kg/m3).

            Real rhoAir=PairM*100.0_rt/(blk_Rgas*TairK*(1.0_rt+0.61_rt*Q));

            //  Kinematic viscosity of dry air (m2/s), Andreas (1989).

            Real VisAir=1.326E-5_rt*(1.0_rt+TairC*(6.542E-3_rt+TairC*
                                     (8.301E-6_rt-4.84E-9_rt*TairC)));
\
            //  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.

            Real Hlv = (2.501_rt-0.00237_rt*cons(i,j,N,Temp_comp))*1.0e6_rt;

            //  Assume that wind is measured relative to sea surface and include
            //  gustiness.

            Real Wgus=0.5_rt;
            Real delW=std::sqrt(wind_mag*wind_mag+Wgus*Wgus);
            Real delQ=Qsea-Q;
            Real delT=cons(i,j,N,Temp_comp)-TairC;

            // Neutral coefficients.
            Real ZoW=0.0001_rt;
            Real u10=delW*std::log(10.0_rt/ZoW)/std::log(blk_ZW/ZoW);
            Real Wstar=0.035_rt * u10;
            Real Zo10=0.011_rt*Wstar*Wstar/g+0.11_rt*VisAir/Wstar;
            Real Cd10 =(vonKar/std::log(10.0_rt/Zo10));
            Cd10 = Cd10 * Cd10;
            Real Ch10 =0.00115_rt;
            Real Ct10 = Ch10/std::sqrt(Cd10);
            Real ZoT10=10.0_rt/std::exp(vonKar/Ct10);
            Real Cd=(vonKar/std::log(blk_ZW/Zo10));
            Cd = Cd * Cd;

            //  Compute Richardson number.
            Real Ct=vonKar/std::log(blk_ZT/ZoT10);  // T transfer coefficient
            Real CC=vonKar*Ct/Cd;

            Real Ribcu = -blk_ZW/(blk_Zabl*0.004_rt*blk_beta*blk_beta*blk_beta);
            Real Ri = -g*blk_ZW*((delT-delTc)+0.61_rt*TairK*delQ)/
                                 (TairK*delW*delW+eps);
            Real Zetu;
            if (Ri < 0.0) {
                Zetu=CC*Ri/(1.0_rt+Ri/Ribcu);       // Unstable
            } else {
                Zetu=CC*Ri/(1.0_rt+3.0_rt*Ri/CC);   // Stable
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
            if (delW > 18.0_rt) {
                charn=0.018_rt;
            } else if ((10.0_rt < delW) and (delW <= 18.0_rt)) {
                charn=0.011_rt+0.125_rt*(0.018_rt-0.011_rt)*(delW-10.0_rt);
            } else {
                charn=0.011_rt;
            }

            //  Iterate until convergence. It usually converges within 3 iterations.
            for (int it=0; it<IterMax; it++) {
                ZoW=charn*Wstar*Wstar/g+0.11_rt*VisAir/(Wstar+eps);
                Real Rr=ZoW*Wstar/VisAir;
                //  Compute Monin-Obukhov stability parameter, Z/L.
                Real ZoQ=std::min(1.15e-4_rt,5.5e-5_rt/std::pow(Rr,0.6_rt));
                Real ZoT=ZoQ;
                Real ZoL=vonKar*g*blk_ZW*(Tstar*(1.0_rt+0.61_rt*Q)+
                             0.61_rt*TairK*Qstar)/
                            (TairK*Wstar*Wstar*(1.0_rt+0.61_rt*Q)+eps);
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
                Real Bf=-g/TairK*Wstar*(Tstar+0.61_rt*TairK*Qstar);
                if (Bf>0.0_rt) {
                    Wgus=blk_beta*std::pow(Bf*blk_Zabl,1.0_rt/3.0_rt);
                } else {
                    Wgus=0.2_rt;
                }
                delW=std::sqrt(wind_mag*wind_mag+Wgus*Wgus);
            }

            // Compute transfer coefficients for momentum (Cd).
            Real Wspeed=std::sqrt(wind_mag*wind_mag+Wgus*Wgus);
            Cd=Wstar*Wstar/(Wspeed*Wspeed+eps);

            // Compute turbulent sensible heat flux (W/m2), Hs.
            Real Hs=-blk_Cpa*rhoAir*Wstar*Tstar;

            //  Compute sensible heat flux (W/m2) due to rainfall (kg/m2/s), Hsr.
            Real diffw=2.11E-5_rt*std::pow(TairK/273.16_rt,1.94_rt);
            Real diffh=0.02411_rt*(1.0_rt+TairC*
                               (3.309E-3_rt-1.44E-6_rt*TairC))/
                               (rhoAir*blk_Cpa+eps);
            cff=Qair*Hlv/(blk_Rgas*TairK*TairK);
            Real wet_bulb=1.0_rt/(1.0_rt+0.622_rt*(cff*Hlv*diffw)/
                                                  (blk_Cpa*diffh));
            Real Hsr=rain(i,j,0)*wet_bulb*blk_Cpw*
                              ((cons(i,j,N,Temp_comp)-TairC)+(Qsea-Q)*Hlv/blk_Cpa);
            SHeat=(Hs+Hsr) * mskr(i,j,0);

            // Compute turbulent latent heat flux (W/m2), Hl.

            Real Hl=-Hlv*rhoAir*Wstar*Qstar;

            // Compute Webb correction (Webb effect) to latent heat flux, Hlw.
            Real upvel=-1.61_rt*Wstar*Qstar-
                        (1.0_rt+1.61_rt*Q)*Wstar*Tstar/TairK;
            Real Hlw=rhoAir*Hlv*upvel*Q;
            LHeat=(Hl+Hlw) * mskr(i,j,0);

            // Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
            Taur=0.85_rt*rain(i,j,0)*wind_mag;

            // Compute wind stress components (N/m2), Tau.
            cff=rhoAir*Cd*Wspeed;
            // amrex::Print() << "rhoAir: " << rhoAir << " Cd: " << Cd << " Wspeed: " << Wspeed << " cff: " << cff << "\n";
            Real sign_u = (uwind(i,j,0) >= 0.0_rt) ? 1 : -1;
            Real sign_v = (vwind(i,j,0) >= 0.0_rt) ? 1 : -1;
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

//            Real one_over_rhow=1.0_rt/rhow;
            lrflx(i,j,0) = LRad*Hscale2;
            lhflx(i,j,0) = -LHeat*Hscale2;
            shflx(i,j,0) = -SHeat*Hscale2;
            // Note: srflx from NetCDF is in W/m², convert to degC m/s by multiplying by Hscale2
            stflux(i,j,0,Temp_comp)=(srflux*Hscale2 + lrflx(i,j,0) + lhflx(i,j,0) + shflx(i,j,0)) * mskr(i,j,0);
            evap(i,j,0) = (LHeat / Hlv+eps) * mskr(i,j,0);
            stflux(i,j,0,Salt_comp) = mskr(i,j,0) * (evap(i,j,0)-rain(i,j,0)) / rhow;
        });

        Real cff_rho = 0.5_rt / solverChoice.rho0;
        ParallelFor(makeSlab(ubx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int ) {
            sustr(i,j,0) = cff_rho*(Taux(i-1,j,0) + Taux(i,j,0)) * msku(i,j,0);
        });
        ParallelFor(makeSlab(vbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int ) {
            svstr(i,j,0) = cff_rho*(Tauy(i,j-1,0) + Tauy(i,j,0)) * mskv(i,j,0);
        });
    }

}
