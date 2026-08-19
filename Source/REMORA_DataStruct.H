#ifndef _REMORA_DATA_STRUCT_H_
#define _REMORA_DATA_STRUCT_H_

#include <string>
#include <iostream>
#include <array>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

#include <REMORA_Constants.H>
#include <REMORA_DateClock.H>
#include "REMORA_IndexDefines.H"

/** \brief Type of coupling between levels in AMR */
enum struct CouplingType {
    one_way, two_way
};

/** \brief Coordinates */
enum class Coord {
    x, y, z
};

/** \brief Horizontal advection schemes */
enum class AdvectionScheme {
    centered2, centered4, upstream3
};

/** \brief Type of initial condition type. Analytic reads from prob.cpp. Netcdf is from file */
enum class IC_Type {
    analytic, netcdf
};

/** \brief Coriolis factor */
enum class Cor_Type {
    analytic, beta_plane, netcdf
};

/** \brief plotfile format */
enum class PlotfileType {
    amrex, netcdf
};

/** \brief vertical mixing type */
enum class VertMixingType {
    analytic, GLS
};

/** \brief horizontal viscosity/diffusion type */
enum class HorizMixingType {
    analytic, constant, scaled_to_grid
};

/** \brief How to scale scaled_to_grid coefficients on AMR levels */
enum class ScaledToGridAMRScaling {
    none, linear
};

/** \brief stability function for GLS */
enum class GLS_StabilityType {
    Canuto_A, Canuto_B, Galperin
};

/** \brief equation of state */
enum class EOSType {
    linear, nonlinear
};

/** \brief bottom stress formulation */
enum class BottomStressType {
    linear, quadratic, logarithmic
};

/** \brief initialization for pm and pn */
enum class GridScaleType {
    constant, analytic
};

/** \brief surface momentum flux */
enum class SMFluxType {
    analytic, netcdf
};

/** \brief surface wind */
enum class WindType {
    analytic, netcdf
};

/** \brief masks */
enum class MaskType {
    none, analytic, netcdf
};

/** \brief harmonic mixing; which surfaces to calculate along */
enum class HarmonicMixingType {
    s, geopotential
};

struct SolverChoice {
  public:
    static BulkForcingType parse_bulk_forcing_type(const std::string& type_string,
                                                   const std::string& param_name,
                                                   bool allow_computed)
    {
        const std::string type = amrex::toLower(type_string);
        if (type == "constant" || type == "const") {
            return BulkForcingType::constant;
        } else if (type == "custom") {
            amrex::Warning((param_name + " uses 'custom'; use 'analytic' instead.").c_str());
            return BulkForcingType::analytic;
        } else if (type == "analytic" || type == "analytical") {
            return BulkForcingType::analytic;
        } else if (type == "netcdf" || type == "file" || type == "nc") {
            return BulkForcingType::netcdf;
        } else if (type == "computed") {
            if (!allow_computed) {
                amrex::Abort((param_name + " does not accept 'computed'; only remora.lwrad_type and remora.eminusp_type do.").c_str());
            }
            return BulkForcingType::computed;
        }
        amrex::Abort(("Don't know this " + param_name).c_str());
        return BulkForcingType::constant;
    }

    /** \brief read in and initialize parameters */
    void init_params(int ncons, int nscalar, const amrex::Vector<std::string>& cons_names)
    {
        amrex::ParmParse pp(pp_prefix);

        pp.queryAdd("flat_bathymetry", flat_bathymetry);

        // Which horizontal advection scheme for tracers
        static std::string tracer_hadv_string = "upstream3";
        pp.queryAdd("tracer_horizontal_advection_scheme",tracer_hadv_string);
        if (tracer_hadv_string == "centered4")
            tracer_Hadv_scheme = AdvectionScheme::centered4;
        else if (tracer_hadv_string == "upstream3")
            tracer_Hadv_scheme = AdvectionScheme::upstream3;
        else
            amrex::Error("Advection scheme unknown.");

        // Which horizontal advection scheme
        static std::string uv_hadv_string = "upstream3";
        pp.queryAdd("uv_horizontal_advection_scheme",uv_hadv_string);
        if (uv_hadv_string == "upstream3")
            uv_Hadv_scheme = AdvectionScheme::upstream3;
        else if (uv_hadv_string == "centered2")
            uv_Hadv_scheme = AdvectionScheme::centered2;
        else
            amrex::Error("UV advection scheme unknown.");

        pp.queryAdd("rdrag", rdrag);
        pp.queryAdd("rdrag2", rdrag2);
        pp.queryAdd("Zos", Zos);
        pp.queryAdd("Zob", Zob);
        pp.queryAdd("Cdb_max", Cdb_max);
        pp.queryAdd("Cdb_min", Cdb_min);

        // Include salinity?
        pp.queryAdd("use_salt", use_salt);

        // Include Coriolis forcing?
        pp.queryAdd("use_coriolis", use_coriolis);

        // Include prestep / lagged predictor / corrections
        pp.queryAdd("use_prestep", use_prestep);

        //This affect forcing and some of the mixing terms for velocity
        pp.queryAdd("use_uv3dmix", use_uv3dmix);

        //This accounts for the main 2d loop but may not include coupling and copying properly
        pp.queryAdd("use_barotropic", use_barotropic);

        pp.queryAdd("use_curvilinear_grid", use_curvilinear_grid);

        // Whether to do rivers. By default, rivers are temp and salt sources. Rivers
        // have to be momentum sources.
        pp.queryAdd("do_rivers", do_rivers);

        pp.queryAdd("do_rivers_temp", do_rivers_temp);
        pp.queryAdd("do_rivers_salt", do_rivers_salt);
        pp.queryAdd("do_rivers_scalar", do_rivers_scalar);

        // Every tracer may take river input, keyed by its own name, so temp and salt keep
        // do_rivers_temp / do_rivers_salt and other tracers use e.g. do_rivers_NO3 or
        // do_rivers_tracer_1. The passive scalars default to do_rivers_scalar; biology
        // tracers default off, since a river concentration for them has to be a deliberate
        // choice. If we aren't doing rivers at all, every flag ends up false.
        do_rivers_cons.assign(ncons, 0);
        for (int icomp = 0; icomp < ncons; ++icomp) {
            bool flag = (icomp == Temp_comp) ? do_rivers_temp
                      : (icomp == Salt_comp) ? do_rivers_salt
                      : (icomp <  Tracer_comp + nscalar) ? do_rivers_scalar
                      : false;
            pp.queryAdd(("do_rivers_" + cons_names[icomp]).c_str(), flag);
            do_rivers_cons[icomp] = (do_rivers && flag) ? 1 : 0;
        }

        pp.queryAdd("init_l1ad_T", init_l1ad_T);

        pp.queryAdd("init_ana_T", init_ana_T);

        pp.queryAdd("init_l0int_T", init_l0int_T);

        static std::string eos_type_string = "linear";
        pp.queryAdd("eos_type",eos_type_string);
        if (eos_type_string == "linear" || eos_type_string == "Linear" ||
              eos_type_string == "lin"    || eos_type_string == "Lin") {
            eos_type = EOSType::linear;
            pp.queryAdd("Tcoef",Tcoef);
            pp.queryAdd("Scoef",Scoef);
        } else if (eos_type_string == "nonlinear" || eos_type_string == "Nonlinear" ||
                   eos_type_string == "non-linear" || eos_type_string == "Non-linear" ||
                   eos_type_string == "nonlin"     || eos_type_string == "Nonlin") {
            eos_type = EOSType::nonlinear;
        } else {
            amrex::Abort("Dont know this eos_type");
        }
        pp.queryAdd("R0",R0);
        pp.queryAdd("S0",S0);
        pp.queryAdd("T0",T0);
        pp.queryAdd("rho0", rho0);

        // Reference date for the model clock, and with it the calendar, as ROMS
        // TIME_REF: a yyyymmdd.dd date, or 0 for 0001-01-01, -1 for the 360-day
        // calendar, -2 for truncated Julian days. See REMORA_DateClock.H.
        pp.queryAdd("time_ref", time_ref);
        if (!remora_time_ref_is_valid(time_ref)) {
            // ROMS's own calendar dispatch has no final else, so a value below
            // -2 leaves the date uninitialized rather than failing.
            amrex::Abort("remora.time_ref = " + std::to_string(time_ref) +
                         " names no calendar. Use a yyyymmdd.dd date, or 0, -1, or -2.");
        }
        if (!remora_time_ref_is_representable(time_ref)) {
            // Fail loudly rather than run on a silently shifted calendar: see
            // remora_time_ref_is_representable.
            amrex::Abort("remora.time_ref = " + std::to_string(time_ref) +
                         " is a yyyymmdd.dd date, which needs eight exact digits and so"
                         " cannot be represented in a single-precision build. Build in"
                         " double precision, or use time_ref = 0, -1, or -2.");
        }

        pp.queryAdd("bulk_fluxes",bulk_fluxes);
        pp.queryAdd("atm2ocn_flux_mode", atm2ocn_flux_mode);
        {
            amrex::ParmParse pp_driver("driver");
            std::string driver_atm2ocn_mode = "state";
            pp_driver.query("atm2ocn_mode", driver_atm2ocn_mode);
            if (amrex::toLower(driver_atm2ocn_mode) == "flux") {
                atm2ocn_flux_mode = true;
            }
        }
        if (bulk_fluxes || atm2ocn_flux_mode) {
            do_salt_flux = true;
            do_temp_flux = true;
        }
        // Outputs forcing variables if true
        pp.queryAdd("output_forcing", output_forcing);
        pp.queryAdd("air_pressure",Pair);
        pp.queryAdd("air_temperature",Tair);
        pp.queryAdd("air_humidity",Hair);
        pp.queryAdd("surface_radiation_flux",srflux);
        pp.queryAdd("uwind", Uwind);
        pp.queryAdd("vwind", Vwind);
        pp.queryAdd("longwave_radiation_flux", longwave_rad);
        pp.query("lwrad", longwave_rad);
        pp.queryAdd("longwave_down", longwave_down);
        pp.queryAdd("longwave_is_net", longwave_is_net);
        pp.queryAdd("longwave_netcdf_varname", longwave_netcdf_varname);
        pp.queryAdd("cloud",cloud);
        pp.queryAdd("rain",rain);
        pp.queryAdd("EminusP", EminusP);
        pp.query("eminusp_value", EminusP);
        pp.queryAdd("blk_ZQ",blk_ZQ);
        pp.queryAdd("blk_ZT",blk_ZT);
        pp.queryAdd("blk_ZW",blk_ZW);
        pp.queryAdd("eminusp",eminusp);
        pp.queryAdd("eminusp_correct_ssh",eminusp_correct_ssh);
        pp.queryAdd("qair_is_percent",qair_is_percent);

        struct BulkTypeInput {
            const char* name;
            int idx;
            const char* default_type;
            bool allow_computed;
        };

        const BulkTypeInput bulk_type_inputs[] = {
            {"uwind_type", BulkFlux::Uwind, "analytic", false},
            {"vwind_type", BulkFlux::Vwind, "analytic", false},
            {"tair_type", BulkFlux::Tair, "constant", false},
            {"qair_type", BulkFlux::Qair, "constant", false},
            {"pair_type", BulkFlux::Pair, "constant", false},
            {"swrad_type", BulkFlux::SWrad, "constant", false},
            {"lwrad_type", BulkFlux::LWrad, "computed", true},
            {"rain_type", BulkFlux::Rain, "constant", false},
            {"cloud_type", BulkFlux::Cloud, "constant", false},
            {"eminusp_type", BulkFlux::EminusP, "computed", true},
        };

        bool uwind_type_specified = false;
        bool vwind_type_specified = false;
        for (const BulkTypeInput& input : bulk_type_inputs) {
            std::string type_string = input.default_type;
            const bool type_specified = pp.queryAdd(input.name, type_string);
            if (type_specified) {
                bulk_flux_type[input.idx] = parse_bulk_forcing_type(type_string,
                                                                     "remora." + std::string(input.name),
                                                                     input.allow_computed);
                if (input.idx == BulkFlux::Uwind) {
                    uwind_type_specified = true;
                } else if (input.idx == BulkFlux::Vwind) {
                    vwind_type_specified = true;
                }
            }
        }

        struct BulkTypeAlias {
            const char* name;
            int idx;
            bool allow_computed;
        };

        const BulkTypeAlias bulk_type_aliases[] = {
            {"srflx_type", BulkFlux::SWrad, false},
            {"longwave_type", BulkFlux::LWrad, true},
            {"longwave_down_type", BulkFlux::LWrad, true},
        };

        for (const BulkTypeAlias& input : bulk_type_aliases) {
            std::string type_string;
            if (pp.query(input.name, type_string)) {
                bulk_flux_type[input.idx] = parse_bulk_forcing_type(type_string,
                                                                     "remora." + std::string(input.name),
                                                                     input.allow_computed);
            }
        }

        const BulkTypeAlias legacy_netcdf_flags[] = {
            {"Tair_from_netcdf", BulkFlux::Tair, false},
            {"qair_from_netcdf", BulkFlux::Qair, false},
            {"Pair_from_netcdf", BulkFlux::Pair, false},
            {"srflx_from_netcdf", BulkFlux::SWrad, false},
            {"longwave_down_from_netcdf", BulkFlux::LWrad, false},
            {"rain_from_netcdf", BulkFlux::Rain, false},
            {"cloud_from_netcdf", BulkFlux::Cloud, false},
            {"EminusP_from_netcdf", BulkFlux::EminusP, false},
        };

        for (const BulkTypeAlias& input : legacy_netcdf_flags) {
            bool from_netcdf = false;
            if (pp.query(input.name, from_netcdf) && from_netcdf) {
                bulk_flux_type[input.idx] = BulkForcingType::netcdf;
            }
        }

        bool legacy_longwave_netcdf_is_net = false;
        if (pp.query("longwave_netcdf_is_net", legacy_longwave_netcdf_is_net) && legacy_longwave_netcdf_is_net) {
            longwave_is_net = true;
            if (bulk_flux_type[BulkFlux::LWrad] == BulkForcingType::computed) {
                bulk_flux_type[BulkFlux::LWrad] = BulkForcingType::netcdf;
            }
        }

        if (bulk_flux_type[BulkFlux::LWrad] != BulkForcingType::computed) {
            longwave_down = !longwave_is_net;
        }

        if (bulk_flux_type[BulkFlux::LWrad] != BulkForcingType::netcdf && longwave_netcdf_varname != "lwrad") {
            amrex::Warning("remora.longwave_netcdf_varname is set but remora.lwrad_type is not netcdf; the value will be ignored.");
        }

        if (eminusp_correct_ssh and !bulk_fluxes) {
            amrex::Abort("If evaporation minus precipitation (E-P) sea surface height correction is on, bulk fluxes must be on as well (remora.bulk_fluxes=true)");
        }
        if (eminusp and !bulk_fluxes) {
            amrex::Abort("Evaporation minus precipitation (E-P) requires bulk flux parametrizations (remora.bulk_fluxes=true)");
        }

        const bool use_external_lwrad =
            bulk_flux_type[BulkFlux::LWrad] != BulkForcingType::computed;

        if (longwave_is_net && !use_external_lwrad) {
            amrex::Abort("remora.longwave_is_net=true requires remora.lwrad_type to be constant, analytic, or netcdf");
        }

        if (longwave_down && !use_external_lwrad) {
            amrex::Abort("remora.longwave_down=true requires remora.lwrad_type to be constant, analytic, or netcdf");
        }

        {
            static bool printed_ep_source = false;
            if (!printed_ep_source) {
                printed_ep_source = true;
                if (bulk_fluxes && eminusp && bulk_flux_type[BulkFlux::EminusP] == BulkForcingType::netcdf) {
                    amrex::Print() << "[REMORA] Active E-P source: NetCDF EminusP (remora.eminusp_type=netcdf).\n";
                } else if (bulk_fluxes && eminusp && bulk_flux_type[BulkFlux::EminusP] == BulkForcingType::analytic) {
                    amrex::Print() << "[REMORA] Active E-P source: analytic EminusP (remora.eminusp_type=analytic).\n";
                } else if (bulk_fluxes && eminusp && bulk_flux_type[BulkFlux::EminusP] == BulkForcingType::constant) {
                    amrex::Print() << "[REMORA] Active E-P source: constant EminusP (remora.eminusp_type=constant).\n";
                } else {
                    amrex::Print() << "[REMORA] Active E-P source: bulk evap-rain diagnostic.\n";
                }
            }
        }

        //Grid stretching
        pp.queryAdd("theta_s",theta_s);
        pp.queryAdd("theta_b",theta_b);
        pp.queryAdd("tcline",tcline);

        //coriolis factor
        pp.queryAdd("coriolis_f0",coriolis_f0);
        pp.queryAdd("coriolis_beta",coriolis_beta);

        pp.queryAdd("Akv_bak",Akv_bak);
        // Minimum/initial vertical diffusivity for both active tracers. Override one of them
        // with remora.Akt_bak_temp or remora.Akt_bak_salt; those keys are resolved below,
        // once cons_names is in scope.
        pp.queryAdd("Akt_bak",Akt_bak_all);


        static std::string grid_scale_type_string = "constant";
        pp.queryAdd("grid_scale_type",grid_scale_type_string);

        if (amrex::toLower(grid_scale_type_string) == "constant") {
            grid_scale_type = GridScaleType::constant;
        } else if (amrex::toLower(grid_scale_type_string) == "custom") {
            amrex::Warning("Initialization of grid scale from prob.cpp is now called 'analytic'. 'custom' will be deprecated");
            grid_scale_type = GridScaleType::analytic;
        } else if (amrex::toLower(grid_scale_type_string) == "analytic") {
            grid_scale_type = GridScaleType::analytic;
        } else {
            amrex::Error("Don't know this grid_scale_type");
        }

        static std::string ic_type_string = "analytic";
        bool found_ic_bc = pp.queryAdd("ic_bc_type", ic_type_string);
        pp.queryAdd("ic_type", ic_type_string);

        if (found_ic_bc) {
            amrex::Warning("remora.ic_bc_type is now called remora.ic_type, and will eventually be deprecated");
        }

        if ( amrex::toLower(ic_type_string) == "custom") {
            amrex::Warning("Problem initialization from prob.cpp is now called 'analytic'. 'custom' will be deprecated");
            ic_type = IC_Type::analytic;
        } else if ( amrex::toLower(ic_type_string) == "analytic") {
            ic_type = IC_Type::analytic;
        } else if ( amrex::toLower(ic_type_string) == "netcdf") {
            ic_type = IC_Type::netcdf;
        } else if ( amrex::toLower(ic_type_string) == "real") {
            amrex::Warning("Problem initialization from NetCDF (remora.ic_type) is now called 'netcdf'. 'real' will be deprecated");
            ic_type = IC_Type::netcdf;
        } else {
            amrex::Error("Don't know this ic_type");
        }

#ifndef REMORA_USE_NETCDF
        // Without this the netcdf branches are compiled out one by one and the run integrates
        // uninitialized state instead of failing: the #ifdef in REMORA::init_only sits inside
        // the ic_type == netcdf branch rather than around it. mask_type below defaults to
        // netcdf for a netcdf run, which would be silently skipped the same way.
        if (ic_type == IC_Type::netcdf) {
            amrex::Abort("remora.ic_type = netcdf requires a NetCDF build. Rebuild with "
                         "USE_PNETCDF=TRUE (GNUmake) or -DREMORA_ENABLE_PNETCDF=ON (CMake)");
        }
#endif

        // Which type of refinement
        static std::string coupling_type_string = "TwoWay";
        pp.queryAdd("coupling_type",coupling_type_string);
        if (amrex::toLower(coupling_type_string) == "twoway" ||
            amrex::toLower(coupling_type_string) == "two_way") {
            coupling_type = CouplingType::two_way;
        } else if (amrex::toLower(coupling_type_string) == "oneway" ||
            amrex::toLower(coupling_type_string) == "one_way") {
            coupling_type = CouplingType::one_way;
        } else {
            amrex::Abort("Dont know this coupling_type");
        }

        // Which type of coriolis forcing
        if (use_coriolis) {
            static std::string coriolis_type_string = "beta_plane";
            pp.queryAdd("coriolis_type",coriolis_type_string);
            if ( amrex::toLower(coriolis_type_string) == "custom") {
                amrex::Warning("Coriolis initialization from prob.cpp is now called 'analytic'. 'custom' will be deprecated");
                coriolis_type = Cor_Type::analytic;
            } else if ( amrex::toLower(coriolis_type_string) == "analytic") {
                coriolis_type = Cor_Type::analytic;
            } else if ((amrex::toLower(coriolis_type_string) == "beta_plane") ||
                       (amrex::toLower(coriolis_type_string) == "betaplane")) {
                coriolis_type = Cor_Type::beta_plane;
            } else if ( (amrex::toLower(coriolis_type_string) == "netcdf")) {
                coriolis_type = Cor_Type::netcdf;
            } else if ( (amrex::toLower(coriolis_type_string) == "real")) {
                amrex::Warning("Coriolis initialization from NetCDF is now called 'netcdf'. 'real' will be deprecated");
                coriolis_type = Cor_Type::netcdf;
            } else {
                amrex::Abort("Don't know this coriolis_type");
            }
        }

        static std::string smflux_type_string = "analytic";
        int smflux_specified = pp.queryAdd("smflux_type",smflux_type_string);
        if ( amrex::toLower(smflux_type_string) == "custom") {
            amrex::Warning("Surface momentum flux initialization from prob.cpp is now called 'analytic'. 'custom' will be deprecated");
            smflux_type = SMFluxType::analytic;
        } else if ( amrex::toLower(smflux_type_string) == "analytic") {
            smflux_type = SMFluxType::analytic;
        } else if ( amrex::toLower(smflux_type_string) == "netcdf") {
            smflux_type = SMFluxType::netcdf;
        } else {
            amrex::Abort("Don't know this smflux_type");
        }

        static std::string wind_type_string = "analytic";
        int wind_specified = pp.queryAdd("wind_type",wind_type_string);
        if ( amrex::toLower(wind_type_string) == "custom") {
            amrex::Warning("Surface wind initialization from prob.cpp is now called 'analytic'. 'custom' will be deprecated");
            wind_type = WindType::analytic;
        } else if ( amrex::toLower(wind_type_string) == "analytic") {
            wind_type = WindType::analytic;
        } else if ( amrex::toLower(wind_type_string) == "netcdf") {
            wind_type = WindType::netcdf;
        } else {
            amrex::Abort("Don't know this smflux_type");
        }

        if (wind_specified) {
            const auto specified_wind_bulk_type = (wind_type == WindType::netcdf) ?
                BulkForcingType::netcdf : BulkForcingType::analytic;
            if (!uwind_type_specified) {
                bulk_flux_type[BulkFlux::Uwind] = specified_wind_bulk_type;
            }
            if (!vwind_type_specified) {
                bulk_flux_type[BulkFlux::Vwind] = specified_wind_bulk_type;
            }
        }

        if (wind_specified && smflux_specified) {
            amrex::Abort("Cannot specify both wind and surface momentum flux");
        }

        static std::string mask_type_string = "none";
        // If initial condition type is netcdf, default to netcdf masks
        if (ic_type == IC_Type::netcdf) {
            mask_type_string = "netcdf";
        }
        pp.queryAdd("mask_type", mask_type_string);
        if (amrex::toLower(mask_type_string) == "none") {
            mask_type = MaskType::none;
        } else if (amrex::toLower(mask_type_string) == "analytic") {
            mask_type = MaskType::analytic;
        } else if (amrex::toLower(mask_type_string) == "netcdf") {
            mask_type = MaskType::netcdf;
        } else {
            amrex::Abort("Don't know this mask_type");
        }

        static std::string bottom_stress_type_string = "linear";
        pp.queryAdd("bottom_stress_type", bottom_stress_type_string);
        if (amrex::toLower(bottom_stress_type_string) == "linear") {
            bottom_stress_type = BottomStressType::linear;
        } else if (amrex::toLower(bottom_stress_type_string) == "quadratic") {
            bottom_stress_type = BottomStressType::quadratic;
        } else if (amrex::toLower(bottom_stress_type_string) == "logarithmic") {
            bottom_stress_type = BottomStressType::logarithmic;
        } else {
            amrex::Abort("Don't know this bottom_stress_type");
        }

        amrex::Real tnu2_salt = zero;
        amrex::Real tnu2_temp = zero;
        amrex::Real tnu2_scalar = zero;
        static std::string horiz_mixing_type_string = "analytic";
        pp.queryAdd("horizontal_mixing_type", horiz_mixing_type_string);
        if (amrex::toLower(horiz_mixing_type_string) == "analytical" ||
            amrex::toLower(horiz_mixing_type_string) == "analytic") {
            horiz_mixing_type = HorizMixingType::analytic;
        } else if (amrex::toLower(horiz_mixing_type_string) == "constant") {
            horiz_mixing_type = HorizMixingType::constant;
        } else if (amrex::toLower(horiz_mixing_type_string) == "scaled_to_grid") {
            horiz_mixing_type = HorizMixingType::scaled_to_grid;
        } else {
            amrex::Abort("Don't know this horizontal mixing type");
        }
        pp.queryAdd("visc2",visc2);
        pp.queryAdd("tnu2_salt",tnu2_salt);
        pp.queryAdd("tnu2_temp",tnu2_temp);
        pp.queryAdd("tnu2_scalar",tnu2_scalar);

        // For scaled_to_grid runs with AMR refinement: optionally scale the coefficients
        // by the horizontal refinement ratio (linear in grid size).
        static std::string scaled_to_grid_amr_scaling_string = "none";
        pp.queryAdd("scaled_to_grid_amr_scaling", scaled_to_grid_amr_scaling_string);
        if (amrex::toLower(scaled_to_grid_amr_scaling_string) == "none") {
            scaled_to_grid_amr_scaling = ScaledToGridAMRScaling::none;
        } else if (amrex::toLower(scaled_to_grid_amr_scaling_string) == "linear") {
            scaled_to_grid_amr_scaling = ScaledToGridAMRScaling::linear;
        } else {
            amrex::Abort("Don't know this scaled_to_grid_amr_scaling option");
        }

        // Horizontal diffusivity is per tracer in ROMS, so mirror that here:
        // remora.tnu2_scalar sets the group default and remora.tnu2_{var} overrides a single
        // tracer by name -- remora.tnu2_NO3, remora.tnu2_tracer_1, and so on.
        //
        // tnu2_temp and tnu2_salt are the pre-existing spelling of that per-name form for the
        // first two components. They keep their own defaults of zero rather than inheriting
        // tnu2_scalar, so every existing input file means exactly what it did before, and the
        // loop below skips them to avoid querying the same key twice.
        tnu2.assign(ncons, tnu2_scalar);
        if (ncons > Temp_comp) {
            tnu2[Temp_comp] = tnu2_temp;
        }
        if (ncons > Salt_comp) {
            tnu2[Salt_comp] = tnu2_salt;
        }
        for (int icomp = Tracer_comp; icomp < ncons; ++icomp) {
            pp.queryAdd(("tnu2_" + cons_names[icomp]).c_str(), tnu2[icomp]);
        }

        // Vertical diffusivity is not. ROMS computes and stores Akt for the active tracers
        // only, and every passive tracer -- dye and biology alike -- mixes with the salinity
        // coefficient, so Akt_bak has NAT entries: remora.Akt_bak sets both and
        // remora.Akt_bak_temp / remora.Akt_bak_salt override one.
        Akt_bak.assign(NAT, Akt_bak_all);
        for (int icomp = 0; icomp < NAT; ++icomp) {
            pp.queryAdd(("Akt_bak_" + cons_names[icomp]).c_str(), Akt_bak[icomp]);
        }
        // A per-tracer key would name a coefficient nothing reads, and reading as though it
        // had been applied is worse than not offering it, so say so rather than ignoring it.
        for (int icomp = Tracer_comp; icomp < ncons; ++icomp) {
            const std::string key = "Akt_bak_" + cons_names[icomp];
            if (pp.contains(key.c_str())) {
                amrex::Abort("remora." + key + " has no effect: vertical diffusivity is carried"
                             " for temperature and salinity only, and every passive tracer mixes"
                             " with the salinity value. Set remora.Akt_bak_salt instead, or"
                             " remora.tnu2_" + cons_names[icomp] + " if you meant the horizontal"
                             " diffusivity.");
            }
        }

        static std::string harmonic_mixing_type_string = "s";
        pp.queryAdd("harmonic_mixing_type", harmonic_mixing_type_string);
        if (amrex::toLower(harmonic_mixing_type_string) == "s") {
            harmonic_mixing_type = HarmonicMixingType::s;
        } else if (amrex::toLower(harmonic_mixing_type_string) == "geopotential" ||
                   amrex::toLower(harmonic_mixing_type_string) == "geo") {
            harmonic_mixing_type = HarmonicMixingType::geopotential;
        } else {
            amrex::Abort("Don't know this harmonic_mixing_type");
        }

        pp.queryAdd("Akk_bak", Akk_bak);
        pp.queryAdd("Akp_bak", Akp_bak);
        static std::string vert_mixing_type_string = "analytic";
        static std::string gls_stability_type_string = "Canuto_A";
        pp.queryAdd("vertical_mixing_type", vert_mixing_type_string);
        pp.queryAdd("gls_stability_type", gls_stability_type_string);
        if (amrex::toLower(vert_mixing_type_string) == "analytical" ||
            amrex::toLower(vert_mixing_type_string) == "analytic") {
            vert_mixing_type = VertMixingType::analytic;
        } else if (amrex::toLower(vert_mixing_type_string) == "gls") {
            vert_mixing_type = VertMixingType::GLS;
            if (amrex::toLower(gls_stability_type_string) == "canuto_a") {
                gls_stability_type = GLS_StabilityType::Canuto_A;
            }
            else if (amrex::toLower(gls_stability_type_string) == "canuto_b") {
                gls_stability_type = GLS_StabilityType::Canuto_B;
            }
            else if (amrex::toLower(gls_stability_type_string) == "galperin") {
                gls_stability_type = GLS_StabilityType::Galperin;
            }
            else {
                amrex::Abort("Don't know this GLS stability type");
            }
        } else {
            amrex::Abort("Don't know this vertical mixing type");
        }
        // Read in GLS params
        if (vert_mixing_type == VertMixingType::GLS) {
            pp.queryAdd("gls_P", gls_p);
            pp.queryAdd("gls_M", gls_m);
            pp.queryAdd("gls_N", gls_n);
            pp.queryAdd("gls_Kmin", gls_Kmin);
            pp.queryAdd("gls_Pmin", gls_Pmin);

            pp.queryAdd("gls_cmu0", gls_cmu0);
            pp.queryAdd("gls_c1", gls_c1);
            pp.queryAdd("gls_c2", gls_c2);
            pp.queryAdd("gls_c3m", gls_c3m);
            pp.queryAdd("gls_c3p", gls_c3p);
            pp.queryAdd("gls_sigk", gls_sigk);
            pp.queryAdd("gls_sigp", gls_sigp);
            if (gls_stability_type == GLS_StabilityType::Canuto_A) {
                gls_Gh0 = amrex::Real(0.0329); // 0.0329 GOTM, 0.0673 Burchard
                gls_Ghcri = amrex::Real(0.03);
                gls_L1 = amrex::Real(0.107);
                gls_L2 = amrex::Real(0.0032);
                gls_L3 = amrex::Real(0.0864);
                gls_L4 = amrex::Real(0.12);
                gls_L5 = amrex::Real(11.9);
                gls_L6 = amrex::Real(0.4);
                gls_L7 = zero;
                gls_L8 = amrex::Real(0.48);
            } else if (gls_stability_type == GLS_StabilityType::Canuto_B) {
                gls_Gh0 = amrex::Real(0.0444); // 0.044 GOTM, 0.0673 Burchard
                gls_Ghcri = amrex::Real(0.0414);
                gls_L1 = amrex::Real(0.127);
                gls_L2 = amrex::Real(0.00336);
                gls_L3 = amrex::Real(0.0906);
                gls_L4 = amrex::Real(0.101);
                gls_L5 = amrex::Real(11.2);
                gls_L6 = amrex::Real(0.4);
                gls_L7 = zero;
                gls_L8 = amrex::Real(0.318);
            } else {
                gls_Gh0 = amrex::Real(0.028);
                gls_Ghcri = amrex::Real(0.02);
            }
        }

        // Read and compute inverse nudging coeffs from inputs given in days,
        // and store in a vector corresponding to BdyVars enum
        amrex::Real tnudg =  zero;
        amrex::Real znudg =  zero;
        amrex::Real m2nudg = zero;
        amrex::Real m3nudg = zero;
        pp.queryAdd("tnudg",tnudg);
        pp.queryAdd("znudg",znudg);
        pp.queryAdd("m2nudg",m2nudg);
        pp.queryAdd("m3nudg",m3nudg);
        pp.queryAdd("obcfac",obcfac);

        nudg_coeff.resize(BdyVars::NumTypes(ncons));
        nudg_coeff[BdyVars::u          ] = (m3nudg > zero) ? one / (m3nudg * amrex::Real(86400.0)) : zero;//BdyVars::u
        nudg_coeff[BdyVars::v          ] = (m3nudg > zero) ? one / (m3nudg * amrex::Real(86400.0)) : zero;//BdyVars::v
        // Every cell-centered tracer -- temp, salt, and any additional passive or
        // biology scalar -- nudges on the single tracer timescale tnudg, as in ROMS
        // where Tnudg defaults to the same value for all tracers.
        for (int icomp = 0; icomp < ncons; ++icomp) {
            nudg_coeff[BdyVars::cons(icomp)] = (tnudg > zero) ? one / (tnudg * amrex::Real(86400.0)) : zero;
        }
        nudg_coeff[BdyVars::ubar(ncons)] = (m2nudg > zero) ? one / (m2nudg * amrex::Real(86400.0)) : zero;//BdyVars::ubar
        nudg_coeff[BdyVars::vbar(ncons)] = (m2nudg > zero) ? one / (m2nudg * amrex::Real(86400.0)) : zero;//BdyVars::vbar
        nudg_coeff[BdyVars::zeta(ncons)] = ( znudg > zero) ? one / ( znudg * amrex::Real(86400.0)) : zero;//BdyVars::zeta

        pp.queryAdd("do_m3_clim_nudg", do_m3_clim_nudg);
        pp.queryAdd("do_m2_clim_nudg", do_m2_clim_nudg);

        // Every cell-centered tracer -- temp, salt, and any additional passive or biology
        // scalar -- can be nudged toward its own climatology. The flag is keyed by the
        // variable's own name, so temp and salt keep do_temp_clim_nudg / do_salt_clim_nudg
        // and a biology tracer uses e.g. do_NO3_clim_nudg.
        do_cons_clim_nudg.assign(ncons, 0);
        for (int icomp = 0; icomp < ncons; ++icomp) {
            bool do_nudg = false;
            pp.queryAdd(("do_" + cons_names[icomp] + "_clim_nudg").c_str(), do_nudg);
            do_cons_clim_nudg[icomp] = do_nudg ? 1 : 0;
            if (do_nudg) { do_any_cons_clim_nudg = true; }
        }

        if (do_m3_clim_nudg || do_m2_clim_nudg || do_any_cons_clim_nudg) {
            do_any_clim_nudg = true;
        }
        // The climatology series are only built along the NetCDF initialization path,
        // so asking for nudging with analytic initial conditions would leave them
        // unallocated. Say so here rather than failing later inside the time step.
        if (do_any_clim_nudg && ic_type != IC_Type::netcdf) {
            amrex::Abort("Climatology nudging requires remora.ic_type = netcdf");
        }
#ifndef REMORA_USE_NETCDF
        if (do_any_clim_nudg) {
            amrex::Abort("Climatology nudging requires building with NetCDF");
        }
#endif
    }

    void display()
    {
        amrex::Print() << "SOLVER CHOICE: " << std::endl;
        amrex::Print() << "use_salt              : " << use_salt     << std::endl;
        amrex::Print() << "use_coriolis          : " << use_coriolis << std::endl;
        amrex::Print() << "use_prestep           : " << use_prestep << std::endl;
        amrex::Print() << "use_uv3dmix           : " << use_uv3dmix << std::endl;
        amrex::Print() << "use_barotropic        : " << use_barotropic << std::endl;
        amrex::Print() << "flat_bathymetry       : " << flat_bathymetry << std::endl;
        amrex::Print() << "spatial_order         : " << spatial_order << std::endl;

        if (ic_type == IC_Type::analytic) {
            amrex::Print() << "Using analytic initial onditions" << std::endl;
        }
        else if (ic_type == IC_Type::netcdf) {
            amrex::Print() << "Using NetCDF initial conditions" << std::endl;
        }

        if (tracer_Hadv_scheme == AdvectionScheme::centered4) {
            amrex::Print() << "Horizontal advection scheme for tracers: " << "Centered 4" << std::endl;
        }
        else if (tracer_Hadv_scheme == AdvectionScheme::upstream3) {
            amrex::Print() << "Horizontal advection scheme for tracers: " << "Upstream 3" << std::endl;
        }
        else {
            amrex::Error("Invalid horizontal advection scheme for tracers.");
        }

        if (uv_Hadv_scheme == AdvectionScheme::centered2) {
            amrex::Print() << "Horizontal advection scheme for momenta: " << "Centered 2" << std::endl;
        }
        else if (uv_Hadv_scheme == AdvectionScheme::upstream3) {
            amrex::Print() << "Horizontal advection scheme for momenta: " << "Upstream 3" << std::endl;
        }
        else {
            amrex::Error("Invalid horizontal advection scheme for momenta.");
        }

        if (coupling_type == CouplingType::two_way) {
            amrex::Print() << "Using two-way coupling " << std::endl;
        } else if (coupling_type == CouplingType::one_way) {
            amrex::Print() << "Using one-way coupling " << std::endl;
        }

        if (use_coriolis) {
            if (coriolis_type == Cor_Type::analytic) {
                amrex::Print() << "Using analytic coriolis forcing " << std::endl;
            } else if (coriolis_type == Cor_Type::beta_plane) {
                amrex::Print() << "Using beta plane coriolis forcing " << std::endl;
            } else if (coriolis_type == Cor_Type::netcdf) {
                amrex::Print() << "Using coriolis forcing loaded from file " << std::endl;
            }
        }
    }

    // Default prefix
    std::string pp_prefix {"remora"};

    bool        flat_bathymetry        = false;

    bool        use_salt               = true;

    // Specify what additional physics/forcing modules we use
    bool        use_coriolis           = false;

    // Specify whether terms are used for debugging purposes
    bool        use_prestep           = true;
    bool        use_uv3dmix           = true;
    bool        use_baroclinic        = true;
    bool        use_barotropic        = true;

    bool        use_curvilinear_grid  = false;

    bool        bulk_fluxes           = false;
    bool        atm2ocn_flux_mode     = false;

    bool        output_forcing        = false;
    bool        do_temp_flux          = false;
    bool        do_salt_flux          = false;
    bool        longwave_down         = false;

    std::array<BulkForcingType, BulkFlux::NumTypes> bulk_flux_type {{
        BulkForcingType::analytic, // Uwind
        BulkForcingType::analytic, // Vwind
        BulkForcingType::constant, // Pair
        BulkForcingType::constant, // Qair
        BulkForcingType::constant, // Tair
        BulkForcingType::constant, // Cloud
        BulkForcingType::constant, // Rain
        BulkForcingType::constant, // SWrad
        BulkForcingType::computed, // LWrad defaults to the internal longwave formula
        BulkForcingType::computed  // EminusP defaults to bulk evap-rain diagnostic
    }};
    bool        longwave_is_net      = false;
    bool        qair_is_percent       = false;
    std::string longwave_netcdf_varname = "lwrad";


    bool        do_rivers             = false;
    bool        do_rivers_temp        = true;
    bool        do_rivers_salt        = true;
    bool        do_rivers_scalar      = false;
    amrex::Vector<int> do_rivers_cons;

    bool init_l1ad_T = false;

    bool init_ana_T = false;

    bool init_l0int_T = true;

    AdvectionScheme tracer_Hadv_scheme;
    AdvectionScheme uv_Hadv_scheme;

    // Coupling options are "OneWay" or "TwoWay"
    CouplingType coupling_type;

    // IC and BC Type: "analytic" or "netcdf"
    IC_Type ic_type;

    // Coriolis forcing type
    Cor_Type coriolis_type;

    // Surface momentum flux type
    SMFluxType smflux_type;

    // Surface wind speed type
    WindType wind_type;

    // EOS type
    EOSType eos_type;

    // Bottom stress type
    BottomStressType bottom_stress_type;

    // Land/sea mask type
    MaskType mask_type;

    // Mixing type and parameters
    VertMixingType vert_mixing_type;
    HorizMixingType horiz_mixing_type;
    ScaledToGridAMRScaling scaled_to_grid_amr_scaling = ScaledToGridAMRScaling::none;
    GLS_StabilityType gls_stability_type;
    HarmonicMixingType harmonic_mixing_type;

    // Type for grid scale (pm and pn)
    GridScaleType grid_scale_type;

    // Stretching and depth parameters which may need to be read from inputs
    amrex::Real theta_s = amrex::Real(3.0);
    amrex::Real theta_b = zero;
    amrex::Real tcline  = amrex::Real(150.0);

    // Linear drag coefficient [m/s]
    amrex::Real rdrag = amrex::Real(3e-4);
    // Quadratic drag coefficient [dimensionless]
    amrex::Real rdrag2 = amrex::Real(3e-3);

    // Momentum stress scales [m]
    amrex::Real Zob = amrex::Real(2e-2);
    amrex::Real Zos = amrex::Real(2e-2);

    amrex::Real Cdb_max = amrex::Real(0.5);
    amrex::Real Cdb_min = amrex::Real(1e-6);

    // Linear equation of state parameters
    amrex::Real R0 = amrex::Real(1028); // background density value (Kg/m3) used in Linear Equation of State
    amrex::Real S0 = amrex::Real(35.0); // background salinity (nondimensional) constant
    amrex::Real T0 = amrex::Real(5.0); // background potential temperature (Celsius) constant
    amrex::Real Tcoef = amrex::Real(1.7e-4); // linear equation of state parameter (1/Celsius)
    amrex::Real Scoef = zero; // linear equation of state parameter (nondimensional)
    amrex::Real rho0 = amrex::Real(1025.0); // Mean density (Kg/m3) used when Boussinesq approx is inferred

    // remora.time_ref: reference date of the model clock (yyyymmdd.dd), or one
    // of the special values 0, -1, -2. Selects the calendar; see
    // Source/Utils/REMORA_DateClock.H.
    amrex::Real time_ref = amrex::Real(0.0);

    // Coriolis forcing
    amrex::Real coriolis_f0 = zero; // f-plane constant (1/s)
    amrex::Real coriolis_beta = zero; // beta-plane constant (1/s/m)

    // Air pressure
    amrex::Real Pair = amrex::Real(1013.48);
    // Air temperature
    amrex::Real Tair = amrex::Real(23.567);
    // Relative humidity (air)
    amrex::Real Hair = amrex::Real(0.776);
    // Cloud cover fraction (0=clear sky, 1=overcast)
    amrex::Real cloud = zero;
    // Precipitation rate (kg/m2/s)
    amrex::Real rain = zero;
    // Height (m) of atmospheric measurements for Bulk fluxes parametrization
    amrex::Real blk_ZQ = amrex::Real(10.0); // air humidity
    amrex::Real blk_ZT = amrex::Real(10.0); // air temperature
    amrex::Real blk_ZW = amrex::Real(10.0); // winds

    bool eminusp = false;
    bool eminusp_correct_ssh = false;

    // Surface radiation flux
    amrex::Real srflux = zero;
    // Surface wind speed
    amrex::Real Uwind = zero;
    amrex::Real Vwind = zero;
    // External longwave radiation flux
    amrex::Real longwave_rad = zero;
    // Prescribed evaporation minus precipitation
    amrex::Real EminusP = zero;

    // Spatial discretization
    int         spatial_order = 2;

    // Horizontal mixing parameters
    amrex::Real visc2 = zero;
    amrex::Vector<amrex::Real> tnu2;

    // GLS params
    amrex::Real gls_p = amrex::Real(3.0);
    amrex::Real gls_m = amrex::Real(1.5);
    amrex::Real gls_n = amrex::Real(-1.0);
    amrex::Real gls_Kmin = amrex::Real(7.6e-6);
    amrex::Real gls_Pmin = amrex::Real(1.0e-12);

    amrex::Real gls_cmu0 = amrex::Real(0.5477);
    amrex::Real gls_c1 = amrex::Real(1.44);
    amrex::Real gls_c2 = amrex::Real(1.92);
    amrex::Real gls_c3m = amrex::Real(-0.4);
    amrex::Real gls_c3p = one;
    amrex::Real gls_sigk = one;
    amrex::Real gls_sigp = amrex::Real(1.3);

    // Turbulence closure
    amrex::Real Akk_bak = amrex::Real(5.0e-6);
    amrex::Real Akp_bak = amrex::Real(5.0e-6);
    amrex::Real Akv_bak = amrex::Real(5.0e-6);
    // remora.Akt_bak: the value both active tracers take unless overridden by name.
    amrex::Real Akt_bak_all = amrex::Real(1.0e-6);
    // NAT entries -- temperature and salinity -- to match the components of vec_Akt.
    // Passive tracers mix with entry Salt_comp; see akt_comp().
    amrex::Vector<amrex::Real> Akt_bak;

    // Params for stability functions.
    amrex::Real gls_Gh0;
    amrex::Real gls_Ghcri;
    amrex::Real gls_Ghmin = amrex::Real(-0.28);
    amrex::Real gls_E2 = amrex::Real(1.33);
    // Params only for Canuto stability
    amrex::Real gls_L1;
    amrex::Real gls_L2;
    amrex::Real gls_L3;
    amrex::Real gls_L4;
    amrex::Real gls_L5;
    amrex::Real gls_L6;
    amrex::Real gls_L7;
    amrex::Real gls_L8;

    // Params for some GLS and also Mellor-Yamada
    amrex::Real my_A1 = amrex::Real(0.92);
    amrex::Real my_A2 = amrex::Real(0.74);
    amrex::Real my_B1 = amrex::Real(16.6);
    amrex::Real my_B2 = amrex::Real(10.1);
    amrex::Real my_C1 = amrex::Real(0.08);
    amrex::Real my_C2 = amrex::Real(0.7);
    amrex::Real my_C3 = amrex::Real(0.2);
    amrex::Real my_E1 = amrex::Real(1.8);
    amrex::Real my_E2 = amrex::Real(1.33);
    amrex::Real my_Gh0 = amrex::Real(0.0233);
    amrex::Real my_Sq = amrex::Real(0.2);
    amrex::Real my_dtfac = amrex::Real(0.05);
    amrex::Real my_lmax = amrex::Real(0.53);
    amrex::Real my_qmin = amrex::Real(1.0E-8);

    // Nudging time scales in 1/s
    amrex::Vector<amrex::Real> nudg_coeff;

    // Factor between passive (outflow) and active (inflow) open boundary
    // conditions.
    amrex::Real obcfac = zero;

    // Whether to do climatoogy nudging
    bool do_m2_clim_nudg = false;
    bool do_m3_clim_nudg = false;
    // Per cell-centered tracer, whether to nudge toward climatology. Indexed by
    // cons component, so entry Temp_comp is the old do_temp_clim_nudg.
    amrex::Vector<int> do_cons_clim_nudg;
    bool do_any_cons_clim_nudg = false;
    bool do_any_clim_nudg = false;

    bool boundary_from_netcdf = false;
};

#endif
