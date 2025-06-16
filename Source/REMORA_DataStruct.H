#ifndef _REMORA_DATA_STRUCT_H_
#define _REMORA_DATA_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

#include <REMORA_Constants.H>
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

/** \brief Type of initial/boundary condition type. Analytic reads from prob.cpp. Netcdf is from file */
enum class IC_BC_Type {
    analytic, netcdf
};

/** \brief Coriolis factor */
enum class Cor_Type {
    analytic, beta_plane, netcdf
};

/** \brief plotfile format */
enum class PlotfileType {
    amrex, netcdf, hdf5
};

/** \brief vertical mixing type */
enum class VertMixingType {
    analytic, GLS
};

/** \brief horizontal viscosity/diffusion type */
enum class HorizMixingType {
    analytic, constant
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

struct SolverChoice {
  public:
    /** \brief read in and initialize parameters */
    void init_params()
    {
        amrex::ParmParse pp(pp_prefix);

        pp.query("flat_bathymetry", flat_bathymetry);

        // Which horizontal advection scheme for tracers
        static std::string tracer_hadv_string = "upstream3";
        pp.query("tracer_horizontal_advection_scheme",tracer_hadv_string);
        if (tracer_hadv_string == "centered4")
            tracer_Hadv_scheme = AdvectionScheme::centered4;
        else if (tracer_hadv_string == "upstream3")
            tracer_Hadv_scheme = AdvectionScheme::upstream3;
        else
            amrex::Error("Advection scheme unknown.");

        // Which horizontal advection scheme
        static std::string uv_hadv_string = "upstream3";
        pp.query("uv_horizontal_advection_scheme",uv_hadv_string);
        if (uv_hadv_string == "upstream3")
            uv_Hadv_scheme = AdvectionScheme::upstream3;
        else if (uv_hadv_string == "centered2")
            uv_Hadv_scheme = AdvectionScheme::centered2;
        else
            amrex::Error("UV advection scheme unknown.");

        pp.query("rdrag", rdrag);
        pp.query("rdrag2", rdrag2);
        pp.query("Zos", Zos);
        pp.query("Zob", Zob);
        pp.query("Cdb_max", Cdb_max);
        pp.query("Cdb_min", Cdb_min);

        // Include salinity?
        pp.query("use_salt", use_salt);

        // Include Coriolis forcing?
        pp.query("use_coriolis", use_coriolis);

        // Include prestep / lagged predictor / corrections
        pp.query("use_prestep", use_prestep);

        //This affect forcing and some of the mixing terms for velocity
        pp.query("use_uv3dmix", use_uv3dmix);

        //This accounts for the main 2d loop but may not include coupling and copying properly
        pp.query("use_barotropic", use_barotropic);

        // Whether to do rivers. By default, rivers are temp and salt sources. Rivers
        // have to be momentum sources.
        pp.query("do_rivers", do_rivers);

        pp.query("do_rivers_temp", do_rivers_temp);
        pp.query("do_rivers_salt", do_rivers_salt);
        pp.query("do_rivers_scalar", do_rivers_scalar);

        do_rivers_cons.resize(NCONS);
        do_rivers_cons[Temp_comp  ] = do_rivers_temp;
        do_rivers_cons[Salt_comp  ] = do_rivers_salt;
        do_rivers_cons[Scalar_comp] = do_rivers_scalar;

        pp.query("init_l1ad_h", init_l1ad_h);
        pp.query("init_l1ad_T", init_l1ad_T);

        pp.query("init_ana_h", init_ana_h);
        pp.query("init_ana_T", init_ana_T);

        pp.query("init_l0int_h", init_l0int_h);
        pp.query("init_l0int_T", init_l0int_T);

        static std::string eos_type_string = "linear";
        pp.query("eos_type",eos_type_string);
        if (eos_type_string == "linear" || eos_type_string == "Linear" ||
              eos_type_string == "lin"    || eos_type_string == "Lin") {
            eos_type = EOSType::linear;
            pp.query("Tcoef",Tcoef);
            pp.query("Scoef",Scoef);
        } else if (eos_type_string == "nonlinear" || eos_type_string == "Nonlinear" ||
                   eos_type_string == "non-linear" || eos_type_string == "Non-linear" ||
                   eos_type_string == "nonlin"     || eos_type_string == "Nonlin") {
            eos_type = EOSType::nonlinear;
        } else {
            amrex::Abort("Dont know this eos_type");
        }
        pp.query("R0",R0);
        pp.query("S0",S0);
        pp.query("T0",T0);
        pp.query("rho0", rho0);

        pp.query("bulk_fluxes",bulk_fluxes);
        if (bulk_fluxes) {
            do_salt_flux = true;
            do_temp_flux = true;
        }
        pp.query("air_pressure",Pair);
        pp.query("air_temperature",Tair);
        pp.query("air_humidity",Hair);
        pp.query("surface_radiation_flux",srflux);
        pp.query("cloud",cloud);
        pp.query("rain",rain);
        pp.query("blk_ZQ",blk_ZQ);
        pp.query("blk_ZT",blk_ZT);
        pp.query("blk_ZW",blk_ZW);
        pp.query("eminusp",eminusp);
        pp.query("eminusp_correct_ssh",eminusp_correct_ssh);

        if (eminusp_correct_ssh and !eminusp) {
            amrex::Abort("If evaporation minus precipitation (E-P) sea surface height correct in is on, E-P must be on as well (remora.eminusp=true)");
        }
        if (eminusp and !bulk_fluxes) {
            amrex::Abort("Evaporation minus precipitation (E-P) requires bulk flux parametrizations (remora.bulk_fluxes=true)");
        }

        //Read in linear eos parameters
        //Grid stretching
        pp.query("theta_s",theta_s);
        pp.query("theta_b",theta_b);
        pp.query("tcline",tcline);

        //coriolis factor
        pp.query("coriolis_f0",coriolis_f0);
        pp.query("coriolis_beta",coriolis_beta);

        pp.query("Akv_bak",Akv_bak);
        pp.query("Akt_bak",Akt_bak);


        static std::string grid_scale_type_string = "constant";
        pp.query("grid_scale_type",grid_scale_type_string);

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

        static std::string ic_bc_type_string = "analytic";
        pp.query("ic_bc_type", ic_bc_type_string);

        if ( amrex::toLower(ic_bc_type_string) == "custom") {
            amrex::Warning("Problem initialization from prob.cpp is now called 'analytic'. 'custom' will be deprecated");
            ic_bc_type = IC_BC_Type::analytic;
        } else if ( amrex::toLower(ic_bc_type_string) == "analytic") {
            ic_bc_type = IC_BC_Type::analytic;
        } else if ( amrex::toLower(ic_bc_type_string) == "netcdf") {
            ic_bc_type = IC_BC_Type::netcdf;
        } else if ( amrex::toLower(ic_bc_type_string) == "real") {
            amrex::Warning("Problem initialization from NetCDF (remora.ic_bc_type) is now called 'netcdf'. 'real' will be deprecated");
            ic_bc_type = IC_BC_Type::netcdf;
        } else {
            amrex::Error("Don't know this ic_bc_type");
        }

        // Which type of refinement
        static std::string coupling_type_string = "TwoWay";
        pp.query("coupling_type",coupling_type_string);
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
            pp.query("coriolis_type",coriolis_type_string);
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
        int smflux_specified = pp.query("smflux_type",smflux_type_string);
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
        int wind_specified = pp.query("wind_type",wind_type_string);
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

        if (wind_specified && smflux_specified) {
            amrex::Abort("Cannot specify both wind and surface momentum flux");
        }

        static std::string bottom_stress_type_string = "linear";
        pp.query("bottom_stress_type", bottom_stress_type_string);
        if (amrex::toLower(bottom_stress_type_string) == "linear") {
            bottom_stress_type = BottomStressType::linear;
        } else if (amrex::toLower(bottom_stress_type_string) == "quadratic") {
            bottom_stress_type = BottomStressType::quadratic;
        } else if (amrex::toLower(bottom_stress_type_string) == "logarithmic") {
            bottom_stress_type = BottomStressType::logarithmic;
        } else {
            amrex::Abort("Don't know this bottom_stress_type");
        }

        amrex::Real tnu2_salt = amrex::Real(0.0);
        amrex::Real tnu2_temp = amrex::Real(0.0);
        amrex::Real tnu2_scalar = amrex::Real(0.0);
        static std::string horiz_mixing_type_string = "analytic";
        pp.query("horizontal_mixing_type", horiz_mixing_type_string);
        if (amrex::toLower(horiz_mixing_type_string) == "analytical" ||
            amrex::toLower(horiz_mixing_type_string) == "analytic") {
            horiz_mixing_type = HorizMixingType::analytic;
        } else if (amrex::toLower(horiz_mixing_type_string) == "constant") {
            horiz_mixing_type = HorizMixingType::constant;
        } else {
            amrex::Abort("Don't know this horizontal mixing type");
        }
        pp.query("visc2",visc2);
        pp.query("tnu2_salt",tnu2_salt);
        pp.query("tnu2_temp",tnu2_temp);
        pp.query("tnu2_scalar",tnu2_scalar);
        tnu2.resize(NCONS);
        if (use_salt) {
            tnu2[0] = tnu2_temp;
            tnu2[1] = tnu2_salt;
            tnu2[2] = tnu2_scalar;
        } else {
            tnu2[0] = tnu2_temp;
            tnu2[1] = tnu2_scalar;
        }

        pp.query("Akk_bak", Akk_bak);
        pp.query("Akp_bak", Akp_bak);
        static std::string vert_mixing_type_string = "analytic";
        static std::string gls_stability_type_string = "Canuto_A";
        pp.query("vertical_mixing_type", vert_mixing_type_string);
        pp.query("gls_stability_type", gls_stability_type_string);
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
            pp.query("gls_P", gls_p);
            pp.query("gls_M", gls_m);
            pp.query("gls_N", gls_n);
            pp.query("gls_Kmin", gls_Kmin);
            pp.query("gls_Pmin", gls_Pmin);

            pp.query("gls_cmu0", gls_cmu0);
            pp.query("gls_c1", gls_c1);
            pp.query("gls_c2", gls_c2);
            pp.query("gls_c3m", gls_c3m);
            pp.query("gls_c3p", gls_c3p);
            pp.query("gls_sigk", gls_sigk);
            pp.query("gls_sigp", gls_sigp);
            if (gls_stability_type == GLS_StabilityType::Canuto_A) {
                gls_Gh0 = amrex::Real(0.0329); // 0.0329 GOTM, 0.0673 Burchard
                gls_Ghcri = amrex::Real(0.03);
                gls_L1 = amrex::Real(0.107);
                gls_L2 = amrex::Real(0.0032);
                gls_L3 = amrex::Real(0.0864);
                gls_L4 = amrex::Real(0.12);
                gls_L5 = amrex::Real(11.9);
                gls_L6 = amrex::Real(0.4);
                gls_L7 = amrex::Real(0.0);
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
                gls_L7 = amrex::Real(0.0);
                gls_L8 = amrex::Real(0.318);
            } else {
                gls_Gh0 = amrex::Real(0.028);
                gls_Ghcri = amrex::Real(0.02);
            }
        }

        // Read and compute inverse nudging coeffs from inputs given in days,
        // and store in a vector corresponding to BdyVars enum
        amrex::Real tnudg =  amrex::Real(0.0);
        amrex::Real znudg =  amrex::Real(0.0);
        amrex::Real m2nudg = amrex::Real(0.0);
        amrex::Real m3nudg = amrex::Real(0.0);
        pp.query("tnudg",tnudg);
        pp.query("znudg",znudg);
        pp.query("m2nudg",m2nudg);
        pp.query("m3nudg",m3nudg);
        pp.query("obcfac",obcfac);

        nudg_coeff.resize(BdyVars::NumTypes);
        nudg_coeff[BdyVars::u    ] = (m3nudg > 0.0) ? amrex::Real(1.0) / (m3nudg * amrex::Real(86400.0)) : 0.0;//BdyVars::u
        nudg_coeff[BdyVars::v    ] = (m3nudg > 0.0) ? amrex::Real(1.0) / (m3nudg * amrex::Real(86400.0)) : 0.0;//BdyVars::v
        nudg_coeff[BdyVars::t    ] = ( tnudg > 0.0) ? amrex::Real(1.0) / ( tnudg * amrex::Real(86400.0)) : 0.0;//BdyVars::t
        nudg_coeff[BdyVars::s    ] = ( tnudg > 0.0) ? amrex::Real(1.0) / ( tnudg * amrex::Real(86400.0)) : 0.0;//BdyVars::s
        nudg_coeff[BdyVars::ubar ] = (m2nudg > 0.0) ? amrex::Real(1.0) / (m2nudg * amrex::Real(86400.0)) : 0.0;//BdyVars::ubar
        nudg_coeff[BdyVars::vbar ] = (m2nudg > 0.0) ? amrex::Real(1.0) / (m2nudg * amrex::Real(86400.0)) : 0.0;//BdyVars::vbar
        nudg_coeff[BdyVars::zeta ] = ( znudg > 0.0) ? amrex::Real(1.0) / ( znudg * amrex::Real(86400.0)) : 0.0;//BdyVars::zeta

        pp.query("do_m3_clim_nudg", do_m3_clim_nudg);
        pp.query("do_m2_clim_nudg", do_m2_clim_nudg);
        pp.query("do_temp_clim_nudg", do_temp_clim_nudg);
        pp.query("do_salt_clim_nudg", do_salt_clim_nudg);

        if (do_m3_clim_nudg || do_m2_clim_nudg || do_temp_clim_nudg || do_salt_clim_nudg) {
            do_any_clim_nudg = true;
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

        if (ic_bc_type == IC_BC_Type::analytic) {
            amrex::Print() << "Using analytic initial and boundary conditions" << std::endl;
        }
        else if (ic_bc_type == IC_BC_Type::netcdf) {
            amrex::Print() << "Using NetCDF initial and boundary conditions" << std::endl;
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

    bool        bulk_fluxes           = false;
    bool        do_temp_flux        = false;
    bool        do_salt_flux        = false;

    bool        do_rivers             = false;
    bool        do_rivers_temp        = true;
    bool        do_rivers_salt        = true;
    bool        do_rivers_scalar      = false;
    amrex::Vector<int> do_rivers_cons;

    bool init_l1ad_h = false;
    bool init_l1ad_T = false;

    bool init_ana_h = false;
    bool init_ana_T = false;

    bool init_l0int_h = true;
    bool init_l0int_T = true;

    AdvectionScheme tracer_Hadv_scheme;
    AdvectionScheme uv_Hadv_scheme;

    // Coupling options are "OneWay" or "TwoWay"
    CouplingType coupling_type;

    // IC and BC Type: "analytic" or "netcdf"
    IC_BC_Type ic_bc_type;

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

    // Mixing type and parameters
    VertMixingType vert_mixing_type;
    HorizMixingType horiz_mixing_type;
    GLS_StabilityType gls_stability_type;

    // Type for grid scale (pm and pn)
    GridScaleType grid_scale_type;

    // Stretching and depth parameters which may need to be read from inputs
    amrex::Real theta_s = amrex::Real(3.0);
    amrex::Real theta_b = amrex::Real(0.0);
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
    amrex::Real Scoef = amrex::Real(0.0); // linear equation of state parameter (nondimensional)
    amrex::Real rho0 = amrex::Real(1025.0); // Mean density (Kg/m3) used when Boussinesq approx is inferred

    // Coriolis forcing
    amrex::Real coriolis_f0 = amrex::Real(0.0); // f-plane constant (1/s)
    amrex::Real coriolis_beta = amrex::Real(0.0); // beta-plane constant (1/s/m)

    // Air pressure
    amrex::Real Pair = amrex::Real(1013.48);
    // Air temperature
    amrex::Real Tair = amrex::Real(23.567);
    // Relative humidity (air)
    amrex::Real Hair = amrex::Real(0.776);
    // Cloud cover fraction (0=clear sky, 1=overcast)
    amrex::Real cloud = amrex::Real(0.0);
    // Precipitation rate (kg/m2/s)
    amrex::Real rain = amrex::Real(0.0);
    // Height (m) of atmospheric measurements for Bulk fluxes parametrization
    amrex::Real blk_ZQ = amrex::Real(10.0); // air humidity
    amrex::Real blk_ZT = amrex::Real(10.0); // air temperature
    amrex::Real blk_ZW = amrex::Real(10.0); // winds

    bool eminusp = false;
    bool eminusp_correct_ssh = false;

    // Surface radiation flux
    amrex::Real srflux = amrex::Real(0.0);

    // Spatial discretization
    int         spatial_order = 2;

    // Horizontal mixing parameters
    amrex::Real visc2 = amrex::Real(0.0);
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
    amrex::Real gls_c3p = amrex::Real(1.0);
    amrex::Real gls_sigk = amrex::Real(1.0);
    amrex::Real gls_sigp = amrex::Real(1.3);

    // Turbulence closure
    amrex::Real Akk_bak = amrex::Real(5.0e-6);
    amrex::Real Akp_bak = amrex::Real(5.0e-6);
    amrex::Real Akv_bak = amrex::Real(5.0e-6);
    amrex::Real Akt_bak = amrex::Real(1.0e-6); // Note: this is a vector with one component per tracer in ROMS

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
    amrex::Real obcfac = amrex::Real(0.0);

    // Whether to do climatoogy nudging
    bool do_m2_clim_nudg = false;
    bool do_m3_clim_nudg = false;
    bool do_temp_clim_nudg = false;
    bool do_salt_clim_nudg = false;
    bool do_any_clim_nudg = false;
};

#endif
