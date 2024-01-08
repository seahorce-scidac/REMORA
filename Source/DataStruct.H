#ifndef _DATA_STRUCT_H_
#define _DATA_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

#include <REMORA_Constants.H>

enum struct CouplingType {
    OneWay, TwoWay
};

enum class Coord {
    x, y, z
};

enum class AdvectionDir {
    x, y, z
};

enum class AdvectionScheme {
    centered4, upstream3
};

enum class IC_BC_Type {
    Custom, Real
};

enum class Cor_Type {
    Custom, Beta_Plane, Real
};

struct SolverChoice {
  public:
    void init_params()
    {
        amrex::ParmParse pp(pp_prefix);

        pp.query("flat_bathymetry", flat_bathymetry);

        // Which horizontal advection scheme
        static std::string hadv_string = "upstream3";
        pp.query("horizontal_advection_scheme",hadv_string);
        if (hadv_string == "centered4")
            Hadv_scheme = AdvectionScheme::centered4;
        else if (hadv_string == "upstream3")
            Hadv_scheme = AdvectionScheme::upstream3;
        else
            amrex::Error("Advection scheme unknown.");

        pp.query("rdrag", rdrag);

        // Order of spatial discretization
        pp.query("spatial_order", spatial_order);

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

        //Read in linear eos parameters
        pp.query("R0",R0);
        pp.query("S0",S0);
        pp.query("T0",T0);
        pp.query("Tcoef",Tcoef);
        pp.query("Scoef",Scoef);
        pp.query("rho0", rho0);

        //Grid stretching
        pp.query("theta_s",theta_s);
        pp.query("theta_b",theta_b);
        pp.query("tcline",tcline);

        //coriolis factor
        pp.query("coriolis_f0",coriolis_f0);
        pp.query("coriolis_beta",coriolis_beta);

        pp.query("ggrav",g);

        static std::string ic_bc_type_string = "custom";
        pp.query("ic_bc_type", ic_bc_type_string);

        if ( !ic_bc_type_string.compare("Custom") ||
             !ic_bc_type_string.compare("custom") ) {
            ic_bc_type = IC_BC_Type::Custom;
        } else if ( !ic_bc_type_string.compare("Real") ||
                    !ic_bc_type_string.compare("real") ) {
            ic_bc_type = IC_BC_Type::Real;
        } else {
            amrex::Error("Don't know this ic_bc_type");
        }

        // Which type of refinement
        static std::string coupling_type_string = "TwoWay";
        pp.query("coupling_type",coupling_type_string);
        if (coupling_type_string == "TwoWay") {
            coupling_type = CouplingType::TwoWay;
        } else if (coupling_type_string == "OneWay") {
            coupling_type = CouplingType::OneWay;
        } else {
            amrex::Abort("Dont know this coupling_type");
        }

        // Which type of coriolis forcing
        if (use_coriolis) {
            static std::string coriolis_type_string = "beta_plane";
            pp.query("coriolis_type",coriolis_type_string);
            if ( (coriolis_type_string == "Custom") || (coriolis_type_string == "custom") ) {
                coriolis_type = Cor_Type::Custom;
            } else if ((coriolis_type_string == "Beta_Plane") ||
                       (coriolis_type_string == "beta_Plane") ||
                       (coriolis_type_string == "beta_plane") ||
                       (coriolis_type_string == "Beta_plane")) {
                coriolis_type = Cor_Type::Beta_Plane;
            } else if ( (coriolis_type_string == "Real") || (coriolis_type_string == "real") ) {
                coriolis_type = Cor_Type::Real;
            } else {
                amrex::Abort("Don't know this coriolis_type");
            }
        }
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

        if (ic_bc_type == IC_BC_Type::Custom) {
            amrex::Print() << "Using custom initial and boundary conditions (No mesoscale forcing!)" << std::endl;
        }
        else if (ic_bc_type == IC_BC_Type::Real) {
            amrex::Print() << "Using REAL initial and boundary conditions (Mesoscale forcing!)" << std::endl;
        }

        if (Hadv_scheme == AdvectionScheme::centered4) {
            amrex::Print() << "Horizontal advection scheme: " << "Centered 4" << std::endl;
        }
        else if (Hadv_scheme == AdvectionScheme::upstream3) {
            amrex::Print() << "Horizontal advection scheme: " << "Upstream 3" << std::endl;
        }
        else {
            amrex::Error("Invalid horizontal advection scheme.");
        }

        if (coupling_type == CouplingType::TwoWay) {
            amrex::Print() << "Using two-way coupling " << std::endl;
        } else if (coupling_type == CouplingType::OneWay) {
            amrex::Print() << "Using one-way coupling " << std::endl;
        }

        if (use_coriolis) {
            if (coriolis_type == Cor_Type::Custom) {
                amrex::Print() << "Using custom coriolis forcing " << std::endl;
            } else if (coriolis_type == Cor_Type::Beta_Plane) {
                amrex::Print() << "Using beta plane coriolis forcing " << std::endl;
            } else if (coriolis_type == Cor_Type::Real) {
                amrex::Print() << "Using coriolis forcing loaded from file " << std::endl;
            }
        }
    }

    // Default prefix
    std::string pp_prefix {"remora"};

    bool        flat_bathymetry        = true;

    bool        use_salt               = true;

    // Specify what additional physics/forcing modules we use
    bool        use_coriolis           = false;

    // Specify whether terms are used for debugging purposes
    bool        use_prestep           = true;
    bool        use_uv3dmix           = true;
    bool        use_baroclinic        = true;
    bool        use_barotropic        = true;

    AdvectionScheme Hadv_scheme;

    // Coupling options are "OneWay" or "TwoWay"
    CouplingType coupling_type;

    // IC and BC Type: "Custom" or "Real"
    IC_BC_Type ic_bc_type;

    // Coriolis forcing type
    Cor_Type coriolis_type;

    // Stretching and depth parameters which may need to be read from inputs
    amrex::Real theta_s = 3.0;
    amrex::Real theta_b = 0.0;
    amrex::Real tcline  = 150.0;

    // Linear drag coefficient
    amrex::Real rdrag = 3e-4;

    // Linear equation of state parameters
    amrex::Real R0 = 1028; // background density value (Kg/m3) used in Linear Equation of State
    amrex::Real S0 = 35.0; // background salinity (nondimensional) constant
    amrex::Real T0 = 5.0; // background potential temperature (Celsius) constant
    amrex::Real Tcoef = 1.7e-4; // linear equation of state parameter (1/Celcius)
    amrex::Real Scoef = 0.0; // linear equation of state parameter (nondimensional)
    amrex::Real rho0 = 1025.0; // Mean density (Kg/m3) used when Boussinesq approx is inferred

    // Coriolis forcing
    amrex::Real coriolis_f0 = 0.0; // f-plane constant (1/s)
    amrex::Real coriolis_beta = 0.0; // beta-plane constant (1/s/m)

    // Gravity
    amrex::Real g = 9.81; // gravitational field strength [kg m/s^2]

    // Spatial discretization
    int         spatial_order = 2;
};

#endif
