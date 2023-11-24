#ifndef _DATA_STRUCT_H_
#define _DATA_STRUCT_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

#include <ROMSX_Constants.H>

enum struct CouplingType {
    OneWay, TwoWay
};

enum class Coord {
    x, y, z
};

enum class Stagger {
    None, x, y, z
};

enum class AdvectedQuantity {
    unity, u, v, w, theta, scalar
};

enum class AdvectingQuantity {
    rho_u, rho_v, rho_w
};

enum class AdvectionDir {
    x, y, z
};

enum class AdvectionScheme {
    centered4, upstream3
};

enum class DiffusionDir {
    x, y, z
};

enum class MomentumEqn {
    x, y, z
};

enum class MolecDiffType {
    None, Constant, ConstantAlpha
};

enum class IC_BC_Type {
    Ideal, Real
};

struct SolverChoice {
  public:
    void init_params()
    {
        amrex::ParmParse pp(pp_prefix);

        pp.query("flat_bathymetry", flat_bathymetry);

        if(!flat_bathymetry)  amrex::Warning("Flat bathymetry needs further testing");

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

        // Include Coriolis forcing?
        pp.query("use_coriolis", use_coriolis);

        // Include prestep / lagged predictor / corrections
        pp.query("use_prestep", use_prestep);

        //This affect forcing and some of the mixing terms for velocity
        pp.query("use_uv3dmix", use_uv3dmix);

        //This accounts for the main 2d loop but may not include coupling and copying properly
        pp.query("use_barotropic", use_barotropic);

        static std::string ic_bc_type_string = "Ideal";
        pp.query("ic_bc_type", ic_bc_type_string);

        if (!ic_bc_type_string.compare("Ideal")) {
            icBcType = IC_BC_Type::Ideal;
        } else if (!ic_bc_type_string.compare("Real")) {
            icBcType = IC_BC_Type::Real;
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
    }

    void display()
    {
        amrex::Print() << "SOLVER CHOICE: " << std::endl;
        amrex::Print() << "use_coriolis          : " << use_coriolis << std::endl;
        amrex::Print() << "spatial_order         : " << spatial_order << std::endl;

        if (icBcType == IC_BC_Type::Ideal) {
            amrex::Print() << "Using idealized initial and boundary conditions (No mesoscale forcing!)" << std::endl;
        }
        else if (icBcType == IC_BC_Type::Real) {
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
    }

    // Default prefix
    std::string pp_prefix {"romsx"};

    bool        flat_bathymetry        = true;

    // Specify what additional physics/forcing modules we use
    bool        use_coriolis           = false;

    // Specify whether terms are used for debugging purposes
    bool        use_prestep           = true;
    bool        use_uv3dmix           = true;
    bool        use_baroclinic        = true;
    bool        use_barotropic        = true;

    AdvectionScheme Hadv_scheme;

    CouplingType coupling_type;

    // Stretching and depth parameters which may need to be read from inputs
    amrex::Real theta_s = 3.0;
    amrex::Real theta_b = 0.0;
    amrex::Real tcline  = 150.0;

    // Linear drag coefficient
    amrex::Real rdrag = 3e-4;

    // Coriolis forcing
    amrex::Real coriolis_factor = 0.0;
    amrex::Real cosphi          = 0.0 ;
    amrex::Real sinphi          = 0.0;

    // Spatial discretization
    int         spatial_order = 2;

    // IC and BC Type: "Ideal" or "Real"
    IC_BC_Type icBcType;
};

#endif
