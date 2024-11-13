#include <AMReX_Vector.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ParmParse.H>

#include <REMORA.H>

using namespace amrex;

void REMORA::init_bcs ()
{
    auto f_set_var_bc = [this] (ParmParse& pp, int bcvar_type, Orientation ori, std::string bc_type_string) {
        if (bc_type_string == "symmetry")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::symmetry;
            domain_bc_type[ori] = "Symmetry";
        }
        else if (bc_type_string == "outflow")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::outflow;
            domain_bc_type[ori] = "Outflow";
        }
        else if (bc_type_string == "inflow")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::inflow;
            domain_bc_type[ori] = "Inflow";

            if (bcvar_type == BCVars::xvel_bc || bcvar_type == BCVars::yvel_bc ||
                bcvar_type == BCVars::zvel_bc) {
                std::vector<Real> v;
                pp.getarr("velocity", v, 0, AMREX_SPACEDIM);
                m_bc_extdir_vals[bcvar_type][ori] = v[bcvar_type - BCVars::xvel_bc];
            } else if (bcvar_type == BCVars::Scalar_bc_comp) {
                Real scalar_in = 0.;
                if (pp.query("scalar", scalar_in))
                m_bc_extdir_vals[BCVars::Scalar_bc_comp][ori] = scalar_in;
            }
        }
        else if (bc_type_string == "noslipwall")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::no_slip_wall;
            domain_bc_type[ori] = "NoSlipWall";

            if (bcvar_type == BCVars::xvel_bc || bcvar_type == BCVars::yvel_bc ||
                bcvar_type == BCVars::zvel_bc) {
                std::vector<Real> v;

                // The values of m_bc_extdir_vals default to 0.
                // But if we find "velocity" in the inputs file, use those values instead.
                if (pp.queryarr("velocity", v, 0, AMREX_SPACEDIM))
                {
                    v[ori.coordDir()] = 0.0_rt;
                    m_bc_extdir_vals[bcvar_type][ori] = v[bcvar_type - BCVars::xvel_bc];
                }
            }
        }
        else if (bc_type_string == "slipwall")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::slip_wall;
            domain_bc_type[ori] = "SlipWall";
        }
        else if (bc_type_string == "clamped")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::clamped;
            domain_bc_type[ori] = "Clamped";
        }
        else if (bc_type_string == "chapman")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::chapman;
            domain_bc_type[ori] = "Chapman";

            if (bcvar_type != BCVars::zeta_bc) {
                amrex::Abort("Chapman BC can only be applied to zeta");
            }
        }
        else if (bc_type_string == "flather")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::flather;
            domain_bc_type[ori] = "Flather";

            if (!(bcvar_type == BCVars::ubar_bc || bcvar_type == BCVars::vbar_bc)) {
                amrex::Abort("Flather BC can only be applied to ubar or vbar");
            }
        }
        else if (bc_type_string == "orlanski_rad")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::orlanski_rad;
            domain_bc_type[ori] = "Orlanski Radiation";
        }
        else if (bc_type_string == "orlanski_rad_nudg")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::orlanski_rad_nudge;
            domain_bc_type[ori] = "Orlanski Radiation with nudging";
        }
        else if (bc_type_string == "periodic")
        {
            if (!geom[0].isPeriodic(ori.coordDir())) {
                amrex::Abort("Periodic boundary specified in a non-periodic direction");
            }
            phys_bc_type[bcvar_type][ori] = REMORA_BC::periodic;
            domain_bc_type[ori] = "Periodic";
        }
        else
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::undefined;
        }

        if (geom[0].isPeriodic(ori.coordDir()))
        {
            domain_bc_type[ori] = "Periodic";
            if (phys_bc_type[bcvar_type][ori] == REMORA_BC::undefined)
            {
                phys_bc_type[bcvar_type][ori] = REMORA_BC::periodic;
            } else if (phys_bc_type[bcvar_type][ori] != REMORA_BC::periodic) {
                amrex::Abort("Wrong BC type for periodic boundary");
            }
        }

        if (phys_bc_type[bcvar_type][ori] == REMORA_BC::undefined)
        {
             amrex::Print() << "BC Type specified for fac is " << bc_type_string << std::endl;
             amrex::Abort("This BC type is unknown");
        }
    };

    auto f_by_side = [this, &f_set_var_bc] (std::string const& bcid, Orientation ori)
    {
        ParmParse pp("bc."+bcid);
        std::string bc_type_in = "null";
        // Default z directions to slipwall
        if (bcid=="zlo" or bcid=="zhi") bc_type_in = "slipwall";
        pp.query("type", bc_type_in);
        std::string bc_type = amrex::toLower(bc_type_in);

        for (int icomp=0; icomp<BCVars::NumTypes; icomp++) {
            f_set_var_bc(pp, icomp, ori, bc_type);
        }
    };

    auto f_by_var = [this, &f_set_var_bc] (std::string const& varname, int bcvar_type)
    {
        amrex::Vector<Orientation> orientations = {Orientation(Direction::x,Orientation::low), Orientation(Direction::y,Orientation::high),Orientation(Direction::x,Orientation::high),Orientation(Direction::y,Orientation::low)}; // west, south, east, north [matches ROMS]
        std::vector<std::string> bc_types = {"null","null","null","null"};
        ParmParse pp("bc."+varname);
        std::string bc_type_in = "null";
        // default zvel to outflow
        if (bcvar_type == BCVars::zvel_bc) {
            bc_types = {"outflow","outflow","outflow","outflow"};
        }
        pp.queryarr("type", bc_types);
        AMREX_ASSERT(bc_types.size() == 4);
        for (int i=0; i<4; i++) {
            std::string bc_type = amrex::toLower(bc_types[i]);
            auto ori = orientations[i];
            f_set_var_bc(pp, bcvar_type, ori, bc_type);
        }
    };

    for (OrientationIter oit; oit; ++oit) {
        Orientation ori = oit();
        // These are simply defaults for Dirichlet faces -- they should be over-written below if needed
        m_bc_extdir_vals[BCVars::Temp_bc_comp  ][ori] = 1.e19_rt;
        m_bc_extdir_vals[BCVars::Salt_bc_comp  ][ori] = 1.e20_rt;
        m_bc_extdir_vals[BCVars::Scalar_bc_comp][ori] = 1.e21_rt;

        m_bc_extdir_vals[BCVars::xvel_bc][ori] = 0.0_rt; // default
        m_bc_extdir_vals[BCVars::yvel_bc][ori] = 0.0_rt;
        m_bc_extdir_vals[BCVars::zvel_bc][ori] = 0.0_rt;

        m_bc_extdir_vals[BCVars::ubar_bc][ori] = 0.0_rt; // default
        m_bc_extdir_vals[BCVars::vbar_bc][ori] = 0.0_rt;
        m_bc_extdir_vals[BCVars::u2d_simple_bc][ori] = 0.0_rt;
        m_bc_extdir_vals[BCVars::v2d_simple_bc][ori] = 0.0_rt;
    }

    // Whether to specify boundary conditions by variable (then side).
    // Alternative is to do it by side by indicating keywords that indicate multiple variables
    set_bcs_by_var = false;

    ParmParse pp("remora");
    pp.query("boundary_per_variable", set_bcs_by_var);
    if (!set_bcs_by_var) {
        f_by_side("xlo", Orientation(Direction::x,Orientation::low));
        f_by_side("xhi", Orientation(Direction::x,Orientation::high));
        f_by_side("ylo", Orientation(Direction::y,Orientation::low));
        f_by_side("yhi", Orientation(Direction::y,Orientation::high));
    } else {
        f_by_var("temp", BCVars::Temp_bc_comp);
        f_by_var("salt", BCVars::Salt_bc_comp);
        f_by_var("scalar", BCVars::Scalar_bc_comp);
        f_by_var("u", BCVars::xvel_bc);
        f_by_var("v", BCVars::yvel_bc);
        f_by_var("w", BCVars::zvel_bc);
        f_by_var("ubar", BCVars::ubar_bc);
        f_by_var("vbar", BCVars::vbar_bc);
        f_by_var("zeta", BCVars::zeta_bc);
        f_by_var("tke", BCVars::tke_bc);
    }

    // Always specify z direction by side keyword
    f_by_side("zlo", Orientation(Direction::z,Orientation::low));
    f_by_side("zhi", Orientation(Direction::z,Orientation::high));

    // *****************************************************************************
    //
    // Here we translate the physical boundary conditions -- one type per face --
    //     into logical boundary conditions for each velocity component
    //
    // *****************************************************************************
    {
        domain_bcs_type.resize(AMREX_SPACEDIM+NCONS+8);
        domain_bcs_type_d.resize(AMREX_SPACEDIM+NCONS+8);

        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            // only do this for xvel and yvel
            for (int i = 0; i < 2; i++) {
                auto const bct = phys_bc_type[BCVars::xvel_bc+i][ori];
                if ( bct == REMORA_BC::symmetry )
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::reflect_even);
                        if (i==1)
                            domain_bcs_type[BCVars::xvel_bc+dir].setLo(dir, REMORABCType::reflect_odd);
                    } else {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::reflect_even);
                        if (i==1)
                            domain_bcs_type[BCVars::xvel_bc+dir].setHi(dir, REMORABCType::reflect_odd);
                    }
                }
                else if (bct == REMORA_BC::outflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::inflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::ext_dir);
                    } else {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
                else if (bct == REMORA_BC::no_slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::ext_dir);
                    } else {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
                else if (bct == REMORA_BC::slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::foextrap);
                        if (i==1) {
                            // Only normal direction has ext_dir
                            domain_bcs_type[BCVars::xvel_bc+dir].setLo(dir, REMORABCType::ext_dir);
                        }

                    } else {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::foextrap);
                        if (i==1) {
                            // Only normal direction has ext_dir
                            domain_bcs_type[BCVars::xvel_bc+dir].setHi(dir, REMORABCType::ext_dir);
                        }
                    }
                }
                else if (bct == REMORA_BC::periodic)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::int_dir);
                    } else {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::int_dir);
                    }
                }
                else if (bct == REMORA_BC::clamped)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::clamped);
                    } else {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::clamped);
                    }
                }
                else if (bct == REMORA_BC::orlanski_rad)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::orlanski_rad);
                    } else {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::orlanski_rad);
                    }
                }
                else if (bct == REMORA_BC::orlanski_rad_nudge)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::orlanski_rad_nudge);
                    } else {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::orlanski_rad_nudge);
                    }
                }
                else
                {
                    amrex::Abort("Velocity boundary condition not validly specified");
                }
            }

            // Always set zvel_bc to foextrap
            if (side == Orientation::low) {
                domain_bcs_type[BCVars::zvel_bc].setLo(dir, REMORABCType::foextrap);
            } else {
                domain_bcs_type[BCVars::zvel_bc].setHi(dir, REMORABCType::foextrap);
            }
        }
    }

    // *****************************************************************************
    //
    // Here we translate the physical boundary conditions -- one type per face --
    //     into logical boundary conditions for each cell-centered variable
    //
    // *****************************************************************************
    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            for (int i = 0; i < NCONS; i++) {
                auto const bct = phys_bc_type[BCVars::cons_bc+i][ori];
                if ( bct == REMORA_BC::symmetry )
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::reflect_even);
                    } else {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::reflect_even);
                    }
                }
                else if ( bct == REMORA_BC::outflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if ( bct == REMORA_BC::no_slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::inflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::ext_dir);
                    } else {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
                else if (bct == REMORA_BC::periodic)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::int_dir);
                    } else {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::int_dir);
                    }
                }
                else if ( bct == REMORA_BC::clamped)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::clamped);
                    } else {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::clamped);
                    }
                }
                else if ( bct == REMORA_BC::orlanski_rad)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::orlanski_rad);
                    } else {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::orlanski_rad);
                    }
                }
                else if ( bct == REMORA_BC::orlanski_rad_nudge)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::orlanski_rad_nudge);
                    } else {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::orlanski_rad_nudge);
                    }
                }
                else
                {
                    amrex::Abort("Scalar/tracer boundary condition not validly specified");
                }
            }
        }
    }

    // *****************************************************************************
    //
    // Here we translate the physical boundary conditions -- one type per face --
    //     into logical boundary conditions for ubar and vbar. Also add simplified
    //     2d boundary condition (corresponding to BCs in bc_2d.F
    //
    // *****************************************************************************
    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            for (int i = 0; i < 2; i++) {
                auto const bct = phys_bc_type[BCVars::ubar_bc+i][ori];
                if ( bct == REMORA_BC::symmetry )
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::ubar_bc+i].setLo(dir, REMORABCType::reflect_even);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setLo(dir, REMORABCType::reflect_even);
                        if (i==1 and dir!=2) {
                            domain_bcs_type[BCVars::ubar_bc+dir].setLo(dir, REMORABCType::reflect_odd);
                            domain_bcs_type[BCVars::u2d_simple_bc+dir].setLo(dir, REMORABCType::reflect_odd);
                        }
                    } else {
                        domain_bcs_type[BCVars::ubar_bc+i].setHi(dir, REMORABCType::reflect_even);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setHi(dir, REMORABCType::reflect_even);
                        if (i==1 and dir!=2) {
                            domain_bcs_type[BCVars::ubar_bc+dir].setHi(dir, REMORABCType::reflect_odd);
                            domain_bcs_type[BCVars::u2d_simple_bc+dir].setHi(dir, REMORABCType::reflect_odd);
                        }
                    }
                }
                else if (bct == REMORA_BC::outflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::ubar_bc+i].setLo(dir, REMORABCType::foextrap);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[BCVars::ubar_bc+i].setHi(dir, REMORABCType::foextrap);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::inflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::ubar_bc+i].setLo(dir, REMORABCType::ext_dir);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setLo(dir, REMORABCType::ext_dir);
                    } else {
                        domain_bcs_type[BCVars::ubar_bc+i].setHi(dir, REMORABCType::ext_dir);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
                else if (bct == REMORA_BC::no_slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::ubar_bc+i].setLo(dir, REMORABCType::ext_dir);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setLo(dir, REMORABCType::ext_dir);
                    } else {
                        domain_bcs_type[BCVars::ubar_bc+i].setHi(dir, REMORABCType::ext_dir);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
                else if (bct == REMORA_BC::slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::ubar_bc+i].setLo(dir, REMORABCType::foextrap);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setLo(dir, REMORABCType::foextrap);
                        if (i==1 and dir!=2) {
                            // Only normal direction has ext_dir
                            domain_bcs_type[BCVars::ubar_bc+dir].setLo(dir, REMORABCType::ext_dir);
                            domain_bcs_type[BCVars::u2d_simple_bc+dir].setLo(dir, REMORABCType::ext_dir);
                        }

                    } else {
                        domain_bcs_type[BCVars::ubar_bc+i].setHi(dir, REMORABCType::foextrap);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setHi(dir, REMORABCType::foextrap);
                        if (i==1 and dir!=2) {
                            // Only normal direction has ext_dir
                            domain_bcs_type[BCVars::ubar_bc+dir].setHi(dir, REMORABCType::ext_dir);
                            domain_bcs_type[BCVars::u2d_simple_bc+dir].setHi(dir, REMORABCType::ext_dir);
                        }
                    }
                }
                else if (bct == REMORA_BC::periodic)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::ubar_bc+i].setLo(dir, REMORABCType::int_dir);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setLo(dir, REMORABCType::int_dir);
                    } else {
                        domain_bcs_type[BCVars::ubar_bc+i].setHi(dir, REMORABCType::int_dir);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setHi(dir, REMORABCType::int_dir);
                    }
                }
                else if (bct == REMORA_BC::clamped)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::ubar_bc+i].setLo(dir, REMORABCType::clamped);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[BCVars::ubar_bc+i].setHi(dir, REMORABCType::clamped);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::flather)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::ubar_bc+i].setLo(dir, REMORABCType::chapman);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setLo(dir, REMORABCType::foextrap);
                        if (i==1 and dir!=2) {
                            // Only normal direction has Flather
                            domain_bcs_type[BCVars::ubar_bc+dir].setLo(dir, REMORABCType::flather);
                            domain_bcs_type[BCVars::u2d_simple_bc+dir].setLo(dir, REMORABCType::foextrap);
                        }

                    } else {
                        domain_bcs_type[BCVars::ubar_bc+i].setHi(dir, REMORABCType::chapman);
                        domain_bcs_type[BCVars::u2d_simple_bc+i].setHi(dir, REMORABCType::foextrap);
                        if (i==1 and dir!=2) {
                            // Only normal direction has Flather
                            domain_bcs_type[BCVars::ubar_bc+dir].setHi(dir, REMORABCType::flather);
                            domain_bcs_type[BCVars::u2d_simple_bc+dir].setHi(dir, REMORABCType::foextrap);
                        }
                    }
                }
                else
                {
                    amrex::Abort("ubar or vbar boundary condition not validly specified");
                }
            }
        }
    }

    // *****************************************************************************
    //
    // Here we translate the physical boundary conditions -- one type per face --
    //     into logical boundary conditions for zeta and tke
    //
    // *****************************************************************************
    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            for (int i = 0; i < 2; i++) {
                auto const bct = phys_bc_type[BCVars::zeta_bc+i][ori];
                if ( bct == REMORA_BC::symmetry )
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::zeta_bc+i].setLo(dir, REMORABCType::reflect_even);
                    } else {
                        domain_bcs_type[BCVars::zeta_bc+i].setHi(dir, REMORABCType::reflect_even);
                    }
                }
                else if ( bct == REMORA_BC::outflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::zeta_bc+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[BCVars::zeta_bc+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if ( bct == REMORA_BC::no_slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::zeta_bc+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[BCVars::zeta_bc+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::zeta_bc+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[BCVars::zeta_bc+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::inflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::zeta_bc+i].setLo(dir, REMORABCType::ext_dir);
                    } else {
                        domain_bcs_type[BCVars::zeta_bc+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
                else if (bct == REMORA_BC::periodic)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::zeta_bc+i].setLo(dir, REMORABCType::int_dir);
                    } else {
                        domain_bcs_type[BCVars::zeta_bc+i].setHi(dir, REMORABCType::int_dir);
                    }
                }
                else if (bct == REMORA_BC::chapman)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::zeta_bc+i].setLo(dir, REMORABCType::chapman);
                    } else {
                        domain_bcs_type[BCVars::zeta_bc+i].setHi(dir, REMORABCType::chapman);
                    }
                }
                else if ( bct == REMORA_BC::clamped)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[BCVars::zeta_bc+i].setLo(dir, REMORABCType::clamped);
                    } else {
                        domain_bcs_type[BCVars::zeta_bc+i].setHi(dir, REMORABCType::clamped);
                    }
                }
                else
                {
                    amrex::Abort("Free surface (zeta) boundary condition not validly specified");
                }
            }
        }
    }

    // *****************************************************************************
    //
    // Here we define a boundary condition that will foextrap while respecting periodicity
    // This is used as a "null BC"
    //
    // *****************************************************************************
    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = phys_bc_type[ori];
            if (side == Orientation::low) {
                domain_bcs_type[BCVars::foextrap_periodic_bc].setLo(dir, REMORABCType::foextrap);
            } else {
                domain_bcs_type[BCVars::foextrap_periodic_bc].setHi(dir, REMORABCType::foextrap);
            }
        }
    }

    // *****************************************************************************
    //
    // Here we define a boundary condition that will unconditionally foextrap
    //
    // *****************************************************************************
    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = phys_bc_type[ori];
            if (side == Orientation::low) {
                domain_bcs_type[BCVars::foextrap_bc].setLo(dir, REMORABCType::foextrap);
            } else {
                domain_bcs_type[BCVars::foextrap_bc].setHi(dir, REMORABCType::foextrap);
            }
        }
    }


#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy
        (domain_bcs_type_d.data(), domain_bcs_type.data(),
         sizeof(amrex::BCRec)*(NCONS+AMREX_SPACEDIM+8));
#else
    std::memcpy
        (domain_bcs_type_d.data(), domain_bcs_type.data(),
         sizeof(amrex::BCRec)*(NCONS+AMREX_SPACEDIM+8));
#endif
}

