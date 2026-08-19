#include <AMReX_Vector.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ParmParse.H>

#include <REMORA.H>

using namespace amrex;

void REMORA::init_bcs ()
{
    const int xvel_bc_idx = xvel_bc();
    const int yvel_bc_idx = yvel_bc();
    const int zvel_bc_idx = zvel_bc();
    const int ubar_bc_idx = ubar_bc();
    const int vbar_bc_idx = vbar_bc();
    const int zeta_bc_idx = zeta_bc();
    const int tke_bc_idx = tke_bc();
    const int foextrap_periodic_bc_idx = foextrap_periodic_bc();
    const int foextrap_bc_idx = foextrap_bc();
    const int u2d_simple_bc_idx = u2d_simple_bc();
    const int v2d_simple_bc_idx = v2d_simple_bc();

    // Every cell-centered tracer is named for itself, so temp, salt, and each additional
    // passive or biology scalar can be configured independently. cons_names is
    // {"temp", "salt", <scalar names>} and is filled by init_scalar_metadata.
    std::vector<std::string> bcvar_names(BCVars::NumTypes(ncons),"");
    for (int icomp = 0; icomp < ncons; ++icomp) {
        bcvar_names[icomp] = cons_names[icomp];
    }
    bcvar_names[xvel_bc_idx] = "u";
    bcvar_names[yvel_bc_idx] = "v";
    bcvar_names[zvel_bc_idx] = "w";
    bcvar_names[ubar_bc_idx] = "ubar";
    bcvar_names[vbar_bc_idx] = "vbar";
    bcvar_names[zeta_bc_idx] = "zeta";
    bcvar_names[tke_bc_idx] = "tke";

    phys_bc_type.assign(num_bc_vars(), {});
    m_bc_extdir_vals.assign(num_bc_vars(), {});

    // One entry per BdyVars slot plus a trailing scratch slot for the variables that
    // have no NetCDF boundary lane at all (see the bdy_index assignment below)
    phys_bc_need_data.assign(num_bdy_vars()+1, {});

    auto uses_velocity_input = [=] (int bcvar_type) noexcept {
        return bcvar_type == xvel_bc_idx || bcvar_type == yvel_bc_idx || bcvar_type == zvel_bc_idx;
    };

    auto uses_scalar_input = [this] (int bcvar_type) noexcept {
        return bcvar_type >= Tracer_comp && bcvar_type < ncons;
    };

    auto f_set_var_bc = [this, uses_velocity_input, uses_scalar_input, xvel_bc_idx, zeta_bc_idx, ubar_bc_idx, vbar_bc_idx, bcvar_names]
        (ParmParse& pp, int bcvar_type, Orientation ori, const std::string& bc_type_string) {
        // A side keyword means the same thing for every variable it covers: no tracer,
        // dye or biology, is quietly given a different condition than the one asked
        // for. A file-driven side therefore expects boundary data for each tracer the
        // run carries, named as ROMS names it (tracer_west, NO3_west, ...), and
        // NCTimeSeriesBoundary aborts by name if the file has none. Runs whose
        // boundary file covers temp and salt alone either drop the extra tracers
        // (remora.nscalar = 0) or move to remora.boundary_per_variable = true, where
        // each tracer's condition is stated outright.

        if (bc_type_string == "symmetry")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::symmetry;
            phys_bc_need_data[bdy_index[bcvar_type]][ori] = false;
            domain_bc_type[ori] = "Symmetry";
        }
        else if (bc_type_string == "outflow")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::outflow;
            phys_bc_need_data[bdy_index[bcvar_type]][ori] = false;
            domain_bc_type[ori] = "Outflow";
        }
        else if (bc_type_string == "inflow")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::inflow;
            phys_bc_need_data[bdy_index[bcvar_type]][ori] = false;
            domain_bc_type[ori] = "Inflow";

            if (uses_velocity_input(bcvar_type)) {
                std::vector<Real> v;
                pp.getarr("velocity", v, 0, AMREX_SPACEDIM);
                m_bc_extdir_vals[bcvar_type][ori] = v[bcvar_type - xvel_bc_idx];
            } else if (uses_scalar_input(bcvar_type)) {
                Real scalar_in = zero;
                if (pp.queryAdd("scalar", scalar_in)) {
                    m_bc_extdir_vals[bcvar_type][ori] = scalar_in;
                }
            }
        }
        else if (bc_type_string == "noslipwall")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::no_slip_wall;
            phys_bc_need_data[bdy_index[bcvar_type]][ori] = false;
            domain_bc_type[ori] = "NoSlipWall";

            if (uses_velocity_input(bcvar_type)) {
                std::vector<Real> v;

                // The values of m_bc_extdir_vals default to 0.
                // But if we find "velocity" in the inputs file, use those values instead.
                if (pp.queryarr("velocity", v, 0, AMREX_SPACEDIM))
                {
                    v[ori.coordDir()] = zero;
                    m_bc_extdir_vals[bcvar_type][ori] = v[bcvar_type - xvel_bc_idx];
                }
            }
        }
        else if (bc_type_string == "slipwall")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::slip_wall;
            phys_bc_need_data[bdy_index[bcvar_type]][ori] = false;
            domain_bc_type[ori] = "SlipWall";
        }
        else if (bc_type_string == "clamped")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::clamped;
            phys_bc_need_data[bdy_index[bcvar_type]][ori] = true;
            domain_bc_type[ori] = "Clamped";
            solverChoice.boundary_from_netcdf = true;
        }
        else if (bc_type_string == "chapman")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::chapman;
            phys_bc_need_data[bdy_index[bcvar_type]][ori] = true;
            domain_bc_type[ori] = "Chapman";
            solverChoice.boundary_from_netcdf = true;

            if (bcvar_type != zeta_bc_idx) {
                amrex::Abort("Chapman BC can only be applied to zeta");
            }
        }
        else if (bc_type_string == "flather")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::flather;
            phys_bc_need_data[bdy_index[bcvar_type]][ori] = true;
            domain_bc_type[ori] = "Flather";
            solverChoice.boundary_from_netcdf = true;

            if (!(bcvar_type == ubar_bc_idx || bcvar_type == vbar_bc_idx)) {
                amrex::Abort("Flather BC can only be applied to ubar or vbar");
            }
        }
        else if (bc_type_string == "orlanski_rad")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::orlanski_rad;
            phys_bc_need_data[bdy_index[bcvar_type]][ori] = false;
            domain_bc_type[ori] = "Orlanski Radiation";
        }
        else if (bc_type_string == "orlanski_rad_nudg")
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::orlanski_rad_nudge;
            phys_bc_need_data[bdy_index[bcvar_type]][ori] = true;
            domain_bc_type[ori] = "Orlanski Radiation with nudging";
            solverChoice.boundary_from_netcdf = true;
        }
        else if (bc_type_string == "periodic")
        {
            if (!geom[0].isPeriodic(ori.coordDir())) {
                amrex::Abort("Periodic boundary specified in a non-periodic direction");
            }
            phys_bc_type[bcvar_type][ori] = REMORA_BC::periodic;
            phys_bc_need_data[bdy_index[bcvar_type]][ori] = false;
            domain_bc_type[ori] = "Periodic";
        }
        else
        {
            phys_bc_type[bcvar_type][ori] = REMORA_BC::undefined;
            phys_bc_need_data[bdy_index[bcvar_type]][ori] = false;
        }

        if (geom[0].isPeriodic(ori.coordDir()))
        {
            domain_bc_type[ori] = "Periodic";
            if (phys_bc_type[bcvar_type][ori] == REMORA_BC::undefined)
            {
                phys_bc_type[bcvar_type][ori] = REMORA_BC::periodic;
                phys_bc_need_data[bdy_index[bcvar_type]][ori] = false;
            } else if (phys_bc_type[bcvar_type][ori] != REMORA_BC::periodic) {
                amrex::Abort("Wrong BC type for periodic boundary");
            }
        }

        if (phys_bc_type[bcvar_type][ori] == REMORA_BC::undefined)
        {
             amrex::Print() << "BC Type specified for variable " << bcvar_names[bcvar_type] << " is " << bc_type_string << std::endl;
             amrex::Abort("This BC type is unknown");
        }
    };

    auto f_by_side = [this, &f_set_var_bc] (std::string const& bcid, Orientation ori)
    {
        ParmParse pp("remora.bc."+bcid);
        std::string bc_type_in = "null";
        // Default z directions to slipwall
        if (bcid=="zlo" or bcid=="zhi") bc_type_in = "slipwall";
        pp.queryAdd("type", bc_type_in);
        std::string bc_type = amrex::toLower(bc_type_in);

        for (int icomp = 0; icomp < num_bc_vars(); ++icomp) {
            f_set_var_bc(pp, icomp, ori, bc_type);
        }
    };

    auto f_by_var = [this, &f_set_var_bc, zvel_bc_idx] (std::string const& varname, int bcvar_type,
                                                        std::string const& fallback_varname = "")
    {
        amrex::Vector<Orientation> orientations = {Orientation(Direction::x,Orientation::low), Orientation(Direction::y,Orientation::high),Orientation(Direction::x,Orientation::high),Orientation(Direction::y,Orientation::low)}; // west, south, east, north [matches ROMS]
        std::vector<std::string> bc_types = {"null","null","null","null"};

        // Additional scalars are addressable by their own name, but fall back to the
        // shared "scalar" keyword when no per-tracer entry is given so that inputs
        // written before tracers were individually addressable keep working
        std::string prefix = "remora.bc." + varname;
        if (!fallback_varname.empty() && !ParmParse(prefix).contains("type")) {
            prefix = "remora.bc." + fallback_varname;
        }
        ParmParse pp(prefix);
        std::string bc_type_in = "null";
        // default zvel to outflow
        if (bcvar_type == zvel_bc_idx) {
            bc_types = {"outflow","outflow","outflow","outflow"};
            for (int i=0; i<4; i++) {
                auto ori = orientations[i];
                if (geom[0].isPeriodic(ori.coordDir())) {
                    bc_types[i] = "periodic";
                }
            }
        }
        pp.queryarr("type", bc_types);
        AMREX_ASSERT(bc_types.size() == 4);
        for (int i=0; i<4; i++) {
            std::string bc_type = amrex::toLower(bc_types[i]);
            auto ori = orientations[i];
            f_set_var_bc(pp, bcvar_type, ori, bc_type);
        }
    };

    // Variables with no entry here point at the trailing scratch slot of
    // phys_bc_need_data, which means "this variable never reads boundary data from file"
    bdy_index.assign(num_bc_vars(), num_bdy_vars());
    for (int icomp = 0; icomp < ncons; ++icomp) {
        bdy_index[BCVars::cons_bc+icomp] = BdyVars::cons(icomp);
    }
    bdy_index[xvel_bc_idx] = BdyVars::u;
    bdy_index[yvel_bc_idx] = BdyVars::v;
    bdy_index[ubar_bc_idx] = bdy_ubar();
    bdy_index[vbar_bc_idx] = bdy_vbar();
    bdy_index[zeta_bc_idx] = bdy_zeta();

    for (OrientationIter oit; oit; ++oit) {
        Orientation ori = oit();
        // These are simply defaults for Dirichlet faces -- they should be over-written below if needed
        m_bc_extdir_vals[BCVars::Temp_bc_comp  ][ori] = Real(1.e19);
        m_bc_extdir_vals[BCVars::Salt_bc_comp  ][ori] = Real(1.e20);
        for (int icomp = Tracer_comp; icomp < ncons; ++icomp) {
            m_bc_extdir_vals[icomp][ori] = Real(1.e21) + static_cast<Real>(icomp - Tracer_comp);
        }

        m_bc_extdir_vals[xvel_bc_idx][ori] = zero; // default
        m_bc_extdir_vals[yvel_bc_idx][ori] = zero;
        m_bc_extdir_vals[zvel_bc_idx][ori] = zero;

        m_bc_extdir_vals[ubar_bc_idx][ori] = zero; // default
        m_bc_extdir_vals[vbar_bc_idx][ori] = zero;
        m_bc_extdir_vals[u2d_simple_bc_idx][ori] = zero;
        m_bc_extdir_vals[v2d_simple_bc_idx][ori] = zero;
    }

    // Whether to specify boundary conditions by variable (then side).
    // Alternative is to do it by side by indicating keywords that indicate multiple variables
    set_bcs_by_var = false;

    ParmParse pp("remora");
    pp.queryAdd("boundary_per_variable", set_bcs_by_var);
    // Any tracer named for itself counts as a per-variable specification too
    bool any_cons_bc_specified = false;
    for (int icomp = 0; icomp < ncons; ++icomp) {
        any_cons_bc_specified = any_cons_bc_specified || pp.contains(("bc."+cons_names[icomp]+".type").c_str());
    }
    // Check whether variable specification matches flag in inputs file
    if (!set_bcs_by_var && (any_cons_bc_specified ||
                pp.contains("bc.scalar.type") ||
                pp.contains("bc.u.type") ||
                pp.contains("bc.v.type") ||
                pp.contains("bc.w.type") ||
                pp.contains("bc.ubar.type") ||
                pp.contains("bc.vbar.type") ||
                pp.contains("bc.zeta.type") ||
                pp.contains("bc.tke.type"))) {
        amrex::Abort("boundary_per_variable set to false, but per-variable boundary conditions are specified. Use bc.{x,y}{lo,hi}.type instead");
    }
    if (set_bcs_by_var && (pp.contains("bc.xlo.type") ||
                pp.contains("bc.xhi.type") ||
                pp.contains("bc.ylo.type") ||
                pp.contains("bc.yhi.type"))) {
        amrex::Abort("boundary_per_variable set to true, but per-side boundary conditions are specified. Use bc.{temp,salt,etc}.type instead");
    }
    if (!set_bcs_by_var) {
        f_by_side("xlo", Orientation(Direction::x,Orientation::low));
        f_by_side("xhi", Orientation(Direction::x,Orientation::high));
        f_by_side("ylo", Orientation(Direction::y,Orientation::low));
        f_by_side("yhi", Orientation(Direction::y,Orientation::high));
    } else {
        f_by_var("temp", BCVars::Temp_bc_comp);
        f_by_var("salt", BCVars::Salt_bc_comp);
        for (int icomp = Tracer_comp; icomp < ncons; ++icomp) {
            f_by_var(cons_names[icomp], icomp, "scalar");
        }
        f_by_var("u", xvel_bc_idx);
        f_by_var("v", yvel_bc_idx);
        f_by_var("w", zvel_bc_idx);
        f_by_var("ubar", ubar_bc_idx);
        f_by_var("vbar", vbar_bc_idx);
        f_by_var("zeta", zeta_bc_idx);
        f_by_var("tke", tke_bc_idx);
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
        domain_bcs_type.resize(num_bc_vars());
        domain_bcs_type_d.resize(num_bc_vars());

        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            // only do this for xvel and yvel
            for (int i = 0; i < 2; i++) {
                auto const bct = phys_bc_type[xvel_bc_idx+i][ori];
                if ( bct == REMORA_BC::symmetry )
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[xvel_bc_idx+i].setLo(dir, REMORABCType::reflect_even);
                        if (i==1)
                            domain_bcs_type[xvel_bc_idx+dir].setLo(dir, REMORABCType::reflect_odd);
                    } else {
                        domain_bcs_type[xvel_bc_idx+i].setHi(dir, REMORABCType::reflect_even);
                        if (i==1)
                            domain_bcs_type[xvel_bc_idx+dir].setHi(dir, REMORABCType::reflect_odd);
                    }
                }
                else if (bct == REMORA_BC::outflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[xvel_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[xvel_bc_idx+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::inflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[xvel_bc_idx+i].setLo(dir, REMORABCType::ext_dir);
                    } else {
                        domain_bcs_type[xvel_bc_idx+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
                else if (bct == REMORA_BC::no_slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[xvel_bc_idx+i].setLo(dir, REMORABCType::ext_dir);
                    } else {
                        domain_bcs_type[xvel_bc_idx+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
                else if (bct == REMORA_BC::slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[xvel_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                        if (i==1) {
                            // Only normal direction has ext_dir
                            domain_bcs_type[xvel_bc_idx+dir].setLo(dir, REMORABCType::ext_dir);
                        }

                    } else {
                        domain_bcs_type[xvel_bc_idx+i].setHi(dir, REMORABCType::foextrap);
                        if (i==1) {
                            // Only normal direction has ext_dir
                            domain_bcs_type[xvel_bc_idx+dir].setHi(dir, REMORABCType::ext_dir);
                        }
                    }
                }
                else if (bct == REMORA_BC::periodic)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[xvel_bc_idx+i].setLo(dir, REMORABCType::int_dir);
                    } else {
                        domain_bcs_type[xvel_bc_idx+i].setHi(dir, REMORABCType::int_dir);
                    }
                }
                else if (bct == REMORA_BC::clamped)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[xvel_bc_idx+i].setLo(dir, REMORABCType::clamped);
                    } else {
                        domain_bcs_type[xvel_bc_idx+i].setHi(dir, REMORABCType::clamped);
                    }
                }
                else if (bct == REMORA_BC::orlanski_rad)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[xvel_bc_idx+i].setLo(dir, REMORABCType::orlanski_rad);
                    } else {
                        domain_bcs_type[xvel_bc_idx+i].setHi(dir, REMORABCType::orlanski_rad);
                    }
                }
                else if (bct == REMORA_BC::orlanski_rad_nudge)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[xvel_bc_idx+i].setLo(dir, REMORABCType::orlanski_rad_nudge);
                    } else {
                        domain_bcs_type[xvel_bc_idx+i].setHi(dir, REMORABCType::orlanski_rad_nudge);
                    }
                }
                else
                {
                    amrex::Abort("Velocity boundary condition not validly specified");
                }
            }

            // Always set zvel_bc to foextrap
            if (side == Orientation::low) {
                domain_bcs_type[zvel_bc_idx].setLo(dir, REMORABCType::foextrap);
            } else {
                domain_bcs_type[zvel_bc_idx].setHi(dir, REMORABCType::foextrap);
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
            for (int i = 0; i < ncons; i++) {
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
                auto const bct = phys_bc_type[ubar_bc_idx+i][ori];
                if ( bct == REMORA_BC::symmetry )
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[ubar_bc_idx+i].setLo(dir, REMORABCType::reflect_even);
                        domain_bcs_type[u2d_simple_bc_idx+i].setLo(dir, REMORABCType::reflect_even);
                        if (i==1 and dir!=2) {
                            domain_bcs_type[ubar_bc_idx+dir].setLo(dir, REMORABCType::reflect_odd);
                            domain_bcs_type[u2d_simple_bc_idx+dir].setLo(dir, REMORABCType::reflect_odd);
                        }
                    } else {
                        domain_bcs_type[ubar_bc_idx+i].setHi(dir, REMORABCType::reflect_even);
                        domain_bcs_type[u2d_simple_bc_idx+i].setHi(dir, REMORABCType::reflect_even);
                        if (i==1 and dir!=2) {
                            domain_bcs_type[ubar_bc_idx+dir].setHi(dir, REMORABCType::reflect_odd);
                            domain_bcs_type[u2d_simple_bc_idx+dir].setHi(dir, REMORABCType::reflect_odd);
                        }
                    }
                }
                else if (bct == REMORA_BC::outflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[ubar_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                        domain_bcs_type[u2d_simple_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[ubar_bc_idx+i].setHi(dir, REMORABCType::foextrap);
                        domain_bcs_type[u2d_simple_bc_idx+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::inflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[ubar_bc_idx+i].setLo(dir, REMORABCType::ext_dir);
                        domain_bcs_type[u2d_simple_bc_idx+i].setLo(dir, REMORABCType::ext_dir);
                    } else {
                        domain_bcs_type[ubar_bc_idx+i].setHi(dir, REMORABCType::ext_dir);
                        domain_bcs_type[u2d_simple_bc_idx+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
                else if (bct == REMORA_BC::no_slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[ubar_bc_idx+i].setLo(dir, REMORABCType::ext_dir);
                        domain_bcs_type[u2d_simple_bc_idx+i].setLo(dir, REMORABCType::ext_dir);
                    } else {
                        domain_bcs_type[ubar_bc_idx+i].setHi(dir, REMORABCType::ext_dir);
                        domain_bcs_type[u2d_simple_bc_idx+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
                else if (bct == REMORA_BC::slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[ubar_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                        domain_bcs_type[u2d_simple_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                        if (i==1 and dir!=2) {
                            // Only normal direction has ext_dir
                            domain_bcs_type[ubar_bc_idx+dir].setLo(dir, REMORABCType::ext_dir);
                            domain_bcs_type[u2d_simple_bc_idx+dir].setLo(dir, REMORABCType::ext_dir);
                        }

                    } else {
                        domain_bcs_type[ubar_bc_idx+i].setHi(dir, REMORABCType::foextrap);
                        domain_bcs_type[u2d_simple_bc_idx+i].setHi(dir, REMORABCType::foextrap);
                        if (i==1 and dir!=2) {
                            // Only normal direction has ext_dir
                            domain_bcs_type[ubar_bc_idx+dir].setHi(dir, REMORABCType::ext_dir);
                            domain_bcs_type[u2d_simple_bc_idx+dir].setHi(dir, REMORABCType::ext_dir);
                        }
                    }
                }
                else if (bct == REMORA_BC::periodic)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[ubar_bc_idx+i].setLo(dir, REMORABCType::int_dir);
                        domain_bcs_type[u2d_simple_bc_idx+i].setLo(dir, REMORABCType::int_dir);
                    } else {
                        domain_bcs_type[ubar_bc_idx+i].setHi(dir, REMORABCType::int_dir);
                        domain_bcs_type[u2d_simple_bc_idx+i].setHi(dir, REMORABCType::int_dir);
                    }
                }
                else if (bct == REMORA_BC::clamped)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[ubar_bc_idx+i].setLo(dir, REMORABCType::clamped);
                        domain_bcs_type[u2d_simple_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[ubar_bc_idx+i].setHi(dir, REMORABCType::clamped);
                        domain_bcs_type[u2d_simple_bc_idx+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::flather)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[ubar_bc_idx+i].setLo(dir, REMORABCType::chapman);
                        domain_bcs_type[u2d_simple_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                        if (i==1 and dir!=2) {
                            // Only normal direction has Flather
                            domain_bcs_type[ubar_bc_idx+dir].setLo(dir, REMORABCType::flather);
                            domain_bcs_type[u2d_simple_bc_idx+dir].setLo(dir, REMORABCType::foextrap);
                        }

                    } else {
                        domain_bcs_type[ubar_bc_idx+i].setHi(dir, REMORABCType::chapman);
                        domain_bcs_type[u2d_simple_bc_idx+i].setHi(dir, REMORABCType::foextrap);
                        if (i==1 and dir!=2) {
                            // Only normal direction has Flather
                            domain_bcs_type[ubar_bc_idx+dir].setHi(dir, REMORABCType::flather);
                            domain_bcs_type[u2d_simple_bc_idx+dir].setHi(dir, REMORABCType::foextrap);
                        }
                    }
                }
                else if (bct == REMORA_BC::orlanski_rad)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[ubar_bc_idx+i].setLo(dir, REMORABCType::orlanski_rad);
                        domain_bcs_type[u2d_simple_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[ubar_bc_idx+i].setHi(dir, REMORABCType::orlanski_rad);
                        domain_bcs_type[u2d_simple_bc_idx+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::orlanski_rad_nudge)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[ubar_bc_idx+i].setLo(dir, REMORABCType::orlanski_rad_nudge);
                        domain_bcs_type[u2d_simple_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[ubar_bc_idx+i].setHi(dir, REMORABCType::orlanski_rad_nudge);
                        domain_bcs_type[u2d_simple_bc_idx+i].setHi(dir, REMORABCType::foextrap);
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
                auto const bct = phys_bc_type[zeta_bc_idx+i][ori];
                if ( bct == REMORA_BC::symmetry )
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[zeta_bc_idx+i].setLo(dir, REMORABCType::reflect_even);
                    } else {
                        domain_bcs_type[zeta_bc_idx+i].setHi(dir, REMORABCType::reflect_even);
                    }
                }
                else if ( bct == REMORA_BC::outflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[zeta_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[zeta_bc_idx+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if ( bct == REMORA_BC::no_slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[zeta_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[zeta_bc_idx+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::slip_wall)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[zeta_bc_idx+i].setLo(dir, REMORABCType::foextrap);
                    } else {
                        domain_bcs_type[zeta_bc_idx+i].setHi(dir, REMORABCType::foextrap);
                    }
                }
                else if (bct == REMORA_BC::inflow)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[zeta_bc_idx+i].setLo(dir, REMORABCType::ext_dir);
                    } else {
                        domain_bcs_type[zeta_bc_idx+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
                else if (bct == REMORA_BC::periodic)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[zeta_bc_idx+i].setLo(dir, REMORABCType::int_dir);
                    } else {
                        domain_bcs_type[zeta_bc_idx+i].setHi(dir, REMORABCType::int_dir);
                    }
                }
                else if (bct == REMORA_BC::chapman)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[zeta_bc_idx+i].setLo(dir, REMORABCType::chapman);
                    } else {
                        domain_bcs_type[zeta_bc_idx+i].setHi(dir, REMORABCType::chapman);
                    }
                }
                else if ( bct == REMORA_BC::clamped)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[zeta_bc_idx+i].setLo(dir, REMORABCType::clamped);
                    } else {
                        domain_bcs_type[zeta_bc_idx+i].setHi(dir, REMORABCType::clamped);
                    }
                }
                else if ( bct == REMORA_BC::orlanski_rad)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[zeta_bc_idx+i].setLo(dir, REMORABCType::orlanski_rad);
                    } else {
                        domain_bcs_type[zeta_bc_idx+i].setHi(dir, REMORABCType::orlanski_rad);
                    }
                }
                else if ( bct == REMORA_BC::orlanski_rad_nudge)
                {
                    if (side == Orientation::low) {
                        domain_bcs_type[zeta_bc_idx+i].setLo(dir, REMORABCType::orlanski_rad_nudge);
                    } else {
                        domain_bcs_type[zeta_bc_idx+i].setHi(dir, REMORABCType::orlanski_rad_nudge);
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
            if (side == Orientation::low) {
                domain_bcs_type[foextrap_periodic_bc_idx].setLo(dir, REMORABCType::foextrap);
            } else {
                domain_bcs_type[foextrap_periodic_bc_idx].setHi(dir, REMORABCType::foextrap);
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
            if (side == Orientation::low) {
                domain_bcs_type[foextrap_bc_idx].setLo(dir, REMORABCType::foextrap);
            } else {
                domain_bcs_type[foextrap_bc_idx].setHi(dir, REMORABCType::foextrap);
            }
        }
    }


#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy
        (domain_bcs_type_d.data(), domain_bcs_type.data(),
         sizeof(amrex::BCRec)*num_bc_vars());
#else
    std::memcpy
        (domain_bcs_type_d.data(), domain_bcs_type.data(),
         sizeof(amrex::BCRec)*num_bc_vars());
#endif
}
