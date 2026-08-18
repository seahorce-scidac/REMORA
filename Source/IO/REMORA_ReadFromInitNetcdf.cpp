#include "REMORA_NCFile.H"
#include "AMReX_FArrayBox.H"
#include "REMORA_DataStruct.H"

using namespace amrex;

#ifdef REMORA_USE_NETCDF
/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_temp_fab     container for temperature data
 * @param NC_salt_fab     container for salinity data
 * @param NC_u_fab        container for u velocity data
 * @param NC_v_fab        container for v velocity data
 */
void
read_data_from_netcdf (int /*lev*/,
                       const Box& domain,
                       const std::string& fname,
                       FArrayBox& NC_temp_fab, FArrayBox& NC_salt_fab,
                       FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab)
{
    amrex::Print() << "Loading initial solution data from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_temp_fab); NC_names.push_back("temp");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 0
    NC_fabs.push_back(&NC_salt_fab); NC_names.push_back("salt");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 1
    NC_fabs.push_back(&NC_xvel_fab); NC_names.push_back("u");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 2
    NC_fabs.push_back(&NC_yvel_fab); NC_names.push_back("v");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 3

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_temp_fab     container for temperature data
 * @param NC_salt_fab     container for salinity data
 * @param NC_u_fab        container for u velocity data
 * @param NC_v_fab        container for v velocity data
 * @param ngrow           number of grow cells to read
 */
void
read_data_full_domain_from_netcdf (int /*lev*/,
                       const Box& domain,
                       const std::string& fname,
                       FArrayBox& NC_temp_fab, FArrayBox& NC_salt_fab,
                       FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                       IntVect ngrow)
{
    amrex::Print() << "Loading initial solution data from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_temp_fab); NC_names.push_back("temp");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 0
    NC_fabs.push_back(&NC_salt_fab); NC_names.push_back("salt");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 1
    NC_fabs.push_back(&NC_xvel_fab); NC_names.push_back("u");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 2
    NC_fabs.push_back(&NC_yvel_fab); NC_names.push_back("v");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 3

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs, false, 0, ngrow);
}

void
read_biology_from_netcdf (int /*lev*/,
                          const Box& domain,
                          const std::string& fname,
                          const Vector<std::string>& biology_names,
                          Vector<FArrayBox>& NC_biology_fab)
{
    if (biology_names.empty()) {
        return;
    }

    amrex::Print() << "Loading initial biology data from NetCDF file " << fname << std::endl;

    NC_biology_fab.resize(biology_names.size());
    Vector<FArrayBox*> NC_fabs;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    for (int ibio = 0; ibio < static_cast<int>(biology_names.size()); ++ibio) {
        NC_fabs.push_back(&NC_biology_fab[ibio]);
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    }

    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, biology_names, NC_dim_types, NC_fabs);
}

void
read_biology_full_domain_from_netcdf (int /*lev*/,
                                      const Box& domain,
                                      const std::string& fname,
                                      const Vector<std::string>& biology_names,
                                      Vector<FArrayBox>& NC_biology_fab,
                                      IntVect ngrow)
{
    if (biology_names.empty()) {
        return;
    }

    amrex::Print() << "Loading high resolution biology data from NetCDF file " << fname << std::endl;

    NC_biology_fab.resize(biology_names.size());
    Vector<FArrayBox*> NC_fabs;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    for (int ibio = 0; ibio < static_cast<int>(biology_names.size()); ++ibio) {
        NC_fabs.push_back(&NC_biology_fab[ibio]);
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    }

    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, biology_names, NC_dim_types, NC_fabs, false, 0, ngrow);
}

/**
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param scalar_names    per passive scalar, name of the variable in the file
 * @param NC_scalar_fab   per passive scalar, container for the dye data
 * @param scalar_in_file  filled on output: per passive scalar, whether the variable was
 *                        found in the file and its FAB was built
 * @param full_domain     whether this is the high resolution full-domain read, which
 *                        needs ngrow grow cells rather than the reader's default
 * @param ngrow           number of grow cells to read, for the full-domain read
 *
 * A passive scalar's initial field is optional. Initial files written for runs that
 * carry no dye -- which is every ROMS initial file predating the dye variables, and
 * most idealized ones -- have no "tracer" variable, and remora.nscalar defaults to 1,
 * so requiring the field would break every existing NetCDF-initialized run. Each name
 * the file does carry is read; the rest keep the zero that init_data_from_netcdf set,
 * which is the behavior those runs had before this read existed.
 *
 * Shared by the per-box and full-domain entry points below: the presence testing is
 * the same for both, and only the grow cells and the wording of the log differ.
 */
namespace {
void
read_scalars_impl (const Box& domain,
                   const std::string& fname,
                   const Vector<std::string>& scalar_names,
                   Vector<FArrayBox>& NC_scalar_fab,
                   Vector<int>& scalar_in_file,
                   bool full_domain,
                   IntVect ngrow)
{
    scalar_in_file.assign(scalar_names.size(), 0);
    if (scalar_names.empty()) {
        return;
    }

    NC_scalar_fab.resize(scalar_names.size());
    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;
    Vector<std::string> missing;

    for (int iscal = 0; iscal < static_cast<int>(scalar_names.size()); ++iscal) {
        if (!QueryNetCDFHasVars(fname, {scalar_names[iscal]})) {
            missing.push_back(scalar_names[iscal]);
            continue;
        }
        scalar_in_file[iscal] = 1;
        NC_fabs.push_back(&NC_scalar_fab[iscal]);
        NC_names.push_back(scalar_names[iscal]);
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    }

    if (!missing.empty()) {
        amrex::Print() << "Initial file " << fname << " does not contain";
        for (const auto& name : missing) { amrex::Print() << " " << name; }
        amrex::Print() << "; those passive scalars start at zero" << std::endl;
    }
    if (NC_names.empty()) {
        return;
    }

    amrex::Print() << "Loading " << (full_domain ? "high resolution " : "initial ")
                   << "passive scalar data from NetCDF file " << fname << std::endl;

    // Read the netcdf file and fill these FABs. The per-box read takes the reader's
    // default grow cells, as the temperature and salinity read does.
    if (full_domain) {
        BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs, false, 0, ngrow);
    } else {
        BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
    }
}
} // namespace

void
read_scalars_from_netcdf (int /*lev*/,
                          const Box& domain,
                          const std::string& fname,
                          const Vector<std::string>& scalar_names,
                          Vector<FArrayBox>& NC_scalar_fab,
                          Vector<int>& scalar_in_file)
{
    read_scalars_impl(domain, fname, scalar_names, NC_scalar_fab, scalar_in_file,
                      false, IntVect(0,0,0));
}

void
read_scalars_full_domain_from_netcdf (int /*lev*/,
                                      const Box& domain,
                                      const std::string& fname,
                                      const Vector<std::string>& scalar_names,
                                      Vector<FArrayBox>& NC_scalar_fab,
                                      Vector<int>& scalar_in_file,
                                      IntVect ngrow)
{
    read_scalars_impl(domain, fname, scalar_names, NC_scalar_fab, scalar_in_file,
                      true, ngrow);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_zeta_fab     container for sea surface height data
 */
void
read_zeta_from_netcdf (int /*lev*/,
                      const Box& domain,
                      const std::string& fname,
                      FArrayBox& NC_zeta_fab)
{
    amrex::Print() << "Loading initial sea surface height from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_zeta_fab )   ; NC_names.push_back("zeta")    ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); // 0

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_zeta_fab     container for sea surface height data
 * @param ngrow           number of grow cells to read in, if not default
 */
void
read_zeta_full_domain_from_netcdf (int /*lev*/,
                      const Box& domain,
                      const std::string& fname,
                      FArrayBox& NC_zeta_fab,
                      IntVect ngrow)
{
    amrex::Print() << "Loading initial sea surface height from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_zeta_fab )   ; NC_names.push_back("zeta")    ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); // 0

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs, false, 0, ngrow);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_h_fab        container for bathymetry data
 */
void
read_bathymetry_from_netcdf (int /*lev*/,
                             const Box& domain,
                             const std::string& fname,
                             FArrayBox& NC_h_fab)
{
    amrex::Print() << "Loading initial bathymetry from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_h_fab )   ; NC_names.push_back("h")     ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 0

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_pm_fab       container for pm data
 * @param NC_pn_fab       container for pn data
 * @param NC_xr_fab       container for x_rho data
 * @param NC_yr_fab       container for y_rho data
 * @param NC_xu_fab       container for x_u data
 * @param NC_yu_fab       container for y_u data
 * @param NC_xv_fab       container for x_v data
 * @param NC_yv_fab       container for y_v data
 * @param NC_xp_fab       container for x_p data
 * @param NC_yp_fab       container for y_p data
 */
void
read_grid_vars_from_netcdf (int /*lev*/,
                             const Box& domain,
                             const std::string& fname,
                             FArrayBox& NC_pm_fab, FArrayBox& NC_pn_fab,
                             FArrayBox& NC_xr_fab, FArrayBox& NC_yr_fab,
                             FArrayBox& NC_xu_fab, FArrayBox& NC_yu_fab,
                             FArrayBox& NC_xv_fab, FArrayBox& NC_yv_fab,
                             FArrayBox& NC_xp_fab, FArrayBox& NC_yp_fab)
{
    amrex::Print() << "Loading grid variables from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_pm_fab)   ; NC_names.push_back("pm")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 1
    NC_fabs.push_back(&NC_pn_fab)   ; NC_names.push_back("pn")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 2
    NC_fabs.push_back(&NC_xr_fab)   ; NC_names.push_back("x_rho") ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 3
    NC_fabs.push_back(&NC_yr_fab)   ; NC_names.push_back("y_rho") ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 4
    NC_fabs.push_back(&NC_xu_fab)   ; NC_names.push_back("x_u")   ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 5
    NC_fabs.push_back(&NC_yu_fab)   ; NC_names.push_back("y_u")   ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 6
    NC_fabs.push_back(&NC_xv_fab)   ; NC_names.push_back("x_v")   ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 7
    NC_fabs.push_back(&NC_yv_fab)   ; NC_names.push_back("y_v")   ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 8
    NC_fabs.push_back(&NC_xp_fab)   ; NC_names.push_back("x_psi") ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 9
    NC_fabs.push_back(&NC_yp_fab)   ; NC_names.push_back("y_psi") ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 10

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_lonp_fab     container for lon_psi data
 * @param NC_latp_fab     container for lat_psi data
 * @returns whether the file carried spherical psi coordinates and they were read
 *
 * ROMS grid files for spherical grids carry lon_psi/lat_psi alongside
 * x_psi/y_psi, and those degree-valued corner arrays are what a coupled driver
 * needs when the two components do not share a projected Cartesian frame
 * (SCRIP/COAWST selects the same way; see Lib/SCRIP_COAWST/read_roms.f).
 * Idealized grid files carry no spherical coordinates, so this read is optional
 * and reports whether it happened rather than aborting.
 */
bool
read_spherical_grid_vars_from_netcdf (int /*lev*/,
                                      const Box& domain,
                                      const std::string& fname,
                                      FArrayBox& NC_lonp_fab, FArrayBox& NC_latp_fab)
{
    Vector<std::string> NC_names;
    NC_names.push_back("lon_psi");
    NC_names.push_back("lat_psi");

    if (!QueryNetCDFHasVars(fname, NC_names)) {
        amrex::Print() << "Grid file " << fname << " has no lon_psi/lat_psi; "
                       << "spherical psi coordinates unavailable" << std::endl;
        return false;
    }

    amrex::Print() << "Loading spherical psi coordinates from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_lonp_fab) ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE);
    NC_fabs.push_back(&NC_latp_fab) ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE);

    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
    return true;
}

/**
 * @param domain          simulation domain at nc_hires_grid_level
 * @param fname           file name to read from
 * @param NC_h_fab        container for bathymetry data
 * @param ngrow           number of grow cells to read in, if not default
 */
void
read_bathymetry_full_domain_from_netcdf (const Box& domain,
                                         const std::string& fname,
                                         FArrayBox& NC_h_fab, IntVect ngrow)
{
    amrex::Print() << "Loading high resolution bathymetry from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_h_fab )   ; NC_names.push_back("h")     ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 0

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs, false, 0, ngrow);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_pm_fab       container for pm data
 * @param NC_pn_fab       container for pn data
 * @param ngrow           number of grow cells to read in, if not default
 */
void
read_grid_vars_full_domain_from_netcdf (const Box& domain,
                             const std::string& fname,
                             FArrayBox& NC_pm_fab, FArrayBox& NC_pn_fab,
                             IntVect ngrow)
{
    amrex::Print() << "Loading high resolution grid variables from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_pm_fab)   ; NC_names.push_back("pm")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 0
    NC_fabs.push_back(&NC_pn_fab)   ; NC_names.push_back("pn")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 1
    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs, false, 0, ngrow);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_fcor_fab     container for Coriolis parameter data
 */
void
read_coriolis_from_netcdf (int /*lev*/,
                           const Box& domain,
                           const std::string& fname,
                           FArrayBox& NC_fcor_fab)
{
    amrex::Print() << "Loading initial coriolis from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_fcor_fab )   ; NC_names.push_back("f")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 0

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_mskr_fab     container for rho-point land/sea mask data
 * @param NC_msku_fab     container for u-point land/sea mask data
 * @param NC_mskv_fab     container for v-point land/sea mask data
 */
void
read_masks_from_netcdf (int /*lev*/,
                        const Box& domain,
                        const std::string& fname,
                        FArrayBox& NC_mskr_fab,
                        FArrayBox& NC_msku_fab,
                        FArrayBox& NC_mskv_fab)
{
    amrex::Print() << "Loading masks from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_mskr_fab )   ; NC_names.push_back("mask_rho")  ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 0
    NC_fabs.push_back(&NC_msku_fab )   ; NC_names.push_back("mask_u")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 1
    NC_fabs.push_back(&NC_mskv_fab )   ; NC_names.push_back("mask_v")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 2

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param lev                 level of data to read
 * @param domain              simulation domain
 * @param fname               file name to read from
 * @param do_m2_clim_nudg     whether to do 2d momentum climatology nudging
 * @param do_m3_clim_nudg     whether to do 3d momentum climatology nudging
 * @param do_cons_clim_nudg   per cons component, whether to do climatology nudging
 * @param cons_names          per cons component, name of the tracer
 * @param NC_M2NC_fab         container for 2d momentum climatology data
 * @param NC_M3NC_fab         container for 3d momentum climatology data
 * @param NC_ConsNC_fab       per cons component, container for tracer climatology data
 * @param cons_coeff_in_file  filled on output: per cons component, whether a spatially
 *                            varying coefficient was found in the file
 */
void
read_clim_nudg_coeff_from_netcdf (int /*lev*/,
                        const Box& domain,
                        const std::string& fname,
                        bool do_m2_clim_nudg,
                        bool do_m3_clim_nudg,
                        const amrex::Vector<int>& do_cons_clim_nudg,
                        const amrex::Vector<std::string>& cons_names,
                        FArrayBox& NC_M2NC_fab,
                        FArrayBox& NC_M3NC_fab,
                        amrex::Vector<FArrayBox>& NC_ConsNC_fab,
                        amrex::Vector<int>& cons_coeff_in_file)
{
    amrex::Print() << "Loading nudging coefficients from NetCDF file " << fname << std::endl;

    const int l_ncons = do_cons_clim_nudg.size();

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    if (do_m3_clim_nudg) {
        NC_fabs.push_back(&NC_M3NC_fab ); NC_names.push_back("M3_NudgeCoef"); NC_dim_types.push_back(NC_Data_Dims_Type::BT_SN_WE);
    }
    if (do_m2_clim_nudg) {
        NC_fabs.push_back(&NC_M2NC_fab ); NC_names.push_back("M2_NudgeCoef"); NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE);
    }
    // Each tracer's spatially varying coefficient is stored under its own name, as in
    // ROMS: temp_NudgeCoef, salt_NudgeCoef, NO3_NudgeCoef, ... A tracer whose field is
    // absent from the file keeps the constant coefficient derived from remora.tnudg, so
    // a file written for temp and salt alone still works when biology tracers are nudged.
    cons_coeff_in_file.assign(l_ncons, 0);
    Vector<std::string> missing;
    for (int icomp = 0; icomp < l_ncons; ++icomp) {
        if (!do_cons_clim_nudg[icomp]) { continue; }
        const std::string coeff_name = cons_names[icomp] + "_NudgeCoef";
        if (!QueryNetCDFHasVars(fname, {coeff_name})) {
            missing.push_back(coeff_name);
            continue;
        }
        cons_coeff_in_file[icomp] = 1;
        NC_fabs.push_back(&NC_ConsNC_fab[icomp]); NC_names.push_back(coeff_name); NC_dim_types.push_back(NC_Data_Dims_Type::BT_SN_WE);
    }
    if (!missing.empty()) {
        amrex::Print() << "Nudging coefficient file " << fname << " does not contain";
        for (const auto& name : missing) { amrex::Print() << " " << name; }
        amrex::Print() << "; those tracers will use the constant coefficient from remora.tnudg"
                       << std::endl;
    }

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param lev            level of data to read
 * @param fnames         file name(s) to read from
 * @param field_name     field name to read
 * @param vec_dat        vector to fill data
 */
//template <typename DType>
void read_vec_from_netcdf (int /*lev*/, const amrex::Vector<std::string>& fnames, const std::string& field_name, amrex::Vector<int>& vec_dat)
{
    AMREX_ALWAYS_ASSERT(!fnames.empty());
    const std::string& fname = fnames[0];

    amrex::Print() << "Reading " << field_name << " from NetCDF file" << std::endl;

    // get x-positions and put in array
    using ARRAY = NDArray<int>;
    amrex::Vector<ARRAY> array_dat(1);
    ReadNetCDFFile(fname, {field_name}, array_dat); // filled only on proc 0
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        int n = array_dat[0].get_vshape()[0];
        for (int i(0); i < n; i++)
        {
            vec_dat.push_back((*(array_dat[0].get_data() + i)));
        }
    }
    int nvals = vec_dat.size();
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&nvals,1,ioproc);
    if (!(amrex::ParallelDescriptor::IOProcessor())) {
        vec_dat.resize(nvals);
    }
    amrex::ParallelDescriptor::Bcast(vec_dat.data(), vec_dat.size(), ioproc);
}

#endif // ROMSX_USE_NETCDF
