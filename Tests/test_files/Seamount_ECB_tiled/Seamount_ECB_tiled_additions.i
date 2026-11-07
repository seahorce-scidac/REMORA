# Added to Tests/test_files/Seamount/Seamount.i by compose_test_input in
# CTestList.cmake, so the Seamount case itself is defined once. Settings here
# override the base's, since ParmParse reads the last occurrence.
#
# Seamount's terrain-following layers sit at different depths over the bump, so
# temperature and layer thickness at each level vary in both x and y, and so does
# the biology after one step. The test runs this input untiled and with small CPU
# tiles and requires the two runs to agree bitwise, so a kernel that mixes up
# columns across tile edges in either direction fails.
remora.max_step = 3
remora.plot_int = 3
# Fill new arrays with signaling NaN and trap invalid operations, so a biology
# kernel that reads scratch it has not written fails instead of using garbage.
amrex.init_snan = 1
amrex.fpe_trap_invalid = 1
remora.bulk_fluxes = true
remora.surface_radiation_flux = 150.0
remora.uwind = 5.0
remora.vwind = 2.0
remora.biology_ic_type = constant
remora.plot_vars_3d = salt temp x_velocity y_velocity z_velocity biology
remora.biology_model = ecb
# Two iterations, so state carried between iterations is checked under tiling too.
remora.ecb.BioIter = 2
remora.biology_initial_NO3 = 8.0
remora.biology_initial_NH4 = 0.5
remora.biology_initial_chlorophyll = 1.0
remora.biology_initial_phytoplankton = 0.3
remora.biology_initial_zooplankton = 0.1
remora.biology_initial_LdetritusN = 0.2
remora.biology_initial_SdetritusN = 0.1
remora.biology_initial_DON = 3.0
remora.biology_initial_PO4 = 0.5
remora.biology_initial_LdetritusC = 1.2
remora.biology_initial_SdetritusC = 0.6
remora.biology_initial_TIC = 2100.0
remora.biology_initial_alkalinity = 2300.0
remora.biology_initial_DOC = 50.0
remora.biology_initial_oxygen = 200.0
remora.biology_initial_Rdon = 1.0
remora.biology_initial_Rdoc = 40.0
remora.biology_initial_calcium = 10000.0
