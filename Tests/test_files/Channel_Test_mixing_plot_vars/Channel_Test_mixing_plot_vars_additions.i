# The vertical mixing fields at initialization, written to a plotfile. tke and gls start uniform at
# gls_Kmin and gls_Pmin. Akv and Akt start at their background values on the interior w-levels
# and zero on the surface and bottom ones, so their cell averages run from half the background,
# in the top and bottom cells, to the full background.
remora.max_step = 0
remora.plot_int = 1
remora.plot_vars_3d = temp Akv Akt tke gls

remora.gls_Kmin = 1.0e-5
remora.gls_Pmin = 2.0e-6
remora.Akv_bak = 2.0e-4
remora.Akt_bak = 1.0e-4
