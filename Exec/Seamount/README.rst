The Seamount problem is a standard test in
`ROMS <https://www.myroms.org/wiki/SEAMOUNT_CASE>`_.
The problem involves an (analytically) stably stratified fluid at rest over a
seamount. In the absence of numerical errors, the fluid will remain at rest.
However, this may not occur due to numerical errors in the calculation of the
horizontal pressure gradient when the vertical coordinates are misaligned with
the geopotential surfaces, as is the case in problems with spatially-varying
bathymetry in ROMS/REMORA.

