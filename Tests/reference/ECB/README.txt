ROMS ECB reference fixtures

Source: https://github.com/myroms/roms
Branch: develop (the default/trunk branch when retrieved 2026-11-03)
Commit: 31cabd318c9379429c365c308b9e0f2ee166a6fc
Unmodified source SHA-256:
  ROMS/Nonlinear/Biology/ecb.h
    d592bc01404fe9b44ba4f83f91994ce33e99cb7a350ce43ba587f33d4a56e38b
  ROMS/Nonlinear/Biology/ecb_mod.h
    6a1af35b587d7ba00298aa29b0a34f216a3f6741baf69f9c1717792d73d946a2

Options: CARBON, OXYGEN, BIO_SEDIMENT, DENITRIFICATION, PO4,
TALK_NONCONSERV, RW14_OXYGEN_SC, RW14_CO2_SC, MASKING.
BULK_FLUXES tested both on and off. All other optional branches are off.
The first four options are implicit requirements of the ECB kernel.

Regenerate with NumPy, cpp, gfortran and a ROMS clone containing that commit:
  python3 generate_reference.py --roms /path/to/roms

The generator reads git objects, never modified working files. It retains the
reference parameter module, ecb_tile and pCO2_water routines. A one-column
adapter supplies the surrounding ROMS scalars, bounds and arrays. It adds the
missing fixed-form continuation marker to the bustr/bvstr argument line,
flattens dual-form continuations, and replaces set_bounds.h with explicit tile
bounds. No equations, process order, clipping or reference constants change.
The unused caldate result is stubbed because atmospheric pCO2 is constant.

The Fortran compiler runs with -O0 -fcheck=all. It emits 24 independent cases
covering 2/3/8/16 layers, 0/1/2/3 iterations, day/night, wind/surface stress,
zero/partial/full resuspension, low oxygen, low nutrients, TIC constraints,
and dry/fractional wet masks. Long steps also exercise characteristics that
cross multiple layers and clear a whole column. Incoming weighted state includes an offset of
7 in every component, so overwriting it with scratch concentrations fails.
Refractory DON/DOC and calcium have zero biological increments in this source.

columns.txt holds parameter names/values followed by case controls, thickness,
interfaces, old concentrations, incoming weighted state, and expected weighted
state. The light control is the kinematic srflx (degC m/s) that ROMS consumes;
REMORA's shortwave is in W/m2, so the native test multiplies it by rho0*Cp
(1025*3985) before calling the model. Matrices are component-major, bottom to top. Synthetic parameters and
profiles exercise the equations; they are not a calibrated scientific setup.
The native test varies dye offsets and checks preservation of physical fields,
dyes, incoming increments, old state and inert reference pools. Expected values
are checked in; production builds and unit tests do not require Fortran/NumPy.
Double-precision tolerance is 128 machine eps times max(1, abs(reference));
single-precision tests use their own machine epsilon with the same scale.
