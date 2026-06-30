# AMReX
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE       = FALSE
TINY_PROFILE  = FALSE
COMM_PROFILE  = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE   = FALSE
USE_GPROF     = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE

USE_CUDA = FALSE
USE_HIP  = FALSE
USE_SYCL = FALSE

# Enable debug flags (rich backtrace)
DEBUG = FALSE
DEBUG = TRUE

TEST = TRUE
USE_ASSERTION = TRUE

# Enable NETCDF support (required for WRFInput or MetGrid initialization)
USE_NETCDF = TRUE

# Enable FFT Poisson solver
USE_FFT = TRUE

# Enable NOAH LSM
USE_NOAHMP = FALSE

# Enable NOAH LSM
USE_NATIVE_SHOC = FALSE

# Enable particles (required for ParticleTests, MultiSpeciesBubble or inputs_SDM in RICO)  
USE_PARTICLES = TRUE

# GNU Make
Bpack := ./Make.package
Blocs := .

REMORA_HOME := ../
REMORA_PROBLEM_DIR = $(REMORA_HOME)/Exec
include $(REMORA_HOME)/Exec/Make.REMORA
