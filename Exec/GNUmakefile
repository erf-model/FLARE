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

TEST = TRUE
USE_ASSERTION = TRUE

# Enable NETCDF support (required for WRFInput or MetGrid initialization)
USE_NETCDF = TRUE

# Enable FFT Poisson solver
USE_FFT = FALSE

# Enable NOAH LSM
#USE_NOAHMP = TRUE

# Enable RRTMGP
#USE_RRTMGP = TRUE

# Enable particles (required for ParticleTests, MultiSpeciesBubble or inputs_SDM in RICO)  
USE_PARTICLES = FALSE

# FLARE: the fire engine (FireEngine/) and, with USE_DUST = TRUE, the dust layer (Dust/)
USE_FIRE = TRUE
USE_DUST = FALSE

# GNU Make: ERF is Submodules/ERF; the problem file is FLARE's own (this directory)
FLARE_HOME := ..
ERF_HOME := $(FLARE_HOME)/Submodules/ERF
ERF_PROBLEM_DIR = $(FLARE_HOME)/Exec
ifeq ($(USE_FIRE), TRUE)
    FIRE_HOME := $(FLARE_HOME)/FireEngine
endif
ifeq ($(USE_DUST), TRUE)
    DUST_HOME := $(FLARE_HOME)/Dust
endif
include $(ERF_HOME)/Exec/Make.ERF
