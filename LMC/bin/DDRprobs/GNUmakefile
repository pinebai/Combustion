#TOP = /scratch2/marc/Valentina
TOP =${HOME}/src/CCSE
#TOP =${HOME}
BOXLIB_DIR = ${TOP}/BoxLib
IAMR_DIR = ${TOP}/IAMR
COMBUSTION_DIR = ${TOP}/Combustion
#COMBUSTION_DIR = /scratch2/marc/Valentina/Combustion

#
# Build configuration
#

PRECISION       = DOUBLE
DEBUG           = TRUE
DEBUG           = FALSE
DIM             = 2
COMP            = g++
FCOMP           = gfortran
USE_MPI         = TRUE
#USE_MPI         = FALSE
USE_OMP         = FALSE
VERBOSE         = FALSE
VERBOSE         = TRUE
PROFILE         = FALSE
COMM_PROFILE    = FALSE
USE_EFIELD      = TRUE
#USE_WBAR      	= TRUE
#USE_WBAR      	= FALSE

# Enable SDC code - add .SDC to executable
# and compile in the #ifdef LMC_SDC
USE_LMC_SDC = FALSE
USE_LMC_SDC = TRUE

# Enable flctin namelist
#USE_FLCTS       = TRUE

#USE_PARTICLES   = FALSE

# Enable ht.velocity_plotfile option
#USE_VELOCITY    = FALSE

# Enable time-dependent forcing term
#USE_TDEPFORCING = FALSE

# Choose model (from list below), and pmf file
#CHEMISTRY_MODEL = DRM19
CHEMISTRY_MODEL = GRI30ION
#CHEMISTRY_MODEL = DODECANE_LU
#CHEMISTRY_MODEL = WANGDODECANE
#fEXE_sources += wang_dodecane_phi_1000.f
#fEXE_sources += drm19Soln_seed_1.00.f
fEXE_sources += gri30ion_phi_1.00.f
#fEXE_sources += gri30ion_phi_1.00_Le1.f
Blocs   := .

#fEXE_sources += pmf.f				# start with charge density from 1D code
#Blocs   := ../../../../Downloads/RHTtoPMF	# start with charge density from 1D code

CEXE_sources +=
FEXE_sources += PROB_$(DIM)D.F PROB_F.F
CEXE_headers +=
FEXE_headers += probdata.H PROB_F.H

include $(COMBUSTION_DIR)/LMC/tools/make/Make.LMC

# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc

#DEFINES += -DDO_JBB_HACK
DEFINES += -DTRIGGER_NEW_J
#DEFINES += -DALWAYS_NEW_J
#DEFINES += -DSOLN_IS_1D

#DEFINES += -DDO_JBB_HACK_TEMP
#DEFINES += -DHACK_TEMP_MIN=273.d0
#DEFINES += -DHACK_TEMP_MAX=1900.d0

#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)
