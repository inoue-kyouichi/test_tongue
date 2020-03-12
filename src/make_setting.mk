##################################################################################
#
# BDIM library

# Copyright (c) 2019- Mechanical and Bioengineering Systems Lab.,
#                     Department of Mechanical Science and Bioengineering,
#                     Graduate School of Engineering Science,
#                     Osaka University.
# All rights reserved.
#
###################################################################################

########################
# Common setting
########################

LIBDIR = /home/totani/lib
INSTALL_BASE = /home/totani/lib
TARGET_DIR  = ../../../

#########################
# Architecture assignment
#########################

# gnu, intel, fx, sx, test (icpc)
ARCH=test

########################
# For each architecture
########################

ifeq ($(ARCH),gnu)
  AR          = ar cru
  RANLIB      = ranlib
  RM          = \rm -f
#  MPI_DIR     = /opt/openmpi
#  MPI_DIR     = /Volumes/Transcend/libs/openmpi-1.10.2
  MPI_DIR = /home/totani/libs
  OMP_FLAGS   = -fopenmp
  UDEF_OPT    =
 CC          = mpicc
 CFLAGS      = -O3
 CXX         = mpicxx
 CXXFLAGS    = -O3 $(OMP_FLAGS)
 FC          = mpif90
 FCFLAGS     = -O3 -ffree-form
 F90         = mpif90
 F90FLAGS    = -O3 -cpp $(OMP_FLAGS) -ffree-form
 #F90FLAGS    = -O3 -cpp  -ffree-form
  LDFLAGS     =
  #LIBS        = -lgfortran
  LIBS        =
  UDEF_LIB_PATH_SPEC =
  UDEF_LIBS_SPEC     =

endif

########################
ifeq ($(ARCH),test)
  AR          = ar cru
  RANLIB      = ranlib
  RM          = \rm -f
  #  MPI_DIR     = /opt/openmpi
  #MPI_DIR     = /Volumes/Transcend/libs/openmpi-1.10.2
  #MPI_DIR = /home/totani/lib/openmpi-2.0.2
  OMP_FLAGS   = -qopenmp
  UDEF_OPT    =
    CC          = icc
    CXX         = icpc
  #CC          = mpicc
  #CXX         = mpicxx
  #CFLAGS      = -O3
  #CXXFLAGS    = -O3 -qopt-report=5 $(OMP_FLAGS)
  CXXFLAGS    = -Wall -Wextra -O3 $(OMP_FLAGS) -std=c++11 -MMD -MP -mkl
  #FC          = ifort
  #FCFLAGS     = -O3
  #F90         = mpif90
  #F90FLAGS    = -O3 -Warn unused -fpp -qopt-report=5 $(OMP_FLAGS)
  #F90FLAGS    = -O3 -Warn unused  -fpp  $(OMP_FLAGS) -traceback -CB -fpe0
  #F90FLAGS    = -O3 -Warn unused -fpp   -traceback -CB -fpe0
  LDFLAGS     = 
  #TEXTPARSER
  TXTP    = $(LIBDIR)/TextParser-1.8.5
  LDFLAGS  += -L$(TXTP)/lib -lTP
  CXXFLAGS  += -I$(TXTP)/include
  #Eigen
  CXXFLAGS  += -I$(LIBDIR)/eigen-3.3.4
  #glog
  GLOG      = $(LIBDIR)/glog
  LDFLAGS += -L$(GLOG)/lib -lglog
  CXXFLAGS += -I$(GLOG)/include

#  LIBS        = -lifport -lifcore
  LIBS        =
  UDEF_LIB_PATH_SPEC =
  UDEF_LIBS_SPEC     =

endif

########################
ifeq ($(ARCH),intel)
  AR          = ar cru
  RANLIB      = ranlib
  RM          = \rm -f
  #  MPI_DIR     = /opt/openmpi
  MPI_DIR = /home/totani/lib/openmpi-3.0.0
  OMP_FLAGS   = -qopenmp
  UDEF_OPT    =
  CC          = mpicc
  CXX         = mpicxx
  #CFLAGS      = -O3
  #CXXFLAGS    = -Wall -Wextra -O3 $(OMP_FLAGS) -std=c++11 -g -MMD -MP -pthread
  CXXFLAGS    = -Wall -Wextra -O3 $(OMP_FLAGS) -std=c++11
  #FC          = ifort
  #FCFLAGS     = -O3
  #F90         = mpif90
  #F90FLAGS    = -O3 -Warn unused -fpp -qopt-report=5 $(OMP_FLAGS)
  #F90FLAGS    = -O3 -Warn unused  -fpp  $(OMP_FLAGS) -traceback -CB -fpe0
  #F90FLAGS    = -O3 -Warn unused -fpp   -traceback -CB -fpe0
  LDFLAGS     = 
#  LIBS        = -lifport -lifcore
  LIBS        =
  UDEF_LIB_PATH_SPEC =
  UDEF_LIBS_SPEC     =

endif

########################
ifeq ($(ARCH),fx)
  AR          = ar cr
  RANLIB      = ranlib
  RM          = \rm -f
  MPI_DIR     =
  OMP_FLAGS   =
  UDEF_OPT    = -D__K_FPCOLL -D__ARCH_FX
  CC          =  mpifccpx
  CFLAGS      = -Kfast,ocl,preex,simd=2,uxsimd,array_private,parallel,openmp -Nsrc
  CXX         = mpiFCCpx
  CXXFLAGS    = -Kfast,ocl,preex,simd=2,uxsimd,array_private,parallel,openmp,optmsg=2 -V -Nsrc
  FC          = mpifrtpx
  FCFLAGS     = -Cpp -Kfast,ocl,preex,simd=2,uxsimd,array_private,auto,parallel,openmp -Qt
  F90         = mpifrtpx
  F90FLAGS    =  -Kfast,parallel,openmp -Qt
  LDFLAGS     =
  LIBS        =
  UDEF_LIB_PATH_SPEC =
  UDEF_LIBS_SPEC     =

## iff double
#CFLAGS     += -D_REAL_IS_DOUBLE_
#CXXFLAGS   += -D_REAL_IS_DOUBLE_
#FCFLAGS    += -CcdRR8
#F90FLAGS   += -CcdRR8

endif

########################
ifeq ($(ARCH),sx)
  AR          = ar cr
  RANLIB      = ranlib
  RM          = \rm -f
  MPI_DIR     =
  OMP_FLAGS   =
  UDEF_OPT    =
  CC          =
  CFLAGS      =
  CXX         =
  CXXFLAGS    =
  FC          =
  FCFLAGS     =
  F90         = sxmpif90
  F90FLAGS    =  -f5 -Cvopt -sxace -R2 -Wf,-pvctl fullmsg -pi -Pauto -f2003
  LDFLAGS     =
  LIBS        =
  UDEF_LIB_PATH_SPEC =
  UDEF_LIBS_SPEC     =

## iff double
#CFLAGS     += -D_REAL_IS_DOUBLE_
#CXXFLAGS   += -D_REAL_IS_DOUBLE_
#FCFLAGS    += -CcdRR8
#F90FLAGS   += -CcdRR8

endif
########################


UDEF_INC_PATH=-I. \
                                -I../../solvers/shapeOptimization/minimumComplianceProblem \
                                -I../../solvers/PDL_analysis/PDL_TOOTH_interaction \
                                -I../../solvers/PDL_analysis/humanPDL \
                                -I../../solvers/PDL_analysis/ratPDL \
                                -I../../solvers/FEM \
                                -I../../solvers/linearSolver \
                                -I../../solvers/RBD \
                                -I../../solvers/base \

UDEF_LIB_PATH= \
                                  -L../../solvers/lib -lMinimumCompliance \
                                  -L../../solvers/lib -lInteraction \
                                  -L../../solvers/lib -lPDL_human \
                                  -L../../solvers/lib -lPDL_rat \
                                  -L../../solvers/lib -lFEM \
                                  -L../../solvers/lib -lRBD \
                                  -L../../solvers/lib -lLIS \
                                  -L../../solvers/lib -lBase \

#for SX-ACE?
# UDEF_LIB_PATH= \
#                                   ./FDM/*.o \
#                                   ./boundary/*.o \
#                                   ./VOF/*.o \
#                                   ./evaluation/*.o \
#                                   ./fileIO/*.o \
#                                   ./linearSolver/*.o \
#                                   ./MPI/*.o \
 UDEF_LIBS = -lmpi $(UDEF_LIBS_SPEC)
 UDEF_LIBS_UTIL = -lmpi $(UDEF_LIBS_SPEC)
