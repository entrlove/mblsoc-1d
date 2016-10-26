# Setting the environment for SLEPc and PETSc
export PETSC_DIR=/usr/local/mathlib/petsc/petsc-3.7.4
export PETSC_ARCH=arch-darwin-cxx-complex-debug
export SLEPC_DIR=/usr/local/mathlib/slepc/slepc-3.7.3

# Setting the MKL
MKL_COMPILER=-m64 -I${MKLROOT}/include
MKL_LINKER=-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl 

# Setting the Eigen and MPIR 
MATH_BASE_PATH=/usr/local/mathlib
EIGEN_COMPILER=-I$(MATH_BASE_PATH)/eigen/eigen-3.2.10
MPIR_COMPILER=-I$(MATH_BASE_PATH)/mpir/mpir-2.7.2/include
MPIR_LINKER=-L$(MATH_BASE_PATH)/mpir/mpir-2.7.2/lib -Wl,-rpath,$(MATH_BASE_PATH)/mpir/mpir-2.7.2/lib -lmpir 

#Setting the project
PROJECT_COMPILER=-std=c++11 -I./inc

VPATH = inc src
sources := $(notdir $(foreach dir, $(VPATH),$(wildcard $(dir)/*.cpp)))
objects := $(patsubst %.cpp, %.o,  $(sources))

CXXFLAGS = $(PROJECT_COMPILER) $(EIGEN_COMPILER) $(MKL_COMPILER) $(MPIR_COMPILER)
MYLIBS   = $(MKL_LINKER) $(MPIR_LINKER) 
MYOBJS   = $(objects) 
EXE      = mbl

all: ${EXE}

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

${EXE}: ${MYOBJS} chkopts
	-${CLINKER} -o ${EXE} ${MYOBJS} $(MYLIBS) ${SLEPC_LIB}
	${RM} ${MYOBJS}
