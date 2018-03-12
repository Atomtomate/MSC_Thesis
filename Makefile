# ========= utility fucntions =========
# Check that given variables are set and all have non-empty values,
# die with an error otherwise.
#
# Params:
#   1. Variable name(s) to test.
#   2. (optional) Error message to print.
check_defined = \
    $(strip $(foreach 1,$1, \
        $(call __check_defined,$1,$(strip $(value 2)))))
__check_defined = \
    $(if $(value $1),, \
      $(error Undefined $1$(if $2, ($2))))


# ========== build variables ========== 
CXX       := mpic++
TARGET    := ImpSolv
ROOT_DIR  := $(shell basename $(CURDIR))
OBJ_DIR   := obj
SRC_DIR   := src
CPP_FILES := $(filter-out $(wildcard $(SRC_DIR)/test_*.cpp), $(wildcard $(SRC_DIR)/*.cpp))
OBJ_FILES  = $(addprefix $(OBJ_DIR)/,$(notdir $(CPP_FILES:.cpp=.o)))
VPATH    := src:src/examples



# ========== compiler flags, includes ========== 
# TODO: remember to use  -guide and  -profile-loops=all 
CXXFLAGS.DEBUG    := -O2 -pg -ggdb -Wextra -D EIGEN_NO_STATIC_ASSERT # -debug parallel  -ftrapuv -traceback -parallel
CXXFLAGS.PROFILE  := -O3 -pg -D EIGEN_NO_STATIC_ASSERT #-xHost -fno-alias -profile-functions -profile-loops=all -profile-loops-report=2 
#-profile-functions -profile-loops=all -profile-loops-report=2 -guide -ipo -no-prec-div
CXXFLAGS.PARALLEL := -O2 -Wall -fno-alias -parallel
CXXFLAGS.FAST     := -O3 -parallel -g -xHost -ansi-alias -ipo -prof-use -D EIGEN_NO_DEBUG
#-mtune=native -xhost -Wall
CXXFLAGS.FP       := -fimf-precision=3 -fp-model=precise -fp-model source -prec-div -prec-sqrt
CXXFLAGS          := -Wno-deprecated -std=c++17 -D DEBUG_MODE -DMKL_ILP64
CLANGFLAGS        := -O2 -std=c++17 -ferror-limit=0 -L/usr/lib64 -lstdc++ -v

#  libraries
# mathgl2 and mathgl-qt needed
# set prefix_dir in .bashrc (on systems with sudo permissions this can be omitted) 
#ifndef PREFIX_DIR
#	$(warning "PREFIX_DIR variable not set")
#	PREFIX_DIR := /usr
#endif
PREFIX_DIR = /usr/lib/x86_64-linux-gnu
BOOST_H := $(PREFIX_DIR)/include/boost
BOOST_L := $(PREFIX_DIR)/lib
EIGEN_H := ./libs
FFTW_H  := ./libs/fftw/include
FFTW_L  := ./libs/fftw/lib
MPI_H   := $(MPI_INCLUDE)
MPI_L   := $(MPI_LIB)

LIBS     := -L/usr/share -L/usr/lib -L/usr/local/lib -L/opt/local/lib -L/usr/lib/x86_64-linux-gnu -L$(FFTW_L) -L$(BOOST_L) -L$(PREFIX_DIR)
#  -I$(MKL_ROOT)/include-L$(MKLROOT)/lib/intel64    -L/home/julian/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/mic
#-L $(SLEPC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib
INCLUDES := -I./includes -I./src/examples -I./examples -I/usr/include -I$(BOOST_H) -I./$(SRC_DIR) -I$(EIGEN_H) -I$(FFTW_H) -I$(PREFIX_DIR)/include
# -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I/home/julian/intel/compilers_and_libraries_2016.3.210/linux/mkl/include 
LDFLAGS  :=  -ltrng4 -l:libfftw3.so -lm -lboost_system -lboost_iostreams -lboost_filesystem -lboost_mpi -lboost_mpi -lboost_serialization -lgsl
#-Wl,-rpath,$(SLEPC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib  -lpetsc -lslepc -lquadmath  -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -ldl 

# ========== rules ==========
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(INCLUDES) -c $< $(CXXFLAGS) -o $@ 
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.tpp
	$(CXX) $(INCLUDES) -c $< $(CXXFLAGS) -o $@
$(OBJ_DIR):
	mkdir $(OBJ_DIR)

clangTest: CXX = clang
clangTest: CXXFLAGS = $(CLANGFLAGS)
clangTest: $(OBJ_DIR) $(OBJ_FILES)
	clang $(CLANGFLAGS) $(filter-out $(OBJ_DIR),$^) $(LIBS) $(LDFLAGS) -o ImpSolv_$@

debug: CXXFLAGS += $(CXXFLAGS.DEBUG)
debug: $(OBJ_DIR) $(OBJ_FILES)
	$(CXX) $(filter-out $(OBJ_DIR),$^) $(CXXFLAGS) $(LIBS) $(LDFLAGS) -o ImpSolv_$@

release: CXXFLAGS += $(CXXFLAGS.FAST)
release: $(OBJ_DIR) $(OBJ_FILES)
	$(CXX) $(filter-out $(OBJ_DIR),$^) $(CXXFLAGS) $(LIBS) $(LDFLAGS) -o ImpSolv_$@

profile: CXXFLAGS += $(CXXFLAGS.PROFILE)
profile: $(OBJ_DIR) $(OBJ_FILES)
	$(CXX) $(filter-out $(OBJ_DIR),$^) $(CXXFLAGS) $(LIBS) $(LDFLAGS) -o ImpSolv_$@

.PHONY : clean
clean:
	rm -f $(OBJ_FILES) $(OBJ_DIR)/$(OBJ_FILES) $(TARGET) ImpSolv ImpSolv_* *.out *.png *.eps
