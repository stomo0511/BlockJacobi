PROJECT_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

CXX = /usr/local/bin/g++-9
CPPFLAGS = -fopenmp -O2
LDFLAGS = -lgomp

MKL_ROOT =  /opt/intel/compilers_and_libraries/mac/mkl
MKL_INC_DIR = $(MKL_ROOT)/include
MKL_LIB_DIR = $(MKL_ROOT)/lib
MKL_LIBS = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

OBJS = BlockJacobi.o Jacobi.o

ifeq ($(BUILD_MODE),debug)
	CFLAGS += -g
else ifeq ($(BUILD_MODE),run)
	CFLAGS += -O2
else
	$(error Build mode $(BUILD_MODE) not supported by this Makefile)
endif

all:	BlockJacobi

BlockJacobi:	$(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS) -L$(MKL_LIB_DIR) $(MKL_LIBS)

%.o:	$(PROJECT_ROOT)%.cpp
	$(CXX) -c $(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) -I$(MKL_INC_DIR) -o $@ $<

clean:
	rm -fr BlockJacobi $(OBJS)
