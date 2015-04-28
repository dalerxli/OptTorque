SHELL=/bin/bash

SCUFFEM = $(HOME)/git/scuff-em-installation
SCUFFEM2 = $(HOME)/git/scuff-em/src/libs/libscuff

CXXFLAGS = -O2 -fopenmp -lm -std=c++11
#CXXFLAGS = -O0 -ggdb -fopenmp -lm -std=c++11

CPPFLAGS += -I$(SCUFFEM)/include/scuff-em
CPPFLAGS += -I$(SCUFFEM2)

LDFLAGS += -L$(SCUFFEM)/lib -Wl,-rpath,${SCUFFEM}/lib

##################################################
# lapack/blas libraries
##################################################
export LB_LIBS=-llapack -lblas -lgfortran -lgomp

##################################################
# HDF5 libraries
##################################################
export HDF5_LIBS=-lhdf5 -lhdf5_hl

##################################################
# HR libraries ###################################
##################################################

 HR_LIBS = -lscuff 
##HR_LIBS = -lscuff -lIncFieldPrime -lMatProp -lMDInterp \
          -lhmat -lSGJC -lTriInt -lhrutil

LIBS = $(HR_LIBS) $(HDF5_LIBS) $(LB_LIBS)

OT = OptTorque.o GHBeam.o OutputModules.o  GetIntegratedIntensity.o

OptTorque12: 	$(OT)
		$(CXX) -Wall $(LDFLAGS) -o $@ $(OT) $(LIBS)

clean:
		rm *.o 
