SHELL=/bin/bash

SCUFFEM = $(HOME)/git/scuff-em-installation
SCUFFEM2 = $(HOME)/git/scuff-em/src/libs/libscuff

NLOPT=$(HOME)/git/nlopt-installation

CXXFLAGS = -O2 -fopenmp -lm -std=c++11
#CXXFLAGS = -O0 -ggdb -lm -std=c++11

CPPFLAGS += -I$(SCUFFEM)/include/scuff-em
CPPFLAGS += -I$(SCUFFEM2)
CPPFLAGS += -I$(NLOPT)/include

LDFLAGS += -L$(SCUFFEM)/lib -Wl,-rpath,${SCUFFEM}/lib
LDFLAGS +=-L$(NLOPT)/lib 

##################################################
# catch nan and inf values 
export SCUFF_ABORT_ON_FPE=1
##################################################
# lapack/blas libraries
##################################################
export LB_LIBS=-llapack -lblas -lgfortran -lgomp

##################################################
# HDF5 libraries
##################################################
export HDF5_LIBS=-lhdf5 -lhdf5_hl

##################################################
# HDF5 libraries
##################################################
export NLOPT_LIBS=-lnlopt

##################################################
# Homer's  libraries #############################
##################################################
HR_LIBS = -lscuff 
##HR_LIBS = -lscuff -lIncFieldPrime -lMatProp -lMDInterp -lhmat -lSGJC -lTriInt -lhrutil

LIBS = $(HR_LIBS) $(HDF5_LIBS) $(LB_LIBS)
LIBS2 = $(HR_LIBS) $(NLOPT_LIBS)


OT = OptTorque.o OTutils.o VBeam.o GHBeam.o OutputModules.o
OptTorque: 	$(OT)
		$(CXX) -Wall $(LDFLAGS) -o $@ $(OT) $(LIBS)

VBT = VBeamTest.o VBeam.o
VBeamTest: 	$(VBT)
		$(CXX) -Wall $(LDFLAGS) -o $@ $(VBT) $(LIBS)

VBTM = VBeamTestMesh.o VBeam.o
VBeamTestMesh:	$(VBTM)
		$(CXX) -Wall $(LDFLAGS) -o $@ $(VBTM) $(LIBS)

GBT = GHBeamTest.o GHBeam.o
GHBeamTest: 	$(GBT)
		$(CXX) -Wall $(LDFLAGS) -o $@ $(GBT) $(LIBS)

OBJ = objective.o OTutils.o VBeam.o
objective: 	$(OBJ)
		$(CXX) -Wall $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

TUT = NLOPT_tutorial.o
tutorial: 	$(TUT)
		$(CXX) $(LDFLAGS) -o $@ $(TUT) $(LIBS2)

clean:
		rm *.o 
