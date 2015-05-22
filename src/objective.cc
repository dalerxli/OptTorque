/*---------------------------------------------------------------
 * objective.cc  
 * compute Torque FOM
 * for a given BEM matrix, OPFT matrix, frequency, and PARMMatrix
 * created 2015.03
 * last updated on 2015.05.22
 *--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <complex>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <libhrutil.h>
#include <libscuff.h>
#include "VBeam.h" 

#ifndef cdouble
  typedef std::complex<double> cdouble;
#endif

#define II cdouble(0.0,1.0)
using namespace scuff;

void ShowPARMMatirx(HMatrix* PARMMatrix); 
double objective(RWGGeometry *G, HMatrix *PARMMatrix,
		 HMatrix *M, HMatrix *Q); 
//---------------------------------------------------------------//
//---------------------------------------------------------------//
///LEFTOVER TASKS
//---------------------------------------------------------------//

int main(int argc, char *argv[])
  /// Objective 
  /// 
  /// READ In:
  ///    parameter matrix PARMMatrix 
  ///    BEM matrix A, LU Factorized
  ///    PFT matrix Q 
{
  /*--------------------------------------------------------------*/
  /*- process options  -------------------------------------------*/
  /*--------------------------------------------------------------*/
  /// Input arguments should be filenames. (pointers) 
  char *GeoFile=0; 
  char *PARMMatrixFile=0;
  char *BEMMatrixFile=0;
  char *PFTMatrixFile=0;
  char *LogLevel=0;

  /* name    type   #args  max_instances  storage  count  description*/
  OptStruct OSArray[]=
   { 
     {"geometry",   PA_STRING,  1, 1, (void *)&GeoFile,   0,  ".scuffgeo file"},
     {"PARMMatrix", PA_STRING,  1, 1, (void *)&PARMMatrixFile, 0, "VParameters file"},
     {"BEMMatrix",  PA_STRING,  1, 1, (void *)&BEMMatrixFile, 0,  "BEM Matrix datafile"},
     {"PFTMatrix",  PA_STRING,  1, 1, (void *)&PFTMatrixFile, 0,  "PFT Matrix datafile"},
     {"LogLevel",   PA_STRING,  1, 1, (void *)&LogLevel,   0, "none | terse | verbose | verbose2"},
     //
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  if (PARMMatrixFile==0)
   OSUsage(argv[0],OSArray,"--PARMMatrix option is mandatory");
  if (BEMMatrixFile==0)
   OSUsage(argv[0],OSArray,"--BEMMatrix option is mandatory");

  RWGGeometry *G = new RWGGeometry(GeoFile);
  HMatrix *PARMMatrix = new HMatrix(PARMMatrixFile, LHM_TEXT);
  ShowPARMMatirx(PARMMatrix); 
  HMatrix *M = new HMatrix(BEMMatrixFile, LHM_TEXT);
  //A->LUFactorize(); //(already done)
  HMatrix *Q = new HMatrix(PFTMatrixFile, LHM_TEXT);
  objective(G, PARMMatrix, M, Q); 

  delete G, PARMMatrix, M, Q; 
}
/***************************************************************/
double objective(unsigned n, const double *x, double *grad, void *my_func_data) 
{
  cdouble Omega = 10.471975511965978; 
  char *FileBase, *GeoFile, *HDF5File; 
  strncpy(FileBase,sprintf("N3_400nm_Mesh60nm")); 
  strncpy(GeoFile,sprintf("%s.scuffgeo",FileBase)); 
  strncpy(HDF5File,sprintf("%s.HDF5",FileBase)); 
  //I want the GeoFile and HDF5File to search for *.scuffgeo and *.HDF5 in the dir. 

  RWGGeometry *G = new RWGGeometry(GeoFile);
  
  // function to compute FOM.  
  double FOM; //=C^T Q C 
  x[0] = ;
  x[1] = ;
  
  
  return FOM; 
}//end objective 
