/*---------------------------------------------------------------
 * GetQ.cc  
 * Compute and store M and OPFT Q matrix 
 * for a given frequency
 * created      2015.05.09
 * last updated 2015.05.19
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
namespace scuff{
  //#define NUMPFT 8
  void GetOPFTMatrices(RWGGeometry *G, int SurfaceIndex, cdouble Omega,
		  HMatrix *QPFT[NUMPFT], bool NeedMatrix[NUMPFT]);
}
double GetQ(RWGGeometry *G, cdouble Omega); 
//---------------------------------------------------------------//
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

  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile=0;

  char *PARMMatrixFile=0;
  char *BEMMatrixFile=0;
  char *PFTMatrixFile=0;
  char *LogLevel=0;

  /* name    type   #args  max_instances  storage  count  description*/
  OptStruct OSArray[]=
   { 
     {"geometry",   PA_STRING,  1, 1, (void *)&GeoFile,   0,  ".scuffgeo file"},
     {"Omega",     PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals, &nOmegaVals,  "(angular) frequency"},
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

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run calculations                        */
  /*******************************************************************/
  HVector *OmegaList1=0, *OmegaList2=0, *OmegaList=0;
  if (OmegaFile) // process --OmegaFile option if present
   { 
     OmegaList1=new HVector(OmegaFile,LHM_TEXT);
     if (OmegaList1->ErrMsg)
      ErrExit(OmegaList1->ErrMsg);
   }
  if (nOmegaVals>0) // process -- Omega options if present
   {
     OmegaList2=new HVector(nOmegaVals, LHM_COMPLEX);
     for(int n=0; n<nOmegaVals; n++)
      OmegaList2->SetEntry(n,OmegaVals[n]);
   }
  if (  OmegaList1 && !OmegaList2 )
   OmegaList=OmegaList1;
  else if ( !OmegaList1 && OmegaList2 )
   OmegaList=OmegaList2;
  else if (  OmegaList1 && OmegaList2  )
   OmegaList=Concat(OmegaList1, OmegaList2);
  else 
   OSUsage(argv[0], OSArray, "you must specify at least one frequency");


  RWGGeometry *G = new RWGGeometry(GeoFile);
  HMatrix *M = new HMatrix(BEMMatrixFile, LHM_TEXT);

  delete G, M, Q; 
}

//this function doesn't seem necessary. you should include it in OptTorque.
HMatrix** GetQ(RWGGeometry *G, IncField *IF, cdouble Omega)
{ 
  // store Q matrices in the HDF5 file opened. 
  printf("GetQ FUNCTION IS CALLED.\n"); 

  //---------------------------------------------------------------//
  int numPFT = 0; // number of PFT matrices returned;  
  PFTOptions *MyPFTOptions=InitPFTOptions();
  HMatrix *QPFT[8]={0,0,0,0,0,0,0,0};
  printf(" Getting PFT matrix Q...\n");
  bool NeedMatrix[8]={false, false, false, false, 
                      false, false, false, false}   
  NeedMatrix[SCUFF_PABS]=true; //which matrices are needed. 
  NeedMatrix[SCUFF_PSCA]=true;
  //NeedMatrix[SCUFF_XFORCE]=true;
  NeedMatrix[SCUFF_ZFORCE]=true;
  NeedMatrix[SCUFF_ZTORQUE]=true;
  numPFT = 4; 
  GetOPFTMatrices(G, 0, Omega, QPFT, NeedMatrix);
  //---------------------------------------------------------------//
  // Store QPFT in HDF5
  return QPFT; 
}//end GetQ
