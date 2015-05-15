/*
/* MatrixFOM.cc  
/* updates for v7 
/*
 */

#include <stdio.h>
#include <math.h>
#include <complex>
#include <fstream>
#include <iostream>
#include <cstdlib>

#include "libhrutil.h"
#include "libscuff.h"
//#include "GHBeam.h"
//#include "OptTorque.h"

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

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main()
{
  /***************************************************************/
  // Import Mat_wvnm.dat 
  // Import the Q matrix. 
  // then compute C_adj by doing LUSolve(M, QbarCbar) 
  /***************************************************************/
  //  char *OmegaFileName=0;
  //  char *FileBase=0;
  bool MatrixOut=true; //

  char HFileInName[100]; 
  snprintf(HFileInName,100,"MAT_599nm.HDF5");

  char HFileOutName[100];
  snprintf(HFileOutName,100,"FOM_599nm.HDF5"); 
  //Read M, RHS, KN from HDF5 input file. 

  cdouble Omega= 10.471975511965978*(1.0+II);
  double wvnm= 2.0*M_PI*1000.0/std::real(Omega);

  RWGGeometry *G = new RWGGeometry("Tricyl_1092.scuffgeo");
  //  RWGGeometry *G = new RWGGeometry("RoundedPoly_N3_5.scuffgeo");

  //  HMatrix *M_orig = new HMatrix(HFileInName, LHM_HDF5, "M_orig");
  // HMatrix *M_lu = new HMatrix(HFileInName, LHM_HDF5, "M_lu");
  HMatrix *M_orig = new HMatrix(HFileInName, LHM_HDF5, "M");
  HMatrix *M_lu = new HMatrix(HFileInName, LHM_HDF5, "M");
  M_lu->LUFactorize(); 
  HVector *RHS  = new HVector(HFileInName, LHM_HDF5, "RHS");
  HVector *KN   = new HVector(HFileInName, LHM_HDF5, "KN");
       if (M_orig->ErrMsg) ErrExit(M_orig->ErrMsg);
       if (M_lu->ErrMsg) ErrExit(M_lu->ErrMsg);
       if (RHS->ErrMsg) ErrExit(RHS->ErrMsg);
       if (KN->ErrMsg)  ErrExit(KN->ErrMsg);
  printf(" Successfully imported M,RHS,KN...\n");
  int NR = M_orig->NR;
  //  printf("NR = %i",NR); 
  void *HDF5Context = HMatrix::OpenHDF5Context(HFileOutName);
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  PFTOptions *MyPFTOptions=InitPFTOptions();
  HMatrix *QPFT[8]={0,0,0,0,0,0,0,0};
  printf(" Getting PFT matrix Q...\n");
  bool NeedMatrix[8]={false, false, false, false, false, false, false, false};

  NeedMatrix[SCUFF_PABS]=true; //which matrices are needed. 
  //NeedMatrix[SCUFF_PSCA]=true;
  NeedMatrix[SCUFF_XFORCE]=true;
  //NeedMatrix[SCUFF_ZFORCE]=true;
  NeedMatrix[SCUFF_ZTORQUE]=true;

  GetOPFTMatrices(G, 0, Omega, QPFT, NeedMatrix);
	 
  QPFT[SCUFF_PABS]->ExportToHDF5(HDF5Context, "QabsOPFT");
  QPFT[SCUFF_XFORCE]->ExportToHDF5(HDF5Context, "QFXOPFT");
  QPFT[SCUFF_ZTORQUE]->ExportToHDF5(HDF5Context, "QTZOPFT");
  printf(" Exported QPFT files to HDF5 format...\n");
  // first calculate fom for absorbed power. 
  // M*C_adj = transpose(QPFT)*conj(C)  
  HMatrix *Qt = new HMatrix(QPFT[SCUFF_PABS]);
  Qt->Transpose(); 
  printf(" Computed Qt...\n");
  HVector *Cconj  = new HVector(NR,LHM_COMPLEX);
  HVector *Cadj = new HVector(NR,LHM_COMPLEX);

  cdouble TEMPVAL = (0.0,0.0);
  cdouble TEMPVAL2 = (0.0,0.0);

  for(int ii=0;ii<NR;ii++)
    {
      TEMPVAL = KN->GetEntry(ii);
      TEMPVAL = std::conj(TEMPVAL); 
      Cconj->SetEntry(ii,TEMPVAL);

      TEMPVAL2 = (0.0,0.0); 
      for(int jj=0;jj<NR;jj++)
        {
          TEMPVAL2 += Qt->GetEntry(ii,jj)*TEMPVAL;
          Cadj->SetEntry(ii,TEMPVAL2);
        }
    }
  //  printf(" Computed Cconj...\n");
  Cconj->ExportToHDF5(HDF5Context, "Cconj"); 
  Qt->ExportToHDF5(HDF5Context, "Qt"); 
  printf("size of Qt is : %i by %i\n",Qt->NR,Qt->NC);   
  printf("size of Cconj is : %i\n",Cconj->N);   
  
  //  Qt->Multiply(Cconj,Cadj); // Cadj =Q t*Cconj
  M_lu->LUSolve(Cadj);
  //printf(" Computed Cadj...\n");
  Cadj->ExportToHDF5(HDF5Context, "Cadj"); //adjoint current solution
     /*-------------------------------------------------------------*/
     /*- export the matrices into separate data files if asked  ----*/
     /*-------------------------------------------------------------*/
     if(MatrixOut)
       {
         char MatFileName[100];
        // store int(lambda[nm]) rather than omega in the filename.
        snprintf(MatFileName, 100, "Mat_QabsOPFT_%inm.dat",int(wvnm));  
        QPFT[SCUFF_PABS]->ExportToText(MatFileName,"--separate,"); 
        snprintf(MatFileName, 100,"Mat_QFXOPFT_%inm.dat",int(wvnm));
        QPFT[SCUFF_XFORCE]->ExportToText(MatFileName,"--separate,"); 
        snprintf(MatFileName, 100,"Mat_QTZOPFT_%inm.dat",int(wvnm));
        QPFT[SCUFF_ZTORQUE]->ExportToText(MatFileName,"--separate,");
        snprintf(MatFileName, 100,"Mat_Cadj_%inm.dat",int(wvnm));
        Cconj->ExportToText(MatFileName,"--separate,");
       }
  HMatrix::CloseHDF5Context(HDF5Context); 
}//end main 


/*
	 void MyOtherFunction(...)
	 {
	   HMatrix *QPAbs = new HMatrix("QPAbs.HDF5");

	   HMatrix *M3 = ...

	     // if QPAbs is MxN and M3 is NxK

	     HMatrix *M7 = new HMatrix(M,K,LHM_COMPLEX);

	   QPAbs->Multiply(M3, M7); // M7 <= QPAbs * M3;
	 }

*/ 




