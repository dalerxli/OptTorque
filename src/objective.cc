/*---------------------------------------------------------------
 * objective.cc  
 * compute Torque FOM
 * for a given BEM matrix, OPFT matrix, frequency, and PARMMatrix
 * created 2015.03
 * last updated on 2015.05.19
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
#define II cdouble(0.0,1.0)
#define MAXSTR 100 
#define MAXFREQ 1 
using namespace scuff;
//---------------------------------------------------------------//
void ShowPARMMatirx(HMatrix* PARMMatrix); 
void ShowPARMMatirx(int numL, HMatrix* PARMMatrix);
double GetIntegratedIntensity(RWGGeometry *G, int SurfaceIndex, HVector *RHSVector);
double objective(RWGGeometry *G, char *HDF5File, HMatrix *PARMMatrix, cdouble Omega);
//---------------------------------------------------------------//
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*- process options  -------------------------------------------*/
  /*--------------------------------------------------------------*/
  /// Input arguments should be filenames. (pointers) 
  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *GeoFile=0; 
  char *HDF5File=0; 
  char *PARMMatrixFile=0;
  char *LogLevel=0;
  /* name    type   #args  max_instances  storage  count  description*/
  OptStruct OSArray[]=
   { 
     {"geometry",   PA_STRING, 1, 1, (void *)&GeoFile,   0,  ".scuffgeo file"},
     {"Omega",     PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals, &nOmegaVals,  "(angular) frequency"},
     {"PARMMatrix", PA_STRING, 1, 1, (void *)&PARMMatrixFile, 0, "VParameters file"},
     {"HDF5File",   PA_STRING, 1, 1, (void *)&HDF5File,  0, "name of HDF5 file for Matrix Data"},
     {"LogLevel",   PA_STRING, 1, 1, (void *)&LogLevel,  0, "none | terse | verbose | verbose2"},
     //
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
    OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  if (PARMMatrixFile==0)
    OSUsage(argv[0],OSArray,"--PARMMatrix option is mandatory");
  if (HDF5File==0)
    OSUsage(argv[0],OSArray,"--HDF5File option is mandatory");

  RWGGeometry *G = new RWGGeometry(GeoFile);
  HMatrix *PARMMatrix = new HMatrix(PARMMatrixFile, LHM_TEXT);
  ShowPARMMatirx(PARMMatrix); 

  HVector *OmegaList=0;
  if (nOmegaVals==1) // process -- Omega options if present
    {
      OmegaList=new HVector(nOmegaVals, LHM_COMPLEX);
      for(int n=0; n<nOmegaVals; n++)
        OmegaList->SetEntry(n,OmegaVals[n]);
      // but currently only use 1 omega value 
    }
  else 
    OSUsage(argv[0], OSArray, "you must specify at least one frequency");

  cdouble Omega; 
  for(int nFreq=0; nFreq<OmegaList->N; nFreq++)
    { 
      Omega = OmegaList->GetEntry(nFreq); 
      objective(G, HDF5File, PARMMatrix, Omega); 
    }

  delete PARMMatrix; 
}
/***************************************************************/
double objective(RWGGeometry *G, char *HDF5File, HMatrix *PARMMatrix, cdouble Omega)
{
  cdouble dFOM =0.0; ///will be the output. 

  char OmegaStr[MAXSTR];
  char WvnmStr[MAXSTR]; 
  double wvnm; 
  z2s(Omega,OmegaStr); // set frequency. 
  wvnm = 2.0*M_PI*1000.0/real(Omega); 
  snprintf(WvnmStr,MAXSTR,"%i",int(wvnm)); 
  Log("Working at frequency %s...",OmegaStr);

  // import matrices. 
  //void *HDF5Context=0;
  //HDF5Context=HMatrix::OpenHDF5Context(HDF5File);
  //Log("Opened HDF5Context ...");
  //HMatrix *M = new HMatrix(HDF5File,LHM_HDF5,"M"); 
  //if (M->ErrMsg) ErrExit(M->ErrMsg);
  HMatrix *MLU = new HMatrix(HDF5File, LHM_HDF5, "MLU");
  if (MLU->ErrMsg) ErrExit(MLU->ErrMsg);
  int NR = MLU->NR; // number of rows 
  printf("NR = %i\n",NR); 
  HMatrix *QabsOPFT = new HMatrix(HDF5File, LHM_HDF5,"QabsOPFT");
  if (QabsOPFT->ErrMsg) ErrExit(QabsOPFT->ErrMsg);
  HMatrix *QFZOPFT = new HMatrix(HDF5File, LHM_HDF5, "QFZOPFT");
  if (QFZOPFT->ErrMsg) ErrExit(QFZOPFT->ErrMsg);
  HMatrix *QTZOPFT = new HMatrix(HDF5File, LHM_HDF5, "QTZOPFT");
  if (QTZOPFT->ErrMsg) ErrExit(QTZOPFT->ErrMsg);

  printf("size of QabsOPFT is : %i by %i\n",QabsOPFT->NR,QabsOPFT->NC);   

  // HVector *RHS  = new HVector(HDF5File, LHM_HDF5, "RHS");
  // if (RHS->ErrMsg) ErrExit(RHS->ErrMsg); 
  // HVector *KN   = new HVector(HDF5File, LHM_HDF5, "KN");
  // if (KN->ErrMsg)  ErrExit(KN->ErrMsg);
  Log(" Successfully imported all matrices...\n");

  // from PARMMatrix, create IF and assemble RHS

  char GeoFileBase[MAXSTR];
  strncpy(GeoFileBase,GetFileBase(G->GeoFileName),MAXSTR);  


  // Where should the normalizing step be ?? 
  Log("  Assembling the RHS vector...");
  IncField *IF=new VBeam(PARMMatrix);
  HVector* RHS=new HVector(NR,LHM_COMPLEX); 
  HVector* KN=new HVector(NR,LHM_COMPLEX); 
  RHS=G->AllocateRHSVector();
  G->AssembleRHSVector(Omega, IF, RHS); //KN and RHS are formed
  double Intensity= GetIntegratedIntensity(G, 0, RHS);
  KN->Copy(RHS);
  printf("size of RHS is : %i \n",RHS->N);   
  MLU->LUFactorize(); 
  MLU->LUSolve(KN);// solved KN. 
  HMatrix* RHSMAT = new HMatrix(NR,1,RHS->RealComplex, LHM_NORMAL, RHS->ZV); 
  printf("size of RHSMAT is : %i by %i\n", RHSMAT->NR, RHSMAT->NC);   
  
  // M*C_adj = transpose(QPFT)*conj(C)  
  printf("KN[1] = %e+%ei \n",real(KN->GetEntry(1)),imag(KN->GetEntry(1)));

  //this is not functioning.    !!
  HMatrix* KNMAT = new HMatrix(NR,1,KN->RealComplex, LHM_NORMAL, KN->ZV); 
  printf("KNMAT[1,1] = %e+%ei \n",real(KNMAT->GetEntry(1,1)),imag(KNMAT->GetEntry(1,1)));   

  KNMAT->Adjoint();
  KNMAT->Transpose(); //KNMAT = conj(KN); 

  HMatrix* Cadj = new HMatrix(NR,1,LHM_COMPLEX); 
  printf("size of KNMAT is : %i x %i \n",KNMAT->NR,KNMAT->NC);   
  printf("size of Cadj is : %i x %i \n",Cadj->NR,Cadj->NC);   
  QabsOPFT->Multiply(KNMAT,Cadj,"--transA C"); 
  printf("Cadj[1,1] = %e+%ei \n",real(Cadj->GetEntry(1,1)),imag(Cadj->GetEntry(1,1)));   

  MLU->LUSolve(Cadj);   //  printf(" Computed Cadj...\n");

  char MatFileName[100];
  snprintf(MatFileName, 100, "Mat_Cadj.dat");  
  Cadj->ExportToText(MatFileName,"--separate,"); 

  HMatrix* dFOMMAT = new HMatrix(1,1,LHM_COMPLEX);  

  //FOM
  Cadj->Transpose();
  Cadj->Multiply(RHSMAT,dFOMMAT);
  dFOM = dFOMMAT->GetEntry(1,1); 

  printf("dFOM = %e+%ei\n",real(dFOM),imag(dFOM));
  //HMatrix::CloseHDF5Context(HDF5Context); 
}//end main 

//--------------------------------------------------------------------//
//--------------------------------------------------------------------//
//--------------------------------------------------------------------//
void ShowPARMMatirx(HMatrix* PARMMatrix)
{
  int iL, L; 
  double aIn, aL, bL; 
  printf("PARMMatrix: \n"); 
  printf("L \t aIn \t\t aL \t\t bL\n"); 
  for(iL=0; iL<PARMMatrix->NR; iL++) //for each row of PARMMatrix
    {
      L   =PARMMatrix->GetEntryD(iL,0);
      aIn =PARMMatrix->GetEntryD(iL,1);
      aL  =PARMMatrix->GetEntryD(iL,2);
      bL  =PARMMatrix->GetEntryD(iL,3);
      printf("%d \t %f \t %f \t %f\n",L,aIn,aL,bL);
  }   
}
//--------------------------------------------------------------------//
void ShowPARMMatirx(int numL, HMatrix* PARMMatrix)
{//length of PARMMatrix specified 
  int iL, L; 
  double aIn, aL, bL; 
  printf("PARMMatrix: \n"); 
  printf("L \t aIn \t\t aL \t\t bL\n"); 
  for(iL=0; iL<numL; iL++) //for each row of PARMMatrix
    {
      L   =PARMMatrix->GetEntryD(iL,0);
      aIn =PARMMatrix->GetEntryD(iL,1);
      aL  =PARMMatrix->GetEntryD(iL,2);
      bL  =PARMMatrix->GetEntryD(iL,3);
      printf("%d \t %f \t %f \t %f\n",L,aIn,aL,bL);
  }   
}
//--------------------------------------------------------------------//
double GetIntegratedIntensity(RWGGeometry *G, int SurfaceIndex, HVector *RHSVector)
{
  RWGSurface *S = G->Surfaces[SurfaceIndex];
  bool IsPEC    = S->IsPEC;
  int BFOffset  = G->BFIndexOffset[SurfaceIndex];
  int NE        = S->NumEdges;
  int NBF       = IsPEC ? NE : 2*NE;

  /***************************************************************/
  /* construct and LU-factorize the overlap matrix               */
  /***************************************************************/
  HMatrix *M=new HMatrix(NBF, NBF, RHSVector->RealComplex);
  M->Zero();
  Log("GetIntegratedIntensity: Assembling S");
  for(int ne=0; ne<S->NumEdges; ne++)
   for(int nep=ne; nep<S->NumEdges; nep++)
    { 
      double OVLP=S->GetOverlap(ne, nep);
      if (IsPEC)
       { M->SetEntry( ne,  nep, OVLP);
         M->SetEntry( nep, ne,  OVLP);
       }
      else
       { 
         M->SetEntry( 2*ne, 2*nep, OVLP);
         M->SetEntry( 2*nep, 2*ne, OVLP);

         M->SetEntry( 2*ne+1, 2*nep+1, OVLP);
         M->SetEntry( 2*nep+1, 2*ne+1, OVLP);
       };
    };
  M->LUFactorize();

  HVector *V=new HVector(NBF, RHSVector->RealComplex);
  HVector *MInvV=new HVector(NBF, RHSVector->RealComplex);
  for(int nbf=0; nbf<NBF; nbf++)
   { V->SetEntry(nbf, RHSVector->GetEntry(BFOffset+nbf));
     MInvV->SetEntry(nbf, RHSVector->GetEntry(BFOffset+nbf));
   };
  M->LUSolve(MInvV);

  double Intensity=0.0;
  for(int nbf=0; nbf<NBF; nbf++)
   Intensity += real( conj(V->GetEntry(nbf)) * MInvV->GetEntry(nbf) );
  
  delete M;
  delete V;
  delete MInvV;

  return ZVAC*ZVAC*Intensity;
/********************************************************************/
/* return 0 if X lies outside the triangle with the given vertices, */
/* or a positive integer otherwise.                                 */
/********************************************************************/
}//end getintegratedintensities
