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
using namespace scuff;
//---------------------------------------------------------------//
void ShowPARMMatirx(HMatrix* PARMMatrix); 
void ShowPARMMatirx(int numL, HMatrix* PARMMatrix);
double objective(char* HDF5File, HMatrix *PARMMatrix,  cdouble Omega); 
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
  //  char *IFile=0;  // intensity is just read automatically ? 
  char *PARMMatrixFile=0;
  char *LogLevel=0;
  /* name    type   #args  max_instances  storage  count  description*/
  OptStruct OSArray[]=
   { 
     {"geometry",   PA_STRING, 1, 1, (void *)&GeoFile,   0,  ".scuffgeo file"},
     {"Omega",     PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals, &nOmegaVals,  "(angular) frequency"},
     {"PARMMatrix", PA_STRING, 1, 1, (void *)&PARMMatrixFile, 0, "VParameters file"},
     //     {"IFile",      PA_STRING, 1, 1, (void *)&IFile,     0, "name of intensity file"},
     {"HDF5File",   PA_STRING, 1, 1, (void *)&HDF5File,  0, "name of HDF5 file for OptTorque Matrix Data"},
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

  HVector *OmegaList=0;
  if (nOmegaVals>0) // process -- Omega options if present
    {
      OmegaList=new HVector(nOmegaVals, LHM_COMPLEX);
      for(int n=0; n<nOmegaVals; n++)
        OmegaList->SetEntry(n,OmegaVals[n]);
    }
  else 
    OSUsage(argv[0], OSArray, "you must specify at least one frequency");

  // RWGGeometry *G = new RWGGeometry(GeoFile);

  HMatrix *PARMMatrix = new HMatrix(PARMMatrixFile, LHM_TEXT);
  ShowPARMMatirx(PARMMatrix); 

  char OmegaStr[MAXSTR];
  char WvnmStr[MAXSTR]; 
  cdouble Omega; 
  double wvnm; 
  for(int nFreq=0; nFreq<OmegaList->N; nFreq++)
    { 
      Omega = OmegaList->GetEntry(nFreq); 
      objective(HDF5File, PARMMatrix, Omega); 
    }

  delete PARMMatrix; 
}
/***************************************************************/
double objective(char *HDF5File, HMatrix *PARMMatrix, cdouble Omega)
//add G
// Intensity also important. 
{
  z2s(Omega,OmegaStr); 
  wvnm = 2.0*M_PI*1000.0/real(Omega); 
  snprintf(WvnmStr,MAXSTR,"%i",int(wvnm)); 
  Log("Working at frequency %s...",OmegaStr);

  char IFilename[MAXSTR]; 
  snprintf(IFilename,MAXSTR,"%s.1freq.Intensity.dat",GeoFileBase);  
  fIntensity=fopen(IFilename,"a");

  fprintf(fIntensity,"%s    %e\n",WvnmStr,Intensity); 
  fclose(fIntensity); 
  

  // import matrices. 
  void *HDF5Context=0;
  HDF5Context=HMatrix::OpenHDF5Context(HDF5File);
  Log("Opened HDF5Context ...");
  HMatrix *M = new HMatrix(HDF5Context, LHM_HDF5, "M_%s",WvnmStr);
  if (M->ErrMsg) ErrExit(M->ErrMsg);
  int NR = M->NR; // number of rows 
  HMatrix *MLU = new HMatrix(HDF5Context, LHM_HDF5, "MLU_%s",WvnmStr);
  if (MLU->ErrMsg) ErrExit(MLU->ErrMsg);
  HMatrix *QabsOPFT = new HMatrix(HDF5Context, LHM_HDF5,"QabsOPFT_%s",WvnmStr);
  if (QabsOPFT->ErrMsg) ErrExit(QabsOPFT->ErrMsg);
  HMatrix *QFZOPFT = new HMatrix(HDF5Context, LHM_HDF5, "QFZOPFT_%s",WvnmStr);
  if (QFZOPFT->ErrMsg) ErrExit(QFZOPFT->ErrMsg);
  HMatrix *QTZOPFT = new HMatrix(HDF5Context, LHM_HDF5, "QTZOPFT_%s",WvnmStr);
  if (QTZOPFT->ErrMsg) ErrExit(QTZOPFT->ErrMsg);
  Log(" Successfully imported all matrices...\n");

  // from PARMMatrix, create IF and assemble RHS 
    

      Log("  Assembling the RHS vector..."); 
      G->AssembleRHSVector(Omega, SSD->kBloch, IFDList, KN);
        SSD->RHS->Copy(SSD->KN); // copy RHS vector for later 


  HVector *RHS  = new HVector(HDF5Context, LHM_HDF5, "RHS_%s",WvnmStr);
  if (RHS->ErrMsg) ErrExit(RHS->ErrMsg); 
  HVector *KN   = new HVector(HDF5Context, LHM_HDF5, "KN_%s",WvnmStr);
  if (KN->ErrMsg)  ErrExit(KN->ErrMsg);

  ///    RHS (sum aM+bN, GetRHSVector )
  double FOM =0.0; ///will be the output. 
  
  
  G->AssembleRHSVector(Omega, IF, RHS); //KN and RHS are formed
  KN->Copy(RHS); // this is the RHS vector assembled. 
  M->LUSolve(KN);// solved KN. 

  
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
  Cadj->ExportToHDF5(HDF5Context, "Cadj"); //adjoint current solution
  
  /*-------------------------------------------------------------*/
  /*- export the matrices into separate data files if asked  ----*/
  /*-------------------------------------------------------------*/
  char DatFile[100];
  snprintf(DatFile, MAXSTR, "%s.%s.FOM.dat",GeoFileName,WvnmStr);  
  f = fopen(DatFile,'a'); 
  THIS->ExportToText(MatFileName,"--separate,"); 
  fclose(f); 
       }
  HMatrix::CloseHDF5Context(HDF5Context); 
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
