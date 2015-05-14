/*
 * VBeamTestMesh.cc   -- Tests VBeam Result for a given X[3] and Omega 
 *
 */
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <complex>
#include <cstdlib> //for abs
#include <string.h>//for memcpy
#include "VBeam.h"
#include <libscuff.h> //to declare RWGSurface 
  #ifndef cdouble
    typedef std::complex<double> cdouble;
  #endif
  #define Eps         cdouble(1.0,0.0)
  #define Mu          cdouble(1.0,0.0)
//  #define Omega       cdouble(5.0,0.0)
  #define ZVAC        376.73031346177
  #define II cdouble(0.0,1.0)
//--------------------------------------------------------------------//
static char *FieldFuncs=const_cast<char *>(
 "|Ex|,|Ey|,|Ez|,"
 "sqrt(|Ex|^2+|Ey|^2+|Ez|^2),"
 "|Hx|,|Hy|,|Hz|,"
 "sqrt(|Hx|^2+|Hy|^2+|Hz|^2)");

static const char *FieldTitles[]=
 {"|Ex|", "|Ey|", "|Ez|", "|E|",
  "|Hx|", "|Hy|", "|Hz|", "|H|",
 };

#define NUMFIELDFUNCS 8
//--------------------------------------------------------------------//
int main(int argc, char *argv[]){
  //    Step 1: import VParameters and construct an IF. 
  //    Step 2: import a mesh file and construct RWGSurface S. 
  //    Step 3: 

  //--------------------------------------------------------------//
  //- process options  -------------------------------------------//
  //--------------------------------------------------------------//
  cdouble OmegaVals[1];        int nOmegaVals;
  char *MeshFile=0; 
  OptStruct OSArray[]=
   { /* name    type   #args  max_instances  storage  count  description*/
     {"Omega",      PA_CDOUBLE, 1, 1, (void *)OmegaVals, &nOmegaVals,  "(angular) frequency"},
     {"FVMesh",    PA_STRING,  1, 1,(void *)&MeshFile, 0,  "field visualization mesh"},
     {0,0,0,0,0,0,0} };
  ProcessOptions(argc, argv, OSArray);

  if (MeshFile==0)
   OSUsage(argv[0],OSArray,"--FVMesh option is mandatory");
  if (nOmegaVals==0)
    OSUsage(argv[0], OSArray, "you must specify at least one frequency");
  //--------------------------------------------------------------//
  char PARMMatrixFile[100],PPFile[100];
  snprintf(PARMMatrixFile,100,"VParameters");  
  HMatrix *PARMMatrix = new HMatrix(PARMMatrixFile, LHM_TEXT);
  snprintf(PPFile,100,"%s.pp",MeshFile);
  //    try to open output file 
  FILE *f=fopen(PPFile,"a");
  if (!f) 
   ErrExit("could not open field visualization file %s",PPFile);
  //    Step 1: import VParameters and construct a VB and IF. 
  cdouble Omega= OmegaVals[1]; 
  VBeam *VB = new VBeam(PARMMatrix); 
  IncField *IF = VB; 
  //    Step 2: import a mesh file and construct RWGSurface S. 
  scuff::RWGSurface *S=new scuff::RWGSurface(MeshFile);
  printf("Creating flux plot for surface %s...\n",MeshFile);
  //    create an Nx3 HMatrix whose columns are the coordinates of 
  //    the flux mesh panel vertices
  HMatrix *XMatrix=new HMatrix(S->NumVertices, 3);
printf("NumVertices of S is:%d \n ",S->NumVertices);
  for(int nv=0; nv<S->NumVertices; nv++)
   { 
     XMatrix->SetEntry(nv, 0, S->Vertices[3*nv + 0]);
     XMatrix->SetEntry(nv, 1, S->Vertices[3*nv + 1]);
     XMatrix->SetEntry(nv, 2, S->Vertices[3*nv + 2]);
   };
  //     get the incident fields at the panel vertices
  int NumFuncs = 6;
  HMatrix *FMatrix=new HMatrix(XMatrix->NR, NumFuncs, LHM_COMPLEX);
  int ii,jj; 
  double X[3];
  cdouble EH[6]; 

  for (jj=0;jj<XMatrix->NR;jj++){
    X[0]=XMatrix->HMatrix::GetEntryD(jj,0); 
    X[1]=XMatrix->HMatrix::GetEntryD(jj,1); 
    X[2]=XMatrix->HMatrix::GetEntryD(jj,2);
    VB->VBeam::GetFields(X,EH); 
    for (ii=0;ii<NumFuncs;ii++){
      FMatrix->SetEntry(jj,ii,EH[ii]);
    }
  }
  /*--------------------------------------------------------------*/
  for(int nff=0; nff<NUMFIELDFUNCS; nff++)
   { 
     fprintf(f,"View \"%s(%s)\" {\n",FieldTitles[nff],z2s(Omega));
     /*--------------------------------------------------------------*/
     for(int np=0; np<S->NumPanels; np++)
      {
        scuff::RWGPanel *P=S->Panels[np];
        int iV1 = P->VI[0];  double *V1 = S->Vertices + 3*iV1;
        int iV2 = P->VI[1];  double *V2 = S->Vertices + 3*iV2;
        int iV3 = P->VI[2];  double *V3 = S->Vertices + 3*iV3;
        fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                   V1[0], V1[1], V1[2],
                   V2[0], V2[1], V2[2],
                   V3[0], V3[1], V3[2],
                   FMatrix->GetEntryD(iV1,nff),
                   FMatrix->GetEntryD(iV2,nff),
                   FMatrix->GetEntryD(iV3,nff));
      };
     /*--------------------------------------------------------------*/
     fprintf(f,"};\n\n");
   };
  fclose(f);

  delete FMatrix;
  delete XMatrix;

  delete S;
  delete VB; 
  delete PARMMatrix;  
}
