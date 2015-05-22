//-------------------------------------------------------------//
//-------------------------------------------------------------//
// elutil.cc 
// file containing utility functions for OptTorque
// using libraries in scuff-em 
// 2015.05.22
//-------------------------------------------------------------//
//-------------------------------------------------------------//
#include <stdio.h>
#include <math.h>
#include <complex>
#include <stdarg.h>
#include <fenv.h>

#include "OptTorque.h"

#define II cdouble(0.0,1.0)
#define MAXSTR   100
//-------------------------------------------------------------//
//-------------------------------------------------------------//
// These static declarations are called in VisualizeIncFields
static char *FieldFuncs=const_cast<char *>(
 "|Ex|,|Ey|,|Ez|,"
 "sqrt(|Ex|^2+|Ey|^2+|Ez|^2),"
 "|Hx|,|Hy|,|Hz|,"
 "sqrt(|Hx|^2+|Hy|^2+|Hz|^2)");
static const char *FieldTitles[]=
 {"|Ex|", "|Ey|", "|Ez|", "|E|",
  "|Hx|", "|Hy|", "|Hz|", "|H|",
 };
#define NUMFIELDFUNCS 8// used in VisualizeFields
//-------------------------------------------------------------//
//-------------------------------------------------------------//
void VisualizeIncField(RWGGeometry *G, IncField *IF, cdouble Omega, 
                       char *MeshFile);
double GetIntegratedIntensity(RWGGeometry *G, int SurfaceIndex, 
                              HVector *RHSVector);
double GetAvgIntensity(RWGGeometry *G, int SurfaceIndex, 
                       HVector *RHSVector);
//double **AllocateByEdgeArray(RWGGeometry *G, int ns);
//void ProcessByEdgeArray(RWGGeometry *G, int ns, cdouble Omega,
//                        double **ByEdge);
//-------------------------------------------------------------//
double GetIntegratedIntensity(RWGGeometry *G, int SurfaceIndex, 
                              HVector *RHSVector)
{
  // function that returns the integrated |E|^2 on the given surface
  // Homer Reid 
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

//-------------------------------------------------------------//
//-------------------------------------------------------------//
void VisualizeIncField(RWGGeometry *G, IncField *IF, cdouble Omega, 
                         char *MeshFile)
{
  // function to visualize a given incfield on a mesh surface
  // without the scattered field 
  char GeoFileBase[100], PPFileName[100];
  strncpy(GeoFileBase,GetFileBase(G->GeoFileName),100);
  snprintf(PPFileName,100,"%s.%s.FV.Incfield.pp",GeoFileBase,z2s(Omega));
  FILE *f=fopen(PPFileName,"a");
  if (!f)
   ErrExit("could not open field visualization file %s",PPFileName);

  scuff::RWGSurface *S=new scuff::RWGSurface(MeshFile);
  printf("Creating flux plot for surface %s...\n",MeshFile);
  //    create an Nx3 XMatrix where each row is x,y,z of FVMesh. 
  HMatrix *XMatrix=new HMatrix(S->NumVertices, 3);
  printf("NumVertices of S is:%d \n ",S->NumVertices);
  for(int nv=0; nv<S->NumVertices; nv++)
    { 
      XMatrix->SetEntry(nv, 0, S->Vertices[3*nv + 0]);
      XMatrix->SetEntry(nv, 1, S->Vertices[3*nv + 1]);
      XMatrix->SetEntry(nv, 2, S->Vertices[3*nv + 2]);
    };
  HMatrix *FMatrix=new HMatrix(XMatrix->NR, NUMFIELDFUNCS, LHM_COMPLEX);
  int ii,jj; 
  double X[3];
  cdouble EH[6]; 
  
  for (jj=0;jj<XMatrix->NR;jj++){
    X[0]=XMatrix->HMatrix::GetEntryD(jj,0); 
    X[1]=XMatrix->HMatrix::GetEntryD(jj,1); 
    X[2]=XMatrix->HMatrix::GetEntryD(jj,2); 
    //     get the incident fields at the panel vertices
    FMatrix=G->GetFields(IF, 0, Omega, 0, XMatrix, 0, FieldFuncs); 
  }
  for(int nff=0; nff<NUMFIELDFUNCS; nff++)
   {
     fprintf(f,"View \"%s(%s)\" {\n",FieldTitles[nff],z2s(Omega));
     for(int np=0; np<S->NumPanels; np++)
      {
        RWGPanel *P=S->Panels[np];
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
     fprintf(f,"};\n\n");
   };
  fclose(f);
  delete FMatrix;
  delete XMatrix;
  delete S;
}
//-------------------------------------------------------------//
//-------------------------------------------------------------//
double **AllocateByEdgeArray(RWGGeometry *G, int ns)
{
  int NE = G->Surfaces[ns]->NumEdges;
  double **ByEdge=(double **)mallocEC(7*sizeof(double *));
  ByEdge[0]=(double *)mallocEC(7*NE*sizeof(double));
  for(int nq=1; nq<7; nq++)
   ByEdge[nq] = ByEdge[nq-1] + NE;

  return ByEdge;
}
//-------------------------------------------------------------//
//-------------------------------------------------------------//
void ProcessByEdgeArray(RWGGeometry *G, int ns, cdouble Omega,
                        double **ByEdge)
{
  static const char *PFTNames[7]
   ={"PAbs","FX","FY","FZ","TX","TY","TZ"};

  char FileName[100];
  snprintf(FileName,100,"%s.pp",GetFileBase(G->GeoFileName));

  for(int nq=0; nq<7; nq++)
   { char Tag[20];
     snprintf(Tag,20,"%s(%s)",PFTNames[nq],z2s(Omega));
     G->Surfaces[ns]->PlotScalarDensity(ByEdge[nq],true,FileName,Tag);
   };
 // free(ByEdge[0]);
 // free(ByEdge);
}
//-------------------------------------------------------------//
//-------------------------------------------------------------//
