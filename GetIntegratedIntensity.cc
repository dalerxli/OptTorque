/*
 *
 */

#include <stdio.h>
#include <math.h>
#include <complex>

#include "libhrutil.h"
#include "libscuff.h"

#define II cdouble(0.0,1.0)

using namespace scuff;

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

}

/********************************************************************/
/* return 0 if X lies outside the triangle with the given vertices, */
/* or a positive integer otherwise.                                 */
