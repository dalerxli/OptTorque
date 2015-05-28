/*---------------------------------------------------------------
 * objective_struct.cc  
 * checo objective function that takes in struct OT_data 
 * Yoonkyung Eunnie Lee 
 * created on 2015.05.28
 * last updated on 2015.05.28
 *--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
//#include <nlopt.h>
#include <complex>
#include <cstdlib>
#include <libhrutil.h>
#include <libscuff.h>
#include "VBeam.h" 

#define II cdouble(0.0,1.0)
#define MAXSTR 100
#define NUML 5 

using namespace scuff;
// using namespace std;
//---------------------------------------------------------------//
// uses objective function and struct OT_data
double GetIntegratedIntensity(RWGGeometry *G, int SurfaceIndex, 
                              HVector *RHSVector);
double objective(unsigned n, double *x, double *grad, void *data);
//---------------------------------------------------------------//
typedef struct {
  cdouble Omega;//frequency (has to be real although in cdouble)
  char FileBase[MAXSTR];//base for geofile&HDF5file
  char QName[MAXSTR]; //selected name for Q (QabsOPFT, QTZOPFT, QFZOPFT)
} OT_data;
//---------------------------------------------------------------//
int main(int argc, char *argv[])
// main function to test objective
{
  unsigned n = NUML*5; 
  //---------------------------------------------------------------//
  // --- create a OT_data *data
  OT_data *data; 
  data->Omega=cdouble(10.471975511965978,0.0); 
  snprintf(data->FileBase,MAXSTR,"N3_400nm_Mesh60nm"); 
  snprintf(data->QName,MAXSTR,"QFZOPFT"); 

  double x[n]; // set it equal to some initial guess. 
  double grad[n]; 

  double FOM; /* the minimum objective value, upon return */
  for (int l = 0; l<NUML; l++)
    {
      x[l*5] = 10.0;  //alpha 
      x[l*5+1] = 1.0; //ar
      x[l*5+2] = 0.0; //br
      x[l*5+3] = 0.0; //ai
      x[l*5+4] = 1.0; //bi
    }

  FOM=objective(n, x, grad, data); 

  for (int l = 0; l<NUML; l++)
    {
      x[l*5] = 10.0;  //alpha 
      x[l*5+1] = 1.01; //ar
      x[l*5+2] = 0.0; //br
      x[l*5+3] = 0.0; //ai
      x[l*5+4] = 1.0; //bi
    }

  FOM=objective(n, x, grad, data); 
}
double objective(unsigned n, double *x, double *grad, void *data)
// objective function for derivative-free optimization 
{
  OT_data *d = (OT_data *)data; 
  cdouble Omega=d->Omega; 

  char GeoFile[MAXSTR], HDF5File[MAXSTR];
  snprintf(GeoFile,MAXSTR,"%s.scuffgeo",d->FileBase);  
  snprintf(HDF5File,MAXSTR,"%s.HDF5",d->FileBase ); 
  
  RWGGeometry *G = new RWGGeometry(GeoFile);
  printf(" Created RWGGeometry G...\n");
  HMatrix *MLU = new HMatrix(HDF5File, LHM_HDF5, "MLU");
  if (MLU->ErrMsg) ErrExit(MLU->ErrMsg);
  
  int NR = MLU->NR; // number of rows 
  printf("NR = %i\n",NR); 
  
  HMatrix* Q = new HMatrix(HDF5File, LHM_HDF5, d->QName);
  if (Q->ErrMsg) ErrExit(Q->ErrMsg);
  printf("size of Q is : %i by %i\n",Q->NR,Q->NC);   
  printf(" Successfully imported all matrices...\n");
  //---------------------------------------------------------------//
  // --- initialize RHS vectors and IF from double *x. 
  printf("  Assembling the RHS vector...");
  HVector* KN1=new HVector(NR,LHM_COMPLEX); 
  HMatrix* PMatrix=new HMatrix(NUML,6,LHM_REAL); 
  
  double cnorm[NUML]; 
  for (int l; l<NUML; l++){
    // for each l, compute cnorm[l]
    // IncField *IF1 = new VBeam(l,x[l*5],x[l*5+1],x[l*5+2],x[l*5+3],x[l*5+4]); 
    IncField *IF1 = new VBeam(l,x[l*5]);     // initialize with VBeam = 1*M; 
    G->AssembleRHSVector(Omega, IF1, KN1);   // KN1 is formed here
    cnorm[l]=GetIntegratedIntensity(G, 0, KN1); // get cnorm for each (l,aIn)
    printf("cnorm[%d]=%e\n",l,cnorm[l]);
    // create a big PMatrix with each row corresponding to one mode 
    // normalize ar,br,ai,bi
    PMatrix->SetEntry(l,0,l); 
    PMatrix->SetEntry(l,1,x[l*5]); 
    PMatrix->SetEntry(l,2,x[l*5+1]/cnorm[l]); 
    PMatrix->SetEntry(l,3,x[l*5+2]/cnorm[l]); 
    PMatrix->SetEntry(l,4,x[l*5+3]/cnorm[l]); 
    PMatrix->SetEntry(l,5,x[l*5+4]/cnorm[l]); 
    delete IF1; 
  }
  // --- allocate necessary vectors and matrices
  HVector* KN=new HVector(NR,LHM_COMPLEX); 
  HVector* RHS=new HVector(NR,LHM_COMPLEX); 

  HMatrix* TEMPMAT = new HMatrix(NR,1,LHM_COMPLEX); 
  HMatrix* FOM = new HMatrix(1,1,LHM_COMPLEX); 
  // --- create IncField from the normalized PMatrix 
  IncField* IF=new VBeam(PMatrix);
  // --- compute RHS, Solve for C=KN
  G->AssembleRHSVector(Omega, IF, KN); 
  RHS->Copy(KN);
  MLU->LUFactorize(); 
  MLU->LUSolve(KN);// solved KN. 
  HMatrix* C    = new HMatrix(NR,1,KN->RealComplex, LHM_NORMAL, KN->ZV); 
  // --- compute FOM = Cconj*Q*C , where C=KN and Q=QTZOPFT
  HMatrix* Cconj = new HMatrix(NR,1,LHM_COMPLEX); 
  Cconj->Copy(C); 
  Cconj->Adjoint(); // row vector for Adjoint(C)
  Q->Multiply(C,TEMPMAT); // TEMPMAT= Q*C
  Cconj->Multiply(TEMPMAT,FOM); // FOM= Cconj^T*Q*C 
  printf("Computed FOM = (%e)+(%e)i.\n",
         real(FOM->GetEntry(0,0)),imag(FOM->GetEntry(0,0)));
  // --- if gradient value is required, use adjoint method 
  // --- M*Cadj = transpose(QPFT)*conj(C)  
  if(grad) { /*NOT DOING ANYTNING NOW*/ }
  delete KN1, PMatrix, KN, RHS, TEMPMAT, FOM, IF, C, Q, Cconj; 
  return FOM->GetEntryD(0,0); 
}
