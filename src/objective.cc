/*---------------------------------------------------------------
 * objective.cc  
 * compute Torque FOM
 * for a given BEM matrix, OPFT matrix, frequency, and PARMMatrix
 * Yoonkyung Eunnie Lee 
 * created on 2015.03.18
 * last updated on 2015.05.27
 *--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <complex>
#include <cstdlib>
#include <libhrutil.h>
#include <libscuff.h>
#include <nlopt.h> 
#include "VBeam.h" 

#define II cdouble(0.0,1.0)
#define NUML 5
#define MAXSTR 100 
#define FileBase "N3_400nm_Mesh60nm" 
#define Omega cdouble(10.471975511965978,0.0)

using namespace scuff;

//---------------------------------------------------------------//
double GetIntegratedIntensity(RWGGeometry *G, int SurfaceIndex, 
                              HVector *RHSVector);
void ShowPARMMatrix(HMatrix* PARMMatrix); 
double objective_Qabs(unsigned n, double *x, double *grad);
double objective_QTZ(unsigned n, double *x, double *grad);
double objective_QFZ(unsigned n, double *x, double *grad);
//---------------------------------------------------------------//
int main()
// main function to test objective_Q
{
  // create a vector x to test. 
  double x[NUML*5]; 
  double grad[NUML*5]; 
  for (int l = 0; l<NUML; l++)
    {
      x[l*5] = 10.0;  //alpha 
      x[l*5+1] = 1.0; //ar
      x[l*5+2] = 0.0; //br
      x[l*5+3] = 0.0; //ai
      x[l*5+4] = 0.0; //bi
    }

  objective_QFZ(NUML*5, x, grad); 

  for (int l = 0; l<NUML; l++)
    {
      x[l*5] = 10.0;  //alpha 
      x[l*5+1] = 1.1; //ar
      x[l*5+2] = 0.0; //br
      x[l*5+3] = 0.0; //ai
      x[l*5+4] = 0.0; //bi
    }

  objective_QFZ(NUML*5, x, grad); 

      // char MatFileName[100];
      // snprintf(MatFileName, 100, "Mat_Cadj.dat");  
      // Cadj->ExportToText(MatFileName,"--separate,"); 

}
//---------------------------------------------------------------//

double objective_Qabs(unsigned n, double *x, double *grad) 
// objective function for power absorption 
{
  // cdouble Omega = 10.471975511965978; 
  //---------------------------------------------------------------//
  // char FileBase[MAXSTR]; 
  char GeoFile[MAXSTR], HDF5File[MAXSTR];
  //snprintf(FileBase,MAXSTR,"N3_400nm_Mesh60nm"); 
  snprintf(GeoFile,MAXSTR,"%s.scuffgeo",FileBase);  
  snprintf(HDF5File,MAXSTR,"%s.HDF5",FileBase ); 
  // I want the GeoFile and HDF5File to search for *.scuffgeo and *.HDF5 in the dir. 
  RWGGeometry *G = new RWGGeometry(GeoFile);
  Log(" Created RWGGeometry G...\n");
  //---------------------------------------------------------------//
  HMatrix *MLU = new HMatrix(HDF5File, LHM_HDF5, "MLU");
  if (MLU->ErrMsg) ErrExit(MLU->ErrMsg);
  int NR = MLU->NR; // number of rows 
  printf("NR = %i\n",NR); 
  HMatrix *QabsOPFT = new HMatrix(HDF5File, LHM_HDF5,"QabsOPFT");
  if (QabsOPFT->ErrMsg) ErrExit(QabsOPFT->ErrMsg);
  // HMatrix *QFZOPFT = new HMatrix(HDF5File, LHM_HDF5, "QFZOPFT");
  // if (QFZOPFT->ErrMsg) ErrExit(QFZOPFT->ErrMsg);
  // HMatrix *QTZOPFT = new HMatrix(HDF5File, LHM_HDF5, "QTZOPFT");
  // if (QTZOPFT->ErrMsg) ErrExit(QTZOPFT->ErrMsg);
  printf("size of Q is : %i by %i\n",QabsOPFT->NR,QabsOPFT->NC);   
  Log(" Successfully imported all matrices...\n");
  //---------------------------------------------------------------//

  // --- initialize RHS vectors and IF from double *x. 
  Log("  Assembling the RHS vector...");
  HVector* KN1=new HVector(NR,LHM_COMPLEX); 
  HMatrix* PMatrix=new HMatrix(NUML,6,LHM_REAL); 
  double cnorm[NUML]; 
  // KN1=G->AllocateRHSVector();
  for (int l; l<NUML; l++){
    // for each l, compute cnorm[l]
    // IncField *IF1 = new VBeam(l,x[l*5],x[l*5+1],x[l*5+2],x[l*5+3],x[l*5+4]); 
    IncField *IF1 = new VBeam(l,x[l*5]);     // initialize with VBeam = 1*M; 
    G->AssembleRHSVector(Omega, IF1, KN1);   // KN1 is formed here
    cnorm[l]=GetIntegratedIntensity(G, 0, KN1); // get cnorm for each (l,aIn)

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
  // --- choose Q 
  HMatrix* Q = new HMatrix(NR,NR,LHM_COMPLEX);
  //  Q->Copy(QTZOPFT); // Choose as torque matrix  
  Q->Copy(QabsOPFT); // Choose as torque matrix 

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
  if(grad)
    { 
      // --- to compute the gradient, we need F_lm which is RHS for a single coeff

      // --- compute adjoint current Cadj
      HMatrix* Cadj = new HMatrix(NR,1,LHM_COMPLEX);// adjoint current
      Cconj->Transpose();                           // column vector for Conj(C)
      Q->Multiply(Cconj,Cadj,"--transA T");           // Cadj = transpose(Q)*conj(C)
      MLU->LUSolve(Cadj);                           // Cadj = M\(transpose(Q)*conj(C))
      printf("Computed Cadj[0,0] = %e+%ei \n",
             real(Cadj->GetEntry(0,0)),imag(Cadj->GetEntry(0,0)));   

      // --- the gradient does not exist for angle aIn. 
      for(int l = 0; l<NUML; l++){
        grad[l*5]=0.0; 
      }
      // --- the gradient exists for ar,br,ai,bi for a fixed aIn and L. 
      for(int l = 0; l<NUML; l++){
        //individual RHS is needed
        HMatrix* RHSMAT = new HMatrix(NR,1,RHS->RealComplex, LHM_NORMAL, RHS->ZV); 
        HMatrix* dFOM = new HMatrix(1,1,LHM_COMPLEX);  // gradient of FOM
        Cadj->Multiply(RHSMAT,dFOM,"--transA T"); // compute dFOM 
        grad[l*5+1]=dFOM->GetEntryD(0,0); // ar
        grad[l*5+2]=0.0; // br
        grad[l*5+3]=0.0; // ai
        grad[l*5+4]=0.0; // bi
        delete RHSMAT, dFOM; 
      }    
      delete Cadj; 
    }
  delete KN1, PMatrix, KN, RHS, TEMPMAT, FOM, IF, C, Q, Cconj; 
  return FOM->GetEntryD(0,0); 
}//end objective 

//---------------------------------------------------------------//
double objective_QTZ(unsigned n, double *x, double *grad) 
{
  // cdouble Omega = 10.471975511965978; 
  //---------------------------------------------------------------//
  // char FileBase[MAXSTR]; 
  char GeoFile[MAXSTR], HDF5File[MAXSTR];
  //snprintf(FileBase,MAXSTR,"N3_400nm_Mesh60nm"); 
  snprintf(GeoFile,MAXSTR,"%s.scuffgeo",FileBase);  
  snprintf(HDF5File,MAXSTR,"%s.HDF5",FileBase ); 
  // I want the GeoFile and HDF5File to search for *.scuffgeo and *.HDF5 in the dir. 
  RWGGeometry *G = new RWGGeometry(GeoFile);
  Log(" Created RWGGeometry G...\n");
  //---------------------------------------------------------------//
  HMatrix *MLU = new HMatrix(HDF5File, LHM_HDF5, "MLU");
  if (MLU->ErrMsg) ErrExit(MLU->ErrMsg);
  int NR = MLU->NR; // number of rows 
  printf("NR = %i\n",NR); 
  // HMatrix *QabsOPFT = new HMatrix(HDF5File, LHM_HDF5,"QabsOPFT");
  // if (QabsOPFT->ErrMsg) ErrExit(QabsOPFT->ErrMsg);
  // HMatrix *QFZOPFT = new HMatrix(HDF5File, LHM_HDF5, "QFZOPFT");
  // if (QFZOPFT->ErrMsg) ErrExit(QFZOPFT->ErrMsg);
  HMatrix *QTZOPFT = new HMatrix(HDF5File, LHM_HDF5, "QTZOPFT");
  if (QTZOPFT->ErrMsg) ErrExit(QTZOPFT->ErrMsg);
  printf("size of Q is : %i by %i\n",QTZOPFT->NR,QTZOPFT->NC);   
  Log(" Successfully imported all matrices...\n");
  //---------------------------------------------------------------//
 

  // --- initialize RHS vectors and IF from double *x. 
  Log("  Assembling the RHS vector...");
  HVector* KN1=new HVector(NR,LHM_COMPLEX); 
  HMatrix* PMatrix=new HMatrix(NUML,6,LHM_REAL); 
  double cnorm[NUML]; 
  // KN1=G->AllocateRHSVector();
  for (int l; l<NUML; l++){
    // for each l, compute cnorm[l]
    // IncField *IF1 = new VBeam(l,x[l*5],x[l*5+1],x[l*5+2],x[l*5+3],x[l*5+4]); 
    IncField *IF1 = new VBeam(l,x[l*5]);     // initialize with VBeam = 1*M; 
    G->AssembleRHSVector(Omega, IF1, KN1);   // KN1 is formed here
    cnorm[l]=GetIntegratedIntensity(G, 0, KN1); // get cnorm for each (l,aIn)

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
  //ShowPARMMatrix(PMatrix);

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
  // --- choose Q 
  HMatrix* Q = new HMatrix(NR,NR,LHM_COMPLEX);
  Q->Copy(QTZOPFT); // Choose as torque matrix  
  //Q->Copy(QabsOPFT); // Choose as torque matrix 

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

  delete KN1, PMatrix, KN, RHS, TEMPMAT, FOM, IF, C, Q, Cconj; 
  return FOM->GetEntryD(0,0); 
}//end objective 

//---------------------------------------------------------------//
double objective_QFZ(unsigned n, double *x, double *grad) 
{
  // cdouble Omega = 10.471975511965978; 
  //---------------------------------------------------------------//
  // char FileBase[MAXSTR]; 
  char GeoFile[MAXSTR], HDF5File[MAXSTR];
  //snprintf(FileBase,MAXSTR,"N3_400nm_Mesh60nm"); 
  snprintf(GeoFile,MAXSTR,"%s.scuffgeo",FileBase);  
  snprintf(HDF5File,MAXSTR,"%s.HDF5",FileBase ); 
  // I want the GeoFile and HDF5File to search for *.scuffgeo and *.HDF5 in the dir. 
  RWGGeometry *G = new RWGGeometry(GeoFile);
  Log(" Created RWGGeometry G...\n");
  //---------------------------------------------------------------//
  HMatrix *MLU = new HMatrix(HDF5File, LHM_HDF5, "MLU");
  if (MLU->ErrMsg) ErrExit(MLU->ErrMsg);
  int NR = MLU->NR; // number of rows 
  printf("NR = %i\n",NR); 
  // HMatrix *QabsOPFT = new HMatrix(HDF5File, LHM_HDF5,"QabsOPFT");
  // if (QabsOPFT->ErrMsg) ErrExit(QabsOPFT->ErrMsg);
  HMatrix *QFZOPFT = new HMatrix(HDF5File, LHM_HDF5, "QFZOPFT");
  if (QFZOPFT->ErrMsg) ErrExit(QFZOPFT->ErrMsg);
  // HMatrix *QTZOPFT = new HMatrix(HDF5File, LHM_HDF5, "QTZOPFT");
  // if (QTZOPFT->ErrMsg) ErrExit(QTZOPFT->ErrMsg);
  printf("size of Q is : %i by %i\n",QFZOPFT->NR,QFZOPFT->NC);   
  Log(" Successfully imported all matrices...\n");
  //---------------------------------------------------------------//
 
  // --- initialize RHS vectors and IF from double *x. 
  Log("  Assembling the RHS vector...");
  HVector* KN1=new HVector(NR,LHM_COMPLEX); 
  HMatrix* PMatrix=new HMatrix(NUML,6,LHM_REAL); 
  double cnorm[NUML]; 
  // KN1=G->AllocateRHSVector();
  for (int l; l<NUML; l++){
    // for each l, compute cnorm[l]
    // IncField *IF1 = new VBeam(l,x[l*5],x[l*5+1],x[l*5+2],x[l*5+3],x[l*5+4]); 
    IncField *IF1 = new VBeam(l,x[l*5]);     // initialize with VBeam = 1*M; 
    G->AssembleRHSVector(Omega, IF1, KN1);   // KN1 is formed here
    cnorm[l]=GetIntegratedIntensity(G, 0, KN1); // get cnorm for each (l,aIn)

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
  //ShowPARMMatrix(PMatrix);

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
  // --- choose Q 
  HMatrix* Q = new HMatrix(NR,NR,LHM_COMPLEX);
  //  Q->Copy(QTZOPFT); // Choose as torque matrix  
  Q->Copy(QFZOPFT); // Choose as torque matrix 

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

  delete KN1, PMatrix, KN, RHS, TEMPMAT, FOM, IF, C, Q, Cconj; 
  return FOM->GetEntryD(0,0); 
}//end objective 
