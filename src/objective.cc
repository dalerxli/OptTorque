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
//#include "OptTorque.h" 

#define II cdouble(0.0,1.0)
#define NUML 3
#define MAXSTR 100 
using namespace scuff;
//---------------------------------------------------------------//
double GetIntegratedIntensity(RWGGeometry *G, int SurfaceIndex, 
                              HVector *RHSVector);

void ShowPARMMatrix(HMatrix* PARMMatrix); 
double objective(unsigned n, double *x, double *grad);
//---------------------------------------------------------------//
int main()
{
  // create a vector x to test. 
  double x[NUML*5]; 
  double grad[2]; 

  for (int l = 0; l<NUML; l++)
    {
      x[l*5] = 10.0;  //alpha 
      x[l*5+1] = 1.0; //ar
      x[l*5+2] = 0.0; //br
      x[l*5+3] = 0.0; //ai
      x[l*5+4] = 0.0; //bi
    }
  objective(NUML*5, x, grad); 
}
//---------------------------------------------------------------//

double objective(unsigned n, double *x, double *grad) 
{
  //---------------------------------------------------------------//
  cdouble Omega = 10.471975511965978; 
  //---------------------------------------------------------------//
  char FileBase[MAXSTR], GeoFile[MAXSTR], HDF5File[MAXSTR];
  snprintf(FileBase,MAXSTR,"N3_400nm_Mesh60nm"); 
  snprintf(GeoFile,MAXSTR,"%s.scuffgeo",FileBase);  
  snprintf(HDF5File,MAXSTR,"%s.HDF5",FileBase ); 
  //I want the GeoFile and HDF5File to search for *.scuffgeo and *.HDF5 in the dir. 
  RWGGeometry *G = new RWGGeometry(GeoFile);
  Log(" Created RWGGeometry G...\n");
  //---------------------------------------------------------------//
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
  Log(" Successfully imported all matrices...\n");
  //---------------------------------------------------------------//
  // initialize RHS vectors and IF from double *x. 
  Log("  Assembling the RHS vector...");
  HVector* KN1=new HVector(NR,LHM_COMPLEX); 
  HMatrix* PMatrix=new HMatrix(NUML,6,LHM_REAL); 
  double cnorm[NUML]; 
  //KN1=G->AllocateRHSVector();
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
  HVector* KN=new HVector(NR,LHM_COMPLEX); 
  HVector* RHS=new HVector(NR,LHM_COMPLEX); 
  IncField* IF=new VBeam(PMatrix);
  G->AssembleRHSVector(Omega, IF, KN); 
  RHS->Copy(KN);
  // printf("size of RHS is : %i \n",RHS->N);
  MLU->LUFactorize(); 
  MLU->LUSolve(KN);// solved KN. 
  HMatrix* KNMAT = new HMatrix(NR,1,KN->RealComplex, LHM_NORMAL, KN->ZV); 
  printf("KNMAT[0,0] = %e+%ei \n",real(KNMAT->GetEntry(0,0)),imag(KNMAT->GetEntry(0,0)));   
  
  HMatrix* TempN = new HMatrix(NR,1,LHM_COMPLEX); 
  HMatrix* Temp1 = new HMatrix(1, 1,LHM_COMPLEX); 
  QTZOPFT->Multiply(KNMAT,TempN);

  KNMAT->Multiply(TempN,Temp1,"--transA C"); //Temp1 = Adj(KN)*Q*KN
  printf("FOM = %e+%ei \n",real(Temp1->GetEntry(0,0)),imag(Temp1->GetEntry(0,0)));
  double FOM = Temp1->GetEntryD(0,0); 
  delete TempN, Temp1; 
  if(grad)
    { // M*C_adj = transpose(QPFT)*conj(C)  
      HMatrix* RHSMAT = new HMatrix(NR,1,RHS->RealComplex, LHM_NORMAL, RHS->ZV); 
      // printf("KN[0] = %e+%ei \n",real(KN->GetEntry(1)),imag(KN->GetEntry(1)));
      KNMAT->Adjoint();
      KNMAT->Transpose(); // KNMAT = conj(KN); 
     
      HMatrix* Cadj = new HMatrix(NR,1,LHM_COMPLEX); 
      printf("size of KNMAT is : %i x %i \n",KNMAT->NR,KNMAT->NC);   
      printf("size of Cadj is : %i x %i \n",Cadj->NR,Cadj->NC);   
      QabsOPFT->Multiply(KNMAT,Cadj,"--transA C"); 
      printf("Cadj[0,0] = %e+%ei \n",real(Cadj->GetEntry(0,0)),imag(Cadj->GetEntry(0,0)));   
      
      MLU->LUSolve(Cadj); 
      Log(" Computed Cadj...\n");
      char MatFileName[100];
      snprintf(MatFileName, 100, "Mat_Cadj.dat");  
      Cadj->ExportToText(MatFileName,"--separate,"); 

      HMatrix* dFOMMAT = new HMatrix(1,1,LHM_COMPLEX);  
      //dFOM is grad. grad[0] and grad[1] are different things. 
      Cadj->Multiply(RHSMAT,dFOMMAT,"transA"); 
      grad[0] = 0.0; 
      grad[1] = dFOMMAT->GetEntryD(0,0); //This is not right . real and imag
      printf("grad[1] = %e+%ei\n",std::real(grad[1]),std::imag(grad[1]));
 
      delete RHSMAT, Cadj, dFOMMAT; 
    }
  delete KNMAT; 
  return FOM; 
}//end objective 
