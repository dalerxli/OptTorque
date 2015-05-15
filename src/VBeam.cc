/*
 * VBeam.cc   -- Vectorbeam implementation of the
 *                IncField class
 * 
 * v15. : changed coefficients for EH summation (-1/(II*Z))
 * v16. : GetField is debugged to take correct HMatrix Input 
 */
#include <stdio.h>
#include <cstdlib> //for abs
#include <string.h>//for memcpy
#include <complex>
#include <math.h>
#include "VBeam.h"
#include <libIncField.h>
#define II cdouble(0.0,1.0)
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
VBeam::VBeam(int NewL, double NewaIn) 
// VBeam initialization for a single mode
// Default is to use radial polarization (bL[L]=1)
{ 
  PMatrix->SetEntry(0,0,NewL);    //L
  PMatrix->SetEntry(0,1,NewaIn);  //aIn
  PMatrix->SetEntry(0,2,0.0);     //aL[L]
  PMatrix->SetEntry(0,3,1.0);     //bL[L]
  numL= 1; 
 } 
/**********************************************************************/
VBeam::VBeam(HMatrix *NewPMatrix) 
{
  if ( NewPMatrix->NC!=4 )
    ErrExit("%s:%i: NewPMatrix should have 4 columns.");

  if ( NewPMatrix->NR>25 )
    ErrExit("%s:%i: NewPMatrix has too many rows >25.");

  numL=NewPMatrix->NR; 
  int iL; 
  for(iL=0;iL<numL;iL++){
    PMatrix->SetEntry(iL,0,NewPMatrix->GetEntryD(iL,0)); 
    PMatrix->SetEntry(iL,1,NewPMatrix->GetEntryD(iL,1)); 
    PMatrix->SetEntry(iL,2,NewPMatrix->GetEntryD(iL,2)); 
    PMatrix->SetEntry(iL,3,NewPMatrix->GetEntryD(iL,3)); 
  }
}
/**********************************************************************/
VBeam::~VBeam(){} // Destructor is not made yet. 
/**********************************************************************/
void VBeam::SetCxyz(double NewCxyz[3])
{ memcpy(Cxyz,NewCxyz,3*sizeof(double)); }
void VBeam::SetnHat(double NewnHat[3]) 
{ memcpy(nHat,NewnHat,3*sizeof(double)); }
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void VBeam::GetFields(const double X[3], cdouble EH[6])
{
  if ( imag(Eps) !=0.0 || imag(Mu) != 0.0 )
    ErrExit("%s:%i: LG beams not implemented for dispersive media");
  if ( imag(Omega) !=0.0 )
    ErrExit("%s:%i: LG beams not implemented for imaginary frequencies");
  if ( LDim != 0 )
    ErrExit("%s:%i: LG beams not implemented for bloch-periodic geometries");
  
  cdouble Z=ZVAC*sqrt(Mu/Eps); //relative wave impedance of exterior medium

  int j;
  for(j=0;j<6;j++)        // initialize EH
    EH[j]=0.0; 
  
  int iL, L; 
  double aIn,aL,bL; 
  for(iL=0; iL<numL; iL++) //for each row of PMatrix
    {
      L   =PMatrix->GetEntryD(iL,0);
      aIn =PMatrix->GetEntryD(iL,1);
      aL  =PMatrix->GetEntryD(iL,2);
      bL  =PMatrix->GetEntryD(iL,3);
      //      printf("aL=%f,bL=%f\n",aL,bL);
      if(aL==0.0&&bL==0.0){}///don't run GetMN if coeff.s are zero
      else{
        VBeam::GetMN(X,L,aIn,M,N);
        EH[0]+=-(aL*M[0]);
        EH[0]+=-(bL*N[0]);
        EH[1]+=-(aL*M[1]);
        EH[1]+=-(bL*N[1]);
        EH[2]+=-(aL*M[2]);
        EH[2]+=-(bL*N[2]);
        EH[3]+=-(bL*M[3])/(II*Z);
        EH[3]+=-(aL*N[3])/(II*Z);
        EH[4]+=-(bL*M[4])/(II*Z);
        EH[4]+=-(aL*N[4])/(II*Z);
        EH[5]+=-(bL*M[5])/(II*Z);
        EH[5]+=-(aL*N[5])/(II*Z);
      }
    }
 }//end GetField
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void VBeam::GetMN(const double X[3], int L, double aIn, cdouble M[3], cdouble N[3])
{
  ////////////////////////////////////////	
  //Constants
  double k = std::real(sqrt(Eps*Mu)*Omega); //wavevector
  double kz = cos(aIn*M_PI/180.0)*k;
  double kt = sin(aIn*M_PI/180.0)*k;
  double kpow2 = pow(k,2);
  double ktpow2 = pow(kt,2);
  double kzpow2 = pow(kz,2);
  ////////////////////////////////////////	
  //Coordinates (x,y,z) and (r, Phi, z)
  double XX[3];
  memcpy(XX, X, 3*sizeof(double));
  double x = XX[0];
  double y = XX[1];
  double z = XX[2];
  double r=sqrt(XX[0]*XX[0] + XX[1]*XX[1]);
  double Phi; 
  if(x==0.0 && y>=0.0)
    Phi=M_PI*0.5; 
  else if(x==0.0 && y<0.0)
    Phi=-M_PI*0.5; 
  else
    Phi=atan2(y,x);
  ////////////////////////////////////////	
  double LL = (double)L; 
  //Bessel Functions Simplified
  double BesP =   jn(L,kt*r);
  double BesPm1 = jn(L-1,kt*r);
  double BesPm2 = jn(L-2,kt*r); 
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  /// initialize M=N=0. 
  for(int j=0;j<=2;j++){///catch errors this way. j should live only here.
    M[j]=0.0; 
    N[j]=0.0;
  }

  /// Set M, N accordingly 
  if(r!=0.0){
    ////////////////////////////////////////	
    /// For Normal Cases 
    M[0]=exp(II*(LL*Phi + kz*z))*((BesPm1*kt*y)/r - (BesP*LL)/(II*x + y));
    M[1]=exp(II*(LL*Phi + kz*z))*(-((BesPm1*kt*x)/r) + (BesP*LL)/(x - II*y));
    M[2]=0.0;
    N[0]=(II*exp(II*(LL*Phi + kz*z))*kz*((BesPm1*kt*x)/r - (BesP*LL)/(x - II*y)))/k;
    N[1]=(exp(II*(LL*Phi + kz*z))*kz*(-((BesP*LL)/(x - II*y)) + (II*BesPm1*kt*y)/r))/k;
    N[2]=(BesP*exp(II*(LL*Phi + kz*z))*ktpow2)/k;
  }else if(r==0.0&&L==1){
    ////////////////////////////////////////	
    /// For L=1, r=0.0 
    M[0]=II*(0.5)*exp(II*kz*z)*kt;
    M[1]=-(exp(II*kz*z)*kt)/2.;
    M[2]=0.0;
    N[0]=(II*(0.5)*exp(II*kz*z)*kt*kz)/k;
    N[1]=-(exp(II*kz*z)*kt*kz)/(2.*k);
    N[2]=0.0;
  }else if(r==0.0&&L==0){
    ////////////////////////////////////////	
    /// For L=0, r=0.0, only Nz is nonzero. 
    N[2]=(exp(II*kz*z)*ktpow2)/k;
  }//end else following if(r==0)
  ////////////////////////////////////////	
  /// For L>1, r=0.0, M=N=0. 
 }//end getMN
