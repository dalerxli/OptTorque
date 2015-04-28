/*
 * VBeam.cc   -- Vectorbeam implementation of the
 *                IncField class
 * v15. : changed coefficients for EH summation (-1/(II*Z))
 * radial polarization
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
// Mode index L is accepted as input
// azimuthal polarization is chosen
{ 
  L=NewL;  aIn=NewaIn;   
  int Li;  
  for(Li=0;Li<LSPAN;Li++) // initialize
    {
      aL[Li]=0.0; 
      bL[Li]=0.0; 
    }
  /*
  if(L>=0)        //set a single component of a or b to 1.0 
    aL[L]=1.0;    // using azimuthal polarization
  else{
    Li = LMAX-L; 
    aL[Li]=1.0;   // using azimuthal polarization
  }
  */ 
  if(L>=0)   /// radial polarization 
       bL[L]=1.0;
   else{
     Li = LMAX-L; 
     bL[Li]=1.0;    
   }
 } 
VBeam::VBeam(double NewaL[LSPAN],double NewbL[LSPAN], double NewaIn)
// VBeam initialization for superposition of multiple modes
// List of coefficients aL and bL are accepted as inputs
{ L=LMAX+1;  aIn=NewaIn; 
  memcpy(aL,NewaL,LSPAN*sizeof(double));
  memcpy(bL,NewbL,LSPAN*sizeof(double));  }
/**********************************************************************/
VBeam::~VBeam(){} // Destructor is not made yet. 
/**********************************************************************/
void VBeam::SetL(int NewL) {   L=NewL; }
void VBeam::Setab(double NewaL[LSPAN],double NewbL[LSPAN])
{ 
  memcpy(aL,NewaL,LSPAN*sizeof(double));
  memcpy(bL,NewbL,LSPAN*sizeof(double));  
}
void VBeam::SetaIn(double NewaIn) {   aIn=NewaIn; }
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

  if(-LMAX<L && L<LMAX) 
    { // run for a single mode if -LMAX<L<LMAX 
      VBeam::GetMN(X,L,aIn,M,N);
      int Li=0; 
      if(L>=0)
        Li=L; 
      else
        Li=LMAX-L;      
      EH[0]=-(aL[Li]*M[0]+bL[Li]*N[0]);
      EH[1]=-(aL[Li]*M[1]+bL[Li]*N[1]);
      EH[2]=-(aL[Li]*M[2]+bL[Li]*N[2]);
      EH[3]=-(bL[Li]*M[0]+aL[Li]*N[0])/(II*Z);
      EH[4]=-(bL[Li]*M[1]+aL[Li]*N[1])/(II*Z);
      EH[5]=-(bL[Li]*M[2]+aL[Li]*N[2])/(II*Z);
    }
  else //run for superposition of modes if L>LMAX or L<-LMAX
    {
      int L0, Li; // input L0 and index Li
      int j;
      for(j=0;j<=6;j++)        // initialize EH
        EH[j]=0.0; 
      
      for(Li=0;Li<LSPAN;Li++)  // sum EH for each entry
        {
          if(Li<=LMAX)
            L0=Li;
          else
            L0=-Li+LMAX; 

          if(aL[Li]==0&&bL[Li]==0){
            //don't run GetMN if coeff.s are zero.
          }else{
            VBeam::GetMN(X,L0,aIn,M,N);
            EH[0]+=-(aL[Li]*M[0]+bL[Li]*N[0]);
            EH[1]+=-(aL[Li]*M[1]+bL[Li]*N[1]);
            EH[2]+=-(aL[Li]*M[2]+bL[Li]*N[2]);
            EH[3]+=-(bL[Li]*M[0]+aL[Li]*N[0])/(II*Z);
            EH[4]+=-(bL[Li]*M[1]+aL[Li]*N[1])/(II*Z);
            EH[5]+=-(bL[Li]*M[2]+aL[Li]*N[2])/(II*Z);
          }
        }
    }//endif(L==LMAX+1)
  
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
  int j;
  for(j=0;j<=2;j++)
    M[j]=0.0; N[j]=0.0;

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
