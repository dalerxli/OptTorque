/*
 * GHBeam.cc   -- gauss-legendre beam implementation of the
 *                IncField class
 * v16. 
 */
#include <stdio.h>
#include <cstdlib> //for abs
#include <string.h>//for memcpy
#include <complex>
#include <math.h>
#include "libIncField.h"
#include "GHBeam.h"
#define II cdouble(0.0,1.0)

double hn_polynomial_value ( int n, double x );
double lm_polynomial (int n, int m, double x);
/**********************************************************************/
/* Gauss- Laguerre                                                    */
/**********************************************************************/
GLBeam::GLBeam(int NewL, int NewP, double Neww0, double NewI0, 
               const char *Newpolname)
{
  L=NewL;   P=NewP;
  w0=Neww0;
  I0=NewI0;  
  polname=Newpolname;
  Cxyz[0] =0.0;
  Cxyz[1] =0.0;
  Cxyz[2] =0.0;
  nHat[0] =0.0; 
  nHat[1] =0.0; 
  nHat[2] =1.0;  
}
/**********************************************************************/
 GLBeam::GLBeam(int NewL, int NewP, double Neww0,double NewI0, 
                const char Newpolname[30],double NewCxyz[3], double NewnHat[3])
{
  L=NewL;   P=NewP;
  w0=Neww0;
  I0=NewI0;  
  polname=Newpolname;
  memcpy(Cxyz,NewCxyz,3*sizeof(double));
  memcpy(nHat, NewnHat, 3*sizeof(double));
}
/**********************************************************************/
GLBeam::~GLBeam(){ }  // Destructor is not made yet. 
/**********************************************************************/
void GLBeam::SetL(int NewL) {   L=NewL; }
void GLBeam::SetP(int NewP) {   P=NewP; }
void GLBeam::Setw0(double Neww0) {   w0=Neww0; }
void GLBeam::SetI0(double NewI0) {   I0=NewI0; }
void GLBeam::Setpol(const char *Newpolname) {   polname=Newpolname; }
void GLBeam::Setpol(const char *Newpolname, cdouble Newalpha, cdouble Newbeta) {
  polname=Newpolname;  alpha=Newalpha;  beta=Newbeta;   }
void GLBeam::SetCxyz(double NewCxyz[3]) {   memcpy(Cxyz,NewCxyz,3*sizeof(double)); }
void GLBeam::SetnHat(double NewnHat[3]) {   memcpy(nHat,NewnHat,3*sizeof(double)); }
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void GLBeam::GetFields(const double X[3], cdouble EH[6])
{
  if ( imag(Eps) !=0.0 || imag(Mu) != 0.0 )
   ErrExit("%s:%i: LG beams not implemented for dispersive media");
  if ( imag(Omega) !=0.0 )
   ErrExit("%s:%i: LG beams not implemented for imaginary frequencies");
  if ( LDim != 0 )
   ErrExit("%s:%i: LG beams not implemented for bloch-periodic geometries");

  GLBeam::GetLG(X,P,L,w0,polname,EH);
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void GLBeam::GetLG(const double X[3],int P,int L,double w0, const char *polname,
           cdouble EH[6])
{
  double rt2=sqrt(2);	
  double w0pow2 = pow(w0,2);
  double w0pow3 = pow(w0,3);
  double w0pow4 = pow(w0,4);
			
  cdouble uG, M[3], N[3];
  //////////////////////////////////////////////////////////////////
  //Coordinates (x,y,z) and (r, Phi, z)
  double XX[3];
  memcpy(XX, X, 3*sizeof(double));
  double x = XX[0];
  double y = XX[1];
  double z = XX[2];
  double r=sqrt(XX[0]*XX[0] + XX[1]*XX[1]);
  double zpow2 = pow(z,2);
  double rpow2 = pow(r,2);
  double rpow3 = pow(r,3);
  double rpow4 = pow(r,4);
  double Phi; 
  if(x==0.0 && y>=0.0)
    Phi=M_PI*0.5; 
  else if(x==0.0 && y<0.0)
    Phi=-M_PI*0.5; 
  else
    Phi=atan2(y,x);
  //////////////////////////////////////////////////////////////////
  //Constants
  double k = std::real(sqrt(Eps*Mu)*Omega); //wavevector
  double kpow2 = pow(k,2);
  cdouble Z=ZVAC*sqrt(Mu/Eps); //relative wave impedance of exterior medium
  //////////////////////////////////////////////////////////////////
  double PP = (double)P; // radial    mode index in (double) 
  double LL = (double)L; // azimuthal mode index in (double) 
  double Labs = std::abs(LL); // |L| in (double)
  double Labspow2 = pow(Labs,2);
  double LLpow2 = pow(LL,2);
  //////////////////////////////////////////////////////////////////
  double Cnorm = 
    sqrt(2.0*std::tgamma(PP+1.0)/(M_PI*w0pow2*std::tgamma(PP+Labs+1.0)));
  if(L==0)
    Cnorm=sqrt(2.0/(M_PI*w0pow2));
  //////////////////////////////////////////////////////////////////
  double zR =k*w0pow2/2.0; // Rayleigh Range
  double zRpow2 = pow(zR,2);
  double zRpow3 = pow(zR,3);
  //////////////////////////////////////////////////////////////////
  double w = w0*sqrt(zpow2/zRpow2+1.0); // beam waist radius w[z]
  double wpow2 = pow(w,2);
  double wpow3 = pow(w,3);
  double wpow4 = pow(w,4);
  double wpow5 = pow(w,5);
  double wpow7 = pow(w,7);
  double RR = z+zRpow2/z;                 // Radius of curvature R[z]
  double RRinv = z/(zpow2 + zRpow2);      // 1/R[z] 
  double RRinvpow2 = pow(RRinv,2);
  double zeta = atan2(z,zR);              // Gouy phase shift
  double Rho = rt2*r/w;                   // Rho[x,y,z]
  double Rhopow2 = pow(Rho,2);
  cdouble Rhopow4 = pow(Rho,4);

  //////////////////////////////////////////////////////////////////
  //Laguerre Functions Simplified
  double LagP, LagPm1, LagPm2; 
  if(P==0){
    LagP=1.0;
    LagPm1=0.0; 
    LagPm2=0.0;
  }else if(P==1){
    LagP= 1.0 + std::abs(LL) - Rhopow2;
    LagPm1=-1.0; 
    LagPm2=0.0;
  }else if(P==2){
    LagP=  lm_polynomial(   P,  std::abs(L), Rhopow2); 
    LagPm1=lm_polynomial(-1+P,1+std::abs(L), Rhopow2); 
    LagPm2=1.0; 
  }else{
    LagP=  lm_polynomial(   P,  std::abs(L), Rhopow2); 
    LagPm1=lm_polynomial(-1+P,1+std::abs(L), Rhopow2); 
    LagPm2=lm_polynomial(-2+P,2+std::abs(L), Rhopow2); 
  }
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  cdouble uG0, M0L1[3], N0L1[3]; // values at r=0 
  uG0=(exp(II*k*z - II*zeta)*w0)/w;
  /// Use M0L1, N0L1 for x=y=0, Abs[L]=1.  
  /// if Abs[L]!=1, M=N=0 
  M0L1[0]=(II*rt2*Cnorm*exp(II*(k*z - 2.0*(1.0 + PP)*zeta))*(1.0 + PP)*w0)/wpow2;
  M0L1[1]=-((rt2*Cnorm*exp(II*(k*z - 2.0*(1.0 + PP)*zeta))*(1.0 + PP)*w0)/wpow2);
  M0L1[2]=0;
  N0L1[0]=(II*rt2*Cnorm*exp(II*(k*z - 2.0*(1.0 + PP)*zeta))*(1.0 + PP)*w0)/wpow2;
  N0L1[1]=-((rt2*Cnorm*exp(II*(k*z - 2.0*(1.0 + PP)*zeta))*(1.0 + PP)*w0)/wpow2);
  N0L1[2]=0;
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  uG=(exp(II*(-0.5)*k*rpow2*RRinv - rpow2/wpow2 + II*k*z - II*zeta)*w0)/w;
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  if(r==0.0) 
  {
    if(Labs==1.0)
      {  /// Use M0L1, N0L1 for x=y=0, Abs[L]=1.  
        memcpy(M, M0L1, 3*sizeof(cdouble));
        memcpy(N, N0L1, 3*sizeof(cdouble));
      }
    else
      {  /// if Abs[L]!=1, M=N=0 
        M[0]=0.0; 
        M[1]=0.0; 
        M[2]=0.0; 
        N[0]=0.0; 
        N[1]=0.0; 
        N[2]=0.0; 
      }
  }
  else //r!=0 
  {
M[0]=(Cnorm*exp(II*(LL*Phi - (Labs + 2.0*PP)*zeta))*pow(Rho,-1.0 + Labs)*uG*(-2.0*rt2*LagPm1*r*Rhopow2*w*y + LagP*(II*LL*Rho*wpow2*x + r*(rt2*Labs*w + r*Rho*(-2.0 - II*k*RRinv*wpow2))*y)))/(rpow2*wpow2);
M[1]=(Cnorm*exp(II*(LL*Phi - (Labs + 2.0*PP)*zeta))*pow(Rho,-1.0 + Labs)*uG*(2.0*rt2*LagPm1*r*Rhopow2*w*x + LagP*(-(rt2*Labs*r*w*x) + rpow2*Rho*(2.0 + II*k*RRinv*wpow2)*x + II*LL*Rho*wpow2*y)))/(rpow2*wpow2);
M[2]=0;
    if(z==0.0)
    {
N[0]=(II*(0.5)*Cnorm*exp(II*LL*Phi - rpow2/w0pow2)*pow(Rho,-1.0 + Labs)*(2.0*(1.0 + Labs + 2.0*PP)*(r*(2.0*LagP*r*Rho - rt2*Labs*LagP*w0 + 2.0*rt2*LagPm1*Rhopow2*w0)*x + II*LagP*LL*Rho*w0pow2*y)*zR + II*k*LagP*LL*Rho*w0pow2*y*(rpow2 - 2.0*zRpow2) + k*r*x*(-(rt2*Labs*LagP*w0*(rpow2 - 2.0*zRpow2)) + 2.0*rt2*LagPm1*Rhopow2*w0*(rpow2 - 2.0*zRpow2) + 2.0*LagP*r*Rho*(rpow2 - w0pow2 - 2.0*zRpow2))))/(k*rpow2*w0pow2*zRpow2);
N[1]=(Cnorm*exp(II*LL*Phi - rpow2/w0pow2)*pow(Rho,-1.0 + Labs)*(2.0*(1.0 + Labs + 2.0*PP)*(LagP*LL*Rho*w0pow2*x + II*r*(2.0*LagP*r*Rho - rt2*Labs*LagP*w0 + 2.0*rt2*LagPm1*Rhopow2*w0)*y)*zR + k*LagP*LL*Rho*w0pow2*x*(rpow2 - 2.0*zRpow2) + II*k*r*y*(-(rt2*Labs*LagP*w0*(rpow2 - 2.0*zRpow2)) + 2.0*rt2*LagPm1*Rhopow2*w0*(rpow2 - 2.0*zRpow2) + 2.0*LagP*r*Rho*(rpow2 - w0pow2 - 2.0*zRpow2))))/(2.*k*rpow2*w0pow2*zRpow2);
N[2]=(Cnorm*exp(II*LL*Phi - rpow2/w0pow2)*pow(Rho,-2.0 + Labs)*(-4.0*LagP*rpow4*Rhopow2 + 4.0*rt2*rpow3*Rho*(Labs*LagP - 2.0*LagPm1*Rhopow2)*w0 + 2.0*rpow2*(-((-1.0 + Labs)*Labs*LagP) + 2.0*(LagP + LagPm1 + 2.0*Labs*LagPm1)*Rhopow2 - 4.0*LagPm2*Rhopow4)*w0pow2 + rt2*r*Rho*(-(Labs*LagP) + 2.0*LagPm1*Rhopow2)*w0pow3 + LagP*LLpow2*Rhopow2*w0pow4))/(k*rpow2*w0pow4);
    }
    else
    {
      N[0]=(-8.0*Cnorm*exp(II*(LL*Phi - (Labs + 2.0*PP)*zeta))*LagPm2*pow(Rho,2.0 + Labs)*uG*w0pow2*x*z)/(k*wpow4*zRpow2) + (Cnorm*exp(II*(LL*Phi - (Labs + 2.0*PP)*zeta))*LagPm1*pow(Rho,Labs)*((8.0*Labs*uG*w0pow2*x*zpow2)/wpow4 + (8.0*(1.0 + Labs)*uG*w0pow2*x*zpow2)/wpow4 + (4.0*rt2*Rho*uG*w0pow2*x*zpow2)/(r*wpow3) - (II*(4.0)*rt2*r*Rho*uG*(II*(-2.0) + k*RRinv*wpow2)*w0pow2*x*zpow2)/wpow5 - (II*(4.0)*rt2*LL*Rho*uG*w0pow2*y*zpow2)/(r*wpow3) + (II*(4.0)*rt2*(Labs + 2.0*PP)*Rho*RRinv*uG*x*zRpow3)/(r*w) + (2.0*rt2*Rho*uG*x*(2.0*wpow2*w0pow2*zpow2 + II*(2.0)*wpow4*zRpow2*(-(k*z) + RRinv*zR) + rpow2*(-4.0*w0pow2*zpow2 + II*k*RRinv*wpow4*(1.0 - 2.0*RRinv*z)*zRpow2)))/(r*wpow5)))/(2.*k*z*zRpow2) + (Cnorm*exp(II*(LL*Phi - (Labs + 2.0*PP)*zeta))*LagP*pow(Rho,Labs)*((-4.0*(-1.0 + Labs)*Labs*uG*w0pow2*x*zpow2)/(Rhopow2*wpow4) - (2.0*rt2*Labs*uG*w0pow2*x*zpow2)/(r*Rho*wpow3) + (2.0*rt2*Labs*r*uG*(2.0 + II*k*RRinv*wpow2)*w0pow2*x*zpow2)/(Rho*wpow5) + (II*(2.0)*rt2*Labs*LL*uG*w0pow2*y*zpow2)/(r*Rho*wpow3) - (II*(2.0)*rt2*Labs*(Labs + 2.0*PP)*RRinv*uG*x*zRpow3)/(r*Rho*w) - (2.0*(Labs + 2.0*PP)*RRinv*uG*(II*(-2.0) + k*RRinv*wpow2)*x*zRpow3)/wpow2 - (2.0*LL*(Labs + 2.0*PP)*RRinv*uG*y*zRpow3)/rpow2 + (LL*uG*y*(II*(2.0)*wpow2*w0pow2*zpow2 + 2.0*wpow4*zRpow2*(k*z - RRinv*zR) + rpow2*(II*(-4.0)*w0pow2*zpow2 + k*RRinv*wpow4*(-1.0 + 2.0*RRinv*z)*zRpow2)))/(rpow2*wpow4) + (rt2*Labs*uG*x*(rpow2*(4.0*w0pow2*zpow2 + II*k*RRinv*wpow4*(-1.0 + 2.0*RRinv*z)*zRpow2) + II*(2.0)*wpow2*(II*w0pow2*zpow2 + wpow2*zRpow2*(k*z - RRinv*zR))))/(r*Rho*wpow5) + (exp(-(rpow2/wpow2) - II*(0.5)*k*(rpow2*RRinv - 2.0*z) - II*zeta)*w0*x*(rpow2*(II*(-2.0) + k*RRinv*wpow2)*(II*(-4.0)*w0pow2*zpow2 + k*RRinv*wpow4*(-1.0 + 2.0*RRinv*z)*zRpow2) + 2.0*wpow2*((6.0 + II*k*RRinv*wpow2)*w0pow2*zpow2 + wpow2*zRpow2*(kpow2*RRinv*wpow2*z + II*k*(-(RRinv*wpow2) - 2.0*z + RRinvpow2*wpow2*(2.0*z + II*zR)) + II*(2.0)*RRinv*zR))))/wpow7))/(2.*k*z*zRpow2);
N[1]=(-8.0*Cnorm*exp(II*(LL*Phi - (Labs + 2.0*PP)*zeta))*LagPm2*pow(Rho,2.0 + Labs)*uG*w0pow2*y*z)/(k*wpow4*zRpow2) + (Cnorm*exp(II*(LL*Phi - (Labs + 2.0*PP)*zeta))*LagPm1*pow(Rho,Labs)*((II*(4.0)*rt2*LL*Rho*uG*w0pow2*x*zpow2)/(r*wpow3) + (8.0*Labs*uG*w0pow2*y*zpow2)/wpow4 + (8.0*(1.0 + Labs)*uG*w0pow2*y*zpow2)/wpow4 + (4.0*rt2*Rho*uG*w0pow2*y*zpow2)/(r*wpow3) - (II*(4.0)*rt2*r*Rho*uG*(II*(-2.0) + k*RRinv*wpow2)*w0pow2*y*zpow2)/wpow5 + (II*(4.0)*rt2*(Labs + 2.0*PP)*Rho*RRinv*uG*y*zRpow3)/(r*w) + (2.0*rt2*Rho*uG*y*(2.0*wpow2*w0pow2*zpow2 + II*(2.0)*wpow4*zRpow2*(-(k*z) + RRinv*zR) + rpow2*(-4.0*w0pow2*zpow2 + II*k*RRinv*wpow4*(1.0 - 2.0*RRinv*z)*zRpow2)))/(r*wpow5)))/(2.*k*z*zRpow2) + (Cnorm*exp(II*(LL*Phi - (Labs + 2.0*PP)*zeta))*LagP*pow(Rho,Labs)*((II*(-2.0)*rt2*Labs*LL*uG*w0pow2*x*zpow2)/(r*Rho*wpow3) - (4.0*(-1.0 + Labs)*Labs*uG*w0pow2*y*zpow2)/(Rhopow2*wpow4) - (2.0*rt2*Labs*uG*w0pow2*y*zpow2)/(r*Rho*wpow3) + (2.0*rt2*Labs*r*uG*(2.0 + II*k*RRinv*wpow2)*w0pow2*y*zpow2)/(Rho*wpow5) + (2.0*LL*(Labs + 2.0*PP)*RRinv*uG*x*zRpow3)/rpow2 - (II*(2.0)*rt2*Labs*(Labs + 2.0*PP)*RRinv*uG*y*zRpow3)/(r*Rho*w) - (2.0*(Labs + 2.0*PP)*RRinv*uG*(II*(-2.0) + k*RRinv*wpow2)*y*zRpow3)/wpow2 + (LL*uG*x*(II*(4.0)*rpow2*w0pow2*zpow2 - II*(2.0)*wpow2*w0pow2*zpow2 + wpow4*zRpow2*(k*(-2.0*z + rpow2*(RRinv - 2.0*RRinvpow2*z)) + 2.0*RRinv*zR)))/(rpow2*wpow4) + (rt2*Labs*uG*y*(rpow2*(4.0*w0pow2*zpow2 + II*k*RRinv*wpow4*(-1.0 + 2.0*RRinv*z)*zRpow2) + II*(2.0)*wpow2*(II*w0pow2*zpow2 + wpow2*zRpow2*(k*z - RRinv*zR))))/(r*Rho*wpow5) + (exp(-(rpow2/wpow2) - II*(0.5)*k*(rpow2*RRinv - 2.0*z) - II*zeta)*w0*y*(rpow2*(II*(-2.0) + k*RRinv*wpow2)*(II*(-4.0)*w0pow2*zpow2 + k*RRinv*wpow4*(-1.0 + 2.0*RRinv*z)*zRpow2) + 2.0*wpow2*((6.0 + II*k*RRinv*wpow2)*w0pow2*zpow2 + wpow2*zRpow2*(kpow2*RRinv*wpow2*z + II*k*(-(RRinv*wpow2) - 2.0*z + RRinvpow2*wpow2*(2.0*z + II*zR)) + II*(2.0)*RRinv*zR))))/wpow7))/(2.*k*z*zRpow2);
N[2]=(-8.0*Cnorm*exp(II*(LL*Phi - (Labs + 2.0*PP)*zeta))*LagPm2*pow(Rho,2.0 + Labs)*uG)/(k*wpow2) + (2.0*Cnorm*exp(II*(LL*Phi - (Labs + 2.0*PP)*zeta))*LagPm1*pow(Rho,Labs)*uG*(rt2*Rho*wpow2 + 2.0*r*(w + 2.0*Labs*w) + 2.0*rt2*rpow2*Rho*(-2.0 - II*k*RRinv*wpow2)))/(k*r*wpow3) + (Cnorm*exp(II*(LL*Phi - (Labs + 2.0*PP)*zeta))*LagP*pow(Rho,-2.0 + Labs)*uG*(-(rt2*Labs*r*Rho*wpow3) + LLpow2*Rhopow2*wpow4 + 2.0*rt2*Labs*rpow3*Rho*w*(2.0 + II*k*RRinv*wpow2) + rpow4*Rhopow2*pow(II*(-2.0) + k*RRinv*wpow2,2) + 2.0*rpow2*wpow2*(Labs - Labspow2 + Rhopow2*(2.0 + II*k*RRinv*wpow2))))/(k*rpow2*wpow4);
    }//end else following if z==0
}///end else following if(r==0.0)

//// For Azimuthal (TE) Polarization, let's use E=aM, H=aN. 
//// EH[0]=sum a[P,L]*M[P,L]
//// EH[1]=sum a[P,L]*N[P,L] 

//// For now, let's use a single mode of M for E and single mode of N for H. 
EH[0]=M[0];
EH[1]=M[1];
EH[2]=M[2];
EH[3]=N[0]/Z;
EH[4]=N[1]/Z;
EH[5]=N[2]/Z;
}//end getLG

/**********************************************************************/
// Gauss-Hermite 
/**********************************************************************/
GHBeam::GHBeam(int NewHM, int NewHN)
{
  HM=NewHM;
  HN=NewHN;
  w0=2.0;
  E0[0]=1.0;
  E0[1]=0.0;
  E0[2]=0.0;
}
GHBeam::GHBeam(int NewHM, int NewHN, double Neww0)
{
  HM=NewHM;
  HN=NewHN;
  w0=Neww0;
}
GHBeam::GHBeam(int NewHM, int NewHN, double Neww0, double NewCxyz[3], cdouble NewE0[3],double NewnHat[3])
{
  HM=NewHM;
  HN=NewHN;
  w0=Neww0;
  memcpy(Cxyz,NewCxyz,3*sizeof(double));
  memcpy(E0, NewE0, 3*sizeof(cdouble));
  memcpy(nHat, NewnHat, 3*sizeof(double));
}
void GHBeam::SetHM(int NewHM) {   HM=NewHM; }
void GHBeam::SetHN(int NewHN) {   HN=NewHN; }
void GHBeam::Setw0(double Neww0) {   w0=Neww0; }
void GHBeam::SetCxyz(double NewCxyz[3]) {   memcpy(Cxyz,NewCxyz,3*sizeof(double)); }
void GHBeam::SetE0(cdouble NewE0[3]) {  memcpy(E0, NewE0, 3*sizeof(cdouble)); }
void GHBeam::SetnHat(double NewnHat[3]) {   memcpy(nHat,NewnHat,3*sizeof(double)); }
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void GHBeam::GetFields(const double X[3], cdouble EH[6])
{
  // get the wavenumber and relative wave impedance of the
  // exterior medium
  cdouble k=sqrt(Eps*Mu)*Omega; //wavevector
  cdouble Z=ZVAC*sqrt(Mu/Eps); //impedance
  cdouble ExpFac = exp(II*k*X[2]);//for now, propagate in z
  //cdouble ExpFac=exp(II*k*(nHat[0]*X[0] + nHat[1]*X[1] + nHat[2]*X[2]));

// convert the evaluation point to cylindrical coordinates
//  double XX[3];
//  memcpy(XX, X, 3*sizeof(double));
//  double r, Phi, z;
//  r=sqrt(XX[0]*XX[0] + XX[1]*XX[1]);
//  Phi=atan2(XX[1], XX[0]);
//  z = XX[2];

  double lambda = 2.0*M_PI*pow(10,-6)/std::real(Omega);
  double zR = std::real(k)*w0*w0/2.0; //Rayleigh Range
  cdouble q = X[2] + II*zR;
  cdouble qconj = X[2] - II*zR;
  double  w = w0*sqrt(pow(X[2]/zR,2.0)+1.0); // waist at z
  // double zeta = atan2(z,zR);

  cdouble HG_ux = pow(pow(2,HM)*std::tgamma(HM+1)*w0,-0.5)
                    *pow(-qconj/q,HM*0.5)
                    *hn_polynomial_value(HM, sqrt(2)*X[0]/w)
                    *exp(-II*k*X[0]*X[0]/2.0/q);

  cdouble HG_uy = pow(pow(2,HN)*std::tgamma(HN+1)*w0,-0.5)
                    *pow(-qconj/q,HN*0.5)
                    *hn_polynomial_value(HN, sqrt(2)*X[1]/w)
                    *exp(-II*k*X[1]*X[1]/2.0/q);

  cdouble HG_u = sqrt(2.0/M_PI)*(II*zR/q)*HG_ux*HG_uy;


/***************************************************************/
  // stuff the field components into the output vector
  // E0 specifies polarization.
  EH[0] = E0[0] *HG_u* ExpFac;
  EH[1] = E0[1] *HG_u* ExpFac;
  EH[2] = E0[2] *HG_u* ExpFac;
    /* H = (nHat \cross E) / Z */
  EH[3] = (nHat[1]*EH[2] - nHat[2]*EH[1]) / Z;
  EH[4] = (nHat[2]*EH[0] - nHat[0]*EH[2]) / Z;
  EH[5] = (nHat[0]*EH[1] - nHat[1]*EH[0]) / Z;
}
/**********************************************************************/
/* Polynomials                            *****************************/
/**********************************************************************/
double hn_polynomial_value ( int n, double x )
//****************************************************************************80
//    HN_POLYNOMIAL_VALUE evaluates Hn(i,x).
//
//    Hn(i,x) is the normalized physicist's Hermite polynomial of degree I.
//    These polynomials satisfy the orthonormality condition:
//      Integral ( -oo < X < +oo )
//        exp ( - X^2 ) * Hn(M,X) Hn(N,X) dX = delta ( N, M )
//  Licensing:
//    This code is distributed under the GNU LGPL license.
//  Modified:
//    26 February 2012
//  Author:
//    John Burkardt
//  Reference:
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//****************************************************************************80
//    Yoonkyung Eunnie Lee, 2015
//    function simplified to compute Hn(x) for a single double value x. 
//    also, r8_pi -> M_Pi
//****************************************************************************80
{
  double pj, pj1, pj2, pnorm, fact, two_power; 
  pnorm = fact = two_power = 1.0;
  if ( n < 0 )
    { printf("polynomial order n cannot be negative, changing n to 1.\n"); 
      return NAN; 
    }else if(n==0)
    {
      return 1.0;
    }else
    {
      pnorm = sqrt(fact*two_power*sqrt(M_PI)); 
      fact = fact*(double)(1); 
      two_power= two_power*2.0; 
      pj2 = 1.0; 
      pj1 = 2.0*x; 
      pnorm = sqrt(fact*two_power*sqrt(M_PI)); 
      fact = fact*(double)(2); 
      two_power= two_power*2.0; 
      for(int j=2;j<=n;j++){
        pj = 2.0*x*pj1 - 2.0*(double)(j-1)*pj2; 
        pj2 = pj1;
        pj1 = pj; 
        pnorm = sqrt(fact*two_power*sqrt(M_PI)); 
        fact = fact*(double)(j+1); 
        two_power= two_power*2.0; 
      }
      return pj1/pnorm; 
    } 
}
double lm_polynomial( int n, int m, double x)
//****************************************************************************80
//    LM_POLYNOMIAL evaluates Laguerre polynomials Lm(n,m,x).
//  Recursion:
//    Lm(0,M,X)   = 1
//    Lm(1,M,X)   = (M+1-X)
//    if 2 <= N:
//      Lm(N,M,X)   = ( (M+2*N-1-X) * Lm(N-1,M,X)
//                   +   (1-M-N)    * Lm(N-2,M,X) ) / N
//  Licensing:
//    This code is distributed under the GNU LGPL license.
//  Modified:
//    10 March 2012
//  Author:
//    John Burkardt
//****************************************************************************80
//    Yoonkyung Eunnie Lee, 2015
//    function simplified to compute Lm(n,x) for a single double value x. 
//****************************************************************************80
{
  int i, j;
  double vj, vj2, vj1 ; 

  if ( n < 0 )
    {
      printf("polynomial order n cannot be negative, changing n to 1.\n"); 
      return NAN; 
    }else if (n==0)
    {
      return x; 
    }else
    {
      vj2 = 1.0; // first term always 1.0; 
      vj1 = (m+1)-x; 
      for(int j=2;j<=n;j++){
        vj = (((m+2*j-1)-x)*vj1 + (-m-j+1)*vj2) / j;
        vj2 = vj1;
        vj1 = vj;         
      }
      return vj1; 
    }
}
