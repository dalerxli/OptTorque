/*
 * GHBeam.cc   -- gauss-legendre beam implementation of the
 *                IncField class
 * v12. 
 */
#include <stdio.h>
#include <cstdlib> //for abs
#include <string.h>//for memcpy
#include <complex>
#include <math.h>
#include "GHBeam.h"

#include <libIncField.h>
#define II cdouble(0.0,1.0)

double *hn_polynomial_value ( int m, int n, double x[] );
double hn_polynomial_value ( int n, double x );
double *lm_polynomial (int mm, int n, int m, double x[]);
double lm_polynomial (int n, int m, double x);

/**********************************************************************/
/* Gauss- Laguerre                                                    */
/**********************************************************************/
GLBeam::GLBeam(int NewL, int NewP, double Neww0, double NewI0, 
               const char *Newpolname)
{
  L=NewL;
  P=NewP;
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
 GLBeam::GLBeam(int NewL, int NewP, double Neww0,double NewI0, 
                const char Newpolname[30], cdouble Newalpha, cdouble Newbeta,
                double NewCxyz[3], double NewnHat[3])
{
  L=NewL;
  P=NewP;
  w0=Neww0;
  I0=NewI0;  
  polname=Newpolname;
  alpha=Newalpha;
  beta=Newbeta; 
  memcpy(Cxyz,NewCxyz,3*sizeof(double));
  memcpy(nHat, NewnHat, 3*sizeof(double));
}
GLBeam::~GLBeam()
{ 
  // Destructor is not made yet. 
}
void GLBeam::SetL(int NewL) {   L=NewL; }
void GLBeam::SetP(int NewP) {   P=NewP; }
void GLBeam::Setw0(double Neww0) {   w0=Neww0; }
void GLBeam::SetI0(double NewI0) {   I0=NewI0; }
void GLBeam::Setpol(const char *Newpolname) {   polname=Newpolname; }
void GLBeam::Setpol(const char *Newpolname, cdouble Newalpha, cdouble Newbeta) {
  polname=Newpolname;
  alpha=Newalpha;
  beta=Newbeta;   
}
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

  // get the wavenumber and relative wave impedance of the exterior medium
  const cdouble IU(0,1);
  double EpsR = real(Eps);
  double MuR  = real(Mu);
  double k    = sqrt(EpsR*MuR)*real(Omega); // wavenumber of medium 
  cdouble Z=ZVAC*sqrt(Mu/Eps);              // impedance
  double ZR   = sqrt(MuR/EpsR);             // relative wave impedance of medium


  // convert the evaluation point to cylindrical coordinates
  double XX[3];
  memcpy(XX, X, 3*sizeof(double));
  double x,y,z, r2, z2, Phi;
  x = XX[0];
  y = XX[1];
  z = XX[2];
  r2 = XX[0]*XX[0] + XX[1]*XX[1];
  z2 = z*z; 
  Phi=atan2(XX[1], XX[0]);

  // Updated for v11
  /***************************************************************/
  double PP = (double)P;
  double LL = (double)L;
  double Cnorm = pow(2.0*std::tgamma(PP+1.0)/(M_PI*w0*w0*std::tgamma(PP+std::abs(LL))),0.5);
  cdouble Exfactor, Eyfactor, Ezfactor, Hxfactor, Hyfactor, Hzfactor ;

  if(strcmp(pollist[0],polname)==0){
    alpha = 1.0;  beta = 0.0;}
  else if(strcmp(pollist[1],polname)==0){
    alpha = 0.0;  beta = 1.0;}
  else if(strcmp(pollist[2],polname)==0){
    alpha = 1.0/sqrt(2); beta = II/sqrt(2); }
  else if(strcmp(pollist[3],polname)==0){
    alpha = 1.0/sqrt(2); beta = -II/sqrt(2); }
  else if(strcmp(pollist[4],polname)==0){
    alpha = x/sqrt(r2);   ///cos(phi)
    beta =  y/sqrt(r2);   ///sin(phi) 
  }
  else if(strcmp(pollist[5],polname)==0){
    alpha = -y/sqrt(r2);  ///-sin(phi)
    beta =  x/sqrt(r2);   ///cos(phi)
  }
  else{ //if none matches, assign radial polarization
    alpha = x/sqrt(r2);   ///cos(phi)
    beta =  y/sqrt(r2);   ///sin(phi) 
  }
  //  printf("polname=%s, alpha, beta = (%e,%e i ),(%e,%e i)\n",polname,real(alpha),imag(alpha),real(beta),imag(beta));
 // Parameters used in Laguerre Gaussian Beams
  double lambda = 2.0*M_PI*pow(10,-6)/std::real(Omega);
  double zR = std::real(k)*w0*w0/2.0; //Rayleigh Range
  //  cdouble q = X[2] + II*zR;
  // cdouble qconj = X[2] - II*zR;
  double  w = w0*sqrt(pow(z/zR,2.0)+1.0); // waist at z
  double zeta = atan2(z,zR);

  cdouble uLG = I0*Cnorm*(w0/w)*pow(2.0,std::abs(LL)/2.0)*exp(II*k*z - (II*0.5*k*r2*z)/(z2 + zR*zR) - r2/w/w + II*LL*atan2(y,x) -II*(1.0+ 2.0*PP + std::abs(LL))*zeta)*
     pow(sqrt(r2)/w,std::abs(LL))*lm_polynomial(P,std::abs(L),2.0*r2/w/w);

//////////////////////////////////////////////////////////////////
/// Ex /uLG, Ey /uLG
//////////////////////////////////////////////////////////////////
Exfactor = alpha; 
Eyfactor = beta; 
Ezfactor = ((-II)*((-II)*(LL*(zR*zR*w*w)*(y*alpha - x*beta) + k*w0*w0*r2*z*(x*alpha + y*beta) - 
          II*2.0*r2*zR*zR*(x*alpha + y*beta)) - (4.0*PP*r2*zR*zR*(x*alpha + y*beta)*std::tgamma(PP)*
          lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w))/
                   (std::tgamma(1.0 + PP)*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)) + (std::abs(LL)*(zR*zR*w*w)*(x*alpha + y*beta))
       ))/(k*w0*w0*r2*(z2+zR*zR));

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
/// Hx /uLG
//////////////////////////////////////////////////////////////////
Hxfactor = (-II*(

       (-4.0*r2*z*beta*lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w))/((zR*zR*w*w)) + 

       (z*beta*lm_polynomial(P,std::abs(L), 2.0*r2/w/w))/zR*zR + 

       (z*beta*std::abs(LL)*lm_polynomial(P,std::abs(L), 2.0*r2/w/w))/zR*zR - 

       (II*0.5*beta*(-2.0*zR*(II*2.0*r2*z*zR + (1.0 + 2.0*PP)*(zR*zR*w*w)) + 
            k*w0*w0*(pow(x,2)*(z2 - zR*zR) + pow(y,2)*(z2 - zR*zR) + 2.0*pow(z2 + zR*zR,2)) - 
            2.0*w0*w0*zR*(z2 + zR*zR)*std::abs(LL))*lm_polynomial(P,std::abs(L), 2.0*r2/w/w))/
        (w0*w0*zR*zR*(z2 + zR*zR)) + 

        (II*4.0*y*lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w)*
          ((zR*zR*w*w)*(x*alpha + y*beta)*std::abs(LL) - (II*
               (II*(-4.0)*r2*zR*zR*(x*alpha + y*beta)*
                  lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w) + 
                 (LL*(zR*zR*w*w)*(y*alpha - x*beta) + k*w0*w0*r2*z*(x*alpha + y*beta) - 
                    II*2.0*r2*zR*zR*(x*alpha + y*beta))*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/
         lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/(k*pow(w0,4.0)*r2*(z2 + zR*zR)) - 

       (II*((II*LL*x)/(r2) + ((-II)*k*w0*w0*y*z - 2.0*y*zR*zR)/(zR*zR*w*w))*
          lm_polynomial(P,std::abs(L), 2.0*r2/w/w)*
          ((zR*zR*w*w)*(x*alpha + y*beta)*std::abs(LL) - (II*
               (II*(-4.0)*r2*zR*zR*(x*alpha + y*beta)*
                  lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w) + 
                 (LL*(zR*zR*w*w)*(y*alpha - x*beta) + k*w0*w0*r2*z*(x*alpha + y*beta) - 
                    II*2.0*r2*zR*zR*(x*alpha + y*beta))*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/
             lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/(k*w0*w0*r2*zR*zR) - 

       (II*y*std::abs(LL)*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)*
          ((zR*zR*w*w)*(x*alpha + y*beta)*std::abs(LL) - (II*
               (II*(-4.0)*r2*zR*zR*(x*alpha + y*beta)*
                  lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w) + 
                 (LL*(zR*zR*w*w)*(y*alpha - x*beta) + k*w0*w0*r2*z*(x*alpha + y*beta) - 
                    II*2.0*r2*zR*zR*(x*alpha + y*beta))*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/
             lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/(k*w0*w0*pow(r2,2)*zR*zR) - 

       (II*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)*
         (beta*(((-II)*(pow(x,2) + 3.0*pow(y,2))*(k*w0*w0*z - II*2.0*zR*zR) + (zR*zR*w*w)*std::abs(LL))/(r2) -              
                (2.0*y*((-II)*k*w0*w0*y*r2*z + II*LL*w0*w0*x*(z2 + zR*zR) + 
                    y*(-2.0*r2*zR*zR + (zR*zR*w*w)*std::abs(LL))))/pow(r2,2) - 
               (16.0*pow(y,2)*pow(zR,4)*pow(lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w),2))/
                ((zR*zR*w*w)*pow(lm_polynomial(P,std::abs(L), 2.0*r2/w/w),2)) + 
               (16.0*pow(y,2)*pow(zR,4)*lm_polynomial(-2+P,2+std::abs(L), 2.0*r2/w/w))/
                ((zR*zR*w*w)*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)) - 
               (4.0*zR*zR*lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w))/
                lm_polynomial(P,std::abs(L), 2.0*r2/w/w)) + 

            (alpha*(-2.0*pow(w0,4)*x*y*pow(z2 + zR*zR,2)*std::abs(LL) + 
                 (-16.0*x*y*pow(r2,2)*pow(zR,4)*pow(lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w),2) - II*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)*
                     (II*16.0*x*y*pow(r2,2)*pow(zR,4)*
                        lm_polynomial(-2+P,2+std::abs(L), 2.0*r2/w/w) + 
                       LL*pow(w0,4)*(pow(x,2) - pow(y,2))*pow(z2 + zR*zR,2)*
                        lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/
                  pow(lm_polynomial(P,std::abs(L), 2.0*r2/w/w),2)))/
             (w0*w0*pow(r2,2)*(z2 + zR*zR))
           ))/(k*w0*w0*zR*zR)

))/((1.0 + z2/zR*zR)*lm_polynomial(P,std::abs(L), 2.0*r2/w/w))/k/Z;


/// Hy /uLG
//////////////////////////////////////////////////////////////////
Hyfactor = (-II*(
       (4.0*r2*z*alpha*lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w))/((zR*zR*w*w)) - 

       (z*alpha*lm_polynomial(P,std::abs(L), 2.0*r2/w/w))/zR*zR - 

       (z*alpha*std::abs(LL)*lm_polynomial(P,std::abs(L), 2.0*r2/w/w))/zR*zR + 

       (II*0.5*alpha*(-2.0*zR*(II*2.0*r2*z*zR + (1.0 + 2.0*PP)*(zR*zR*w*w)) + 
            k*w0*w0*(pow(x,2)*(z2 - zR*zR) + pow(y,2)*(z2 - zR*zR) + 2.0*pow(z2 + zR*zR,2)) - 
            2.0*w0*w0*zR*(z2 + zR*zR)*std::abs(LL))*lm_polynomial(P,std::abs(L), 2.0*r2/w/w))/(w0*w0*zR*zR*(z2 + zR*zR)) - 

        (II*4.0*x*
          lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w)*
          ((zR*zR*w*w)*(x*alpha + y*beta)*std::abs(LL) - (II*
               (II*(-4.0)*r2*zR*zR*(x*alpha + y*beta)*
                  lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w) + 
                 (LL*(zR*zR*w*w)*(y*alpha - x*beta) + k*w0*w0*r2*z*(x*alpha + y*beta) - 
                    II*2.0*r2*zR*zR*(x*alpha + y*beta))*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/
             lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/(k*pow(w0,4)*r2*(z2 + zR*zR)) + 

       (II*(((-II)*LL*y)/(r2) + ((-II)*k*w0*w0*x*z - 2.0*x*zR*zR)/((zR*zR*w*w)))*
          lm_polynomial(P,std::abs(L), 2.0*r2/w/w)*
          ((zR*zR*w*w)*(x*alpha + y*beta)*std::abs(LL) - (II*
               (II*(-4.0)*r2*zR*zR*(x*alpha + y*beta)*
                  lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w) + 
                 (LL*(zR*zR*w*w)*(y*alpha - x*beta) + k*w0*w0*r2*z*(x*alpha + y*beta) - 
                    II*2.0*r2*zR*zR*(x*alpha + y*beta))*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/
             lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/(k*w0*w0*r2*zR*zR) + 

       (II*x*std::abs(LL)*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)*
          ((zR*zR*w*w)*(x*alpha + y*beta)*std::abs(LL) - (II*
               (II*(-4.0)*r2*zR*zR*(x*alpha + y*beta)*
                  lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w) + 
                 (LL*(zR*zR*w*w)*(y*alpha - x*beta) + k*w0*w0*r2*z*(x*alpha + y*beta) - 
                    II*2.0*r2*zR*zR*(x*alpha + y*beta))*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/
             lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/(k*w0*w0*pow(r2,2)*zR*zR) + 

       (II*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)*
          (alpha*(((-II)*(3.0*r2)*(k*w0*w0*z - II*2.0*zR*zR) + (zR*zR*w*w)*std::abs(LL))/(r2) - 
               (2.0*x*((-II)*k*w0*w0*x*r2*z - 2.0*x*r2*zR*zR - II*LL*w0*w0*y*(z2 + zR*zR) + 
                    w0*w0*x*(z2 + zR*zR)*std::abs(LL)))/pow(r2,2) - 
               (16.0*pow(x,2)*pow(zR,4)*pow(lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w),2))/
                ((zR*zR*w*w)*pow(lm_polynomial(P,std::abs(L), 2.0*r2/w/w),2)) + 
               (16.0*pow(x,2)*pow(zR,4)*lm_polynomial(-2+P,2+std::abs(L),2.0*r2/w/w))/
                ((zR*zR*w*w)*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)) - 
               (4.0*zR*zR*lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w))/
                lm_polynomial(P,std::abs(L), 2.0*r2/w/w)) + 
            (beta*(-2.0*pow(w0,4)*x*y*pow(z2 + zR*zR,2)*std::abs(LL) + 
                 (-16.0*x*y*pow(r2,2)*pow(zR,4)*pow(lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w),2) - II*lm_polynomial(P,std::abs(L), 2.0*r2/w/w)*
                     (II*16.0*x*y*pow(r2,2)*pow(zR,4)*
                        lm_polynomial(-2+P,2+std::abs(L),2.0*r2/w/w) + 
                       LL*pow(w0,4)*(pow(x,2) - pow(y,2))*pow(z2 + zR*zR,2)*
                        lm_polynomial(P,std::abs(L), 2.0*r2/w/w)))/
                  pow(lm_polynomial(P,std::abs(L), 2.0*r2/w/w),2)))/
             (w0*w0*pow(r2,2)*(z2 + zR*zR))))/(k*w0*w0*zR*zR)
))    
           /((1.0 + z2/zR*zR)*lm_polynomial(P,std::abs(L), 2.0*r2/w/w))/k/Z;

/// Hz /uLG
//////////////////////////////////////////////////////////////////
Hzfactor = (II*4.0*r2*zR*zR*(-(y*alpha) + x*beta)*lm_polynomial(-1+P,1+std::abs(L), 2.0*r2/w/w)-(k*w0*w0*r2*z*(-(y*alpha) + x*beta) - II*2.0*r2*zR*zR*(-(y*alpha) + x*beta) + LL*(zR*zR*w*w)*(x*alpha + y*beta))*lm_polynomial(P,std::abs(L), 2.0*r2/w/w) + II*(zR*zR*w*w)*(y*alpha - x*beta)*std::abs(LL)*lm_polynomial(P,std::abs(L), 2.0*r2/w/w))/(w0*w0*r2*(z2 + zR*zR)*lm_polynomial(P,std::abs(L), 2.0*r2/w/w))/k/Z;
  EH[0] = alpha*uLG;
  EH[1] = beta*uLG;
  EH[2] = uLG*Ezfactor;  
  EH[3] = uLG*Hxfactor;
  EH[4] = uLG*Hyfactor;
  EH[5] = uLG*Hzfactor;
}

/**********************************************************************/
/* Gauss-Hermite                          *****************************/
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
//  Purpose:
//    HN_POLYNOMIAL_VALUE evaluates Hn(i,x).
//
//  Discussion:
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
//  Parameters:
//    Input, int M, the number of evaluation points.
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//    Input, double X[M], the evaluation points.
//    Output, double HN_POLYNOMIAL_VALUE[M*(N+1)], the values of the first
//    N+1 Hermite polynomials at the evaluation points.
//****************************************************************************80
//    r8_pi -> M_Pi, Yoonkyung Eunnie Lee, 2015
//  This function is modified to return the highest order polynomial value itself
//  given a single location input (double x).
//****************************************************************************80
{
  double fact, two_power, pval;
  int i, j;
  double* p;
  if ( n < 0 )
  { return NAN; }
  p = new double[n+1];
  p[0] = 1.0;
  if ( n == 0 )
  {  return p[0];}
  p[1] = 2.0*x;
  for (j=2; j<=n;j++)
  { p[j] = 2.0*x*p[j-1] - 2.0*(double)(j-1)*p[j-2];  }
//  Normalize.
  fact = 1.0;
  two_power = 1.0;
  for ( j = 0; j <= n; j++ )
  {  p[j] = p[j] / sqrt(fact*two_power*sqrt(M_PI));
     fact = fact*(double)(j+1);
     two_power = two_power * 2.0;
  }
  pval = p[n];
  delete [] p;
  return pval;
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
//  Modified by Y. Eunnie Lee 2015, to input and output real numbers instead of arrays
//  Parameters:
//    Input, int N, the highest order polynomial to compute.
//    Input, int M, the parameter.  M must be nonnegative.
//    Input, double X, the evaluation points.
//    Output, double LM_POLYNOMIAL[MM*(N+1)], the function values.
{
  int i, j;
  double *v, vval;
  v = new double[n+1];
  if ( n < 0 )
  {    printf("polynomial order n cannot be negative, changing n to 1.\n"); 
       n=0;  
  }
  v = new double[1*(n+1)];
  for ( j = 0; j <= n; j++ )
    {      v[j] = 0.0;  }//initialize 
  v[0] = 1.0; //first term always 1.0 ; 
  if ( n == 0 )
  {    return v[0];  }

  v[1] = (double)(m+1)-x;

  for ( j = 2; j <= n; j++ )
  {   v[j] = (((double)(m+2*j-1)-x)*v[j-1]
             + (double)(-m-j+1)*v[j-2])/(double)(j);
  }
  vval = v[n];
  delete [] v;
  return vval;
}

/**********************************************************************/
/**********************************************************************/
double *hn_polynomial_value ( int m, int n, double x[] )
//****************************************************************************80
//  Purpose:
//    HN_POLYNOMIAL_VALUE evaluates Hn(i,x).
//  Discussion:
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
//  Parameters:
//    Input, int M, the number of evaluation points.
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//    Input, double X[M], the evaluation points.
//    Output, double HN_POLYNOMIAL_VALUE[M*(N+1)], the values of the first
//    N+1 Hermite polynomials at the evaluation points.
//    r8_pi -> M_Pi, Yoonkyung Eunnie Lee, 2015
{
  double fact;
  int i, j;
  double *p;
  double two_power;

  if ( n < 0 )
  { return NULL; }
  p = new double[m*(n+1)];
  for ( i = 0; i < m; i++ )
  {     p[i+0*m] = 1.0;  }
  if ( n == 0 )
  {    return p;  }
  for ( i = 0; i < m; i++ )
  {    p[i+1*m] = 2.0 * x[i];  }
  for ( j = 2; j <= n; j++ )
  {
      for ( i = 0; i < m; i++ )
    {
      p[i+j*m] = 2.0 * x[i] * p[i+(j-1)*m]
        - 2.0 * ( double ) ( j - 1 ) * p[i+(j-2)*m];
    }
  }
//  Normalize.
  fact = 1.0;
  two_power = 1.0;
  for ( j = 0; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      p[i+j*m] = p[i+j*m] / sqrt ( fact * two_power * sqrt ( M_PI ) );
    }
    fact = fact * ( double ) ( j + 1 );
    two_power = two_power * 2.0;
  }
  return p; //this also has memory leaks.
}
double *lm_polynomial ( int mm, int n, int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LM_POLYNOMIAL evaluates Laguerre polynomials Lm(n,m,x).
//
//  First terms:
//
//    M = 0
//
//    Lm(0,0,X) =   1
//    Lm(1,0,X) =  -X   +  1
//    Lm(2,0,X) =   X^2 -  4 X   +  2
//    Lm(3,0,X) =  -X^3 +  9 X^2 -  18 X   +    6
//    Lm(4,0,X) =   X^4 - 16 X^3 +  72 X^2 -   96 X +     24
//    Lm(5,0,X) =  -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x   +  120
//    Lm(6,0,X) =   X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720
//
//    M = 1
//
//    Lm(0,1,X) =    0
//    Lm(1,1,X) =   -1,
//    Lm(2,1,X) =    2 X - 4,
//    Lm(3,1,X) =   -3 X^2 + 18 X - 18,
//    Lm(4,1,X) =    4 X^3 - 48 X^2 + 144 X - 96
//
//    M = 2
//
//    Lm(0,2,X) =    0
//    Lm(1,2,X) =    0,
//    Lm(2,2,X) =    2,
//    Lm(3,2,X) =   -6 X + 18,
//    Lm(4,2,X) =   12 X^2 - 96 X + 144
//
//    M = 3
//
//    Lm(0,3,X) =    0
//    Lm(1,3,X) =    0,
//    Lm(2,3,X) =    0,
//    Lm(3,3,X) =   -6,
//    Lm(4,3,X) =   24 X - 96
//
//    M = 4
//
//    Lm(0,4,X) =    0
//    Lm(1,4,X) =    0
//    Lm(2,4,X) =    0
//    Lm(3,4,X) =    0
//    Lm(4,4,X) =   24
//
//  Recursion:
//
//    Lm(0,M,X)   = 1
//    Lm(1,M,X)   = (M+1-X)
//
//    if 2 <= N:
//
//      Lm(N,M,X)   = ( (M+2*N-1-X) * Lm(N-1,M,X)
//                   +   (1-M-N)    * Lm(N-2,M,X) ) / N
//
//  Special values:
//
//    For M = 0, the associated Laguerre polynomials Lm(N,M,X) are equal
//    to the Laguerre polynomials L(N,X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int MM, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.
//
//    Input, int M, the parameter.  M must be nonnegative.
//
//    Input, double X[MM], the evaluation points.
//
//    Output, double LM_POLYNOMIAL[MM*(N+1)], the function values.
//
{
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    printf("polynomial order n cannot be negative, changing n to 1.\n"); 
    n = 0  ; 
  }

  v = new double[mm*(n+1)];

  for ( j = 0; j <= n; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = 0.0;
    }
  }

  for ( i = 0; i < mm; i++ )
  {
    v[i+0*mm] = 1.0;
  }

  if ( n == 0 )
  {
    return v;
  }

  for ( i = 0; i < mm; i++ )
  {
    v[i+1*mm] = ( double ) ( m + 1 ) - x[i];
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = ( ( ( double ) (   m + 2 * j - 1 ) - x[i] ) * v[i+(j-1)*mm]
                    + ( double ) ( - m     - j + 1 )          * v[i+(j-2)*mm] )
                    / ( double ) (           j     );
    }
  }

  return v;
}
