/*
 * GHBeam.h   -- definition for the GHBeam implementation
 *            -- of the IncField class
 * v13
 */

#ifndef GHBEAM_H
#define GHBEAM_H
#include <libhrutil.h>
#include <libIncField.h>
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
 class GLBeam: public IncField
 {
 public:
   int L, P;          // Gauss-Hermite beam indices
   double w0;         // waist radius
   double I0; 
   cdouble alpha; 
   cdouble beta; 
   const char *polname;     // polarization name. maximum size 30
   const char *pollist[6] = {"xp","yp","rcp","lcp","radial","azimuthal"};
   double Cxyz[3];    // beam center position [um]
   double nHat[3];    // propagation direction

   ///////////////////////////////////////////////////////////
   GLBeam(int NewL, int NewP, double Neww0, double NewI0, const char *Newpolname);
   GLBeam(int NewL, int NewP, double Neww0, double NewI0, const char *Newpolname,
          double NewCxyz[3], double NewnHat[3]);
   ~GLBeam();
   ///////////////////////////////////////////////////////////
   void SetL(int NewL);
   void SetP(int NewP);
   void Setw0(double Neww0);
   void SetI0(double NewI0); 
   void Setpol(const char *Newpolname);
   void Setpol(const char *Newpolname, cdouble Newalpha, cdouble Newbeta); 
   void SetCxyz(double NewCxyz[3]);
   void SetnHat(double NewnHat[3]);
   ///////////////////////////////////////////////////////////
   void GetFields(const double X[3], cdouble EH[6]);
   void GetLG(const double X[3],int P,int L,double w0, 
              const char *polname, cdouble EH[6]);         
 };
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
class GHBeam: public IncField
 {
 public:
   int HM, HN;        // Gauss-Hermite beam indices
   double w0;         // waist radius
   double Cxyz[3];    // center coordinate [um]
   cdouble E0[3];     // complex polarization vector
   double nHat[3];    // propagation direction

   GHBeam(int NewHM, int NewHN);
   GHBeam(int NewHM, int NewHN, double Neww0);
   GHBeam(int NewHM, int NewHN, double Neww0, double NewCxyz[3],cdouble NewE0[3], double NewnHat[3]);

   void SetHM(int NewHM);
   void SetHN(int NewHN);
   void Setw0(double Neww0);
   void SetCxyz(double NewCxyz[3]);
   void SetE0(cdouble NewE0[3]);
   void SetnHat(double NewnHat[3]);
   void GetFields(const double X[3], cdouble EH[6]);
 };
#endif // #ifndef GHBEAM_H
