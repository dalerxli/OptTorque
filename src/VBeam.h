/*
 * VBeam.h   -- definition for the VectorBeam implementation
 *            -- of the IncField class
 * v14
 */
#ifndef VBEAM_H
#define VBEAM_H
#define LMAX 5
#define LSPAN 2*LMAX+1

#include <libhrutil.h>
#include <libIncField.h>

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
 class VBeam: public IncField
 {
   ///////////////////////////
   // class for VectorBeam construction ( allows superposed beam input )
   // VBeam::GetMN returns the orthogonal vector fields M[L] and N[L]  
   // VBeam::GetFields constructs EH for one or more modes with index L
   // where EH is constructed as a weighted sum of M[L] and N[L]. 
   ///////////////////////////
   // Bessel Beams: GetMN
   // L modes with -LMAX<L<LMAX 
   // Default aperture angle is 0.5degree
   ///////////////////////////
 public:
   int L; 
   double aL[LSPAN], bL[LSPAN]; // array of coefficients 
   double aIn;        // beam opening angle in degrees (Numerical Aperture)
   double Cxyz[3];    // beam center position [um]
   double nHat[3];    // propagation direction
   cdouble M[3],N[3]; // Orthogonal Vector Fields
   ///////////////////////////////////////////////////////////
   VBeam(int NewL, double NewaIn);
   VBeam(double NewaL[LSPAN], double NewbL[LSPAN],double NewaIn);
   ~VBeam();
   ///////////////////////////////////////////////////////////
   void SetL(int NewL);
   void SetaIn(double NewaIn);
   void Setab(double NewaL[LSPAN], double NewbL[LSPAN]);
   void SetCxyz(double NewCxyz[3]);
   void SetnHat(double NewnHat[3]);
   ///////////////////////////////////////////////////////////
   void GetFields(const double X[3], cdouble EH[6]);
   void GetMN(const double X[3], int L, double aIn, 
              cdouble M[3], cdouble N[3]);         
 };
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
#endif // #ifndef VBEAM_H
