/*
 * VBeam.h   -- definition for the VectorBeam implementation
 *            -- of the IncField class
 * v15
 */
#ifndef VBEAM_H
#define VBEAM_H
#define LMAX 25
#include <libhrutil.h>
#include <libhmat.h>
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
 public:
   HMatrix *PMatrix=new HMatrix(LMAX,4,LHM_REAL);  // Parameter matrix 
   int numL; //number of modes summed in PMatrix (nonzero mode number);

   double Cxyz[3];    // beam center position [um]
   double nHat[3];    // propagation direction
   cdouble M[3],N[3]; // Orthogonal Vector Fields

   VBeam(int NewL, double NewaIn);
   VBeam(HMatrix *NewPMatrix);
   ~VBeam();

   void SetCxyz(double NewCxyz[3]);
   void SetnHat(double NewnHat[3]);

   void GetFields(const double X[3], cdouble EH[6]);
   void GetMN(const double X[3], int L, double aIn, 
              cdouble M[3], cdouble N[3]);         
 };
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
#endif // #ifndef VBEAM_H
