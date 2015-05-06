//--------------------------------------------------------------------//
// VBeamTest_v15.cc 
//--------------------------------------------------------------------//
/// Test PARMMatrix input 
/// Verify VBeam Profile. 
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <complex>
#include <cstdlib> //for abs
#include <string.h>//for memcpy
#include "VBeam.h"
  #ifndef cdouble
    typedef std::complex<double> cdouble;
  #endif
  #define Eps         cdouble(1.0,0.0)
  #define Mu          cdouble(1.0,0.0)
  #define Omega       cdouble(5.0,0.0)
  #define ZVAC        376.73031346177
  #define II cdouble(0.0,1.0)

//--------------------------------------------------------------------//
int main(int argc, char *argv[]){
  printf("Testing Bessel Beam Implementation for v15\n");
  char *PARMMatrixFile=0;
  OptStruct OSArray[]=
   { 
     /* name    type   #args  max_instances  storage  count  description*/
     {"PARMMatrix", PA_STRING,  1, 1, (void *)&PARMMatrixFile, 0, "VParameters file"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (PARMMatrixFile==0)
   OSUsage(argv[0],OSArray,"--PARMMatrix option is mandatory");

  HMatrix *PARMMatrix = new HMatrix(PARMMatrixFile, LHM_TEXT);
  VBeam *VB=new VBeam(PARMMatrix); 

////////////////////////////////////////////////////////////////////////
printf("input three real values for x y z [microns]\n");  
double X[3];
std::cin >>  X[0];
std::cin >>  X[1];
std::cin >>  X[2];
////////////////////////////////////////////////////////////////////////
printf("Test Bessel Function of C\n");
double besP= jn(0,X[0]); 
printf("jn(0,x)=%e\n",besP);
////////////////////////////////////////////////////////////////////////
cdouble EH[6];
 VB->VBeam::GetFields(X,EH); 

printf("new     EX= (%3.3f,%3.3f),  EY= (%3.3f,%3.3f),  EZ= (%3.3f,%3.3f)\n",std::real(EH[0]),std::imag(EH[0]),std::real(EH[1]),std::imag(EH[1]),std::real(EH[2]),std::imag(EH[2]));
printf("new     HX= (%3.3e,%3.3e),  HY= (%3.3e,%3.3e),  HZ= (%3.3e,%3.3e)\n",std::real(EH[3]),std::imag(EH[3]),std::real(EH[4]),std::imag(EH[4]),std::real(EH[5]),std::imag(EH[5]));

printf(" |E|=%e, |H|=%e \n",sqrt(pow(std::abs(EH[0]),2)+pow(std::abs(EH[1]),2)+pow(std::abs(EH[2]),2)),sqrt(pow(std::abs(EH[3]),2)+pow(std::abs(EH[4]),2)+pow(std::abs(EH[5]),2))); 
return 0;
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
