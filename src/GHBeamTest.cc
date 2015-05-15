/*
 * GHBeamTest.cc   -- Tests GHBeam Result for a given X[3] and Omega 
 *
 */
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <complex>
#include <cstdlib> //for abs
#include <string.h>//for memcpy
#include "GHBeam.h"
  #define Eps         cdouble(1.0,0.0)
  #define Mu          cdouble(1.0,0.0)
  #define Omega       cdouble(10.471975511965978,0.0)
  #define ZVAC        376.73031346177
  #define II cdouble(0.0,1.0)
//--------------------------------------------------------------------//
  /// Test modified polynomial fcns in GHBeam.cc against old functions 
  /// added 2015.05.14
double hn_polynomial_value ( int n, double x );
double lm_polynomial (int n, int m, double x);

double hn_polynomial_value_old( int n, double x );
double lm_polynomial_old(int n, int m, double x);
//--------------------------------------------------------------------//
int main()
{
  printf("input three real values for x y z [microns]\n");  
  double X[3];
  std::cin >>  X[0] >> X[1] >>  X[2];
  //--------------------------------------------------------------------//
  printf("Test Hermite-Gaussian Function\n");
  for(int jj=0; jj<5; jj++){
    int n = jj; 
    double x = X[0]; 
    printf("known = %e, computed = %e\n",hn_polynomial_value_old(n,x),
           hn_polynomial_value(n,x));  
  }
  //--------------------------------------------------------------------//
  printf("Test Laguerre-Gaussian Function\n");
  for(int jj=0; jj<5; jj++){
    int n = jj; 
    int m = jj; 
    double x = X[0]; 
    printf("known = %e, computed = %e\n",lm_polynomial_old(n,m,x),
           lm_polynomial(n,m,x));  
  }
  //--------------------------------------------------------------------//
  printf("GLBeam test\n");
  int P,L;
  printf("input two integers for P and L\n");  
  std::cin >> P >> L; 
  double w0 = 2.0;                    // beam waist radius at z=0
  const char *polname = "radial";

  cdouble EH[6];
  GLBeam *GLB=new GLBeam(P,L,w0,1.0,polname); 
 
  GLB->GLBeam::GetFields(X,EH); 
  printf("Ex Ey Ez: %e %e %e %e %e %e\n",
         std::real(EH[0]),std::imag(EH[0]),
         std::real(EH[1]),std::imag(EH[1]),
         std::real(EH[2]),std::imag(EH[2]));
  printf("Hx Hy Hz: %e %e %e %e %e %e\n",
         std::real(EH[3]),std::imag(EH[3]),
         std::real(EH[4]),std::imag(EH[4]),
         std::real(EH[5]),std::imag(EH[5]));
  printf(" |E|=%e, |H|=%e \n",
           sqrt(pow(std::abs(EH[0]),2)+pow(std::abs(EH[1]),2)+pow(std::abs(EH[2]),2)),
           sqrt(pow(std::abs(EH[3]),2)+pow(std::abs(EH[4]),2)+pow(std::abs(EH[5]),2))
           ); 
  return 0;
}

/**********************************************************************/
// Old Polynomials for comparison 
/**********************************************************************/
double hn_polynomial_value_old ( int n, double x )
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

double lm_polynomial_old( int n, int m, double x)
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
