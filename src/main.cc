/*---------------------------------------------------------------
 * main.cc  
 * runs the nlopt routine using objective function 
 * created 2015.05.21
 * last updated on 2015.05.22
 *--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <nlopt.h>
#include <complex>
#include <fstream>
#include <iostream>
#include <cstdlib>
#define II cdouble(0.0,1.0)
#define MAXSTR 100 
using namespace std;
// strategy 
// main function will be calling all other functions. 
// objective function gets a vector of size n. 
double objective(unsigned n, const double *x, double *grad) 
//---------------------------------------------------------------//
double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
  //Owen comment on const "const correctness"  for input. notforscuff
    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}
//---------------------------------------------------------------//
typedef struct {
    double a, b;
} my_constraint_data;
//---------------------------------------------------------------//
double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
    my_constraint_data *d = (my_constraint_data *) data;
    double a = d->a, b = d->b;
    if (grad) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
 }
//---------------------------------------------------------------//
int main(int argc, char *argv[])
  /// Objective 
  /// 
  /// READ In:
  ///    parameter matrix PARMMatrix 
  ///    BEM matrix A, LU Factorized
  ///    PFT matrix Q 
{
  /*--------------------------------------------------------------*/
  /*- process options  -------------------------------------------*/
  /*--------------------------------------------------------------*/
  /// Input arguments should be filenames. (pointers) 
  char *GeoFile=0; 
  char *PARMMatrixFile=0;
  char *BEMMatrixFile=0;
  char *PFTMatrixFile=0;
  char *LogLevel=0;
  /* name    type   #args  max_instances  storage  count  description*/
  OptStruct OSArray[]=
    { 
      {"geometry",   PA_STRING,  1, 1, (void *)&GeoFile, 0,  ".scuffgeo file"},
      {"PARMMatrix", PA_STRING,  1, 1, (void *)&PARMMatrixFile, 0, "VParameters file"},
      {"BEMMatrix",  PA_STRING,  1, 1, (void *)&BEMMatrixFile, 0,  "BEM Matrix datafile"},
      {"PFTMatrix",  PA_STRING,  1, 1, (void *)&PFTMatrixFile, 0,  "PFT Matrix datafile"},
      {"LogLevel",   PA_STRING,  1, 1, (void *)&LogLevel, 0, "none | terse | verbose | verbose2"},
      {0,0,0,0,0,0,0}
    };
  ProcessOptions(argc, argv, OSArray);
  // if (GeoFile==0)
  //   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  // if (PARMMatrixFile==0)
  //   OSUsage(argv[0],OSArray,"--PARMMatrix option is mandatory");
  // if (BEMMatrixFile==0)
  //  OSUsage(argv[0],OSArray,"--BEMMatrix option is mandatory");
  RWGGeometry *G = new RWGGeometry(GeoFile);
  HMatrix *PARMMatrix = new HMatrix(PARMMatrixFile, LHM_TEXT);
  ShowPARMMatirx(PARMMatrix); 
  HMatrix *M = new HMatrix(BEMMatrixFile, LHM_TEXT);
  //A->LUFactorize(); //(already done)
  HMatrix *Q = new HMatrix(PFTMatrixFile, LHM_TEXT);
  objective(G, PARMMatrix, M, Q); 
  
  delete G, PARMMatrix, M, Q; 
}
 //  nlopt_set_upper_bounds; 
 //  my_constraint_data data[2] = {{},{}}; 
 //  nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-8);
 //  nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-8);
