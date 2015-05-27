/*---------------------------------------------------------------
 * main.cc  
 * contains examples of NLOPT in C 
 * created 2015.05.21
 * last updated on 2015.05.27
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
{

 //  nlopt_set_upper_bounds; 
 //  my_constraint_data data[2] = {{},{}}; 
 //  nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-8);
 //  nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-8);

}
