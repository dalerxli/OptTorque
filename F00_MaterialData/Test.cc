#include <stdio.h>
#include <cmath>

int main(int argc, char *argv[])
{
  int i=0; 
  printf("\n cmdline args count=%d", argc); 
  /* First argument is executable name only */
  printf("\nexe name=%s", argv[0]);

  for (i=1; i<argc; i++)
    {
      printf("\n arg%d=%s",i,argv[i]); 
    }
  printf("\n"); 
  return 0; 
}

/*
double* MakeOmegaFile(lambdamin, lambdamax, numpts, OmegaFileName)

    % Omegafile : omega is in the order of 2*pi*10^(-6)/lambda.
    lambda = linspace(lambdamin, lambdamax, numpts) 
    omega = 2*pi*1000./lambda;

  dlmsave(OmegaFileName, omega); 
end
*/
