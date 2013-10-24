#include <stdio.h>
#include <math.h>
#include "magnus.h"
#include "cubature.h"


int main(){
  double wo,we,wp;
  unsigned ndim=2;
  unsigned fdim=1;

  double ws[2];
  double x[ndim];
  //  double ff[fdim];
  ws[0]=1.0;
  ws[1]=2.0;
  x[0]=0.5;
  x[1]=0.1;
  
  wo=1.0;
  we=1.0;
  wp=2.0;

  double xmin[ndim];
  xmin[0]=-1.0;
  xmin[1]=0;
  double xmax[ndim];
  xmax[0]=1.0;
  xmax[1]=1.0;
  int maxEval=0;
  double reqAbsError=0.001;
  double reqRelError=0.001;
  error_norm norm=ERROR_L2;
  double val[fdim];
  double error[fdim];

  error[0]=0.0;    
  val[0]=0.0;

  hcubature(fdim,magnus2a,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,norm,val,error);

  fprintf(stdout,"val %.16e +/- %.16e\n",val[0],error[0]);
  return 0;
}
