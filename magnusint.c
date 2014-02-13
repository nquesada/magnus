#include <stdio.h>
#include <math.h>
#include "magnus.h"
#include "cubature.h"


#define MAX_EVAL_INT 10000000
#define REQ_ABS_ERROR 1e-8
#define REQ_REL_ERROR 1e-4



int magnus2aint(double wa, double waa, double *res){
  unsigned ndim=2;
  unsigned fdim=1;
  double xmin[ndim];
  xmin[0]=0.0;
  xmin[1]=-1.0;
  double xmax[ndim];
  xmax[0]=1.0;
  xmax[1]=1.0;
 
  double val[fdim];
  double error[fdim];
  double ws[2];
  ws[0]=wa;
  ws[1]=waa;
  double maxEval=MAX_EVAL_INT;
  double reqAbsError=REQ_ABS_ERROR;
  double reqRelError=REQ_REL_ERROR;
  double enorm=ERROR_L2;
  

  hcubature(fdim,magnus2a,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,enorm,val,error);
  res[0]=val[0];
  res[1]=error[0];
  return 0;
}



int magnus2bint(double wb, double wbb, double *res){
  unsigned ndim=2;
  unsigned fdim=1;
  double xmin[ndim];
  xmin[0]=0.0;
  xmin[1]=-1.0;
  double xmax[ndim];
  xmax[0]=1.0;
  xmax[1]=1.0;
 
  double val[fdim];
  double error[fdim];
  double ws[2];
  ws[0]=wb;
  ws[1]=wbb;
  double maxEval=MAX_EVAL_INT;
  double reqAbsError=REQ_ABS_ERROR;
  double reqRelError=REQ_REL_ERROR;
  double enorm=ERROR_L2;
  

  hcubature(fdim,magnus2b,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,enorm,val,error);
  res[0]=val[0];
  res[1]=error[0];
  return 0;
}






int magnus3int(double wa, double wb, double *res){
  unsigned ndim=4;
  unsigned fdim=1;
  double xmin[ndim];
  xmin[0]=0.0;
  xmin[1]=0.0;
  xmin[2]=-1.0;
  xmin[3]=-1.0;
  double xmax[ndim];
  xmax[0]=1.0;
  xmax[1]=1.0;
  xmax[2]=1.0;
  xmax[3]=1.0;

  double val[fdim];
  double error[fdim];
  double ws[2];
  ws[0]=wa;
  ws[1]=wb;
  double maxEval=MAX_EVAL_INT;
  double reqAbsError=REQ_ABS_ERROR;
  double reqRelError=REQ_REL_ERROR;
  double enorm=ERROR_L2;
  

  hcubature(fdim,magnus3,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,enorm,val,error);
  //HERE IMPORTANT this is the prefactor of 3 in the definition of the integrals.
  res[0]=val[0]/3.0;
  res[1]=error[0]/3.0;
  return 0;
}
