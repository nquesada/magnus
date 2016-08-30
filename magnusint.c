#include <stdio.h>
#include <math.h>
#include "magnus.h"
#include "cubature.h"


#define MAX_EVAL_INT 10000000
#define REQ_ABS_ERROR 1e-8
#define REQ_REL_ERROR 1e-7



int magnus2aint(double wa, double waa, double *res){
  const unsigned ndim=2;
  const unsigned fdim=1;
  const double xmin[2]={0.0,-1.0};
  double xmax[2]={1,1};
 
  double val[fdim];
  double error[fdim];
  double ws[2]={wa,waa};

  const double maxEval=MAX_EVAL_INT;
  const double reqAbsError=REQ_ABS_ERROR;
  const double reqRelError=REQ_REL_ERROR;
 

  hcubature(fdim,magnus2a,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,ERROR_L2,val,error);
  res[0]=val[0];
  res[1]=error[0];
  return 0;
}



int magnus2bint(double wb, double wbb, double *res){
  const unsigned ndim=2;
  const unsigned fdim=1;
  const double xmin[2]={0.0,-1.0};
  double xmax[2]={1.0,1.0};
 
  double val[fdim];
  double error[fdim];
  double ws[2]={wb,wbb};

  const double maxEval=MAX_EVAL_INT;
  const double reqAbsError=REQ_ABS_ERROR;
  const double reqRelError=REQ_REL_ERROR;
  

  hcubature(fdim,magnus2b,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,ERROR_L2,val,error);
  res[0]=val[0];
  res[1]=error[0];
  return 0;
}






int magnus3int(double wa, double wb, double *res){
  const unsigned ndim=4;
  const unsigned fdim=1;
  const double xmin[4]={0.0,0.0,-1.0,-1.0};
  const double xmax[4]={1.0,1.0,1.0,1.0};

  double val[fdim];
  double error[fdim];
  double ws[2]={wa,wb};

  const double maxEval=MAX_EVAL_INT;
  const double reqAbsError=REQ_ABS_ERROR;
  const double reqRelError=REQ_REL_ERROR;
  

  hcubature(fdim,magnus3,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,ERROR_L2,val,error);
  //HERE IMPORTANT this is the prefactor of 3 in the definition of the integrals.
  res[0]=val[0];
  res[1]=error[0];
  return 0;
}


int magnus3sint(double wa, double wb, double *res){
  const unsigned ndim=2;
  const unsigned fdim=1;
  const double xmin[2]={-1.0,-1.0};
  double xmax[2]={1.0,1.0};
 
  double val[fdim];
  double error[fdim];
  double ws[2]={wa,wb};

  const double maxEval=MAX_EVAL_INT;
  const double reqAbsError=REQ_ABS_ERROR;
  const double reqRelError=REQ_REL_ERROR;
  

  hcubature(fdim,magnus3s,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,ERROR_L2,val,error);
  res[0]=val[0];
  res[1]=error[0];
  return 0;
}


int magnus3wint(double wa, double wb, double *res){
  const unsigned ndim=3;
  const unsigned fdim=1;
  double xmin[3]={0.0,-1.0,-1.0};
  double xmax[3]={1.0,1.0,1.0};

  double val[fdim];
  double error[fdim];
  double ws[2]={wa,wb};

  const double maxEval=MAX_EVAL_INT;
  const double reqAbsError=REQ_ABS_ERROR;
  const double reqRelError=REQ_REL_ERROR;
  

  hcubature(fdim,magnus3w,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,ERROR_L2,val,error);
  //HERE IMPORTANT this is the prefactor of 3 in the definition of the integrals.
  res[0]=val[0];
  res[1]=error[0];
  return 0;
}
