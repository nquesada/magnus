#include <stdio.h>
#include <math.h>
#include "c_magnus.h"
#include "cubature.h"


#define MAX_EVAL_INT 10000000
#define REQ_ABS_ERROR 1e-8
#define REQ_REL_ERROR 1e-7



int c_magnus2aint(double wa, double waa, double *res){
  const static unsigned ndim=2;
  const static unsigned fdim=2;
  const static double xmin[2]={0.0,-1.0};
  //  xmin[0]=0.0;
  //  xmin[1]=-1.0;
  double xmax[2]={1.0,1.0};
  //  xmax[0]=1.0;
  //  xmax[1]=1.0;
 
  double val[fdim];
  double error[fdim];
  double ws[2]={wa,waa};
  //  ws[0]=wa;
  //  ws[1]=waa;
  const static double maxEval=MAX_EVAL_INT;
  const static double reqAbsError=REQ_ABS_ERROR;
  const static double reqRelError=REQ_REL_ERROR;
  const static double enorm=ERROR_L2;
 

  hcubature(fdim,c_magnus2a,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,enorm,val,error);
  res[0]=val[0];
  res[1]=error[0];
  res[2]=val[1];
  res[3]=error[1];
  return 0;
}



int c_magnus2bint(double wb, double wbb, double *res){
  unsigned ndim=2;
  unsigned fdim=2;
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
  

  hcubature(fdim,c_magnus2b,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,enorm,val,error);
  res[0]=val[0];
  res[1]=error[0];
  res[2]=val[1];
  res[3]=error[1];
  return 0;
}






int c_magnus3int(double wa, double wb, double *res){
  unsigned ndim=4;
  unsigned fdim=2;
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
  

  hcubature(fdim,c_magnus3,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,enorm,val,error);
  //HERE IMPORTANT this is the prefactor of 3 in the definition of the integrals.
  res[0]=val[0];
  res[1]=error[0];
  res[2]=val[1];
  res[3]=error[1];
  return 0;
}


int c_magnus3sint(double wa, double wb, double *res){
  unsigned ndim=2;
  unsigned fdim=3;
  double xmin[ndim];
  xmin[0]=-1.0;
  xmin[1]=-1.0;
  double xmax[ndim];
  xmax[0]=1.0;
  xmax[1]=1.0;
 
  double val[fdim];
  double error[fdim];
  double ws[2];
  ws[0]=wa;
  ws[1]=wb;
  double maxEval=MAX_EVAL_INT;
  double reqAbsError=REQ_ABS_ERROR;
  double reqRelError=REQ_REL_ERROR;
  double enorm=ERROR_L2;
  

  hcubature(fdim,c_magnus3s,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,enorm,val,error);
  res[0]=val[0];
  res[1]=error[0];
  res[2]=val[1];
  res[3]=error[1];
  return 0;
}


int c_magnus3wint(double wa, double wb, double *res){
  unsigned ndim=3;
  unsigned fdim=2;
  double xmin[ndim];
  xmin[0]=0.0;
  xmin[1]=-1.0;
  xmin[2]=-1.0;
  double xmax[ndim];
  xmax[0]=1.0;
  xmax[1]=1.0;
  xmax[2]=1.0;

  double val[fdim];
  double error[fdim];
  double ws[2];
  ws[0]=wa;
  ws[1]=wb;
  double maxEval=MAX_EVAL_INT;
  double reqAbsError=REQ_ABS_ERROR;
  double reqRelError=REQ_REL_ERROR;
  double enorm=ERROR_L2;
  

  hcubature(fdim,c_magnus3w,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,enorm,val,error);
  //HERE IMPORTANT this is the prefactor of 3 in the definition of the integrals.
  res[0]=val[0];
  res[1]=error[0];
  res[2]=val[1];
  res[3]=error[1];
  return 0;
}

