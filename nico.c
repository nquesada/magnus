#include <stdio.h>
#include <math.h>
#include "cubature.h"

#define s0 1.1
#define se 0.9
#define sp 1.0
#define d 1
#define sigma 0.1
#define nu 2


double G(double w0, double we, double wp){
  double tmp1,tmp2;
  tmp1=s0*w0+se*we-sp*wp+d;
  tmp2=sigma*(wp-nu);
  return exp(-tmp1*tmp1-tmp2*tmp2);
}


int integranduv(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
  double u,v,we1,we2;
  u=x[0];
  v=x[1];
  //  fprintf(stdout,"from the function u=%lf v=%lf \n",u,v);
  we1=((double *) fdata)[0];
  we2=((double *) fdata)[1];
  //  fprintf(stdout,"from the function we1=%lf we2=%lf \n",we1,we2);
  fval[0]=G(u+v,we1,u-v+we1)*G(u+v,we2,u-v+we2)/v;
  return 0;
}


int integrandts(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){

  
  if(x[0]>1 || x[0]<-1 || x[1]>1 || x[1]<0){
    fval[0]=0;
    fprintf(stdout,"from the function t=%lf s=%lf \n",x[0],x[1]);
  }
  else{
  
  double u,v,we1,we2;
  u=x[0]/(1-x[0]*x[0]);
  v=x[1]/(1-x[1]);
  //jacobians
  double tmp1,tmp2;
  tmp1=1-x[0]*x[0];
  tmp1=(1+x[0]*x[0])/(tmp1*tmp1);
  tmp2=1-x[1];
  tmp2=1.0/(tmp2*tmp2);


  we1=((double *) fdata)[0];
  we2=((double *) fdata)[1];
  //  fprintf(stdout,"from the function we1=%lf we2=%lf \n",we1,we2);
  //  fprintf(stdout,"from the function jac1=%lf jac2=%lf \n",tmp1,tmp2);
  fval[0]=tmp1*tmp2*(G(u+v,we1,u-v+we1)*G(u+v,we2,u-v+we2)-G(u-v,we1,u+v+we1)*G(u-v,we2,u+v+we2))/v;

  //  fprintf(stdout,"from the function u=%lf v=%lf fval=%lf\n",u,v,fval[0]);
  }
  return 0;
}

int main(){
  double wo,we,wp;
  unsigned ndim=2;
  unsigned fdim=1;

  double ws[2];
  double x[ndim];
  double ff[fdim];
  ws[0]=1.0;
  ws[1]=2.0;
  x[0]=0.5;
  x[1]=0.1;
  
  wo=1.0;
  we=1.0;
  wp=2.0;
  //  fprintf(stdout,"%lf \n",F(wo,we,wp));
  integrandts(ndim,x,ws,fdim,ff);
  //  fprintf(stdout,"%lf \n",ff[0]);

  int res;
  res=integrandts(ndim,x,ws,fdim,ff);

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

  hcubature(fdim,integrandts,ws,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,norm,val,error);

  fprintf(stdout,"val %.16e +/- %.16e\n",val[0],error[0]);
  return 0;
}
