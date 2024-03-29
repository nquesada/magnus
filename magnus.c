#include"functionF.h"
#include<stdio.h>

int magnus2a(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
  double a,aa,p,q;
  p=x[0]/(1-x[0]);
  q=x[1]/(1-x[1]*x[1]);
  //jacobians
  double tmp0,tmp1;
  tmp0=1-x[0];
  tmp0=1.0/(tmp0*tmp0);
  tmp1=1-x[1]*x[1];
  tmp1=(1+x[1]*x[1])/(tmp1*tmp1);
  a=((double *) fdata)[0];
  aa=((double *) fdata)[1];

  fval[0]=tmp0*tmp1*(-(F(a,q,a - p + q)*F(aa,q,aa - p + q)) + 
		    F(a,q,a + p + q)*F(aa,q,aa + p + q))/p;

  return 0;
}




int magnus2b(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
  double b,bb,p,q;
  p=x[0]/(1-x[0]);
  q=x[1]/(1-x[1]*x[1]);
  //jacobians
  double tmp0,tmp1;
  tmp0=1-x[0];
  tmp0=1.0/(tmp0*tmp0);
  tmp1=1-x[1]*x[1];
  tmp1=(1+x[1]*x[1])/(tmp1*tmp1);
  b=((double *) fdata)[0];
  bb=((double *) fdata)[1];

  fval[0]=tmp0*tmp1*(-(F(q,b,b - p + q)*F(q,bb,bb - p + q)) + 
		     F(q,b,b + p + q)*F(q,bb,bb + p + q))/p;

  return 0;
}




int magnus3(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
  double a,b,c,d,p,q;
  p=x[0]/(1-x[0]);
  q=x[1]/(1-x[1]);
  c=x[2]/(1-x[2]*x[2]);
  d=x[3]/(1-x[3]*x[3]);
  //jacobians
  double tmp1,tmp2,tmp3,tmp4;
  tmp1=1-x[0];
  tmp1=1.0/(tmp1*tmp1);
  tmp2=1-x[1];
  tmp2=1.0/(tmp2*tmp2);
  tmp3=1-x[2]*x[2];
  tmp3=(1+x[2]*x[2])/(tmp3*tmp3);
  tmp4=1-x[3]*x[3];
  tmp4=(1+x[3]*x[3])/(tmp4*tmp4);


  a=((double *) fdata)[0];
  b=((double *) fdata)[1];
  //The next line was generated automatically from mathematica, do not touch!
  fval[0]=tmp1*tmp2*tmp3*tmp4*((F(a,d,a + d - p)*(F(c,b,b + c - q)*F(c,d,c + d - p - q) - 
						  F(c,b,b + c + q)*F(c,d,c + d - p + q)) + 
				F(a,d,a + d + p)*(-(F(c,b,b + c - q)*F(c,d,c + d + p - q)) + 
						  F(c,b,b + c + q)*F(c,d,c + d + p + q)))/(p*q));
  return 0;
}



int magnus3s(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
  double a,b,p,q;
  p=x[0]/(1-x[0]*x[0]);
  q=x[1]/(1-x[1]*x[1]);
  //jacobians
  double tmp0,tmp1;
  
  tmp0=1-x[0]*x[0];
  tmp0=(1+x[0]*x[0])/(tmp0*tmp0);

  tmp1=1-x[1]*x[1];
  tmp1=(1+x[1]*x[1])/(tmp1*tmp1);
  
  a=((double *) fdata)[0];
  b=((double *) fdata)[1];

  fval[0]=tmp0*tmp1*(F(p,q,p+q)*F(a,q,a+q)*F(p,b,p+b));

  return 0;
}



int magnus3w(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
  double a,b,w,p,q;
  p=x[0]/(1-x[0]);
  q=x[1]/(1-x[1]*x[1]);
  w=x[2]/(1-x[2]*x[2]);
  //jacobians
  double tmp1,tmp2,tmp3;
  tmp1=1-x[0];
  tmp1=1.0/(tmp1*tmp1);
  tmp2=1-x[1]*x[1];
  tmp2=(1+x[1]*x[1])/(tmp2*tmp2);
  tmp3=1-x[2]*x[2];
  tmp3=(1+x[2]*x[2])/(tmp3*tmp3);


  a=((double *) fdata)[0];
  b=((double *) fdata)[1];
  //The next line was generated automatically from mathematica, do not touch!
  fval[0]=tmp1*tmp2*tmp3*(-(F(a,w,a + w)*F(q,b,b - p + q)*
			    F(q,w,-p + q + w)) + 
			  F(a,w,a + w)*F(q,b,b + p + q)*
			  F(q,w,p + q + w) - 
			  F(a,q,a - p + q)*F(w,b,b + w)*
			  F(w,q,-p + q + w) + 
			  F(a,q,a + p + q)*F(w,b,b + w)*
			  F(w,q,p + q + w))/p;
  return 0;
}

