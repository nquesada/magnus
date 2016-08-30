#include"functionF.h"
#include<stdio.h>


int c_magnus2a(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
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

  double f1[2], f1b[2], f2[2], f2b[2];

  c_F(a,q,a + q + p,f1);
  c_F(a,q,a + q - p,f1b);
  c_F(aa,q,aa + q + p,f2);
  c_F(aa,q,aa + q - p,f2b);
  
  tmp0=tmp0*tmp1/p;

  fval[0]=tmp0*(f1[0]*f2[0] + f1[1]*f2[1] - f1b[0]*f2b[0] - f1b[1]*f2b[1]);
  fval[1]=tmp0*(f1[1]*f2[0] - f1[0]*f2[1] - f1b[1]*f2b[0] + f1b[0]*f2b[1]);
  return 0;
}




int c_magnus2b(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
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

  double f1[2], f1b[2], f2[2], f2b[2];
  
  c_F(q,b,b + p + q,f1);
  c_F(q,b,b - p + q,f1b);
  c_F(q,bb,bb + p + q,f2);
  c_F(q,bb,bb - p + q,f2b);

  tmp0=tmp0*tmp1/p;

  fval[0]=tmp0*(f1[0]*f2[0] + f1[1]*f2[1] - f1b[0]*f2b[0] - f1b[1]*f2b[1]);
  fval[1]=tmp0*(f1[1]*f2[0] - f1[0]*f2[1] - f1b[1]*f2b[0] + f1b[0]*f2b[1]);

  return 0;
}



int c_magnus3(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
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

  double f1[2],f1p[2],f1q[2],f1pq[2],f2[2],f2p[2],f3[2],f3q[2];
  
  c_F(c,d,c + d + p + q,f1);
  c_F(c,d,c + d - p + q,f1p);
  c_F(c,d,c + d + p - q,f1q);
  c_F(c,d,c + d - p - q,f1pq);
  c_F(a,d,a + d + p,f2);
  c_F(a,d,a + d - p,f2p);
  c_F(c,b,b + c + q,f3);
  c_F(c,b,b + c - q,f3q);

  tmp1=tmp1*tmp2*tmp3*tmp4/(p*q);

  //The next line was generated automatically from mathematica, do not touch!
  fval[0]=tmp1*(f1[0]*f2[0]*f3[0] + f1[1]*f2[1]*f3[0] - f1p[0]*f2p[0]*f3[0] - 
	      f1p[1]*f2p[1]*f3[0] + f1[1]*f2[0]*f3[1] - f1[0]*f2[1]*f3[1] - 
	      f1p[1]*f2p[0]*f3[1] + f1p[0]*f2p[1]*f3[1] - f1q[0]*f2[0]*f3q[0] - 
	      f1q[1]*f2[1]*f3q[0] + f1pq[0]*f2p[0]*f3q[0] + f1pq[1]*f2p[1]*f3q[0] - 
	      f1q[1]*f2[0]*f3q[1] + f1q[0]*f2[1]*f3q[1] + f1pq[1]*f2p[0]*f3q[1] - 
	      f1pq[0]*f2p[1]*f3q[1]);
  fval[1]=tmp1*(-f1[1]*f2[0]*f3[0] + f1[0]*f2[1]*f3[0] + 
	      f1p[1]*f2p[0]*f3[0] - f1p[0]*f2p[1]*f3[0] + f1[0]*f2[0]*f3[1] + 
	      f1[1]*f2[1]*f3[1] - f1p[0]*f2p[0]*f3[1] - f1p[1]*f2p[1]*f3[1] + 
	      f1q[1]*f2[0]*f3q[0] - f1q[0]*f2[1]*f3q[0] - f1pq[1]*f2p[0]*f3q[0] + 
	      f1pq[0]*f2p[1]*f3q[0] - f1q[0]*f2[0]*f3q[1] - f1q[1]*f2[1]*f3q[1] + 
	      f1pq[0]*f2p[0]*f3q[1] + f1pq[1]*f2p[1]*f3q[1]);

  return 0;
}


int c_magnus3s(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
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

  double f1[2],f2[2],f3[2];
  c_F(p,q,p + q,f1);
  c_F(a,q,a + q,f2);
  c_F(p,b,b + p,f3);

  tmp0=tmp0*tmp1;


  fval[0]=tmp0*(f1[0]*f2[0]*f3[0] + f1[1]*f2[1]*f3[0] + f1[1]*f2[0]*f3[1] - f1[0]*f2[1]*f3[1]);
  fval[1]=tmp0*(-f1[1]*f2[0]*f3[0] + f1[0]*f2[1]*f3[0] + f1[0]*f2[0]*f3[1] + f1[1]*f2[1]*f3[1]);
  return 0;
}

/*

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
  //  fprintf(stdout,"%lf %lf %lf %lf %lf\n",x[0],x[1],x[2],x[3],fval[0]);
  return 0;
}

*/
