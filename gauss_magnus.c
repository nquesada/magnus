#include<stdio.h>
#include<math.h>
#include <gsl/gsl_math.h>
#include "cubature.h"


/* Values for cubature */
#define MAX_EVAL_INT 10000000 //Maximum number of evaluations
#define REQ_ABS_ERROR 1e-10 //Required absolute error
#define REQ_REL_ERROR 1e-8 //Required relative error

/* Value used to truncate J_3 using the bound for J_3
obtained in the paper. If the bound on J_3 is less
eps the J_3 is not calculated explicitly and it is 
assumed to be zero. eps should be at least as small
as REQ_ABS_ERROR */
#define eps 1e-8


/* This is simply 4/sqrt(3) to 20 decimal places*/
#define M_4_S_3 2.3094010767585030580


/* Values that define the problem */
#define sa  34.0
#define sb  36.0
#define sc  40.0
#define s   1.0 //\sigma

/* derived quantitities, do not modify! */

#define ea  sc-sa //\eta_a
#define eb  sc-sb //\eta_b


#define eab sb-sa  //\eta_{ab}
#define eba sa-sb  //\eta_{ba}


#define c2  (s*s+(sc-sa)*(sc-sb)) //\chi^2
#define ca2 (s*s+(sc-sa)*(sc-sa)) //\chi_a^2
#define cb2 (s*s+(sc-sb)*(sc-sb)) //\chi_b^2
#define R4  (4*(s*s+(sc-sa)*(sc-sa))*(s*s+(sc-sb)*(sc-sb))-(s*s+(sc-sa)*(sc-sb))*(s*s+(sc-sa)*(sc-sb))) //R^4
#define M4  (4*(s*s+(sc-sa)*(sc-sa))*(s*s+(sc-sb)*(sc-sb))-3*(s*s+(sc-sa)*(sc-sb))*(s*s+(sc-sa)*(sc-sb))) //M^4

//This are the elements of the different matrices in the paper
// Xa indicates the element (1,1), Xb twice the elemen (1,2)=(2,1)
// Xc indicate the element 2,2

#define Na (ca2)
#define Nb (2*c2)
#define Nc (cb2)
#define Q1a (M4*ca2/R4)
#define Q1b (2*c2*c2*c2/R4)
#define Q1c (M4*cb2/R4)
#define Q3a (2*ca2)
#define Q3b (2*c2)
#define Q3c (2*cb2)

/* This quantity defines the prefactor inside the cosine */
#define v M_4_S_3*s*eab





int gaussianpvarg(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
  double p,q;
  p=x[0]/(1-x[0]);
  q=x[1]/(1-x[1]);
  //jacobians
  double tmp0,tmp1;
  tmp0=1-x[0];
  tmp0=1.0/(tmp0*tmp0);
  tmp1=1-x[1];
  tmp1=1.0/(tmp1*tmp1);

  double wa,wb;
  wa=((double *) fdata)[0];
  wb=((double *) fdata)[1];
  double quad,lin;



  quad=-(Q3a*p*p+Q3b*p*q+Q3c*q*q);
  lin=v*(wa*q+wb*p);
  fval[0]=tmp0*tmp1*exp(quad)*cos(lin);

  //  fprintf(stdout,"from gaussianpvarg %lf %lf %.16e\n",p,q,fval[0]);
  return 0;
}






int gaussianpvval(double fdata[2], double *res){
  unsigned ndim=2;
  unsigned fdim=1;
  double xmin[ndim];
  xmin[0]=0.0;
  xmin[1]=0.0;
  double xmax[ndim];
  xmax[0]=1.0;
  xmax[1]=1.0;
 
  double val[fdim];
  double error[fdim];
  double maxEval=MAX_EVAL_INT;
  double reqAbsError=REQ_ABS_ERROR;
  double reqRelError=REQ_REL_ERROR;
  double enorm=ERROR_L2;
  fprintf(stdout,"from gaussianpvval\n");

  hcubature(fdim,gaussianpvarg,fdata,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,enorm,val,error);
  res[0]=val[0];
  res[1]=error[0];
  return 0;
}




double J3(double wa,double wb){
  double res[2];
  double ws[2];

  double W=M_PI*M_PI*M_2_SQRTPI*s*s*s/(3*sqrt(R4))*exp(-(wa*wa*Q1a+wa*wb*Q1b+wb*wb*Q1c));
  double Z=4*sqrt(M_PI)*s*s*s*exp(-(wa*wa*Na+wb*wa*Nb+wb*wb*Nc)/3.0);

  double bound=W+Z*M_PI/(2*sqrt(R4));
  if(bound<eps){
    return eps;
  }
  else{
  ws[0]=wa;
  ws[1]=wb;
  gaussianpvval(ws,res);
  double V=res[0];
  return -W+V*Z;
  }







}







int main(){
  //  fprintf(stdout,"%lf %lf\n",(double)a,(double)b);
  double res[2];
  double ws[2];
  double wa,wb;
  wa=1.0/2;
  wb=-1.0/3;
  ws[0]=wa;
  ws[1]=wb;
  fprintf(stdout,"sa=%.16e\n",sa);
  fprintf(stdout,"sb=%.16e\n",sb);
  fprintf(stdout,"sc=%.16e\n",sc);


  fprintf(stdout,"ca=%.16e\n",ca2);
  fprintf(stdout,"cb=%.16e\n",cb2);
  fprintf(stdout,"cc=%.16e\n",c2);


  fprintf(stdout,"ea=%.16e\n",ea);
  fprintf(stdout,"eb=%.16e\n",eb);
  fprintf(stdout,"eab=%.16e\n",eab);


  fprintf(stdout,"Q1a=%.16e\n",Q1a);
  fprintf(stdout,"Q1b=%.16e\n",Q1b);
  fprintf(stdout,"Q1c=%.16e\n",Q1c);

  fprintf(stdout,"v=%.16e\n",v);

  fprintf(stdout,"Na=%.16e\n",Na);
  fprintf(stdout,"Nb=%.16e\n",Nb);
  fprintf(stdout,"Nc=%.16e\n",Nc);


  fprintf(stdout,"M4=%.16e\n",M4);
  fprintf(stdout,"R4=%.16e\n",R4);


  fprintf(stdout,"Q3a=%.16e\n",Q3a);
  fprintf(stdout,"Q3b=%.16e\n",Q3b);
  fprintf(stdout,"Q3c=%.16e\n",Q3c);


  gaussianpvval(ws,res);
  fprintf(stdout,"V=%lf +/- %lf\n",res[0],res[1]);
  double V=res[0];
  double W=M_PI*M_PI*M_2_SQRTPI*s*s*s/(3*sqrt(R4))*exp(-(wa*wa*Q1a+wa*wb*Q1b+wb*wb*Q1c));
  double Z=4*sqrt(M_PI)*s*s*s*exp(-(wa*wa*Na+wb*wa*Nb+wb*wb*Nc)/3.0);


  fprintf(stdout,"W=%lf\n",W);
  fprintf(stdout,"Z=%lf\n",Z);
  fprintf(stdout,"V*Z=%lf\n",V*Z);
  fprintf(stdout,"\n\nJ_3=%lf\n",-2*M_PI*(-W+V*Z));
  fprintf(stdout,"\n\nJ_3=%lf\n",-2*M_PI*J3(wa,wb));

  return 0;
}
