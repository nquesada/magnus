#include<stdio.h>
#include<math.h>
#include <gsl/gsl_math.h>
#include "cubature.h"


/* Values for cubature */
#define MAX_EVAL_INT 10000000 //Maximum number of evaluations
#define REQ_ABS_ERROR 1e-7 //Required absolute error
#define REQ_REL_ERROR 1e-5 //Required relative error

/* Value used to truncate J_3 using the bound for J_3
obtained in the paper. If the bound on J_3 is less
eps the J_3 is not calculated explicitly and it is 
assumed to be zero. eps should be at least as small
as REQ_ABS_ERROR */
#define eps 1e-7


/* This is simply 4/sqrt(3) to 20 decimal places*/
#define M_4_S_3 2.3094010767585030580


/* Values that define the problem */
/*#define sa  320.049
#define sb  317.973
#define sc  319.008
#define s   1.0 //\sigma
*/
/* 
#define sa  588.751
#define sb  590.769
#define sc  589.745
*/

#define sa  1.0
#define sb  -1.0
#define sc  0.0

#define s   1.0 //\this is the \tau of the text

/* derived quantitities, do not modify! */

#define ea  sc-sa //\eta_a
#define eb  sc-sb //\eta_b




#define c2  (s*s+(sc-sa)*(sc-sb)) //\chi^2
#define ca2 (s*s+(sc-sa)*(sc-sa)) //\chi_a^2
#define cb2 (s*s+(sc-sb)*(sc-sb)) //\chi_b^2
#define R4  (4*(s*s+(sc-sa)*(sc-sa))*(s*s+(sc-sb)*(sc-sb))-(s*s+(sc-sa)*(sc-sb))*(s*s+(sc-sa)*(sc-sb))) //R^4
#define M4  (4*(s*s+(sc-sa)*(sc-sa))*(s*s+(sc-sb)*(sc-sb))-3*(s*s+(sc-sa)*(sc-sb))*(s*s+(sc-sa)*(sc-sb))) //M^4

// For the variables defined up to this point a,b,c are used to label properties of the different mode 



//This are the elements of the different matrices in the paper
// Xa indicates the element (1,1), Xb twice the element (1,2)=(2,1)
// Xc indicates the element 2,2

#define Na (ca2)
#define Nb (2*c2)
#define Nc (cb2)
#define Qa (M4*ca2/R4)
#define Qb (2*c2*c2*c2/R4)
#define Qc (M4*cb2/R4)
#define Ma (2*ca2)
#define Mb (2*c2)
#define Mc (2*cb2)

/* This quantity defines the prefactor inside the cosine */
#define vv M_4_S_3*s*(sb-sa)
#define Wpref M_PI*M_PI*M_2_SQRTPI*s*s*s/(3.0*sqrt(R4))
#define Zpref 4*sqrt(M_PI)*s*s*s
int gaussianarg(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
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
  double quad,lin,pr;

  quad=-(Ma*p*p+Mb*p*q+Mc*q*q);

  
  pr=(wa*q)+(wb*p);
  lin=vv*pr;
  fval[0]=tmp0*tmp1*exp(quad)*cos(lin);

  return 0;
}

int gaussianval(double fdata[2], double *res){
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
  hcubature(fdim,gaussianarg,fdata,ndim,xmin,xmax,(size_t)maxEval,reqAbsError,reqRelError,enorm,val,error);
  res[0]=val[0];
  res[1]=error[0];
  return 0;
}

double J3(double wa,double wb){
  double res[2];
  double ws[2];
  double W=Wpref*exp(-(wa*wa*Qa+wa*wb*Qb+wb*wb*Qc));
  double Zp=exp(-(wa*wa*Na+wb*wa*Nb+wb*wb*Nc)/3.0);

  double bound=W+Wpref*2*Zp;
  if(bound<eps){
    return eps;
  }
  else{
    ws[0]=wa;
    ws[1]=wb;
    gaussianval(ws,res);
    double V=res[0];
    return W-V*Zp*Zpref;
  }
}


int main(){


  double res[2];
  double ws[2];
  /*
  double wa,wb;
  wa=1.0/2;
  wb=-1.0/3;
  ws[0]=wa;
  ws[1]=wb;
  */


  
  fprintf(stdout,"sa=%.16e\n",sa);
  fprintf(stdout,"sb=%.16e\n",sb);
  fprintf(stdout,"sc=%.16e\n",sc);


  fprintf(stdout,"ca2=%.16e\n",ca2);
  fprintf(stdout,"cb2=%.16e\n",cb2);
  fprintf(stdout,"cc2=%.16e\n",c2);


  fprintf(stdout,"ea=%.16e\n",ea);
  fprintf(stdout,"eb=%.16e\n",eb);


  fprintf(stdout,"Qa=%.16e\n",Qa);
  fprintf(stdout,"Qb=%.16e\n",Qb);
  fprintf(stdout,"Qc=%.16e\n",Qc);

  fprintf(stdout,"v=%.16e\n",vv);

  fprintf(stdout,"Na=%.16e\n",Na);
  fprintf(stdout,"Nb=%.16e\n",Nb);
  fprintf(stdout,"Nc=%.16e\n",Nc);


  fprintf(stdout,"M4=%.16e\n",M4);
  fprintf(stdout,"R4=%.16e\n",R4);


  fprintf(stdout,"Ma=%.16e\n",Ma);
  fprintf(stdout,"Mb=%.16e\n",Mb);
  fprintf(stdout,"Mc=%.16e\n",Mc);

  
  /*  gaussianval(ws,res);
  fprintf(stdout,"V=%lf +/- %lf\n",res[0],res[1]);
  double V=res[0];
  double W=M_PI*M_PI*M_2_SQRTPI*s*s*s/(3.0*sqrt(R4))*exp(-(wa*wa*Qa+wa*wb*Qb+wb*wb*Qc));
  double Z=4*sqrt(M_PI)*s*s*s*exp(-(wa*wa*Na+wb*wa*Nb+wb*wb*Nc)/3.0);


  fprintf(stdout,"W=%lf\n",W);
  fprintf(stdout,"Z=%lf\n",Z);
  fprintf(stdout,"V*Z=%lf\n",V*Z);
  fprintf(stdout,"J_3(%lf,%lf)=%lf\n",wa,wb,-W+V*Z);
  fprintf(stdout,"J_3(%lf,%lf)=%lf\n",wa,wb,J3(wa,wb));
  */

  /*
  int ndim=2;
  double x[2];
  x[0]=0.5;
  x[1]=0.5;
  double fdata[2];
  fdata[0]=fdata[1]=0.0;
  int fdim=1;
  double fval[1];
  gaussianarg(ndim,x,fdata,fdim,fval);
  */


  FILE *pf;
  pf=fopen("dataEnt.dat","w");
  double wa,wb;
  double ll=1.5;
  double dl=0.01;
  wa=wb=1.0;//1.5265566588595902e-15+1;
  ws[0]=wa;
  ws[1]=wb;
  gaussianval(ws,res);

  fprintf(stdout,"At the origin %lf +/- %lf\n",res[0],res[1]);
  wa=0.0;
  wb=0.0;
  fprintf(stdout,"J3(%lf,%lf)=%lf\n",wa,wb,J3(wa,wb));

  wa=1.0;
  wb=1.0;
  fprintf(stdout,"J3(%lf,%lf)=%lf\n",wa,wb,J3(wa,wb));

  wa=1.0;
  wb=2.0;
  fprintf(stdout,"J3(%lf,%lf)=%lf\n",wa,wb,J3(wa,wb));
  /*
  for(wa=-ll;wa<=ll;wa+=dl){
    for(wb=-ll;wb<=ll;wb+=dl){
      fprintf(pf,"%.16e ",J3(wa,wb));
    }
    fprintf(pf,"\n");
    fflush(pf);
    }*/
  fclose(pf);

  return 0;
}
