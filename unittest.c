#include <stdio.h>
#include <math.h>
#include "magnus.h"
#include "functionF.h"


int main(){
  double x[4];
  double ws[2];
  double res[1];
  x[0]=0.1;
  x[1]=0.2;
  x[2]=0.3;
  x[3]=0.4;
  ws[0]=1;
  ws[1]=2;

  fprintf(stdout,"F(%lf,%lf,%lf)=%lf\n",ws[0],ws[1],ws[0]+ws[1],F(ws[0],ws[1],ws[0]+ws[1]));
  magnus2a(2,x,ws,1,res);
  fprintf(stdout,"magnus2a(%lf,%lf;%lf,%lf)=%.16e\n",x[0],x[1],ws[0],ws[1],res[0]);
  magnus2b(2,x,ws,1,res);
  fprintf(stdout,"magnus2b(%lf,%lf;%lf,%lf)=%.16e\n",x[0],x[1],ws[0],ws[1],res[0]);
  magnus3x(3,x,ws,1,res);
  fprintf(stdout,"magnus3x(%lf,%lf,%lf;%lf,%lf)=%.16e\n",x[0],x[1],x[2],ws[0],ws[1],res[0]);
  magnus3r(3,x,ws,1,res);
  fprintf(stdout,"magnus3r(%lf,%lf,%lf;%lf,%lf)=%.16e\n",x[0],x[1],x[2],ws[0],ws[1],res[0]);
  magnus3i(3,x,ws,1,res);
  fprintf(stdout,"magnus3r(%lf,%lf,%lf %lf ;%lf,%lf)=%.16e\n",x[0],x[1],x[2],x[3],ws[0],ws[1],res[0]);



  return 0;
}
