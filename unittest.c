#include <stdio.h>
#include <math.h>
#include "magnus.h"
#include "functionF.h"


int main(){
  double x[4];
  double ws[2];
  double res[1];
  x[0]=0.618034;
  x[1]=0.618034;
  x[2]=0.0;
  x[3]=0.0;
  ws[0]=1.0;
  ws[1]=1.0;

  fprintf(stdout,"This little program performs a unit test of the routines in magnus.c \n");
  magnus2a(2,x,ws,1,res);
  fprintf(stdout,"magnus2a(%lf,%lf ;%lf,%lf)=%.16e\n",x[0],x[1],ws[0],ws[1],res[0]);
  magnus2b(2,x,ws,1,res);
  fprintf(stdout,"magnus2b(%lf,%lf ;%lf,%lf)=%.16e\n",x[0],x[1],ws[0],ws[1],res[0]);
  magnus3(4,x,ws,1,res);
  fprintf(stdout,"magnus3(%lf,%lf,%lf %lf ;%lf,%lf)=%.16e\n",x[0],x[1],x[2],x[3],ws[0],ws[1],res[0]);

  magnus3s(2,x,ws,1,res);
  fprintf(stdout,"magnus3s(%lf,%lf ;%lf,%lf)=%.16e\n",x[0],x[1],ws[0],ws[1],res[0]);



  return 0;
}
