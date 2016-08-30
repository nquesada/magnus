#include "functionF.h"
#include "magnusint.h"
#include "c_magnusint.h"
#include <stdio.h>
#include <math.h>

#define twopi 6.28319



int main(){
  double wa=1/2.0;
  double wb=-1/3.0;
  double waa=1.0;
  double wbb=1.0;
  double res[2];
  double res4[4];
  magnus2aint(wa,waa,res);
  fprintf(stdout,"magnus2aint(%lf, %lf ):     %.16e +/- %.16e\n",wa,waa,res[0],res[1]);
  c_magnus2aint(wa,waa,res4);
  fprintf(stdout,"c_magnus2aint(%lf, %lf ):   %.16e +/- %.16e + i  %.16e +/- %.16e \n",wa,waa,res4[0],res4[1],res4[2],res4[3]);
  fprintf(stdout,"\n");
  magnus2bint(wb,wbb,res);
  fprintf(stdout,"magnus2bint(%lf, %lf ):     %.16e +/- %.16e\n",wb,wbb,res[0],res[1]);
  c_magnus2bint(wb,wbb,res4);
  fprintf(stdout,"c_magnus2bint(%lf, %lf ):   %.16e +/- %.16e + i  %.16e +/- %.16e \n",wb,wbb,res4[0],res4[1],res4[2],res4[3]);
  fprintf(stdout,"\n");
  magnus3int(wa,wb,res);
  fprintf(stdout,"magnus3int(%lf, %lf ):     %.16e +/- %.16e\n",wa,wb,res[0],res[1]);
  c_magnus3int(wa,wb,res4);
  fprintf(stdout,"c_magnus3int(%lf, %lf ):   %.16e +/- %.16e + i  %.16e +/- %.16e \n",wa,wb,res4[0],res4[1],res4[2],res4[3]);
  fprintf(stdout,"|c_magnus3int(%lf, %lf )|:  %.16e \n",wa,wb,sqrt(res4[0]*res4[0]+res4[2]*res4[2]));
  fprintf(stdout,"\n");
  magnus3sint(wa,wb,res);
  fprintf(stdout,"magnus3sint(%lf, %lf ):     %.16e +/- %.16e\n",wa,wb,res[0],res[1]);
  c_magnus3sint(wa,wb,res4);
  fprintf(stdout,"c_magnus3sint(%lf, %lf ):   %.16e +/- %.16e + i  %.16e +/- %.16e \n",wa,wb,res4[0],res4[1],res4[2],res4[3]);
  fprintf(stdout,"|c_magnus3sint(%lf, %lf )|: %.16e \n",wa,wb,sqrt(res4[0]*res4[0]+res4[2]*res4[2]));
  fprintf(stdout,"\n");
  magnus3wint(wa,wb,res);
  fprintf(stdout,"magnus3wint(%lf, %lf ):    %.16e +/- %.16e\n",wa,wb,res[0],res[1]);
  c_magnus3wint(wa,wb,res4);
  fprintf(stdout,"c_magnus3wint(%lf, %lf ):  %.16e +/- %.16e + i  %.16e +/- %.16e \n",wa,wb,res4[0],res4[1],res4[2],res4[3]);
  fprintf(stdout,"|c_magnus3wint(%lf, %lf )|: %.16e \n",wa,wb,sqrt(res4[0]*res4[0]+res4[2]*res4[2]));
  fprintf(stdout,"\n");
  return 0;
}
