#include "functionF.h"
#include "magnusint.h"
#include "c_magnusint.h"
#include <stdio.h>
#include <math.h>
#include <time.h>

#define twopi 6.28319



int main(){
  double wa=1/2.0;
  double wb=-1/3.0;
  double waa=1.0;
  double wbb=1.0;
  double res[2];
  double res4[4];
  double res5[4];
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
  clock_t launch = clock();
  c_magnus3wint(wa,wb,res4);
  clock_t done = clock();
  fprintf(stdout,"c_magnus3wint(%lf, %lf ):  %.16e +/- %.16e + i  %.16e +/- %.16e \n",wa,wb,res4[0],res4[1],res4[2],res4[3]);
  fprintf(stdout,"|c_magnus3wint(%lf, %lf )|: %.16e \n",wa,wb,sqrt(res4[0]*res4[0]+res4[2]*res4[2]));
  double diff = (done - launch);
  fprintf(stdout,"Time in c_magnus3wint %lf\n",diff);
  fprintf(stdout,"\n");

  launch = clock();
  c_magnus3waint(wa,wb,res4);
  c_magnus3wbint(wa,wb,res5);
  done = clock();
  res4[0]+=res5[0];
  res4[1]+=res5[1];
  res4[2]+=res5[2];
  res4[3]+=res5[3];
  fprintf(stdout,"c_magnus3w(a+b)int(%lf, %lf ):  %.16e +/- %.16e + i  %.16e +/- %.16e \n",wa,wb,res4[0],res4[1],res4[2],res4[3]);
  fprintf(stdout,"|c_magnus3(a+b)wint(%lf, %lf )|: %.16e \n",wa,wb,sqrt(res4[0]*res4[0]+res4[2]*res4[2]));
  diff = (done - launch);
  fprintf(stdout,"Time in c_magnus3wint %lf\n",diff);
  fprintf(stdout,"\n");
  return 0;
}
