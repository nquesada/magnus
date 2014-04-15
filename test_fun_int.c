#include "functionF.h"
#include "magnusint.h"
#include <stdio.h>

#define twopi 6.28319



int main(){
  double wa=1/2.0;
  double wb=-1/3.0;
  double waa=1.0;
  double wbb=1.0;
  double res[2];

  magnus2aint(waa,waa,res);
  fprintf(stdout,"magnus2aint(%lf, %lf ): %.16e +/- %.16e\n",waa,waa,res[0],res[1]);
  magnus2bint(waa,waa,res);
  fprintf(stdout,"magnus2bint(%lf, %lf ): %.16e +/- %.16e\n",waa,waa,res[0],res[1]);
  magnus3sint(wa,wb,res);
  fprintf(stdout,"magnus3sint(%lf, %lf ): %.16e +/- %.16e\n",wa,wb,res[0],res[1]);
  magnus3int(wa,wb,res);
  fprintf(stdout,"magnus3int(%lf, %lf ): %.16e +/- %.16e\n",wa,wb,res[0],res[1]);
  return 0;
}
