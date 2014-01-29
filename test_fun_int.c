#include "functionF.h"
#include "magnusint.h"
#include <stdio.h>

#define twopi 6.28319



int main(){
  double wa=1/3.0;
  double wb=-1/3.0;
  double waa=wa;
  double wbb=wb;
  double res[2];

  magnus2aint(wa,waa,res);
  fprintf(stdout,"magnus2aint(%lf, %lf ): %.16e +/- %.16e\n",wa,waa,res[0],res[1]);
  magnus2bint(wb,wbb,res);
  fprintf(stdout,"magnus2bint(%lf, %lf ): %.16e +/- %.16e\n",wb,wbb,res[0],res[1]);
  magnus3int(wa,wb,res);
  fprintf(stdout,"magnus3int(%lf, %lf ): %.16e +/- %.16e\n",wa,wb,res[0],res[1]);
  

  return 0;
}
