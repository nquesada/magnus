/* \file magnusint.h This function internally perform the integrals necessary
 * to obtain the magnus correction G_2 (magnus2aint), H_2 (magnus2bint) and 
 * magnus3int (J_3) the results of the integration are returned in the
 * dimensions array *res, whose first component is the value of the integral
 * and the second one is the value of the error.
 * MAX_EXAL_INT specifies the number of times the function in magnus.c
 * can be called by the integrator
 * REQ_ABS_ERROR specifies the absolute error that one desires from the 
 * integration
 * REQ_REL_ERROR specifies the relative error.
 */


#ifndef _c_magnusint_
#define _c_magnusint_

int c_magnus2aint(double wa, double waa, double *res);
int c_magnus2bint(double wb, double wbb, double *res);
int c_magnus3int(double wa, double wb, double *res);
int c_magnus3sint(double wa, double wb, double *res);
int c_magnus3wint(double wa, double wb, double *res);
int c_magnus3waint(double wa, double wb, double *res);
int c_magnus3wbint(double wa, double wb, double *res);

#endif
