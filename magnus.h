/*! \file magnus.h
 * \brief Defines the functions that are required by the cubature library
 * to perform the numerical integration. magnus2a, magnus2b and magnus3
 * calculate the second order frequency correction terms for modes a and b
 * and magnus 3 calculates the third order sequeezing correction.
 * Input: ndim, the number of dimensions in which the integral will be performed
 * which is 2 for magnus2a and magnus 2b and 3 for magnus 3.
 * *x which is the vector containing the values of the variables over which 
 * integration will be perfomed. Using the notation of the paper we have 
 * for magnus 2a:
 * x[0]=\omega_p
 * x[1]=\omega_d
 * for magnus 2b
 * x[0]=\omega_p
 * x[1]=\omega_c
 * and for magnus 3a:
 * x[0]=\omega_p
 * x[1]=\omega_q
 * x[2]=\omega_c
 * x[3]=\omega_d
 * Note that variables c,d would be in the range (-oo,oo) are internally 
 * transformed to the range (-1,1) via x->x/(1-x^2) using a change of variables, 
 * the functions return the function evaluated in the transformed variables
 * times the jacobian of the transformation. For p and q that range in (0,oo)
 * the change of variables is x->x/(1-x) and they now live in [0,1)
 * void *fdata contains the two frequencies omega_a, omega_a' for magnus_2a
 * \omega_b, \omega_b' for magnus_2b and \omega_a and \omega_b of which the
 * integration is a function, i.e., they are not integration variables
 * fdim is the dimension of the integration which is one in all cases.
 * Output: *fval is a one dimensional array containing the value of the integrand 
 * at the specified sampling point *x and parameter *fdata.
 */



#ifndef _magnus_
#define _magnus_

int magnus2a(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int magnus2b(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int magnus3(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int magnus3s(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
#endif
