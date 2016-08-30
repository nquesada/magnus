/*! \file functionF.h
 * \brief
 * Provides the function that returns the product of the phase matching 
 * function and the pump function. For this particular example the function 
 * is gaussian and the exponential is evaluated only if the argument is 
 * bigger that -36. If the number if smaller than \exp(-36) is rounded to zero.
 */


#include<math.h>
#include<stdio.h>

#define so 3.4
#define se 3.6
#define sp 4.0

double F(double wo, double we, double wp){
  double momentum,energy;
  momentum=so*wo + se*we - sp*wp;
  energy=wp;
  double tmp=momentum*momentum+energy*energy;
  double res;
  if(tmp>36){
    res=0.0;
  }
  else{
    res=exp(-tmp);
  }
  return res;

}

void c_F(double wo, double we, double wp, double res[2]){
  double momentum,energy;
  momentum=so*wo + se*we - sp*wp;
  energy=wp;
  double tmp=momentum*momentum+energy*energy;
  if(tmp>36){
    res[0]=0.0;
    res[1]=0.0;
  }
  else{
    res[0]=M_SQRT1_2*exp(-tmp);
    res[1]=M_SQRT1_2*exp(-tmp);
  }

}
