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

