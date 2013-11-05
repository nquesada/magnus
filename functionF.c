#include<math.h>
#include<stdio.h>

#define so 77.08923800880757
#define se 81.58043418887748
#define sp -91.12584870545336
#define d 0.0
#define sigma 1
#define nu 0.0

#define eps 1.0e-15

double F(double wo, double we, double wp){
  double momentum,energy;
  momentum=so*wo + se*we + sp*wp;
  energy=(wp);//Note that I am assuming sigma=1 and that the axes have been shifted by nu/2,nu/2 and nu
  energy=energy*energy;
  if(momentum<eps && momentum>-eps){
    momentum=1.0;
  }
  else{
    momentum=sin(momentum)/momentum;
  }
  energy=exp(-energy);

  return momentum*energy;
}
