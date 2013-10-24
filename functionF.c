#include<math.h>

#define s0 1.1
#define se 0.9
#define sp 1.0
#define d 1
#define sigma 0.1
#define nu 2


double F(double w0, double we, double wp){
  double tmp1,tmp2;
  tmp1=s0*w0+se*we-sp*wp+d;
  tmp2=sigma*(wp-nu);
  return exp(-tmp1*tmp1-tmp2*tmp2);
}
