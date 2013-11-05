#include"functionF.h"
#include<stdio.h>


//This function need a prefactor of 2*pi and corresponds to G in the manuscript
int magnus2a(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){

  double u,v,wa1,wa2;
  u=x[0]/(1-x[0]*x[0]);
  v=x[1]/(1-x[1]);
  //jacobians
  double tmp1,tmp2;
  tmp1=1-x[0]*x[0];
  tmp1=(1+x[0]*x[0])/(tmp1*tmp1);
  tmp2=1-x[1];
  tmp2=1.0/(tmp2*tmp2);


  wa1=((double *) fdata)[0];
  wa2=((double *) fdata)[1];

  //The next line was generated automatically from mathematica, do not touch!
  fval[0]=tmp1*tmp2*(-(F(wa1,u - v,u + v + wa1)*F(wa2,u - v,u + v + wa2)) + 
		     F(wa1,u + v,u - v + wa1)*F(wa2,u + v,u - v + wa2))/v;
  return 0;
}


//This function need a prefactor of 2*pi and corresponds to H in the manuscript
int magnus2b(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
  
  double u,v,wb1,wb2;
  u=x[0]/(1-x[0]*x[0]);
  v=x[1]/(1-x[1]);
  //jacobians
  double tmp1,tmp2;
  tmp1=1-x[0]*x[0];
  tmp1=(1+x[0]*x[0])/(tmp1*tmp1);
  tmp2=1-x[1];
  tmp2=1.0/(tmp2*tmp2);


  wb1=((double *) fdata)[0];
  wb2=((double *) fdata)[1];

  //The next line was generated automatically from mathematica, do not touch!
  fval[0]=tmp1*tmp2*(-(F(u - v,wb1,u + v + wb1)*F(u - v,wb2,u + v + wb2)) + 
		     F(u + v,wb1,u - v + wb1)*F(u + v,wb2,u - v + wb2))/v;

  return 0;
}


//This function need a prefactor of 2*pi^2 and corresponds to J in the manuscript
int magnus3x(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){

  double u,v,w,wa,wb;
  u=x[0]/(1-x[0]*x[0]);
  v=x[1]/(1-x[1]);
  w=x[2]/(1-x[2]*x[2]);
  //jacobians
  double tmp1,tmp2,tmp3;
  tmp1=1-x[0]*x[0];
  tmp1=(1+x[0]*x[0])/(tmp1*tmp1);
  tmp2=1-x[1];
  tmp2=1.0/(tmp2*tmp2);
  tmp3=1-x[2]*x[2];
  tmp3=(1+x[2]*x[2])/(tmp3*tmp3);



  wa=((double *) fdata)[0];
  wb=((double *) fdata)[1];

  //The next line was generated automatically from mathematica, do not touch!
  fval[0]=tmp1*tmp2*tmp3*(-(F(w,u - v,u + v + w)*F(w,wb,w + wb)*F(wa,u - v,u + v + wa)) + 
			  F(w,u + v,u - v + w)*F(w,wb,w + wb)*F(wa,u + v,u - v + wa) + 
			  (-(F(u - v,w,u + v + w)*F(u - v,wb,u + v + wb)) + 
			   F(u + v,w,u - v + w)*F(u + v,wb,u - v + wb))*F(wa,w,w + wa))/v;
  return 0;
}



//This function need a prefactor of pi^2/3 and corresponds to J in the manuscript
int magnus3r(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
  double u,v,w,wa,wb;
  u=x[0]/(1-x[0]*x[0]);
  v=x[1]/(1-x[1]);
  w=x[2]/(1-x[2]*x[2]);
  //jacobians
  double tmp1,tmp2,tmp3;
  tmp1=1-x[0]*x[0];
  tmp1=(1+x[0]*x[0])/(tmp1*tmp1);
  tmp2=1-x[1];
  tmp2=1.0/(tmp2*tmp2);
  tmp3=1-x[2]*x[2];
  tmp3=(1+x[2]*x[2])/(tmp3*tmp3);


  wa=((double *) fdata)[0];
  wb=((double *) fdata)[1];
  //The next line was generated automatically from mathematica, do not touch!
  fval[0]=tmp1*tmp2*tmp3*(-((F(w,u - v,2*u + wa - wb)*F(w,wb,u - v + wa) - 
			     F(w,u - v,-w + 2*(u + v + wa) - wb)*F(w,wb,u + v - w + 2*wa) + 
			     F(w,u - v,u - v + wa)*F(w,wb,2*u + wa - wb) - 
			     F(w,u - v,u + v - w + 2*wa)*F(w,wb,-w + 2*(u + v + wa) - wb))*F(wa,u - v,u + v + wa))
			  + (F(w,u + v,2*u + wa - wb)*F(w,wb,u + v + wa) - 
			     F(w,u + v,-w + 2*(u - v + wa) - wb)*F(w,wb,u - v - w + 2*wa) + 
			     F(w,u + v,u + v + wa)*F(w,wb,2*u + wa - wb) - 
			     F(w,u + v,u - v - w + 2*wa)*F(w,wb,-w + 2*(u - v + wa) - wb))*F(wa,u + v,u - v + wa) + 
			  F(u + v,wb,u - v + wb)*(F(u + v,w,w + wa)*F(wa,w,u - v + w) + 
						  F(u + v,w,u - v + w)*F(wa,w,w + wa) - F(u + v,w,-2*v + wa + wb)*F(wa,w,u - 3*v + wb) - 
						  F(u + v,w,u - 3*v + wb)*F(wa,w,-2*v + wa + wb)) - 
			  F(u - v,wb,u + v + wb)*(F(u - v,w,w + wa)*F(wa,w,u + v + w) + 
						  F(u - v,w,u + v + w)*F(wa,w,w + wa) - F(u - v,w,2*v + wa + wb)*F(wa,w,u + 3*v + wb) - 
						  F(u - v,w,u + 3*v + wb)*F(wa,w,2*v + wa + wb)) - 
			  12*F(u - 2*v - 3*w,u - 2*v + 3*w,2*(u + v))*
			  (F(u - 2*v - 3*w,wb,3*(u + w) - wb)*F(wa,u - 2*v + 3*w,u - 2*v + 3*w + wa) - 
			   F(u - 2*v - 3*w,wb,3*(u + 2*v + w) - wb)*F(wa,u - 2*v + 3*w,u + 4*v + 3*w + wa) + 
			   F(u - 2*v - 3*w,wb,u - 2*v + 3*w + wa)*F(wa,u - 2*v + 3*w,3*(u + w) - wb) - 
			   F(u - 2*v - 3*w,wb,u + 4*v + 3*w + wa)*F(wa,u - 2*v + 3*w,3*(u + 2*v + w) - wb)) + 
			  12*F(u + 2*v - 3*w,u + 2*v + 3*w,2*u - 2*v)*
			  (-(F(u + 2*v - 3*w,wb,3*(u - 2*v + w) - wb)*F(wa,u + 2*v + 3*w,u - 4*v + 3*w + wa)) + 
			   F(u + 2*v - 3*w,wb,3*(u + w) - wb)*F(wa,u + 2*v + 3*w,u + 2*v + 3*w + wa) + 
			   F(u + 2*v - 3*w,wb,u + 2*v + 3*w + wa)*F(wa,u + 2*v + 3*w,3*(u + w) - wb) - 
			   F(u + 2*v - 3*w,wb,u - 4*v + 3*w + wa)*F(wa,u + 2*v + 3*w,3*(u - 2*v + w) - wb)))/v;

  return 0;
}

//This function needs a prefactor of i*pi/3
int magnus3i(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
  double v,y,z,u,wa,wb;
  v=x[0]/(1-x[0]);//I was using x but x is the name of the argument of the function
  y=x[1]/(1-x[1]);
  z=x[2]/(1-x[2]*x[2]);
  u=x[3]/(1-x[3]*x[3]);
  //jacobians
  double tmp1,tmp2,tmp3,tmp4;
  tmp1=1-x[0];
  tmp1=1.0/(tmp1*tmp1);
  tmp2=1-x[1];
  tmp2=1.0/(tmp2*tmp2);
  tmp3=1-x[2]*x[2];
  tmp3=(1+x[2]*x[2])/(tmp3*tmp3);
  tmp4=1-x[3]*x[3];
  tmp4=(1+x[3]*x[3])/(tmp4*tmp4);


  wa=((double *) fdata)[0];
  wb=((double *) fdata)[1];
  //The next line was generated automatically from mathematica, do not touch!
  fval[0]=tmp1*tmp2*tmp3*tmp4*(2*F(wa,2*(v + y + z),wa + 2*z)*F(u - 2*v - y - 3*z,wb,u - wb - y + z)*
			       F(u - 2*v - y - 3*z,2*(v + y + z),u - y - z) - 
			       2*F(wa,2*(-v + y + z),wa + 2*z)*F(u + 2*v - y - 3*z,wb,u - wb - y + z)*
			       F(u + 2*v - y - 3*z,2*(-v + y + z),u - y - z) - 
			       2*F(wa,2*(v - y + z),wa + 2*z)*F(u - 2*v + y - 3*z,wb,u - wb + y + z)*
			       F(u - 2*v + y - 3*z,2*(v - y + z),u + y - z) + 
			       2*F(wa,-2*(v + y - z),wa + 2*z)*F(u + 2*v + y - 3*z,wb,u - wb + y + z)*
			       F(u + 2*v + y - 3*z,-2*(v + y - z),u + y - z) - 
			       F(wa,v - y + z,-v + wa + y + z)*F(u - v - z,wb,u + v + wb - z)*
			       F(u - v - z,v - y + z,u + y) + F(wa,v + y + z,-v + wa - y + z)*
			       F(u - v - z,wb,u + v + wb - z)*F(u - v - z,v + y + z,u - y) + 
			       F(wa,-v - y + z,v + wa + y + z)*F(u + v - z,wb,u - v + wb - z)*
			       F(u + v - z,-v - y + z,u + y) - F(wa,-v + y + z,v + wa - y + z)*
			       F(u + v - z,wb,u - v + wb - z)*F(u + v - z,-v + y + z,u - y) - 
			       8*F(wa,-u - 2*v - y + 3*z,3*u + 2*v - wa + 3*y - z)*
			       F(2*u - 2*y - 2*z,-u - 2*v - y + 3*z,u + 2*v + y + z)*
			       F(2*(u - y - z),wb,wb + 2*(u + y - z)) + 
			       8*F(wa,-u + 2*v - y + 3*z,3*u - 2*v - wa + 3*y - z)*
			       F(2*u - 2*y - 2*z,-u + 2*v - y + 3*z,u - 2*v + y + z)*
			       F(2*(u - y - z),wb,wb + 2*(u + y - z)) + 
			       8*F(wa,-u - 2*v + y + 3*z,3*u + 2*v - wa - 3*y - z)*
			       F(2*(u + y - z),wb,2*u + wb - 2*y - 2*z)*
			       F(2*(u + y - z),-u - 2*v + y + 3*z,u + 2*v - y + z) - 
			       8*F(wa,-u + 2*v + y + 3*z,3*u - 2*v - wa - 3*y - z)*
			       F(2*(u + y - z),wb,2*u + wb - 2*y - 2*z)*
			       F(2*(u + y - z),-u + 2*v + y + 3*z,u - 2*v - y + z) + 
			       2*F(wa,u + 2*v + y - 3*z,u - wa + y + z)*F(-2*(v + y - z),wb,wb + 2*z)*
			       F(-2*(v + y - z),u + 2*v + y - 3*z,u + y - z) + 
			       F(wa,u + v - z,u - v + wa - z)*F(-v - y + z,wb,v + wb + y + z)*
			       F(-v - y + z,u + v - z,u + y) - F(wa,u - v - z,u + v + wa - z)*
			       F(v - y + z,wb,-v + wb + y + z)*F(v - y + z,u - v - z,u + y) - 
			       2*F(wa,u - 2*v + y - 3*z,u - wa + y + z)*F(2*(v - y + z),wb,wb + 2*z)*
			       F(2*(v - y + z),u - 2*v + y - 3*z,u + y - z) - 
			       F(wa,u + v - z,u - v + wa - z)*F(-v + y + z,wb,v + wb - y + z)*
			       F(-v + y + z,u + v - z,u - y) - 2*F(wa,u + 2*v - y - 3*z,u - wa - y + z)*
			       F(2*(-v + y + z),wb,wb + 2*z)*F(2*(-v + y + z),u + 2*v - y - 3*z,u - y - z) + 
			       F(wa,u - v - z,u + v + wa - z)*F(v + y + z,wb,-v + wb - y + z)*
			       F(v + y + z,u - v - z,u - y) + 2*F(wa,u - 2*v - y - 3*z,u - wa - y + z)*
			       F(2*(v + y + z),wb,wb + 2*z)*F(2*(v + y + z),u - 2*v - y - 3*z,u - y - z) - 
			       8*F(wa,2*(u - v - z),wa + 2*(u + v - z))*
			       F(-u - v - 2*y + 3*z,wb,3*u + 3*v - wb + 2*y - z)*
			       F(-u - v - 2*y + 3*z,2*u - 2*v - 2*z,u + v + 2*y + z) + 
			       8*F(wa,2*(u + v - z),2*u - 2*v + wa - 2*z)*
			       F(-u + v - 2*y + 3*z,wb,3*u - 3*v - wb + 2*y - z)*
			       F(-u + v - 2*y + 3*z,2*(u + v - z),u - v + 2*y + z) + 
			       8*F(wa,2*(u - v - z),wa + 2*(u + v - z))*
			       F(-u - v + 2*y + 3*z,wb,3*u + 3*v - wb - 2*y - z)*
			       F(-u - v + 2*y + 3*z,2*u - 2*v - 2*z,u + v - 2*y + z) - 
			       8*F(wa,2*(u + v - z),2*u - 2*v + wa - 2*z)*
			       F(-u + v + 2*y + 3*z,wb,3*u - 3*v - wb - 2*y - z)*
			       F(-u + v + 2*y + 3*z,2*(u + v - z),u - v - 2*y + z))/(v*y);
  return 0;
}

