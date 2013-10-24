#include"functionF.h"

int magnus2a(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
  
  double u,v,we1,we2;
  u=x[0]/(1-x[0]*x[0]);
  v=x[1]/(1-x[1]);
  //jacobians
  double tmp1,tmp2;
  tmp1=1-x[0]*x[0];
  tmp1=(1+x[0]*x[0])/(tmp1*tmp1);
  tmp2=1-x[1];
  tmp2=1.0/(tmp2*tmp2);


  we1=((double *) fdata)[0];
  we2=((double *) fdata)[1];

  fval[0]=tmp1*tmp2*(F(u+v,we1,u-v+we1)*F(u+v,we2,u-v+we2)-F(u-v,we1,u+v+we1)*F(u-v,we2,u+v+we2))/v;

  return 0;
}

int magnus2b(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
  
  double u,v,wo1,wo2;
  u=x[0]/(1-x[0]*x[0]);
  v=x[1]/(1-x[1]);
  //jacobians
  double tmp1,tmp2;
  tmp1=1-x[0]*x[0];
  tmp1=(1+x[0]*x[0])/(tmp1*tmp1);
  tmp2=1-x[1];
  tmp2=1.0/(tmp2*tmp2);


  wo1=((double *) fdata)[0];
  wo2=((double *) fdata)[1];
  fval[0]=tmp1*tmp2*(F(wo1,u+v,u-v+wo1)*F(wo2,u+v,u-v+wo2)-F(wo1,u-v,u+v+wo1)*F(wo2,u-v,u+v+wo2))/v;

  return 0;
}

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
  fval[0]=tmp1*tmp2*tmp3*(-(F(u - v,w,u + v + w)*(F(u - v,wa,u + v + wa)*F(w,wb,w + wb) +
						  F(u - v,wb,u + v + wb)*F(wa,w,w + wa))) + 
			  F(u + v,w,u - v + w)*(F(u + v,wa,u - v + wa)*F(w,wb,w + wb) + 
						F(u + v,wb,u - v + wb)*F(wa,w,w + wa)))/v;

  return 0;
}


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
  fval[0]=tmp1*tmp2*tmp3*(-((F(w,u - v,w + wa)*F(w,wa,u + v + w) + F(w,u - v,u + v + w)*F(w,wa,w + wa) - 
			     F(w,u - v,2*v + wa + wb)*F(w,wa,u + 3*v + wb) - 
			     F(w,u - v,u + 3*v + wb)*F(w,wa,2*v + wa + wb))*F(wb,u - v,u + v + wb)) + 
			  (F(w,u + v,w + wa)*F(w,wa,u - v + w) + F(w,u + v,u - v + w)*F(w,wa,w + wa) - 
			   F(w,u + v,-2*v + wa + wb)*F(w,wa,u - 3*v + wb) - 
			   F(w,u + v,u - 3*v + wb)*F(w,wa,-2*v + wa + wb))*F(wb,u + v,u - v + wb) + 
			  F(u + v,wa,u - v + wa)*(F(u + v,w,2*u + wa - wb)*F(wb,w,u + v + wa) - 
						  F(u + v,w,-w + 2*(u - v + wa) - wb)*F(wb,w,u - v - w + 2*wa) + 
						  F(u + v,w,u + v + wa)*F(wb,w,2*u + wa - wb) - 
						  F(u + v,w,u - v - w + 2*wa)*F(wb,w,-w + 2*(u - v + wa) - wb)) - 
			  F(u - v,wa,u + v + wa)*(F(u - v,w,2*u + wa - wb)*F(wb,w,u - v + wa) - 
						  F(u - v,w,-w + 2*(u + v + wa) - wb)*F(wb,w,u + v - w + 2*wa) + 
						  F(u - v,w,u - v + wa)*F(wb,w,2*u + wa - wb) - 
						  F(u - v,w,u + v - w + 2*wa)*F(wb,w,-w + 2*(u + v + wa) - wb)) - 
			  12*F(u - 2*v - 3*w,u - 2*v + 3*w,2*(u + v))*
			  (F(u - 2*v - 3*w,wa,3*u - 3*w - wb)*F(wb,u - 2*v + 3*w,u - 2*v - 3*w + wa) - 
			   F(u - 2*v - 3*w,wa,3*u + 6*v - 3*w - wb)*F(wb,u - 2*v + 3*w,u + 4*v - 3*w + wa) + 
			   F(u - 2*v - 3*w,wa,u - 2*v - 3*w + wa)*F(wb,u - 2*v + 3*w,3*u - 3*w - wb) - 
			   F(u - 2*v - 3*w,wa,u + 4*v - 3*w + wa)*F(wb,u - 2*v + 3*w,3*u + 6*v - 3*w - wb)) + 
			  12*F(u + 2*v - 3*w,u + 2*v + 3*w,2*u - 2*v)*
			  (-(F(u + 2*v - 3*w,wa,3*u - 6*v - 3*w - wb)*F(wb,u + 2*v + 3*w,u - 4*v - 3*w + wa)) + 
			   F(u + 2*v - 3*w,wa,3*u - 3*w - wb)*F(wb,u + 2*v + 3*w,u + 2*v - 3*w + wa) + 
			   F(u + 2*v - 3*w,wa,u + 2*v - 3*w + wa)*F(wb,u + 2*v + 3*w,3*u - 3*w - wb) - 
			   F(u + 2*v - 3*w,wa,u - 4*v - 3*w + wa)*F(wb,u + 2*v + 3*w,3*u - 6*v - 3*w - wb)))/v;

  return 0;
}

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
  fval[0]=tmp1*tmp2*tmp3*tmp4*(2*F(2*u,wa,wa + 2*(u + v + y))*F(2*u,-3*u - v - 2*y + z,-u - v + z)*
			       F(wb,-3*u - v - 2*y + z,u - wb + v + 2*y + z) - 
			       2*F(2*u,wa,wa + 2*(u - v + y))*F(2*u,-3*u + v - 2*y + z,-u + v + z)*
			       F(wb,-3*u + v - 2*y + z,u - wb - v + 2*y + z) - 
			       2*F(2*u,wa,wa + 2*(u + v - y))*F(2*u,-3*u - v + 2*y + z,-u - v + z)*
			       F(wb,-3*u - v + 2*y + z,u - wb + v - 2*y + z) + 
			       2*F(2*u,wa,2*u + wa - 2*v - 2*y)*F(2*u,-3*u + v + 2*y + z,-u + v + z)*
			       F(wb,-3*u + v + 2*y + z,u - wb - v - 2*y + z) + 
			       F(wb,-u + v + z,-u + wb - v + z)*F(u - v - y,wa,u + wa + v + y)*
			       F(u - v - y,-u + v + z,y + z) - 
			       F(wb,-u - v + z,-u + wb + v + z)*F(u + v - y,wa,u + wa - v + y)*
			       F(u + v - y,-u - v + z,y + z) - 
			       F(wb,-u + v + z,-u + wb - v + z)*F(u - v + y,wa,u + wa + v - y)*
			       F(u - v + y,-u + v + z,-y + z) + 
			       F(wb,-u - v + z,-u + wb + v + z)*F(u + v + y,wa,u + wa - v - y)*
			       F(u + v + y,-u - v + z,-y + z) - 
			       8*F(wb,3*u - v - 2*y - z,-u - wb + 3*v + 2*y + 3*z)*
			       F(-2*(u + v - z),wa,wa + 2*(-u + v + z))*
			       F(-2*(u + v - z),3*u - v - 2*y - z,u + v + 2*y + z) + 
			       8*F(wb,3*u - v + 2*y - z,-u - wb + 3*v - 2*y + 3*z)*
			       F(-2*(u + v - z),wa,wa + 2*(-u + v + z))*
			       F(-2*(u + v - z),3*u - v + 2*y - z,u + v - 2*y + z) - 
			       F(wb,u + v - y,u + wb - v + y)*F(-u - v + z,wa,-u + wa + v + z)*
			       F(-u - v + z,u + v - y,y + z) + 
			       F(wb,u + v + y,u + wb - v - y)*F(-u - v + z,wa,-u + wa + v + z)*
			       F(-u - v + z,u + v + y,-y + z) + 
			       F(wb,u - v - y,u + wb + v + y)*F(-u + v + z,wa,-u + wa - v + z)*
			       F(-u + v + z,u - v - y,y + z) - 
			       F(wb,u - v + y,u + wb + v - y)*F(-u + v + z,wa,-u + wa - v + z)*
			       F(-u + v + z,u - v + y,-y + z) + 
			       8*F(wb,3*u + v - 2*y - z,-u - wb - 3*v + 2*y + 3*z)*
			       F(2*(-u + v + z),wa,wa - 2*(u + v - z))*
			       F(2*(-u + v + z),3*u + v - 2*y - z,u - v + 2*y + z) - 
			       8*F(wb,3*u + v + 2*y - z,-u - wb - 3*v - 2*y + 3*z)*
			       F(2*(-u + v + z),wa,wa - 2*(u + v - z))*
			       F(2*(-u + v + z),3*u + v + 2*y - z,u - v - 2*y + z) + 
			       2*F(wb,2*u,wb + 2*(u + v + y))*F(-3*u - v - 2*y + z,2*u,-u - v + z)*
			       F(-3*u - v - 2*y + z,wa,u - wa + v + 2*y + z) - 
			       2*F(wb,2*u,wb + 2*(u - v + y))*F(-3*u + v - 2*y + z,2*u,-u + v + z)*
			       F(-3*u + v - 2*y + z,wa,u - wa - v + 2*y + z) - 
			       2*F(wb,2*u,wb + 2*(u + v - y))*F(-3*u - v + 2*y + z,2*u,-u - v + z)*
			       F(-3*u - v + 2*y + z,wa,u - wa + v - 2*y + z) + 
			       2*F(wb,2*u,2*u + wb - 2*v - 2*y)*F(-3*u + v + 2*y + z,2*u,-u + v + z)*
			       F(-3*u + v + 2*y + z,wa,u - wa - v - 2*y + z) - 
			       8*F(wb,2*(u - v - z),wb + 2*(u + v - z))*
			       F(-u - v - 2*y + 3*z,wa,3*u - wa + 3*v + 2*y - z)*
			       F(-u - v - 2*y + 3*z,2*u - 2*v - 2*z,u + v + 2*y + z) + 
			       8*F(wb,2*(u + v - z),2*u + wb - 2*v - 2*z)*
			       F(-u + v - 2*y + 3*z,wa,3*u - wa - 3*v + 2*y - z)*
			       F(-u + v - 2*y + 3*z,2*(u + v - z),u - v + 2*y + z) + 
			       8*F(wb,2*(u - v - z),wb + 2*(u + v - z))*
			       F(-u - v + 2*y + 3*z,wa,3*u - wa + 3*v - 2*y - z)*
			       F(-u - v + 2*y + 3*z,2*u - 2*v - 2*z,u + v - 2*y + z) - 
			       8*F(wb,2*(u + v - z),2*u + wb - 2*v - 2*z)*
			       F(-u + v + 2*y + 3*z,wa,3*u - wa - 3*v - 2*y - z)*
			       F(-u + v + 2*y + 3*z,2*(u + v - z),u - v - 2*y + z))/(v*y);
  return 0;
}

