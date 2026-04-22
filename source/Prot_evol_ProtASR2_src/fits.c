#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include "fits.h"

double Logistic_fit(double *tau, double *r, float *yy, int step, int n)
{
  // y(t)=y(0)+(y_inf-y(0))(1-exp(-t/tau)) =>
  // (y(t+s)-y(0)) = (y(t)-y(0))*exp(-s/tau)+(y_inf-y(0))(1-exp(-s/tau))
  int m=n-2, i; 
  double xy=0, x1=0, y1=0, x2=0, y2=0;
  float x=yy[1]-yy[0], y;
  for(i=2; i<n; i++){
    y=yy[i]-yy[0];
    x1+=x; x2+=x*x; xy+=x*y;
    y1+=y; y2+=y*y;
    x=y;
  }
  x2=m*x2-x1*x1; y2=m*y2-y1*y1; xy=m*xy-x1*y1;
  if((x2<=0)||(y2<=0)){
    printf("WARNING in logistic fit, all data are zero\n");
    *tau=0; *r=0; return(0);
  }
  *r=xy/sqrt(x2*y2);
  float slope=xy/x2, offset=(y1-(slope)*x1)/m;
  if((slope<=0)||(slope>=1)){
    printf("WARNING in logistic fit, slope= %.3f not in ]0,1[\n", slope);
    *tau=0; 
  }else{
    (*tau)=-step/log(slope);
  }
  double y_inf=yy[0]+offset/(1.-slope);
  return(y_inf);
}
