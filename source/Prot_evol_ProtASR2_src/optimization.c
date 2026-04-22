#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "optimization.h"

float Find_max_quad(float x1, float x2, float x3,
		    float y1, float y2, float y3,
		    float MIN, float MAX)
{
  //printf("x1= %.2f x2= %.2f x3= %.2f\n", x1, x2, x3);
  //printf("y1= %.2f y2= %.2f y3= %.2f\n", y1, y2, y3);
  /*if((y2<=y1)&&(y2<=y3)){ // Minimum
    //printf("WARNING, minimum\n");
    if(y1>=y3){
      return((x1+MIN)/2);
    }else{
      return((x3+MAX)/2);
    }
    }*/

  /*if(x1==x2){
    printf("WARNING, singular values x1= %.6g x2= %.6g x3= %.3g\n",x1,x2,x3);
    return((x2+x3)/2);
  }else if(x1==x3){
    printf("WARNING, singular values x1= %.6g x2= %.6g x3= %.3g\n",x1,x2,x3);
    return((x1+x2)/2);
  }else if(x2==x3){
    printf("WARNING, singular values x1= %.6g x2= %.6g x3= %.3g\n",x1,x2,x3);
    return((x1+x3)/2);
    }*/

  float yprime_2=((y1-y2)/(x1-x2));
  float yprime_3=((y1-y3)/(x1-x3));
  if(yprime_2 == yprime_3){
    // Singular data, return maximum value
    printf("WARNING in Find_max_quad, equal derivatives\n");
    float xmax=x2, ymax=y2;
    if(y1>ymax){ymax=y1; xmax=x1;}
    if(y3>ymax){xmax=x3;}
    return(xmax);
  }
  float a=(yprime_2-yprime_3)/(x2-x3);
  if(a>=0 || isnan(a)){
    // Second derivative is positive, no maximum!
    //printf("WARNING, second derivative >0\n");
    printf("WARNING in Find_max_quad, second derivative %.2g >0\n", a);
    float xmax=x2, ymax=y2;
    if(y1>ymax){ymax=y1; xmax=x1*0.99;}
    if(y3>ymax){xmax=x3*1.01;}
    return(xmax);
  }
  float x0=0.5*(x1+x2)-yprime_2/(2*a);
  //printf("a= %.2f x0= %.2f\n", a, x0);
  if(x0<MIN){x0=(x1+MIN)/2;}
  //else if(x0>MAX){x0=(x3+MAX)/2;}
  return(x0);
}
    
void Rearrange_points(float *x1, float *x2, float *x3,
		      float *y1, float *y2, float *y3,
		      float x0, float y0)
{
  if(y0<*y1){
    *x3=*x2; *y3=*y2;
    *x2=*x1; *y2=*y1;
    *x1=x0;  *y1=y0;
  }else if(y0<*y2){
    *x3=*x2; *y3=*y2;
    *x2=x0;  *y2=y0;
  }else if(y0<*y3){
    *x1=*x2; *y1=*y2;
    *x2=x0; *y2=y0;
  }else{
    *x1=*x2; *y1=*y2;
    *x2=*x3; *y2=*y3;
    *x3=x0;  *y3=y0;
  }
}


