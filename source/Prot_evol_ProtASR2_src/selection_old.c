#include <math.h>
#include <stdio.h>
#include "selection.h"
#include "random3.h"
#define LOG_2  0.69314718

int ini;

float Fitness_pow_E_a(float x_E, float x_alpha, float sel_coeff)
{
  // if((E>=0)||(alpha<=0)){fitness=0;}else{
  // fitness = 1/[1+(E/E0)^-S+(a/a0)^-S]

  float f;
  if((x_alpha<=0)||(x_E<=0))return(0);
  f=1+pow(x_E,-sel_coeff)+pow(x_alpha,-sel_coeff);
  return(1/f);
}


float Fitness_pow_min(float x_E, float x_alpha, float sel_coeff)
{
  // if((E>=0)||(alpha<=0)){fitness=0;}else{
  // fitness = 1/[1+min(E/E0,a/a0)^-S]

  float f, x;
  if((x_alpha<=0)||(x_E<=0))return(0);
  if(x_alpha<x_E){x=x_alpha;}else{x=x_E;}
  f=1+pow(x,-sel_coeff);
  return(1/f);
}


float Fitness(float E, float E_thr, float alpha, float alpha_thr,
	      float sel_coeff, char *fit_def){
  // if((E>=0)||(alpha<=0)){fitness=0;}else{
  // fitness = [1+exp(-s)]/[1+exp(s*max(E0/E-1,a0/a-1)]
  float x, y;
  if(ini==0){
    sprintf(fit_def,
	    "0 if(E>=0)||(a<=0); [1+exp(-s)]/[1+exp(s*max(E0/E-1,a0/a-1)]\0");
    ini=1;
  }

  if((E>=0)||(alpha<=0))return(0);
  x=E_thr/E-1;  y=alpha_thr/alpha-1;
  if(y>x)x=y;
  return((1.+exp(-sel_coeff))/(1+exp(sel_coeff*x)));
}



float Fitness_Andreas(float E, float E_thr, float alpha, float alpha_thr,
		      float sel_coeff, char *fit_def){
  // fitness = 1/[1+exp(s(E-E0)/|E0|)+exp(s(a0-a)/a0)]
  float x=1.-E/E_thr, y=1.-alpha/alpha_thr;
  if(ini==0){
    sprintf(fit_def, "1/[1+exp(s(E-E0)/|E0|)+exp(s(a0-a)/a0)]\0");
    ini=1;
  }
  return(1./(1+exp(sel_coeff*x)+exp(sel_coeff*y)));
}

float Fitness_Ugo(float E, float E_thr, float alpha, float alpha_thr,
		  float sel_coeff, char *fit_def){
  // fitness = 1/[1+exp(s*max((E-E0)/|E0|,(a0-a)/a0)]
  float x=1.-E/E_thr, y=1.-alpha/alpha_thr;
  if(ini==0){
    sprintf(fit_def, "1/[1+exp(s*max((E-E0)/|E0|,(a0-a)/a0)]\0");
    ini=1;
  }
  if(y>x)x=y;
  return(1./(1+exp(sel_coeff*x)));
}

int Selection(float fitness, float fitness_old, int N_pop)
{
  // Moran's process:
  /* P_fix = (1-exp(-Df))/(1-exp(-N*Df)) */
  double f_ratio, x;
  if((int)N_pop==1)return(1);
  if(fitness<=0)return(0);
  f_ratio= fitness_old / fitness;
  x= (1.-f_ratio)/(1.-pow(f_ratio, N_pop));
  if(RandomFloating() < x)return(1); return(0);
}


int Selection_old(float E, float E_thr, float alpha, float alpha_thr,
		  float zz, float zz_thr)
{
  if((E<=E_thr)&&(alpha>=alpha_thr)&&(zz<=zz_thr))return(1);
  return(0);
}
