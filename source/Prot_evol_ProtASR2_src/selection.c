//#include "selection.h"
#include <stdlib.h>
#include <math.h>

float Selection_neutral(float *fitness, float DG,
			float DG_thr, float one_over_N){
  // Returns P_fix
  if(DG<DG_thr){*fitness=1; return(one_over_N);}
  else{*fitness=0; return(0);}
}

float Selection_N(float *fitness, float DG, int N_pop,
		  float one_over_N, float fitness_wt)
{
  // If fitness_wt<0 computes only fitness, otherwise returns P_fix
  if(DG > 100){*fitness=0; return(0);}
  *fitness=1./(1+exp(DG));
  if(fitness_wt<0)return(1);
  double f_ratio= fitness_wt/(*fitness);
  if((f_ratio>0.99999999)&&(f_ratio<1.00000001)){
    return(one_over_N);
  }else if(f_ratio > 10000000){
    return(0);
  }else{
    float P_fix= (1.-f_ratio)/(1.-pow(f_ratio, N_pop));
    if(isnan(P_fix)){
      printf("ERROR P_fix= %.2g f_wt= %.2g f= %.2g ratio= %.2g\n",
	     P_fix, fitness_wt, *fitness, f_ratio);
      exit(8);
    }
    return(P_fix);
  }
}
