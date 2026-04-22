#include "wag.h"
#include "jtt.h"
#include "gen_code.h"
#include "diagonalize.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

static void Analyze_exchange(float *eigen_value, float **eigen_vector,
			     char *MODEL, float *mut_par);
static void Store_exchange(float **exch, char *MODEL);

// Genetic code
// char *codon[64], coded_aa[64]

int main(int argc, char **argv){

  float mut_par[7];

  float eigen_value[20], *eigen_vector[20]; int a;
  for(a=0; a<20; a++)eigen_vector[a]=malloc(20*sizeof(float));

  Analyze_exchange(eigen_value, eigen_vector, "JTT", mut_par);
  Analyze_exchange(eigen_value, eigen_vector, "WAG", mut_par);
  Analyze_exchange(eigen_value, eigen_vector, "MUT", mut_par);


}

void Analyze_exchange(float *eigen_value, float **eigen_vector,
		      char *MODEL, float *mut_par)
{
  float *exch[20]; int a;
  for(a=0; a<20; a++)exch[a]=malloc(20*sizeof(float));
  if(strncmp(MODEL, "MUT", 3)==0){
    int GET_FREQ=0, CpG=0;
    float  P_mut_a[20], *Q_cod[64], P_cod[64];
    for(a=0; a<64; a++)Q_cod[a]=malloc(64*sizeof(float));
    Get_mut_par(mut_par, P_mut_a, P_cod, Q_cod, GET_FREQ,
		0, NULL, 0, NULL, CpG);
    Compute_exchange_mut(exch, P_cod, Q_cod);
  }else{
    Store_exchange(exch, MODEL);
  }
  f_Diagonalize(20, exch, eigen_value, eigen_vector, 1);
  for(a=0; a<20; a++)free(exch[a]);
}

void Store_exchange(float **exch, char *MODEL)
{
  float *exch1[20]; int a;
  if(strncmp(MODEL, "WAG", 3)==0){
    for(a=0; a<20; a++)exch1[a]=WAG[a];
  }else{
    strcpy(MODEL, "JTT");
    for(a=0; a<20; a++)exch1[a]=JTT[a];
  }
  for(a=0; a<20; a++){
    for(int b=0; b<20; b++){
      if(b<a){exch[a][b]=exch1[a][b];}
      else{exch[a][b]=exch1[b][a];}
    }
  }
}
