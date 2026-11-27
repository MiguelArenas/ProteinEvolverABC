#include "codes.h"
#include "random3.h" 
#include <stdio.h>
#include <stdlib.h>
#include "mutation.h"

char Mutate_nuc(char nuc, float *freq_nuc, float *rate, float tt_ratio);

void Ini_count(char *dna_seq, int len_dna, int *count){
  int i;
  for(i=0; i<4; i++)count[i]=0;
  for(i=0; i<len_dna; i++){
    count[Code_nuc(dna_seq[i])]++;
  }
} 

int Mutate_seq(char *dna_seq, int len_dna, char **codon, char *coded_aa,
	       short *aa_seq, int len_nat,
	       float *freq_nuc, float tt_ratio, int *count, float *rate, 
	       int *nuc_mut, char *nuc_new, int *res_mut, int *aa_new)
{
  int i; float ran, norm=0, sum=0;
  char nuc_old;

  for(i=0; i<4; i++)norm+=count[i]*rate[i];
  ran=RandomFloating();  //ran=RandomVar0();
  ran*=norm;

  for(i=0; i<len_dna; i++){
    sum+=rate[Code_nuc(dna_seq[i])];
    if(sum>=ran){
      *nuc_new=Mutate_nuc(dna_seq[i], freq_nuc, rate, tt_ratio);
      *nuc_mut=i; break;
    }
  }

  *res_mut=(*nuc_mut)/3; i=(*res_mut)*3; nuc_old=dna_seq[*nuc_mut];
  dna_seq[*nuc_mut]=*nuc_new;
  *aa_new=Coded_aa(dna_seq+i, codon, coded_aa);
  dna_seq[*nuc_mut]=nuc_old;

  if(*aa_new==-1)return(-1);
  if((*aa_new)==aa_seq[*res_mut])return(1);

  return(0);
}

char Mutate_nuc(char nuc, float *freq_nuc, float *rate, float tt_ratio){
  int i_nuc=Code_nuc(nuc), j_nuc;
  float ran, p, sum=0;
  ran= RandomFloating(); if(ran >=1)ran=0.999999;
  //ran= RandomVar0();
  ran*=rate[i_nuc];
  for(j_nuc=0; j_nuc<4; j_nuc++){
    if(j_nuc==i_nuc)continue;
    p=freq_nuc[j_nuc];
    if(j_nuc==Transition(NUC_CODE[i_nuc]))p*=tt_ratio;
    sum+=p; if(ran<=sum)return(NUC_CODE[j_nuc]);
  }
  printf("Error in mutate, nuc=%c, sum=%.4f ran=%.4f\n", nuc, sum, ran);
  exit(8);
}


int Compute_rates(float *rate, float *freq_nuc, float trans_ratio){
  float k=trans_ratio-1; int i;

  printf("Stationary frequencies: ");
  for(i=0; i<4; i++)printf("%c= %.3f ", NUC_CODE[i], freq_nuc[i]);
  printf("\n");
  printf("Mutation rates: ");
  for(i=0; i<4; i++){
    rate[i]=1-freq_nuc[i]+k*freq_nuc[Transition(NUC_CODE[i])];
    printf("%c= %.3f ", NUC_CODE[i], rate[i]);
  }
  printf("\n");
  return(0);
}
