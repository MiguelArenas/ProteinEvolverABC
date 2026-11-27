#include "REM.h"
#include "subst_exhaustive.h"
#include "codes.h"
#include "gen_code.h"
#include "random3.h"           /* Generating random numbers */
//#include "selection.h"
#include "externals.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Global variables //
float P_large=0.99999;
int VBT=1; // Verbatim
double *P_subst_sum=NULL;
int *knuc_mut=NULL;
int *knuc_new=NULL;
int *kres_mut=NULL;
int *kres_new=NULL;
int *neut=NULL;
int num_mut;
int ini_subst;
double *E1i_mut, *E2i_mut;
double *E1i_wt, *E2i_wt;

int Set_codon_fam(int *cod_family, char **codon, char *coded_aa);
float **Set_trans_prob(int N_Fam, int *cod_family,
		       float *freq_nuc, float trans_ratio);
int Set_trans_CpG(float ***Mut_mat_fam_CpG_00, float ***Mut_mat_fam_CpG_01,
		  float ***Mut_mat_fam_CpG_10, float ***Mut_mat_fam_CpG_11,
		  int N_Fam, int *cod_family, float *freq_nuc,
		  float trans_ratio, float CpG);
long Waiting_time(float *P_nofix);
int Binomial_var(int n, double p);
int Extract_codon(int *codon_new, int j_fam, int *cod_family, float *freq_nuc);
void Record_loads(struct load *load);

int Substitute_exhaustive(int *Res_mut, int *AA_new,
			  int *Nuc_mut, int *Nuc_new,
			  struct load *mut_load, struct load *trans_load,
			  int it_load, double *fitness_mut,
			  struct REM *E_mut, struct REM *E_wt,
			  long *naa_mut, long *nsyn_mut, long *syn_subst,
			  int NEUTRAL, float DG_thr, int N_pop,
			  float **C_nat, int *i_sec, char *c_sec,
			  short *aa_seq,
			  char *dna_seq, short *nuc_seq, int len_dna, 
			  char **codon, char *coded_aa,
			  float tt_ratio, float *mut_par)
{
  // WARNING: E_wt and E_mut must be initialized before calling
  // this routine

  // Allocate
  if(ini_subst==0){
    num_mut=19*E_wt->L;
    P_subst_sum=malloc(num_mut*sizeof(double));
    knuc_mut=malloc(num_mut*sizeof(int));
    knuc_new=malloc(num_mut*sizeof(int));
    kres_mut=malloc(num_mut*sizeof(int));
    kres_new=malloc(num_mut*sizeof(int));
    ini_subst=1;
  }

  // Wild-type stability
  E_wt->DeltaG=
    Compute_DG_overT_contfreq(E_wt, aa_seq, C_nat, i_sec, c_sec, 0);
  // Fitness
  float fitness_wt;
  if(NEUTRAL){
    if(E_wt->DeltaG<DG_thr){fitness_wt=1;}
    else{fitness_wt=0;}
  }else{
    fitness_wt=1./(1+exp(E_wt->DeltaG));
  }

  /*int nuc_C=Code_nuc('C'), nuc_G=Code_nuc('G'), CpG_true=0, l1=len_dna-1;
    float CpG_trans=1.00;*/

  // Store
  int n_aa=0, aa_store[19];
  float P_fix_store[19], fitness_store[19], DG_store[19];

  // Probabilities
  double P_mut_sum=0;
  double P_syn_sum=0; // Sum of synonymous fixations
  double one_over_N=1./N_pop; // Fixation prob. of neutral mutants
  P_subst_sum[0]=0; // sum of P_fix

  // Load calculation
  double f_mut=0, dG_mut=0;
  double f_trans=0, dG_trans=0;
  long w_trans=0;

  // Mutations
  int k_mut=0, res_old=-1, aa_old, i;
  char codon_new[4]; int synonymous;

  for(int i_nuc=0; i_nuc<len_dna; i_nuc++){
    int res_mut=i_nuc/3;
    if(res_mut!=res_old){
      aa_old=aa_seq[res_mut];
      res_old=res_mut; n_aa=0;
    }
    
    // Mutation
    int ibase=nuc_seq[i_nuc];
    int itrans=Transition(dna_seq[i_nuc]);
    //if((ibase==nuc_C)&&(i_nuc<l1)&&(nuc_seq[i_nuc+1]==nuc_G))CpG_true=1;
    for(int base=0; base<4; base++){
      if(base==ibase)continue;
      float P_mut=mut_par[base];
      if(base==itrans)P_mut*=tt_ratio; //if(CpG_true)P_mut*=CpG;
      
      // Amino acid
      int i=res_mut*3;
      for(int j=0; j<3; j++){
	if(i!=i_nuc){codon_new[j]=dna_seq[i];}
	else{codon_new[j]=Nuc_code(base);}
	i++;
      }
      //if(VBT)printf("%d %d %d %s\n", i_nuc, base, res_mut, codon_new);

      float DG, fitness, P_fix;
      int aa_new=Coded_aa(codon_new, codon, coded_aa);
      if(aa_new==aa_old){
	P_fix= one_over_N; fitness=fitness_wt; DG=E_wt->DeltaG;
	synonymous=1;
	goto p_sum;
      }
      synonymous=0;
      if(aa_new<0){
	P_fix=0; fitness=0; DG=0; goto p_sum; // Stop codon
      }
      for(int i=0; i<n_aa; i++){
	if(aa_new==aa_store[i]){
	  P_fix=P_fix_store[i]; fitness=fitness_store[i]; DG=DG_store[i];
	  goto p_sum;
	}
      }
      
      // Folding thermodynamics and fitness
      Copy_E_REM(E_mut, E_wt);
      DG=Mutate_DG_overT_contfreq(E_mut, aa_seq, C_nat, i_sec, c_sec,
				  res_mut, aa_new);
      
      // Fitness
      if(NEUTRAL){
	if(DG<DG_thr){fitness=1; P_fix=one_over_N;}
	else{fitness=0; P_fix=0;}
      }else{
	if(DG > 100){fitness=0; P_fix=0; goto store;}
	fitness=1./(1+exp(DG));
	double f_ratio= fitness_wt/fitness;
	if((f_ratio>0.99999999)&&(f_ratio<1.00000001)){
	  P_fix=one_over_N;
	}else if(f_ratio > 10000000){
	  P_fix=0;
	}else{
	  P_fix= (f_ratio-1.)/(pow(f_ratio, N_pop)-1.);
	}
	if(isnan(P_fix)){
	  printf("ERROR P_fix= %.2g f_wt= %.2g f= %.2g ratio= %.2g t-t0=%d\n",
		 P_fix, fitness_wt, fitness, f_ratio, it_load);
	  exit(8);
	}
      }
      // Store new mutation
    store:
      if(n_aa>19){
	printf("ERROR, too many amino acid changes: %d\n", n_aa);
	exit(8);
      }
      aa_store[n_aa]=aa_new;
      P_fix_store[n_aa]=P_fix;
      fitness_store[n_aa]=fitness;
      DG_store[n_aa]=DG;
      n_aa++;
      
      //if(CpG_true && (ibase == nuc_G))CpG_true=0;
    p_sum: 
      if(synonymous){
	P_syn_sum+=P_mut;
      }else{
	P_mut_sum+=P_mut;
	if(aa_new<0)continue;
	float P_subst=P_mut*P_fix;
	if(k_mut==0){P_subst_sum[0]=P_subst;}
	else{P_subst_sum[k_mut]=P_subst_sum[k_mut-1]+P_subst;}
	knuc_mut[k_mut]=i_nuc; knuc_new[k_mut]=base;
	kres_mut[k_mut]=res_mut; kres_new[k_mut]=aa_new;
	k_mut++;
      }

      // Sample loads
      if(it_load>0){
	if(aa_new>=0){
	  /* Translation load
	     Stop codons do not contribute to the translation load, since
	     the release factor affecting ribosome termination have a
	     smaller error rate than other mistranslations */
	  f_trans+=fitness;
	  dG_trans+=DG;
	  w_trans++;
	}
	// Mutation load
	f_mut+=fitness*P_mut;
	dG_mut+=DG*P_mut;
      }
    }
  }

  if(it_load>0){
    double w_mut=P_mut_sum+P_syn_sum;
    mut_load->df=fitness_wt-f_mut/w_mut;
    mut_load->dG=dG_mut/w_mut-E_wt->DeltaG;
    Record_loads(mut_load);

    trans_load->df=fitness_wt-f_trans/w_trans;
    trans_load->dG=dG_trans/w_trans-E_wt->DeltaG;
    Record_loads(trans_load);
  }
  
  // Extract substitution
  if(k_mut>=num_mut){
    printf("ERROR, unexpected number of mutations %d instead of %d\n",
	   k_mut, num_mut); exit(8);
  }
  double P_subst_aa=P_subst_sum[k_mut-1];
  if(isnan(P_subst_aa)){
    for(i=0; i<k_mut; i++)printf("%.2g", P_subst_sum[i]);
    printf("\nERROR k_mut= %d Sum_P_subst= %.2g t-t0=%d\n",
	   k_mut, P_subst_aa, it_load); exit(8);
  }else if(P_subst_aa==0){
    printf("ERROR k_mut= %d Sum_P_subst= %.2g subst-trans=%d\n",
	   k_mut, P_subst_aa, it_load); i=k_mut/2; exit(8);
  }
  double ran= P_subst_aa*RandomFloating();
  for(i=0; i<k_mut; i++){
    if(P_subst_sum[i]>=ran)break;
  }
  if(i>=k_mut){
    if(ran>P_subst_aa){
      printf("ERROR, substitution not found k_mut=%d ran=%.3f norm=%.3f\n",
	     k_mut, ran, P_subst_aa);
    }
    i=k_mut-1;
  }  
  *Nuc_mut=knuc_mut[i];
  *Nuc_new=knuc_new[i];
  *Res_mut=kres_mut[i];
  *AA_new=kres_new[i];
  
  // Update thermodynamic parameters
  Copy_E_REM(E_mut, E_wt);
  E_mut->DeltaG=
    Mutate_DG_overT_contfreq(E_mut, aa_seq, C_nat, i_sec, c_sec,
			     *Res_mut, *AA_new);
  if(NEUTRAL){
    if(E_mut->DeltaG<DG_thr){*fitness_mut=1;}else{*fitness_mut=0;}
  }else{
    *fitness_mut=1./(1+exp(E_mut->DeltaG));
  }

  // Debugging
  if(0){
    printf("%d\t%d", k_mut, i);
    if(i>0){
      printf("\t%.3g", (P_subst_sum[i]-P_subst_sum[i-1])/P_subst_aa);
    }else{
      printf("\t%.3g", P_subst_sum[0]/P_subst_aa);
    }
    printf("\t%.4f\t%.2f", fitness_wt, E_wt->DeltaG);
    printf("\t%.4f\t%.2f\n", *fitness_mut, E_mut->DeltaG);
    printf("%.1f\t%.2f\t%.2f\t%.2f\t%.2f\n", E_mut->E2,
	   E_mut->E2cont1, E_mut->E2cont2*E_mut->E1,
	   E_mut->E2site1, E_mut->E2site2);
    aa_old=aa_seq[*Res_mut]; aa_seq[*Res_mut]=*AA_new;
    E_mut->DeltaG=
      Compute_DG_overT_contfreq(E_mut, aa_seq, C_nat, i_sec, c_sec, 0);
    aa_seq[*Res_mut]=aa_old;
    printf("DG: %.2f\n", E_mut->DeltaG);
    printf("%.1f\t%.2f\t%.2f\t%.2f\t%.2f\n", E_mut->E2,
	   E_mut->E2cont1, E_mut->E2cont2*E_mut->E1,
	   E_mut->E2site1, E_mut->E2site2);
    printf("\n");
  }

  // Extract failed and synonymous mutations 
  P_syn_sum*=one_over_N;
  float P_nofix= P_syn_sum/(P_subst_aa+P_syn_sum);
  *syn_subst=Waiting_time(&P_nofix);
  P_nofix=1.-one_over_N;
  *nsyn_mut=0;
  for(i=0; i< *syn_subst; i++){
    *nsyn_mut+=Waiting_time(&P_nofix);
  }
  P_nofix= 1.-  P_subst_aa/P_mut_sum;
  *naa_mut=Waiting_time(&P_nofix);

  return(0);
}

long Waiting_time(float *P_nofix){
  // P(n)=(1-P_nofix)P_nofix^n
  // n < log(1-ran)/log(P_nofix) < n+1
  double ran=RandomFloating();
  if(ran > P_large)ran=P_large;
  if(*P_nofix> P_large)*P_nofix=P_large;
  long n= log(1.-ran)/log(*P_nofix);
  //if(n<0){printf("n= %ld P_fix= %.5g\n", n, *P_nofix);}
  return(n);
}

int Binomial_var(int n, double p){
  double q=1.-p, ran=RandomFloating(), pp=1, pq, P_sum=0;
  int k;
  if(n==0)return(0);
  pp=pow(q, n); pq=p/q;
  for(k=0; k<=n; k++){
    P_sum+=pp;
    if(ran <= P_sum)break;
    pp*=(pq*(float)(n-k)/(float)(k+1));
  }
  if(k>n){
    k=n;
    if((ran-P_sum)>0.0001){ 
      printf("ERROR in Binomial_var! n=%d p=%.3f k=%d\n", n, p, k+1);
      exit(8);
    }
  }
  return(k);
}

void Record_loads(struct load *load){
    load->df_ave+=load->df;
    load->df_dev+=load->df*(load->df);
    load->dG_ave+=load->dG;
    load->dG_dev+=load->dG*(load->dG);
    load->num++;
}
