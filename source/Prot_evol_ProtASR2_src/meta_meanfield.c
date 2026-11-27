#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "REM.h"
#include "allocate.h"
#include "meanfield.h"
#include "optimization.h"
#include "codes.h"
#include "gen_code.h"
#include "externals.h"
#define GET_DG_OPT 0

int Make_meanfield(struct MF_results *opt_res,
		   float **P_MF_ia, float *P_mut_a,
		   float **C_nat, int *i_sec, char *c_sec,
		   float **f_reg_ia, float **f_msa_ia, float *wi,
		   short *aa_seq, struct REM E_wt,
		   FILE *file_out, FILE *file_summary,
		   int L, int Naa, int type, char *name);
void Print_test(char *name_file, struct MF_results res,
		struct REM E_wt, int L, char *what, int open);
void Print_fitness(char *name_file, float *mut_par,
		   struct MF_results res, struct REM E_wt,
		   int L, char *model);
//float Find_max(float *y, float *x, int n, float MIN, float MAX);
//static void Copy_vars(struct MF_results *store, struct MF_results *tmp);
/*static void Optimize_DG(struct MF_results *DG_opt,
			float **P_DG_opt_ia, float *dDG_opt,
			float *P_mut_a, float **P_MF_ia, int L,
			float *Lambda_k, float *y, int N_Lambda,
			float DG_OPT, int **C_nat, int *i_sec,
			float **f_reg_ia, float *wi, struct REM E_wt);*/

int Fixed_Lambda(struct MF_results *MF_res, float **P_MF_ia, float *P_mut_a,
		 float **C_nat, int *i_sec, char *c_sec,
		 float **f_reg_ia, float **f_msa_ia, float *wi,
		 struct REM E_wt, float Lambda, char *name_file,
		 int L, int Naa)
{
  struct MF_results nat_opt, all_opt; 

  printf("### Mean-field computations\n");
  printf("Mean field distribution, native state Lambda=%.3f\n", Lambda);
  meanfield(P_MF_ia, &nat_opt, P_mut_a, NULL, C_nat, i_sec, c_sec,
	    Lambda, E_wt, L, Naa, 0); // Only unfolding
  //nat_opt.lik=Compute_lik(P_MF_ia, f_reg_ia, f_msa_ia, L, Naa);
  Compute_score(&nat_opt, P_MF_ia, f_reg_ia, f_msa_ia, wi, L, Naa);
  Print_test(name_file, nat_opt, E_wt, L, "nat", 0);

  printf("Mean field distribution, native and misfolding Lambda=%.3f\n",
	 Lambda);
  meanfield(P_MF_ia, &all_opt, P_mut_a, NULL, C_nat, i_sec, c_sec,
	    Lambda, E_wt, L, Naa, 1); // Also misfolded state
  //all_opt.lik=Compute_lik(P_MF_ia, L, f_reg_ia);
  Compute_score(&all_opt, P_MF_ia, f_reg_ia, f_msa_ia, wi, L, Naa);
  Print_test(name_file, all_opt, E_wt, L, "MF", 0);
  Print_fitness(name_file, mut_par, all_opt, E_wt, L, "MF");
  *MF_res=all_opt;
  return(0);
}

int Make_meanfield(struct MF_results *opt_res,
		   float **P_MF_ia, float *P_mut_a,
		   float **C_nat, int *i_sec, char *c_sec,
		   float **f_reg_ia, float **f_msa_ia, float *wi,
		   short *aa_seq, struct REM E_wt,
		   FILE *file_out, FILE *file_summary,
		   int L, int Naa, int type, char *name)
{
  int end=0, k;
  float Lambda_ini=0.75, Lambda_step=0.10, Lambda=Lambda_ini;
  float Lambda_k[3], score_k[3], Entropy_min=0;
  struct MF_results res; 
  float eps=0.0001;

  // Meanfield
  float **P_opt_ia=Allocate_mat2_f(L, Naa);
  printf("\n### Mean-field computation type %s\n", name);
  fprintf(file_out, "### Mean-field computation type %s\n", name);
  fprintf(file_out, "# 1=Lambda 2=DG "
	  "3=score 4=KL_model 5=KL_reg 6=lik_reg 7=h 8=Tf 9=conv\n");
  // Initialize
  res.conv=meanfield(P_MF_ia, &res, P_mut_a, NULL, //P_MF_ia,
		     C_nat, i_sec, c_sec, Lambda, E_wt, L, Naa, type);
  Compute_score(&res, P_MF_ia, f_reg_ia, f_msa_ia, wi, L, Naa);
  fprintf(file_out, "%.2f %6.1f %.4f %.4f %.4f %.4f %.3f %.2f %d\n",
	  Lambda, res.DG,
	  res.score, res.KL_mod, res.KL_reg, res.lik_reg,
	  res.h, res.Tf, res.conv);
  score_k[1]=res.score; Lambda_k[1]=Lambda;
  *opt_res=res; opt_res->Lambda[0]=Lambda; opt_res->Lambda[1]=0;
  Copy_P(P_opt_ia, P_MF_ia, L, Naa);
  int n_down=0;
  //
  printf("Mean field %s distribution\n", name);
  for(int iter=0; iter<100; iter++){
    if(n_down==2)break;
    Lambda_k[0]=Lambda_k[1]-Lambda_step;
    Lambda_k[2]=Lambda_k[1]+Lambda_step; Lambda_step*=0.7;
    for(k=0; k<=2; k+=2){
      Lambda=Lambda_k[k];
      res.conv=meanfield(P_MF_ia, &res, P_mut_a, NULL, //P_MF_ia,
			 C_nat, i_sec, c_sec, Lambda, E_wt, L, Naa, type);
      Compute_score(&res, P_MF_ia, f_reg_ia, f_msa_ia, wi, L, Naa);
      fprintf(file_out, "%.2f %6.1f %.4f %.4f %.4f %.4f %.3f %.2f %d\n",
	      Lambda, res.DG,
	      res.score, res.KL_mod, res.KL_reg, res.lik_reg,
	      res.h, res.Tf, res.conv);

      score_k[k]=res.score;
      if((res.entropy<=Entropy_min)||(isnan(res.score)))end=1;
      if(res.score > opt_res->score){
	*opt_res=res; opt_res->Lambda[0]=Lambda;
	Copy_P(P_opt_ia, P_MF_ia, L, Naa);
      }
    }
    if(end)break;
    Lambda=Find_max_quad(Lambda_k[0],Lambda_k[1],Lambda_k[2],
			 score_k[0], score_k[1],score_k[2], 0, 100);
    res.conv=meanfield(P_MF_ia, &res, P_mut_a, NULL, // P_MF_ia,
		       C_nat, i_sec, c_sec, Lambda, E_wt, L, Naa, type);
    Compute_score(&res, P_MF_ia, f_reg_ia, f_msa_ia, wi, L, Naa);
    fprintf(file_out, "%.2f %6.1f %.4f %.4f %.4f %.4f %.3f %.2f %d\n",
	    Lambda, res.DG,
	    res.score, res.KL_mod, res.KL_reg, res.lik_reg,
	    res.h, res.Tf, res.conv);
    Lambda_k[1]=Lambda;
    score_k[1]=res.score;
    if((res.entropy<=Entropy_min)||(isnan(res.score)))break;
    if(res.score > opt_res->score){
      float score_old=opt_res->score;
      *opt_res=res; opt_res->Lambda[0]=Lambda;
      Copy_P(P_opt_ia, P_MF_ia, L, Naa); 
      if(fabs(res.score-score_old)<eps)break;
      n_down=0;
    }else{
      n_down++;
    }
  }

  opt_res->Lambda[1]=0;
  printf("Maximum score %s model: %.2f optimal Lambda=%.3f\n",
	 name, opt_res->score, opt_res->Lambda[0]);
  char label[10];
  if(type==0){strcpy(label,"nat");}else{strcpy(label,"MF");}
  fprintf(file_out, "# Max_score_meanfield_%s: ",label);
  fprintf(file_out, "Lambda= %.3f lik= %.4f DG= %.1f h=%.3f\n",
	  opt_res->Lambda[0], opt_res->score, opt_res->DG, opt_res->h);
  Print_results(*opt_res, label, file_summary);
  Copy_P(P_MF_ia, P_opt_ia, L, Naa);
  Empty_matrix_f(P_opt_ia, L);
  return(0);
}

int Optimize_Lambda(struct MF_results *MF_res,
		    float **P_MF_ia, float *P_mut_a,
		    float **C_nat, int *i_sec, char *c_sec,
		    float **f_reg_ia, float **f_msa_ia, float *wi,
		    short *aa_seq, struct REM E_wt, float DG_OPT,
		    int GET_FREQ, char *MODEL, char *name_file,
		    FILE *file_summary, int L, int Naa)
{

  // Open output file
  char name[500]; sprintf(name, "%s_Lambda.dat", name_file);
  FILE *file_out=fopen(name, "w");
  fprintf(file_out, "# Optimization of Lambda for the mean-field model\n");
  fprintf(file_out, "# T=%.3f SU= %.2f SC= %.2f\n",
	  E_wt.T,  E_wt.S_U,  E_wt.S_C);
  if(SEC_STR){fprintf(file_out, "# Sec.str. used, factor= %.3f\n",SEC_STR);}
  else{fprintf(file_out, "# Sec.str. not used\n");}
  char sfreq[200];
  if(GET_FREQ==0){
    sprintf(sfreq, "# Mutation model (nucs) obtained from input\n");
  }else if (GET_FREQ==1){
    sprintf(sfreq, "# Mutation model (nucs) obtained from fit\n");
  }else if (GET_FREQ==2){
    sprintf(sfreq, "# Mutation model (nucs) obtained from fit + aa freq\n");
  }else if (GET_FREQ==3){
    sprintf(sfreq, "# Global amino acid frequencies obtained from MSA\n");
  }
  fprintf(file_out, "%s", sfreq);
  // Wild type sequence
  double h_wt=0;
  for(int i=0; i<L; i++){if(aa_seq[i]>=0)h_wt+=hydro[aa_seq[i]];} h_wt/=L;
  fprintf(file_out, "# DG(PDB_sequence)= %.2f h_wt= %.3f\n",
	  E_wt.DeltaG, h_wt);
  fprintf(file_out, "# Optimizing Lambda L=%d\n", L);
  /*fprintf(file_out, "# Mutation model: DG= %.2f lik= %.4f h= %.3f Tf= %.2f\n",
    mut_res.DG, mut_res.lik, mut_res.h, mut_res.Tf);*/
  fprintf(file_out, "# Lambda DG_ score h Tf_nat Conv_all\n");

  // Native meanfield
  if(1){
    Make_meanfield(MF_res,P_MF_ia,P_mut_a,C_nat,i_sec, c_sec,
		   f_reg_ia, f_msa_ia, wi, aa_seq, E_wt,
		   file_out, file_summary, L, Naa, 0, "unfolding");
  }

  // Native+misfolding meanfield
  Make_meanfield(MF_res,P_MF_ia,P_mut_a, C_nat, i_sec, c_sec,
		 f_reg_ia, f_msa_ia, wi, aa_seq, E_wt,
		 file_out, file_summary, L, Naa, 1, "misfolding");

  Print_fitness(name_file, mut_par, *MF_res, E_wt, L, "MF");

  /*if(GET_DG_OPT){
    float y[N_Lambda];
    for(k=0; k<N_Lambda; k++)y[k]=-fabs(all_res[k].DG-DG_OPT); 
    Lambda=Find_max(y, Lambda_k, N_Lambda, Lambda_ini/2, Lambda_end);
    printf("Complete model with DG ~ %.3f optimal Lambda=%.3f\n",
	   DG_OPT, Lambda);
    meanfield(P_MF_ia, &tmp, P_mut_a, P_MF_ia,
	      C_nat, i_sec, Lambda, E_wt, L, Naa, 1);
    if(fabs(all_opt.DG-DG_OPT) < dDG_opt){
      dDG_opt=fabs(all_opt.DG-DG_OPT);
      Copy_P(P_DG_opt_ia, P_MF_ia, L, Naa);
      //tmp.lik=Compute_lik(P_MF_ia, L, f_reg_ia);
      Compute_score(&tmp, P_MF_ia, L, f_reg_ia, f_msa_ia, wi);
      tmp.h=Hydro_ave(P_MF_ia, hydro, L);
      DG_opt=tmp; //Copy_vars(&DG_opt, &tmp);
    }
    // Interpolate quadratically to identify DG_OPT
    Optimize_DG(&DG_opt, P_DG_opt_ia, &dDG_opt, P_mut_a, P_MF_ia, L,
		Lambda_k, y, N_Lambda, DG_OPT, C_nat, i_sec,f_reg_ia,wi,E_wt);
    fprintf(file_out, "# Optimal_DG_meanfield_all: ");
    fprintf(file_out, "Lambda= %.3f lik= %.4f DG= %.1f h=%.3f\n",
	    DG_opt.Lambda[0], DG_opt.lik, DG_opt.DG, DG_opt.h);
	    }*/

  printf("### End of mean-field computations\n");
  printf("Writing %s\n", name);
  fclose(file_out);

  // Output distribution
  /*printf("Meanfield model, type %s\n", MODEL);
  if(strcmp(MODEL, "NAT")==0){
    Copy_P(P_MF_ia, P_nat_opt_ia, L, Naa); *MF_res=nat_opt;
  }else if(strcmp(MODEL, "ALL")==0){
    Copy_P(P_MF_ia, P_all_opt_ia, L, Naa); *MF_res=all_opt;
  }else if(strcmp(MODEL, "DG")==0){
    Copy_P(P_MF_ia, P_DG_opt_ia, L, Naa);  *MF_res=DG_opt;
  }else{
    printf("ERROR, %s is not an allowed optimization criterion.\n", MODEL);
    printf("Allowed criteria are NAT ALL and DG\n");
    exit(8);
    }*/

  return(0);
}

void Print_results(struct MF_results res, char *what, FILE *file_out)
{
  //Model Lambda likelihood KL DG Tf h
  // Previously it was -Entr_ave-res.lik_reg
  fprintf(file_out, "%s\t"
	  "%.4f\t%.4f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.1f\t%.2f\t%.3f\n",
	  what, res.Lambda[0], res.Lambda[1], res.lik_MSA, res.KL_mod,
	  res.KL_reg, -res.score, res.KL_mut, res.entropy,
	  res.DG, res.Tf, res.h);
}

/*void Optimize_DG(struct MF_results *DG_opt, float **P_DG_opt_ia,
		 float *dDG_opt, float *P_mut_a, float **P_MF_ia,
		 int L, float *Lambda_k, float *y, int N_Lambda, float DG_OPT,
		 int **C_nat, int *i_sec, float **f_reg_ia, float *wi,
		 struct REM E_wt)
{
  // Iterations to find DG_OPT
  int IT_MAX=4, iter=0, k;
  float L0, L1, L2;
  float y0, y1, y2;
  float Lambda=DG_opt->Lambda[0],
    Lambda_ini=Lambda_k[0], Lambda_end=Lambda_k[N_Lambda-1]+0.2;
  struct MF_results tmp;
  while( (((DG_opt->DG)>0)||((DG_opt->DG)<2*DG_OPT))&&(iter<IT_MAX)){
    printf("Lambda= %.3f DG= %.2f iter=%d\n", Lambda, DG_opt->DG, iter);
    iter++;
    float yy=-fabs(DG_opt->DG-DG_OPT);
    if(iter==1){
      if(Lambda < Lambda_k[0]){
	L0=Lambda, L1=Lambda_k[0], L2=Lambda_k[1];
	y0=yy; y1=y[0]; y2=y[1];
      }else if(Lambda > Lambda_k[N_Lambda-1]){
	L0=Lambda_k[N_Lambda-2]; L1=Lambda_k[N_Lambda-1]; L2=Lambda;
	y0=y[N_Lambda-2]; y1=y[N_Lambda-1]; y2=yy;
      }else{
	for(k=1; k<N_Lambda; k++)if(Lambda_k[k]>Lambda)break;
	if(k==N_Lambda)k--;
	if((Lambda < Lambda_k[k-1])||(Lambda > Lambda_k[k])){
	  printf("WARNING, badly ordered series of Lambda: %.3f %.3f %.3f\n",
		 Lambda_k[k-1], Lambda, Lambda_k[k]); break;
	}
	L0=Lambda_k[k-1]; L1=Lambda; L2=Lambda_k[k];
	y0=y[k-1]; y1=yy; y2=y[k];
      }
    }else{
      if(Lambda < L0){
	L2=L1; y2=y1; L1=L0; y1=y0; L0=Lambda; y0=yy;
      }else if(Lambda < L1){
	L2=L1; y2=y1; L1=Lambda; y1=yy;
      }else if(Lambda < L2){
	L0=L1; y0=y1; L1=Lambda; y1=yy;
      }else{
	L0=L1; y0=y1; L1=L2; y1=y2; L2=Lambda; y2=yy;
      }
    }
    printf("L0= %.2f L1= %.2f L2= %.2f\n", L0, L1, L2);
    printf("y0= %.2f y1= %.2f y2= %.2f\n", y0, y1, y2);
    
    Lambda =Find_max_quad(L0,L1,L2,y0,y1,y2,Lambda_ini/4,Lambda_end*10);
    printf("Performing attempt %d to optimize Lambda, Lambda=%.3f\n",
	   iter+1, Lambda);
    meanfield(P_MF_ia, &tmp, P_mut_a, P_MF_ia,
	      C_nat, i_sec, Lambda, E_wt, L, Naa, 1);
    float ddG=-fabs(tmp.DG-DG_OPT);
    if(ddG > *dDG_opt){
      *dDG_opt=ddG;
      Copy_P(P_DG_opt_ia, P_MF_ia, L, Naa);
      //tmp.lik=Compute_lik(P_MF_ia, L, f_reg_ia);
      Compute_score(&tmp, P_MF_ia, L, f_reg_ia, wi);
      tmp.h=Hydro_ave(P_MF_ia, hydro, L);
      *DG_opt=tmp; //Copy_vars(DG_opt, &tmp);
    }
  }
}
*/

void Print_fitness(char *name_file, float *mut_par,
		   struct MF_results res, struct REM E_wt,
		   int L, char *model)
{
  char name[500]; sprintf(name, "%s_fitness.dat", name_file);
  FILE *file_out=fopen(name, "w");
  double f=1./(1+exp(res.DG)), LL;
  printf("Writing %s\n", name);
  fprintf(file_out, "#Model %s DG_WT= %.3f\n", model, E_wt.DeltaG);
  if(SEC_STR)
    fprintf(file_out, "#Local interactions used, factor= %.3f\n",SEC_STR);
  fprintf(file_out, "# tt_ratio kCpG mu  TEMP Lambda_eff = ");
  fprintf(file_out, " %.2f %.1f %.4f", mut_par[4], mut_par[5], mut_par[6]);
  fprintf(file_out, "   %.2f ", E_wt.T);
  fprintf(file_out, "   %.3f ", res.Lambda[0]);
  fprintf(file_out, "\n");
  fprintf(file_out, "#len   G C A T  ");
  fprintf(file_out, " DG fitness likelihood PDB Lambda\n");
  fprintf(file_out, "%d   %.3f %.3f %.3f %.3f", L,
	  mut_par[Code_nuc('G')], mut_par[Code_nuc('C')],
	  mut_par[Code_nuc('A')], mut_par[Code_nuc('T')]);
  fprintf(file_out, "   %.2f ", res.DG*E_wt.T);
  fprintf(file_out, "   %.6f ", f);
  fprintf(file_out, "   %.3f ", res.lik_reg);
  name[5]='\0';
  fprintf(file_out, " \"%s\"", name);
  LL=log(res.Lambda[0]*E_wt.T);
  if(f >= 0.5){LL-=res.DG;} // log(Lambda*T/(f(1-f));
  else{LL-=log(0.25);}
  fprintf(file_out, "   %.4g\n", LL);
  fclose(file_out);
}

void Print_test(char *name_file, struct MF_results res,
		struct REM E_wt, int L, char *what, int open)
{
  char name[500]; FILE *file_out;
  sprintf(name, "%s_likelihood.dat", name_file);
  if(open){
    file_out=fopen(name, "w");
    printf("Writing %s\n", name);
    fprintf(file_out, "# Third moment of energy used: ");
    if(E_wt.REM==3){fprintf(file_out,"YES\n");}
    else{fprintf(file_out,"NO\n");}
    fprintf(file_out, "#Model DG/T log_lik/L Lambda= %.2f L=%d\n",
	    res.Lambda[0], L);
    fprintf(file_out, "DeltaG/T_%s Pairwise= %7.4f\n", what, E_wt.DeltaG);
  }else{
    file_out=fopen(name, "a");
  }
  fprintf(file_out, "%s %7.4f %7.4f\n", what, res.DG, res.lik_reg);
  fclose(file_out);
}

/*float Find_max(float *y, float *x, int n, float MIN, float MAX)
{
  int i, k=0; float ymax=y[0];
  for(i=1; i<n; i++){
    if(isnan(ymax)==0||isfinite(ymax))break; ymax=y[i]; k=i;
  }
  for(i=1; i<n; i++)if(y[i]>ymax){ymax=y[i]; k=i;}
  if(k==0){
    printf("WARNING, maximum attained at lowest x (x= %.3g max= %.3g)\n",
	   x[k], y[k]); k=1;
    //for(i=0; i<n; i++)printf("%.3g %.3g\n", x[i], y[i]);
  }else if(k==(n-1)){
    printf("WARNING, maximum attained at highest x (x= %.3g max= %.3g)\n",
	   x[k], y[k]); k=n-2;
    //for(i=0; i<n; i++)printf("%.3g %.3g\n", x[i], y[i]);
  }
  
  float x0=Find_max_quad(x[k-1], x[k], x[k+1], y[k-1], y[k], y[k+1],MIN,MAX);
  return(x0);
  }*/

float Hydro_ave(float **P_ia, float *hydro, int L){
  int i, a; double sum=0, h;
  for(i=0; i<L; i++){
    float *P_a=P_ia[i];
    h=0; for(a=0; a<20; a++)h+=P_a[a]*hydro[a];
    sum+=h;
  }
  return(sum/L);
}

/*void Copy_vars(struct MF_results *store, struct MF_results *tmp)
{
  store->DG=tmp->DG; store->lik_reg=tmp->lik_reg; store->Lambda=tmp->Lambda; 
  store->Tf=tmp->Tf;  store->KL_mut=tmp->KL_mut;  store->h=tmp->h;
  store->KL_reg=tmp->KL_reg; store->entropy=tmp->entropy;
  store->score=tmp->score;
  }*/

void Compute_P_sel(float *P_sel, float **P_MF_ia, float *P_mut_a,
		   int L, int Naa){
  int i, a;
  for(a=0; a<Naa; a++){
    double sum=0;
    for(i=0; i<L; i++)sum+= P_MF_ia[i][a];
    P_sel[a]=sum/(L*P_mut_a[a]);
  }
}

