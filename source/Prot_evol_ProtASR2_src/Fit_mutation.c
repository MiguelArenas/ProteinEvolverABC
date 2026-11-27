#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "codes.h"
#include "input.h"
#include "gen_code.h"
#include "allocate.h"
#include "optimization.h"
#include "externals.h"

// mut_par: 0=A 1=T 2=G 3=C 4=CpG 5=tt_ratio 6=TWONUC 7=(1-2*CpG) 8=CpG*[GC]

int **neighbor=NULL;
int n_neigh=0;
int Stop_codon[64];
float twonuc;
int VBS=1;

// Multiplicity of nucleotides in codons
int INI_MULT=0;
int Multi[64][4];
int Mut_CpG[64];

// Computation of factorials
#define Nf 11
float lfact[Nf];
float logfact_aa[21];
float logfact_dna[5];
float log2pi_2;
int INI_FACT_aa=0;
int INI_FACT_dna=0;
int INI_FACT=0;

void Minimize_AIC(float *x, int Npar, float *lmax,
		  float *P_mut, float *P_cod, float **Q_cod,
		  char **codon, char *coded_aa,
		  float *num_aa, float *num_dna, int L_aa,
		  FILE *file_out);
float Compute_AIC(float lik, int np, int no);
float Steepest_ascent(float *x, int *freezed, int N,
		      float *P_mut, float *P_cod, float **Q_cod,
		      char **codon, char *coded_aa,
		      float *num_aa, float *num_dna);
float Quadratic_line_search(float *dF, float *x, float F0, int N,
			    float *P_mut, float *P_cod, float **Q_cod,
			    char **codon, char *coded_aa,
			    float *num_aa, float *num_dna);
int Initialize_codon_neighbors(int **neighbor, float twonuc,
			       char **codon, char *coded_aa);

// Auxilliary
void Print_mut_fit(FILE *file_out, float *mut_par, float *P_mut, float lmax,
		   int GET_FREQ, float *num_aa, float *num_dna);
void Update_freq(float *fnew, float *fold, float *num_aa,
		 char **codon, char *coded_aa);
void Get_freq_dna(float *P_dna, float *P_cod);
void Vcopy(float *v1, float *v2, int n);
float Compare_distr(float *f1, float *f2, int n);
float Logfact(int n);
float Log_lik(float *num_aa, float *p_aa, int Naa);
void AA_distr(float *P_mut, float *P_cod, char *coded_aa);
int Combine_P_mut(float *P_mut_a, float scale, float *num_aa);
void Print_mutation(FILE *file_out, int np, float lik, float AIC,
		    float *x, int Npar);
void Normalize_rates(float **Q, float *P, int n);

// Main
void Codon_Transition(float **Q_cod, float *mut_par,
		      char **codon, char *coded_aa);
void Codon_Transition_1mut(float **Q_cod, float *mut_par,
			   char **codon, char *coded_aa);
float Fit_nuc_freq_1(float *mut_par, float *num_aa,
		     char **codon, char *coded_aa);
float Compute_P_mut(float *P_mut, float *P_cod, float **Q_cod,
		    float *mut_par, char **codon, char *coded_aa,
		    float *num_aa, float *num_dna);
float Limit_distr(float *Pi1, float **Q, int n, int *Stop_codon);
void Update_x_pos(float *x_new, float *x, float *dir, float *COEFF, int N);
void Update_x_pos_log_constraint(float *x_new, float *x, float *dir,
				 float *COEFF, int N);
int Step(double *F, float *x, float *x_new, float DELTA, int j,
	 float *P_mut, float *P_cod, float **Q_cod,
	 float *num_aa, float *num_dna);

void Get_mut_par(float *mut_par, float *P_mut, float *P_cod, float **Q_cod,
		 int GET_FREQ, float *num_aa, int L_aa, float *num_dna,
		 char *name_file, int iter)
{
  if(mut_par[6]<0.00001 || mut_par[6]>=1)mut_par[6]=0;
  if(mut_par[6]){NPAR=7;}else{NPAR=6;} // Consider Twonuc
  if(neighbor==NULL){
    twonuc=mut_par[6];
    neighbor=Allocate_mat2_i(64,64);
    n_neigh=Initialize_codon_neighbors(neighbor,twonuc,codon, coded_aa);
    for(int c=0; c<64; c++){
      if(coded_aa[c]=='*'){Stop_codon[c]=1;}
      else{Stop_codon[c]=0;}
    }
  }
  int a; float p=1./61;
  for(a=0; a<64; a++)if(Stop_codon[a]){P_cod[a]=0;}else{P_cod[a]=p;}
  printf("Computing limit mutation distribution\n");
  float lmax=
    Compute_P_mut(P_mut,P_cod,Q_cod,mut_par,codon,coded_aa,num_aa,num_dna);

  char name[100]; sprintf(name, "%s_mutation.dat", name_file);
  FILE *file_out; printf("\nWriting %s\n", name);
  if(iter==0){file_out=fopen(name, "w");}
  else{file_out=fopen(name, "a");}

  if(GET_FREQ){ // Fit mutation parameters
    // Initialize mutation parameters
    //lmax=Fit_nuc_freq_1(mut_par, num_aa, codon, coded_aa);
    printf("Initial Mut_par: ");
    for(a=0;a<NPAR; a++)printf(" %.3g", mut_par[a]);
    printf("\nLikelihood= %.2f Npar=%d \n\n", lmax, NPAR);
    Minimize_AIC(mut_par, NPAR, &lmax, P_mut, P_cod, Q_cod,
		 codon, coded_aa, num_aa, num_dna, L_aa, file_out);
    if(GET_FREQ==2){Combine_P_mut(P_mut, 400, num_aa);}
    if(GET_FREQ==3){
      float sum=0; int a;
      for(a=0; a<20; a++){P_mut[a]=num_aa[a]; sum+=num_aa[a];}
      for(a=0; a<20; a++)P_mut[a]/=sum;
    }
  }

  Print_mut_fit(file_out, mut_par, P_mut, lmax, GET_FREQ, num_aa, num_dna);
  fclose(file_out);
}


float Compute_P_mut(float *P_mut, float *P_cod, float **Q_cod,
		    float *mut_par, char **codon, char *coded_aa,
		    float *num_aa, float *num_dna)
{
  Normalize_distr(mut_par, 4);
  if(NPAR>6){ // twonuc!
    Codon_Transition(Q_cod, mut_par, codon, coded_aa);
  }else{
    Codon_Transition_1mut(Q_cod, mut_par, codon, coded_aa);
  }
  /*float p=1./61;
    for(int i=0; i<64; i++)if(Stop_codon[i]){P_cod[i]=0;}else{P_cod[i]=p;}*/
  float x=Limit_distr(P_cod, Q_cod, 64, Stop_codon);
  if(x<0.1)AA_distr(P_mut, P_cod, coded_aa);
  //for(a=0; a<20; a++)printf(" %.3f", P_mut[a]); printf("\n");
  float lik=0;
  if(num_aa)lik+=Log_lik(num_aa, P_mut, 20);
  if(num_dna){
    float P_dna[4]; Get_freq_dna(P_dna, P_cod);
    lik+=Log_lik(num_dna, P_mut, 4);
  }
  return(lik);
}

void Get_freq_dna(float *P_dna, float *P_cod){
  int i, j, c;
  for(i=0; i<4; i++)P_dna[i]=0;
  for(c=0; c<64; c++){
    char *cod=codon[c];
    for(j=0; j<3; j++)P_dna[Code_nuc(cod[j])]+=P_cod[c];
  }
  Normalize_distr(P_dna, 4);
}

float Fit_nuc_freq_1(float *mut_par, float *num_aa,
		     char **codon, char *coded_aa)
{
  // Get nucleotide frequencies from maximum likelihood
  int n, iter, NPAR=4;
  for(n=0; n<4; n++)mut_par[n]=0.25;
  printf("# Fast fit of nucleotide frequencies from amino acid sequence\n");
  printf(" observed freq: ");
  for(n=0; n<20; n++)printf(" %.0f", num_aa[n]);
  printf("\n");

  // Maximum likelihood
  float EPS=0.00000000001, l0=-10000;
  float *freq_new=malloc(NPAR*sizeof(float));
  float *freq_opt=malloc(NPAR*sizeof(float));
  float p_aa[20], P_cod[64];
  for(int c=0; c<64; c++)P_cod[c]=0;
  for(iter=0; iter<50; iter++){
    for(int c=0; c<64; c++){
      if(coded_aa[c]!='*')
	P_cod[c]=Weight_codon(codon[c], mut_par);
    }
    AA_distr(p_aa, P_cod, coded_aa);
    float l=Log_lik(num_aa, p_aa, 20);
    if(l>l0){l0=l; Vcopy(freq_opt, mut_par, NPAR);}
    Update_freq(freq_new, mut_par, num_aa, codon, coded_aa);
    float d=Compare_distr(freq_new, mut_par, 4);
    printf("%d %.2f %.4f\n", iter+1, l, sqrt(d));
    Vcopy(mut_par, freq_new, NPAR);
    if(d<EPS)break;
  }
  Vcopy(mut_par, freq_opt, NPAR);
  free(freq_new); free(freq_opt);
  // Printing
  printf("Fitted nucleotide frequencies:  ");
  for(n=0; n<4; n++)printf("%c %.3f ", Nuc_code(n), mut_par[n]);
  printf("\n");
  float sum=0; for(n=0; n<4; n++)sum+=mut_par[n];
  if(sum<=0){
    for(n=0; n<4; n++)mut_par[n]=0.25;
    printf("WARNING, convergence problem, setting equal nucl. frequencies\n");
  }else{
    for(n=0; n<4; n++)mut_par[n]/=sum;
  }
  return(l0);
}

float Log_lik(float *num, float *p, int N)
{
  int a, *INI; float *logfact;
  if(N==20){INI=&INI_FACT_aa; logfact=logfact_aa;}
  else{INI=&INI_FACT_dna; logfact=logfact_dna;}
  if(*INI==0){
    for(a=0; a<N; a++)logfact[a]=Logfact((int)num[a]);
    int L=0; for(a=0; a<N; a++)L+=num[a];
    logfact[N]=Logfact(L);
    *INI=1;
  }
  double sum=logfact[N];
  for(a=0; a<N; a++){
    if(p[a])sum+=num[a]*log(p[a])-logfact[a];
  }
  return(sum);
}

float Logfact(int n){
  if(INI_FACT==0){
    INI_FACT=1; log2pi_2=log(6.283)/2; int i;
    lfact[0]=0; for(i=1; i<Nf; i++)lfact[i]=lfact[i-1]+log(i); 
  }
  if(n<Nf)return(lfact[n]);
  return(log(n)*(n+0.5)-n+log2pi_2);
}

void Update_freq(float *fnew, float *fold, float *num_aa,
		 char **codon, char *coded_aa)
{
  int c, i, j, a;
  if(INI_MULT==0){
    for(c=0; c<64; c++){
      char *cod=codon[c], cod2[3]; int *g=Multi[c];
      for(i=0; i<4; i++)g[i]=0;
      for(j=0; j<3; j++)g[Code_nuc(cod[j])]++;
      if((strncmp(cod, "CG",2)==0)||(strncmp(cod+1, "CG",2)==0)){
	Mut_CpG[c]=c;
      }else if((strncmp(cod,  "CA",2)==0)||(strncmp(cod,  "TG",2)==0)){
	cod2[0]='C'; cod2[1]='G'; cod2[2]=cod[2];
	Mut_CpG[c]=Code_codon(cod2, codon);
      }else if((strncmp(cod+1,"CA",2)==0)||(strncmp(cod+1,"TG",2)==0)){
	cod2[0]=cod[0]; cod2[1]='C'; cod2[2]='G';
	Mut_CpG[c]=Code_codon(cod2, codon);
      }else{
	Mut_CpG[c]=-1;
      }
    }
    INI_MULT=1;
  }

  double norm[20], mult[20][4];
  for(a=0; a<20; a++){
    norm[a]=0; for(i=0; i<4; i++)mult[a][i]=0;
  }
  for(c=0; c<64; c++){
    if(coded_aa[c]=='*')continue;
    a=Code_AA(coded_aa[c]);
    float w=Weight_codon(codon[c], fold);
    norm[a]+=w;
    for(i=0; i<4; i++)mult[a][i]+=w*Multi[c][i];
  }

  // Compute new frequencies
  for(i=0; i<4; i++){
    double f=0;
    for(a=0; a<20; a++)f+=(num_aa[a]*mult[a][i]/norm[a]);
    fnew[i]=f;
  }
  Normalize_distr(fnew, 4);
}

float Compare_distr(float *f1, float *f2, int n){
  double sum=0; float d; int i;
  for(i=0; i<n; i++){d=f1[i]-f2[i]; sum+=d*d;}
  return(sum/n);
}

void Vcopy(float *v1, float *v2, int n){
  int i; for(i=0; i<n; i++)v1[i]=v2[i];
}

void Codon_Transition(float **Q_cod, float *mut_par,
		      char **codon, char *coded_aa)
{
  int c, d, j, cc[3], isCpG, isCpG_b;
  float mu2=mut_par[6], mu3=mu2*mu2;
  for(c=0; c<64; c++)for(d=0; d<64; d++)Q_cod[c][d]=0;
  for(c=0; c<64; c++){
    if(coded_aa[c]=='*')continue;  // Stop codon
    char *cod_c=codon[c];
    for(j=0; j<3; j++)cc[j]=Code_nuc(cod_c[j]);
    // CpG forwards and backwards
    if(strncmp(cod_c, "CG", 2)==0){isCpG=1;}
    else if(strncmp(cod_c+1, "CG", 2)==0){isCpG=2;}
    else{isCpG=0;}
    if(     strncmp(cod_c,  "TG", 2)==0){isCpG_b=0;}
    else if(strncmp(cod_c,  "CA", 2)==0){isCpG_b=1;}
    else if(strncmp(cod_c+1,"TG", 2)==0){isCpG_b=1;}
    else if(strncmp(cod_c+1,"CA", 2)==0){isCpG_b=2;}
    else{isCpG_b=-1;}

    for(d=0; d<c; d++){
      if(coded_aa[d]=='*')continue; // Stop codon
      float Rf=1, Rb=1; // forward and backward
      int n=0;
      for(j=0; j<3; j++){
	if(cod_c[j]!=codon[d][j]){
	  n++;
	  Rf*=mut_par[Code_nuc(codon[d][j])];
	  Rb*=mut_par[cc[j]];
	  if(Transition(codon[d][j])==cc[j]){
	    if(((isCpG==1)&&((j==0)||(j==1)))||
	       ((isCpG==2)&&((j==1)||(j==2)))){
	      Rf*=mut_par[4]; // k_CpG
	      Rb*=mut_par[5]; // tt_ratio
	    }else if(isCpG_b==j){
	      Rf*=mut_par[5]; // tt_ratio
	      Rb*=mut_par[4]; // k_CpG
	    }else{
	      Rf*=mut_par[5]; // tt_ratio
	      Rb*=mut_par[5]; // tt_ratio
	    }
	  }
	}
      } // end mutated codon
      if(n==2){Rf*=mu2; Rb*=mu2;}
      else if(n==3){Rf*=mu3; Rb*=mu3;}
      Q_cod[c][d]=Rf;
      Q_cod[d][c]=Rb;
    }
  }

  // Diagonal elements
  for(c=0; c<64; c++){
    double sum=0;
    for(d=0; d<64; d++)if(d!=c)sum+=Q_cod[c][d];
    Q_cod[c][c]=-sum;
  }

  return;
}

void Codon_Transition_1mut(float **Q_cod, float *mut_par,
			   char **codon, char *coded_aa)
{
  int c, d, isCpG=0;
  for(c=0; c<64; c++)for(d=0; d<64; d++)Q_cod[c][d]=0;
  for(c=0; c<64; c++){
    if(coded_aa[c]=='*')continue; // Stop codon
    int j, n; char cod2[3], *cod=codon[c];
    for(j=0; j<3; j++)cod2[j]=cod[j];
    double rate=0;
    for(j=0; j<3; j++){
      if(((strncmp(cod,  "CG", 2)==0)&&(j!=2))||
	 ((strncmp(cod+1,"CG", 2)==0)&&(j!=0))){isCpG=1;}
      else{isCpG=0;}
      int nj=Code_nuc(cod[j]), nt=Transition(cod[j]);
      for(n=0; n<4; n++){
	if(n==nj)continue;
	cod2[j]=Nuc_code(n);
	int c2=Code_codon(cod2, codon);
	if(coded_aa[c2]=='*')continue; // Stop codon
	float q=mut_par[n];
	// Transition
	if(n==nt){
	  if(isCpG){q*=mut_par[4];}else{q*=mut_par[5];}
	}
	Q_cod[c][c2]=q; rate+=q;
      } // End n
      cod2[j]=cod[j];
    } // End j
    Q_cod[c][c]=-rate;
  } // End c loop
  return;
}

int Initialize_codon_neighbors(int **neighbor, float twonuc,
			       char **codon, char *coded_aa)
{
  int c, nmax=0;
  for(c=0; c<64; c++){
    int nc=0;
    if(coded_aa[c]=='*')goto end_codon; // Stop codon
    if(twonuc){
      for(int d=0; d<64; d++){
	if((d!=c)&&(coded_aa[d]!='*')){neighbor[c][nc]=d; nc++;}
      }
    }else{
      int j1, n1; char cod2[3], *cod=codon[c];
      for(j1=0; j1<3; j1++)cod2[j1]=cod[j1];
      for(j1=0; j1<3; j1++){
	int nj=Code_nuc(cod[j1]);
	for(n1=0; n1<4; n1++){
	  if(n1==nj)continue;
	  cod2[j1]=Nuc_code(n1);
	  int c2=Code_codon(cod2, codon);
	  if(coded_aa[c2]=='*')continue; // Stop codon
	  neighbor[c][nc]=c2; nc++;
	} // End n1
	cod2[j1]=cod[j1];
      } // End j1
    } // end if
  end_codon:
    neighbor[c][nc]=-1;
    if(nc>nmax)nmax=nc;
  } // End c loop
  return(nmax);
}

void AA_distr(float *P_mut, float *P_cod, char *coded_aa)
{
  int a, c; 
  for(a=0; a<20; a++)P_mut[a]=0;
  for(c=0; c<64; c++){
    if(coded_aa[c]!='*')
      P_mut[Code_AA(coded_aa[c])]+=P_cod[c];
  }
}

float Limit_distr(float *Pi1, float **Q, int n, int *Forbidden)
{
  // Stationary distribution: PQ=0 
  int i, k, ik=0;
  float Pi2[n], *p1, *p2;

  // Normalize
  double sum=0, norm, d;
  for(i=0; i<n; i++)sum+=Pi1[i];
  if(sum==0){ // Not initialized!
    printf("WARNING, P is not initialized\n");
    for(i=0; i<n; i++){
      if(Forbidden[i]==0){Pi1[i]=1; sum++;}
      else{Pi1[i]=0;}
    }
  }
  for(i=0; i<n; i++)Pi1[i]/=sum;

  // Normalize Q 
  for(i=0; i<n; i++){
    double R=0; float *Qi=Q[i];
    for(k=0; k<n; k++){if(k!=i)R+=*Qi; Qi++;}
    Q[i][i]=-R;
  }

  int Kmax=40;
  for(k=0; k<Kmax; k++){
    if(ik==0){p1=Pi1; p2=Pi2; ik=1;}else{p1=Pi2; p2=Pi1; ik=0;}
    sum=0; norm=0;
    for(i=0; i<n; i++){
      if((Forbidden[i])||(Q[i][i]>=0)){
	p2[i]=0;
      }else{
	double p=0; int *l=neighbor[i];
	while(*l>=0){
	  p+=p1[*l]*Q[*l][i]; l++;
	}
	p/=(-Q[i][i]);
	sum+=p; p2[i]=p;
      }
    }
    if(sum<=0){
      printf("ERROR in Limit_distr, sum= %.3g\n", sum); goto end;
    }
    for(i=0; i<n; i++){
      p2[i]/=sum; d=p1[i]-p2[i]; norm+=d*d;
    }
    if(norm<0.000001)goto end;
  }
  printf("WARNING, limit mutation distribution not reached norm=%.2g\n",
	 norm);
  if(0){
    printf("P= "); for(i=0; i<n; i++){printf("%.2f ",Pi1[i]);} printf("\n");
    printf("Q_ii= "); for(i=0; i<n; i++){printf("%.2f ",Q[i][i]);}
    printf("\n");
  }

 end:
  //printf("Convergence of P_cod: %d %.1g\n", k, norm);
  return(norm);
}

void Minimize_AIC(float *x, int Npar, float *lmax,
		  float *P_mut, float *P_cod, float **Q_cod,
		  char **codon, char *coded_aa,
		  float *num_aa, float *num_dna, int L_aa,
		  FILE *file_out)
{
  x[4]=1; // tt_ratio
  x[5]=1; // CpG
  x[6]=0.; // twonuc
  fprintf(file_out, "Mutation parameters (DNA): A T G C tt_ratio kCpG");
  if(Npar==7)fprintf(file_out, " twonuc");
  fprintf(file_out, "\n");

  int N=4, i,a, No=19, np=0, npopt=0; //L_aa-1;
  float x_new[Npar]; for(i=0; i<Npar; i++)x_new[i]=x[i];
  float AIC_min=Compute_AIC(*lmax, np, No), AIC, l; 
  Print_mutation(file_out, np, *lmax/L_aa, AIC_min, x, Npar);

  int freezed[Npar];
  for(i=0; i<Npar; i++){freezed[i]=0;} freezed[3]=1;
  np=3;
  for(N=4; N<=Npar; N++){
    // N is the number of parameters, but free parameters are N-1
    for(i=0; i<N; i++)x_new[i]=x[i];
    /*for(i=0; i<4; i++)x_new[i]=0.25;
      x_new[4]=1.0; x_new[5]=1.0; x_new[6]=0;*/
    if(N==5){x_new[4]=2;}
    else if(N==6){x_new[5]=2;}
    else if(N==7){x_new[6]=0.05;}  // twonuc
    l=Steepest_ascent(x_new, freezed, N, P_mut, P_cod, Q_cod,
		      codon, coded_aa, num_aa, num_dna);
    AIC=Compute_AIC(l, np, No);
    Print_mutation(file_out, np, l/L_aa, AIC, x_new, Npar);
    if(AIC < AIC_min){
      *lmax=l; AIC_min=AIC; npopt=np; np++;
      for(i=0; i<N; i++)x[i]=x_new[i];
    }else{
      fprintf(file_out,
	      "WARNING, AIC increases, freezing the parameter\n"); 
      printf("WARNING, AIC increases, freezing the parameter\n");
      freezed[N-1]=1;
    }
  }
  fprintf(file_out, "Optimal mutation parameters (DNA):\n");
  Print_mutation(file_out, npopt, *lmax/L_aa, AIC_min, x, Npar);

  if(num_aa){
    int sum=0; for(a=0; a<20; a++)sum+=num_aa[a];
    float Pa[20]; for(a=0; a<20; a++)Pa[a]=num_aa[a]/(float)sum;
    l=Log_lik(num_aa, Pa, 20);
    AIC=Compute_AIC(l, 19, No);
    Print_mutation(file_out, 19, l/L_aa, AIC, Pa, 20);
  }

  *lmax=
    Compute_P_mut(P_mut, P_cod, Q_cod, x, codon, coded_aa, num_aa, num_dna);
  Normalize_rates(Q_cod, P_cod, 64);
}

float Steepest_ascent(float *x, int *freezed, int N,
		      float *P_mut, float *P_cod, float **Q_cod,
		      char **codon, char *coded_aa,
		      float *num_aa, float *num_dna)
{
  int KMAX=100;      // Iterations
  float DELTA=0.012;  // Increment to compute the gradient
  float EPS=0.00000000001;
  float x_new[NPAR]; double F;
  int j, k, np=0;

  printf("Maximize ML by steepest ascent, ");
  //printf("#\n starting: ");
  //for(j=0; j<NPAR; j++)printf(" %.3g", x[j]);
  for(j=0; j<N; j++)if(freezed[j]==0)np++;
  printf(" Free par.= %d\n", np);
  printf("#iter\tlik");
  for(j=0; j<4; j++)printf("\tf%c", Nuc_code(j));
  printf("\tCpG"); printf("\ttt_ratio");
  if(NPAR==7)printf("\tTwonuc");
  printf("\n");

  /* Set limits
  float x_min[NPAR], x_max[NPAR];
  for(j=0; j<4; j++){x_min[j]=0.05; x_max[j]=1-3*0.05;}
  x_min[4]=0.1; x_max[4]=100; x_min[5]=0.1; x_max[5]=100;
  x_min[6]=0; x_max[6]=0.5;*/

  for(k=0; k<KMAX; k++){
    F=Compute_P_mut(P_mut, P_cod, Q_cod, x, codon, coded_aa, num_aa, num_dna);
    float F_old=F;
    if(VBS){
      printf("%d\t%.4g", k, F);
      for(j=0; j<NPAR; j++)printf("\t%.3g", x[j]);
      printf("\n");
      for(j=0; j<NPAR; j++)x_new[j]=x[j];
    }
    if(isnan(F))return(F);

    // Gradient computation
    for(j=0; j<N; j++){
      if(j==3)continue; // 4th nucleotide
      if(freezed[j])continue;
      double xopt=x[j], x3; if(j<3)x3=x[3];
      if(Step(&F, x, x_new, DELTA, j, P_mut,P_cod,Q_cod,num_aa,num_dna)){
	xopt=x_new[j]; if(j<3)x3=x_new[3]; goto update;
      }
      if(Step(&F, x, x_new,-DELTA, j, P_mut,P_cod,Q_cod,num_aa,num_dna)){
	xopt=x_new[j]; if(j<3)x3=x_new[3];
      }
      // Update
    update:
      x[j]=xopt; x_new[j]=xopt; if(j<3){x[3]=x3; x_new[3]=x3;}
    }
    if((k)&&((F-F_old)<EPS*fabs(F_old)))break;
  }
  F=Compute_P_mut(P_mut, P_cod, Q_cod, x, codon, coded_aa, num_aa, num_dna);
  printf("%d\t%.4g", k, F);
  for(j=0; j<NPAR; j++)printf("\t%.3g", x[j]);
  printf("\n");
  printf("#Pmut: ");
  for(j=0; j<20; j++){printf(" %.3f", P_mut[j]);} printf("\n");
  float sum=0; for(j=0; j<20; j++){sum+=num_aa[j];}
  printf("#Pobs: ");
  for(j=0; j<20; j++){printf(" %.3f", num_aa[j]/sum);} printf("\n#\n");

  return(F);
}

int Step(double *F, float *x, float *x_new, float DELTA, int j,
	 float *P_mut, float *P_cod, float **Q_cod,
	 float *num_aa, float *num_dna)
{
  float d=DELTA*x[j];
  if(j<3){
    while(1){
      x_new[j]=x[j]+d; x_new[3]=x[3]-d;
      if((x_new[3]>0)&&(x_new[j]<1)&&(x_new[j]>0))break;
      d/=1.5;
    }
  }else{  
    while(1){
      x_new[j]=x[j]+d; if(x_new[j]>0)break; d/=1.5;
    }
  }
  double Fj=
    Compute_P_mut(P_mut,P_cod,Q_cod,x_new,codon,coded_aa,num_aa,num_dna);
  if(Fj > *F){*F=Fj; return(1);}
  return(0);
}
	
void Steepest_ascent_old(float *x, int N,
			 float *P_mut, float *P_cod, float **Q_cod,
			 char **codon, char *coded_aa,
			 float *num_aa, float *num_dna)
{
  int KMAX=40;      // Iterations
  float ALPHA=0.0;   // Integration step
  float DELTA=0.02;  // Increment to compute the gradient
  float EPS=0.00001; // Convergence
  float *x_new=malloc(N*sizeof(float));
  float *dF=malloc(N*sizeof(float)), F, Fj;
  float *x_opt=malloc(N*sizeof(float)), Fopt;
  int k, j; double norm=0;
  printf("#iter\tscore\t|grad|\talpha");
  for(j=0; j<4; j++)printf("\tf%c", Nuc_code(j));
  printf("\tCpG"); printf("\ttt_ratio");
  if(NPAR==7)printf("\tTwonuc");
  printf("\n");

  for(k=0; k<KMAX; k++){
    F=Compute_P_mut(P_mut, P_cod, Q_cod, x, codon, coded_aa, num_aa, num_dna);
    if((k==0)||(F>Fopt)){
      Fopt=F; for(j=0; j<N; j++)x_opt[j]=x[j];
    }
    printf("%2d\t%.4g\t%.2g\t%.3g", k, F, norm, ALPHA);
    for(j=0; j<N; j++)printf("\t%.3g", x[j]);
    printf("\n");

    // Gradient computation
    for(j=0; j<N; j++)x_new[j]=x[j];
    norm=0;
    for(j=0; j<N; j++){
      if(j==3)continue;
      x_new[j]+=DELTA; if(j<3)x_new[3]-=DELTA;
      Fj=Compute_P_mut(P_mut, P_cod, Q_cod, x_new, codon, coded_aa,
		       num_aa,num_dna);
      dF[j]=(Fj-F); norm+=dF[j]*dF[j];
      x_new[j]=x[j]; if(j<3)x_new[3]=x[3];
    }
    if(norm<EPS)return;
    norm=sqrt(norm); for(j=0; j<N; j++)dF[j]/=norm;
    ALPHA=Quadratic_line_search(dF, x, F, N, P_mut, P_cod, Q_cod,
    				codon, coded_aa, num_aa, num_dna);
    Update_x_pos_log_constraint(x, x, dF, &ALPHA, N);
    //if(ALPHA<0.00001){k==KMAX; break;}
  }

  if(k==KMAX)
    printf("WARNING, steepest ascent optimization did not converge\n");
  printf("%d\t%.4g\t%.2g\t%.3g", k, Fopt, norm, ALPHA);
  for(j=0; j<N; j++)printf("\t%.3g", x_opt[j]);
  printf("\n");
  for(j=0; j<N; j++)x[j]=x_opt[j];
  free(x_new); free(x_opt); free(dF);
}

void Update_x_pos(float *x_new, float *x, float *dir, float *COEFF, int N){
  int j; float *y=malloc(N*sizeof(float));
  while(1){
    for(j=0; j<N; j++){
      y[j]=x[j]+(*COEFF)*dir[j]; if(y[j]<0)break;
    }
    if(j==N)break;
    (*COEFF)/=1.5;
    printf("WARNING,negative x\n");
  }
  for(j=0; j<N; j++)x_new[j]=y[j];
  free(y);
}

void Update_x_pos_log_constraint(float *x_new, float *x, float *dir,
				 float *COEFF, int N)
{
  int j; 
  while(1){
    float sum=0;
    for(j=0; j<N; j++){
      if(j==3)continue;
      x_new[j]=x[j]*exp((*COEFF)*dir[j]*x[j]);
      if(j<3)sum+=x_new[j];
    }
    x_new[3]=1.-sum;
    if((sum>0)&&(x_new[3]>0))break;
    (*COEFF)/=1.5;
  }
}

float Quadratic_line_search(float *dF, float *x, float F0, int N,
			    float *P_mut, float *P_cod, float **Q_cod,
			    char **codon, char *coded_aa,
			    float *num_aa, float *num_dna)
{
  float h0=0, h1=0.01, h2=0.1, F1, F2=F0, F3;
  float *x_new=malloc(N*sizeof(float));
  int i2=0;

  Update_x_pos_log_constraint(x_new, x, dF, &h1, N);
  F1=Compute_P_mut(P_mut, P_cod, Q_cod, x_new, codon,coded_aa,num_aa,num_dna);
  while(F1<F0){
    printf("WARNING, descending in the gradient direction, %.5g %.5g\n",F0,F1);
    i2++;
    if(i2==2){
      printf("ERROR in line search\n"); return(-h1);
    }
    if(F1>F2){h2=h1; F2=F1;} h1=h1/2;
    Update_x_pos_log_constraint(x_new, x, dF, &h1, N);
    F1=Compute_P_mut(P_mut, P_cod, Q_cod, x_new, codon,coded_aa,num_aa,num_dna);
  }
  if(i2==0){
    Update_x_pos_log_constraint(x_new, x, dF, &h2, N);
    F2=Compute_P_mut(P_mut, P_cod, Q_cod, x_new, codon,coded_aa,num_aa,num_dna);
    while(F2>F1){
      i2++;
      if(i2==3){
	printf("No correct bounds in line search\n"); break;
      }
      h1=h2; F1=F2; h2=h2*2.5;
      Update_x_pos_log_constraint(x_new, x, dF, &h2, N);
      F2=Compute_P_mut(P_mut, P_cod, Q_cod, x_new,
		       codon,coded_aa, num_aa, num_dna);
    }
  }
  float h=Find_max_quad(h0,h1,h2,F0,F1,F2,0,h2*4);
  /*printf("Find_max: h=%.2g  h123= %.4g %.4g %.4g  F123= %.4g %.4g %.4g\n",
    h,h0,h1,h2,F0,F1,F2);*/
  Update_x_pos_log_constraint(x_new, x, dF, &h, N);
  F3=Compute_P_mut(P_mut, P_cod, Q_cod, x_new, codon, coded_aa, num_aa,num_dna);
  free(x_new);
  if(F3>F1){return(h);}
  else{return(h1);}
}

int Combine_P_mut(float *P_mut_a, float scale, float *num_aa)
{
  // num_aa = observed counts
  int i; double sum=0;
  for(i=0; i<20; i++){
    P_mut_a[i]=scale*P_mut_a[i]+num_aa[i];
    sum+=P_mut_a[i];
  }
  for(i=0; i<20; i++)P_mut_a[i]/=sum;
  return(0);
}

float Compute_AIC(float lik, int np, int no){
  float AIC=np-lik;
  //AIC+=(np+1)*(np+2)/(no-np-2);
  return(2*AIC);
}

void Print_mutation(FILE *file_out, int np, float lik, float AIC,
		    float *x, int Npar)
{
  fprintf(file_out, "%d mutation parameters, l= %.4f AIC= %.1f   Mut_par: ",
	  np, lik, AIC);
  for(int i=0; i<Npar; i++)fprintf(file_out, " %.3f", x[i]);
  fprintf(file_out, "\n");
}
 
void Normalize_rates(float **Q, float *P, int n){
  double R=0; int i, j;
  for(i=0; i<n; i++){R+=P[i]*Q[i][i];} R=-R;
  for(i=0; i<n; i++)for(j=0; j<n; j++)Q[i][j]/=R;
}

void Print_mut_fit(FILE *file_out, float *mut_par, float *P_mut, float lmax,
		   int GET_FREQ, float *num_aa, float *num_dna)
{

  char buf[50000], tmp[80]; int a;
  sprintf(buf,"Likelihood computed from ");
  if(num_aa)strcat(buf,"amino acid sequence ");
  if(num_dna)strcat(buf,"and nucleotide sequence");
  strcat(buf, "\nObtaining mutation distribution from ");
  if(GET_FREQ==0){
    strcat(buf, "mutation model, input parameters\n");
  }else if(GET_FREQ==3){
    strcat(buf, "amino acid frequencies (no fit)\n");
  }else if(GET_FREQ==2){
    strcat(buf, "mutation model, ML fit plus amino acid sequence\n");
  }else{
    strcat(buf,"mutation model, ML fit\n");
  }
  strcat(buf, "#AA:   ");
  for(a=0; a<20; a++){
    sprintf(tmp,"\t%c", Amin_code(a)); strcat(buf, tmp);
  }
  strcat(buf,"\n#Pobs: ");
  double sum_aa=0; for(a=0; a<20; a++)sum_aa+=num_aa[a];
  for(a=0; a<20; a++){
    sprintf(tmp,"\t%.3f", num_aa[a]/sum_aa); strcat(buf, tmp);
  }
  strcat(buf,"\n#Pmut: ");
  for(a=0; a<20; a++){
    sprintf(tmp,"\t%.3f", P_mut[a]); strcat(buf, tmp);
  }
  strcat(buf,"\n#DP:   ");
  for(a=0; a<20; a++){
    sprintf(tmp,"\t%.3f", num_aa[a]/sum_aa-P_mut[a]); strcat(buf, tmp);
  }
  strcat(buf, "\n");
  if(GET_FREQ<3){
    strcat(buf, "Mut_par: ");
    for(a=0;a<NPAR; a++){
      sprintf(tmp, "\t%.3g", mut_par[a]); strcat(buf, tmp);
    }
    sprintf(tmp, "\nLikelihood= %.2f\n", lmax); strcat(buf, tmp);
  }
  printf("%s", buf);
  fprintf(file_out, "%s", buf);
  printf("###\n### End of computation of background distributions\n###\n\n");
}
