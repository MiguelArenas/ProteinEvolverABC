#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "Print_meanfield.h"
#include "codes.h"
#include "output.h"
#include "gen_code.h"
#include "allocate.h"
#include "wag.h"
float RATE=0.1; // Mutation rate at nucleotide level, needed to compute
                // double nucleotide mutations (not implemented yet)
char AA_string[]="ACDEFGHIKLMNPQRSTVWY";
//#define AMIN_CODE "AEQDNLGKSVRTPIMFYCWHX"
char AA_PAML[]="ARNDCQEGHILKMFPSTWYV";

float hscale[]={0.1366,-0.0484,0.0325,-0.1233,-0.0345,0.4251,-0.0464,-0.0101,-0.0432,0.4083,0.0363,0.0589,0.0019,0.4172,0.1747,0.4076,0.3167,0.2745,0.2362,0.0549};

static void Change_AA_order(int *iaa, char *AA);
static void Compute_Q_sel(float **Q_sel, float *P_MF);
static void Compute_rates_mut(float **Rates_mut, float *P_mut_a,
			      float *freq_nuc, float tt_ratio,
			      char **codon, char *coded_aa,
			      float TWONUC);
static void Update_Q(float *Q_mut, float *Q, char *cod2,
		     int j, int n, int nt,
		     float *rate_nuc, float tt_ratio,
		     char **codon, char *coded_aa);
float Entropy(float *P, int n);
double Mean_hydro(float *P, float *hydro, int n);
void Sum_matrix(float *ncont, int **Cnat, int L);
float Corr_coeff(float *xx, float *yy, int n);

/******************* Codes ****************************/
int Print_profiles(float **P_MF_ia, float *P_mut_a, double DG_ave,
		   float Lambda, short *aa_seq, int L, char *name_file)
{
  FILE *file_out=Output_file(name_file, "AA_profiles", "txt");
  int *iaa=malloc(20*sizeof(int)), i, a;
  Change_AA_order(iaa, AA_string);

  fprintf(file_out, "# Mean field distribution, Lambda= %.3f\n", Lambda);
  fprintf(file_out, "# ave(DeltaG/T)= %.2f corresponding to Npop= %.3g\n",
	  DG_ave, Lambda*exp(-DG_ave));
  fprintf(file_out, "#pos AA");
  for(a=0; a<20; a++)fprintf(file_out, "    %c", AA_string[a]);
  fprintf(file_out, " Entropy\n");

  // Print mutational distribution
  fprintf(file_out, "#P_MUT");
  double S=0;
  for(a=0; a<20; a++){
    int ia=iaa[a]; float p=P_mut_a[ia]; S+=p*log(p);
    fprintf(file_out, " %.3f", p);
  }
  fprintf(file_out, "   %.2f\n", -S);

  // Print mean-field distribution
  for(i=0; i<L; i++){
    S=0;
    fprintf(file_out, "%3d %c ", i+1, Amin_code(aa_seq[i]));
    for(a=0; a<20; a++){
      int ia=iaa[a]; float p=P_MF_ia[i][ia]; S+=p*log(p);
      fprintf(file_out, " %.3f", p);
    }
    fprintf(file_out, "   %.2f\n", -S);
  }
  fclose(file_out);
  free(iaa);
  return(0);
}

int Print_profile_evo(char *name, double **P_ia, short *aa_seq, int L,
		      double DG_ave, long it_sum)
{
  FILE *file_out=Output_file(name, "AA_profiles_evo", "txt");
  int *iaa=malloc(20*sizeof(int)), i, a;
  Change_AA_order(iaa, AA_string);

  fprintf(file_out, "# Evolutionary distribution, %ld substitutions\n",
	  it_sum);
  fprintf(file_out, "# ave(DeltaG/T)= %.2f\n", DG_ave/it_sum);
  fprintf(file_out, "#pos AA");
  for(a=0; a<20; a++)fprintf(file_out, "    %c", AA_string[a]);
  fprintf(file_out, " Entropy\n");

  // Print distribution
  for(i=0; i<L; i++){
    fprintf(file_out, "%3d %c ", i+1, Amin_code(aa_seq[i]));
    double S=0, norm=0;
    for(a=0; a<20; a++)norm+=P_ia[i][a];
    for(a=0; a<20; a++){
      int ia=iaa[a]; float p=P_ia[i][ia]/norm; if(p)S+=p*log(p);
      fprintf(file_out, " %.3f", p);
    }
    fprintf(file_out, "   %.2f\n", -S);
  }
  fclose(file_out);
  free(iaa);
  return(0);
}

int Print_subst_rate(float **P_MF_ia, float *P_mut_a,
		     float *freq_nuc, float tt_ratio,
		     char **codon, char *coded_aa, short *aa_seq,
		     int **Cnat, int L, char *name_file,
		     int PRINT_MAT, float TWONUC)
{
  int *iaa=malloc(20*sizeof(int)), i, a, b;
  Change_AA_order(iaa, AA_string);
  float *P_sel=malloc(20*sizeof(float));
  float **Q_sel=Allocate_mat2_f(20, 20);
  float **Rates_mut=Allocate_mat2_f(20, 20);
  float **Rates_MF=Allocate_mat2_f(20, 20);
  float *Rates=malloc(20*sizeof(float));
  float *ncont=malloc(L*sizeof(float));
  Sum_matrix(ncont, Cnat, L);
  float *entropy=malloc(L*sizeof(float));
  float *Rsum=malloc(L*sizeof(float));
  float *hydro_i=malloc(L*sizeof(float));
  float *hydro_ave=malloc(L*sizeof(float));

  FILE *file_mat;
  if(PRINT_MAT){
    char what[100];
    sprintf(what, "rate_matrix_tt%.2f_2nuc%.2f", tt_ratio, TWONUC);
    file_mat=Output_file(name_file, what, "txt");
    fprintf(file_mat, "AA ");
    for(a=0; a<20; a++)fprintf(file_mat, "  %c  ", AA_string[a]);
    fprintf(file_mat, "\n");
  }
  FILE *file_out=Output_file(name_file, "rate_profile", "dat");
  fprintf(file_out, "#SITE AA ncont rate entropy hydro ave_hydro\n"); 
  Compute_rates_mut(Rates_mut, P_mut_a,freq_nuc,tt_ratio,codon,coded_aa,TWONUC);
  for(i=0; i<L; i++){
    float *p=P_MF_ia[i];
    for(a=0; a<20; a++)P_sel[a]=p[a]/P_mut_a[a];
    Compute_Q_sel(Q_sel, P_sel);
    Rsum[i]=0;
    for(a=0; a<20; a++){
      float *Qm_a=Rates_mut[a];
      float *Qs_a=Q_sel[a];
      double sum=0;
      for(b=0; b<20; b++){
	if(b!=a){
	  float R=Qm_a[b]*Qs_a[b]; sum+=R;
	  if(PRINT_MAT)Rates_MF[a][b]=R;
	}
      }
      Rates_MF[a][a]=-sum;
      Rsum[i] += p[a]*sum;
    }
    if(PRINT_MAT){
      fprintf(file_mat, "SITE %3d %c\n", i+1, Amin_code(aa_seq[i]));
      for(a=0; a<20; a++){
	float *R=Rates_MF[iaa[a]];
	fprintf(file_mat, "%c", AA_string[a]);
	for(b=0; b<20; b++)fprintf(file_mat, "\t%.4f", R[iaa[b]]); //Rsum[i]
	fprintf(file_mat, "\n");
      }
    }
    entropy[i]=Entropy(P_MF_ia[i], 20);
    hydro_i[i]=hscale[aa_seq[i]];
    hydro_ave[i]=Mean_hydro(P_MF_ia[i], hscale, 20);
    fprintf(file_out, "%3d %c  %2.0f  %5.3f  %5.2f %5.2f %5.2f\n",
	    i+1, Amin_code(aa_seq[i]), ncont[i], entropy[i], Rsum[i],
	    hydro_i[i], hydro_ave[i]);
  }
  if(PRINT_MAT)fclose(file_mat);
  float r=Corr_coeff(ncont, hydro_i, L);
  fprintf(file_out, "# Corr(ncont, hydro)=     %6.3f %d sites\n", r, L);
  r=Corr_coeff(ncont, hydro_ave, L);
  fprintf(file_out, "# Corr(ncont, ave_hydro)= %6.3f %d sites\n", r, L);
  r=Corr_coeff(ncont, entropy, L);
  fprintf(file_out, "# Corr(ncont, entropy)=   %6.3f %d sites\n", r, L);
  r=Corr_coeff(ncont, Rsum, L);
  fprintf(file_out, "# Corr(ncont, rate)=      %6.3f %d sites\n", r, L);
  r=Corr_coeff(entropy, Rsum, L);
  fprintf(file_out, "# Corr(rate, entropy)=    %6.3f %d sites\n", r, L);

  fclose(file_out);
  free(iaa);
  free(ncont);
  free(entropy);
  free(Rsum);
  free(hydro_i);
  free(hydro_ave);
  Empty_matrix_f(Rates_mut, 20);
  Empty_matrix_f(Rates_MF, 20);
  Empty_matrix_f(Q_sel, 20);
  free(P_sel);
  return(0);
}

int Print_exchange(float **P_MF_ia, float *P_mut_a,
		   float *freq_nuc, float tt_ratio,
		   char **codon, char *coded_aa, short *aa_seq,
		   int **Cnat, int L, char *name_file,
		   int PRINT_MAT, float TWONUC,
		   int FORMAT, int EXCHANGE)
{
  int *iaa=malloc(20*sizeof(int)), i, a, b;
  if(FORMAT==0){Change_AA_order(iaa, AA_string);}
  else{Change_AA_order(iaa, AA_PAML);}
  float *P_sel=malloc(20*sizeof(float));
  float **Q_sel=Allocate_mat2_f(20, 20);
  float **Rates_mut=Allocate_mat2_f(20, 20);
  float **Rates_MF=Allocate_mat2_f(20, 20);
  float *Rates=malloc(20*sizeof(float));
  float *ncont=malloc(L*sizeof(float));
  Sum_matrix(ncont, Cnat, L);
  float *entropy=malloc(L*sizeof(float));
  float *Rsum=malloc(L*sizeof(float));
  float *hydro_i=malloc(L*sizeof(float));
  float *hydro_ave=malloc(L*sizeof(float));
  double **Fsum=Allocate_mat2_d(20, 20);
  double *Psum=malloc(20*sizeof(double));
  for(a=0; a<20; a++)Psum[a]=0;

  FILE *file_mat; char name[100];
  if(PRINT_MAT){
    if(EXCHANGE==0){
      sprintf(name, "exchangeability_tt%.2f_2nuc%.2f", tt_ratio, TWONUC);
      file_mat=Output_file(name_file, name, "txt");
      fprintf(file_mat,
	      "# Exchangeability: genetic code, tt_ratio=%.2f TWONUC=%.3f\n",
	      tt_ratio, TWONUC);
    }else{
      sprintf(name, "exchangeability_WAG");
      file_mat=Output_file(name_file, name, "txt");
      fprintf(file_mat,"# Exchangeability: WAG matrix\n");
    }
    if(FORMAT==0){
      fprintf(file_mat, "AA ");
      for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_string[a]);
    }else{
      for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
    }
    fprintf(file_mat, "\n");
  }

  float exchange[20][20];
  if(EXCHANGE==1){
    // Use WAG empirical exchangeability matrix
    int *iwag=malloc(20*sizeof(int));
    Change_AA_order(iwag, AA_WAG);
    for(a=0; a<20; a++){
      for(b=0; b<20; b++){
	if(b<a){exchange[iwag[a]][iwag[b]]=WAG[a][b];}
	else{exchange[iwag[a]][iwag[b]]=WAG[b][a];}
      }
    }
    free(iwag);
  }
  
  FILE *file_out=Output_file(name_file, "rate_profile", "dat");
  fprintf(file_out, "#SITE AA ncont rate entropy hydro ave_hydro\n");
  Compute_rates_mut(Rates_mut,P_mut_a,freq_nuc,tt_ratio,codon,coded_aa,TWONUC);
  float norm=1;
  for(i=0; i<L; i++){
    float *p=P_MF_ia[i];
    for(a=0; a<20; a++)P_sel[a]=p[a]/P_mut_a[a];
    Compute_Q_sel(Q_sel, P_sel);
    Rsum[i]=0; 
    for(a=0; a<20; a++){
      float *Qm_a=Rates_mut[a];
      float *Qs_a=Q_sel[a], R;
      double sum=0;
      Psum[a]+=p[a];
      for(b=0; b<20; b++){
	if(b!=a){
	  if(EXCHANGE==0){R=Qm_a[b]*Qs_a[b];}
	  else{R=exchange[a][b]*p[b];}
	  sum+=R;
	  if(PRINT_MAT){
	    Rates_MF[a][b]=R;
	    if(b>a)Fsum[a][b]+=p[a]*R;
	  }
	}
      }
      Rates_MF[a][a]=-sum;
      Rsum[i] += p[a]*sum;
    }
    if((i==0)&&(EXCHANGE==0))norm=Rsum[0];
    if(PRINT_MAT){
      // Normalize by the total rate of the first site Rsum[0]
      fprintf(file_mat, "SITE %3d %c\n", i+1, Amin_code(aa_seq[i]));
      if(FORMAT==0){
	for(a=0; a<20; a++){
	  float *R=Rates_MF[iaa[a]];
	  fprintf(file_mat, "%c", AA_string[a]);
	  for(b=0; b<20; b++){
	    fprintf(file_mat, "\t%.4f", R[iaa[b]]/(p[iaa[b]]*norm));
	  }
	  fprintf(file_mat, "\n");
	}
      }else{
	for(a=1; a<20; a++){
	  float *R=Rates_MF[iaa[a]];
	  for(b=0; b<a; b++)
	    fprintf(file_mat, "%.4f\t", R[iaa[b]]/(p[iaa[b]]*norm));
	  fprintf(file_mat, "\n");
	}
	fprintf(file_mat, "\n");
	for(a=0; a<20; a++)fprintf(file_mat, "%.4f\t", p[iaa[a]]);
	fprintf(file_mat, "\n");
      }
    }
    entropy[i]=Entropy(P_MF_ia[i], 20);
    hydro_i[i]=hscale[aa_seq[i]];
    hydro_ave[i]=Mean_hydro(P_MF_ia[i], hscale, 20);
    fprintf(file_out, "%3d %c  %2.0f  %5.3f  %5.2f %5.2f %5.2f\n",
	    i+1, Amin_code(aa_seq[i]), ncont[i], entropy[i], Rsum[i],
	    hydro_i[i], hydro_ave[i]);
  }
  if(PRINT_MAT && FORMAT){
    fprintf(file_mat, "\n// this is the end of the file. The rest are notes\n");
    for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
    fprintf(file_mat, "\n");
  }

  // Print site averaged matrix
  if(PRINT_MAT ){
    fclose(file_mat);
    if(EXCHANGE==0){
      sprintf(name, "exchangeability_tt%.2f_2nuc_allsites%.2f",tt_ratio,TWONUC);
      file_mat=Output_file(name_file, name, "txt");
      fprintf(file_mat,
	      "# Exchangeability: genetic code, tt_ratio=%.2f TWONUC=%.3f\n",
	      tt_ratio, TWONUC);
    }else{
      sprintf(name, "exchangeability_WAG_allsites");
      file_mat=Output_file(name_file, name, "txt");
      fprintf(file_mat,"# Exchangeability: WAG matrix\n");
    }
    double sum=0; for(a=0; a<20; a++)sum+=Psum[a];
    for(a=0; a<20; a++){
      Psum[a]/=sum;
      for(b=0; b<a; b++)Fsum[a][b]=Fsum[b][a];
    }
    if(FORMAT==0){
      fprintf(file_mat, "AA ");
      for(a=0; a<20; a++)fprintf(file_mat, "\t%c", AA_string[a]);
      fprintf(file_mat, "\n");
      for(a=0; a<20; a++){
	double *F=Fsum[iaa[a]];
	double Pa=Psum[iaa[a]];
	fprintf(file_mat, "%c", AA_string[a]);
	for(b=0; b<20; b++){
	  fprintf(file_mat, "\t%.4f", F[iaa[b]]/(L*Rsum[0]*Pa*Psum[iaa[b]]));
	}
	fprintf(file_mat, "\n");
      }
    }else{
      for(a=1; a<20; a++){
	double *F=Fsum[iaa[a]];
	double Pa=Psum[iaa[a]];
	for(b=0; b<a; b++)
	  fprintf(file_mat, "%.4f\t", F[iaa[b]]/(L*norm*Pa*Psum[iaa[b]]));
	fprintf(file_mat, "\n");
      }
      fprintf(file_mat, "\n");
      for(a=0; a<20; a++)fprintf(file_mat, "%.4f\t", Psum[iaa[a]]);
      fprintf(file_mat, "\n\n");
      fprintf(file_mat,"// this is the end of the file. The rest are notes\n");
      for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
      fprintf(file_mat, "\nP_mut:\n");
      for(a=0; a<20; a++)fprintf(file_mat, "%.4f\t", P_mut_a[iaa[a]]);
      fprintf(file_mat, "\n");
    }
    fclose(file_mat);
  }


  float r=Corr_coeff(ncont, hydro_i, L);
  fprintf(file_out, "# Corr(ncont, hydro)=     %6.3f %d sites\n", r, L);
  r=Corr_coeff(ncont, hydro_ave, L);
  fprintf(file_out, "# Corr(ncont, ave_hydro)= %6.3f %d sites\n", r, L);
  r=Corr_coeff(ncont, entropy, L);
  fprintf(file_out, "# Corr(ncont, entropy)=   %6.3f %d sites\n", r, L);
  r=Corr_coeff(ncont, Rsum, L);
  fprintf(file_out, "# Corr(ncont, rate)=      %6.3f %d sites\n", r, L);
  r=Corr_coeff(entropy, Rsum, L);
  fprintf(file_out, "# Corr(rate, entropy)=    %6.3f %d sites\n", r, L);

  fclose(file_out);
  free(iaa);
  free(ncont);
  free(entropy);
  free(Rsum);
  free(hydro_i);
  free(hydro_ave);
  Empty_matrix_f(Rates_MF, 20);
  Empty_matrix_f(Rates_mut, 20);
  Empty_matrix_f(Q_sel, 20);
  Empty_matrix_d(Fsum, 20);
  free(P_sel);
  free(Psum);
  return(0);
}

void Change_AA_order(int *iaa, char *AA){
  int a; for(a=0; a<20; a++)iaa[a]=Code_AA(AA[a]);
}

void Compute_rates_mut(float **Q_mut, float *P_mut_a,
		       float *freq_nuc, float tt_ratio,
		       char **codon, char *coded_aa, float TWONUC)
{
  // Q_mut(a,b)=sum_{c\in C(a) c'\in C(b) d(c,c')=1} P(c|a)Q(c,c')
  int a, b, c;
  for(a=0; a<20; a++)for(b=0; b<20; b++)Q_mut[a][b]=0;
  double *norm=malloc(20*sizeof(double));
  for(a=0; a<20; a++)norm[a]=0;
  float rate_nuc[4];
  for(a=0; a<4; a++)rate_nuc[a]=freq_nuc[a];
  if(TWONUC>0)for(a=0; a<4; a++)rate_nuc[a]*=TWONUC;

  for(c=0; c<64; c++){
    if(coded_aa[c]=='*')continue; // Stop codon
    float w=Weight_codon(codon[c], freq_nuc);
    a=Code_AA(coded_aa[c]); norm[a]+=w;
    float *Qa=Q_mut[a];
    int j, n; char cod2[3], *cod=codon[c];
    for(j=0; j<3; j++)cod2[j]=cod[j];
    for(j=0; j<3; j++){
      int nt=Transition(cod[j]);
      for(n=0; n<4; n++){
	if(Nuc_code(n)==cod[j])continue;
	float Q=w;
	Update_Q(Qa,&Q,cod2, j,n,nt,rate_nuc,tt_ratio,codon,coded_aa);
	if(TWONUC<=0)continue;
	int j2, n2, j3, n3;
	for(j2=j+1; j2<3; j2++){
	  int nt2=Transition(cod[j2]);
	  for(n2=0; n2<4; n2++){
	    if(Nuc_code(n2)==cod[j2])continue;
	    float Q2=Q;
	    Update_Q(Qa,&Q2,cod2, j2,n2,nt2,rate_nuc,tt_ratio,codon,coded_aa);
	    for(j3=j2+1; j3<3; j3++){
	      int nt3=Transition(cod[j3]);
	      for(n3=0; n3<4; n3++){
		if(Nuc_code(n3)==cod[j3])continue;
		float Q3=Q2;
		Update_Q(Qa,&Q3,cod2,
			 j3,n3,nt3,rate_nuc,tt_ratio,codon,coded_aa);
		cod2[j3]=cod[j3];
	      }
	    }
	    cod2[j2]=cod[j2];
	  }
	}
	// End TWONUC
	cod2[j]=cod[j];
      }
    }
  }
  
  // Impose detailed balance
  for(a=0; a<20; a++){
    for(b=0; b<a; b++){
      Q_mut[a][b]/=norm[a];
      Q_mut[b][a]=Q_mut[a][b]*P_mut_a[a]/P_mut_a[b];
    }
  }
  // Impose normalization
  for(a=0; a<20; a++){
    double sum=0;
    for(b=0; b<20; b++){
      if(a!=b)sum+=Q_mut[a][b];
    }
    Q_mut[a][a]=-sum;
  }
}

void Update_Q(float *Q_mut, float *Q, char *cod2, int j, int n, int nt,
	      float *rate_nuc, float tt_ratio, char **codon, char *coded_aa)
{
  cod2[j]=Nuc_code(n);
  int c2=Code_codon(cod2, codon);
  if(coded_aa[c2]=='*')return; // Stop codon
  int b=Code_AA(coded_aa[c2]);
  (*Q)*=rate_nuc[n];
  if(n==nt)(*Q)*=tt_ratio;
  Q_mut[b]+=*Q;
}
  
void Compute_Q_sel(float **Q_sel, float *P_MF)
{
  int a, b;
  for(a=0; a<20; a++){
    double sum=0;
    for(b=0; b<20; b++){
      if(b==a)continue;
      if(P_MF[b]>=P_MF[a]){Q_sel[a][b]=1;}
      else{Q_sel[a][b]=P_MF[b]/P_MF[a];}
      sum+=Q_sel[a][b];
    }
    Q_sel[a][a]=-sum;
  }
}

void Sum_matrix(float *ncont, int **Cnat, int L){
  int i, j; float sum;
  for(i=0; i<L; i++){
    sum=0; for(j=0; j<L; j++)sum+=Cnat[i][j];
    ncont[i]=sum;
  }
}

float Entropy(float *PP, int n){
  double S=0, norm=0; float *p=PP; int i; 
  for(i=0; i<n; i++){
    if(*p)S+=(*p)*log(*p); norm+=(*p); p++;
  }
  S = S/norm -log(norm);
  return(-S);
}

float Corr_coeff(float *xx, float *yy, int n){
  double x1=0, x2=0, y1=0, y2=0, xy=0;
  int i; float *x=xx, *y=yy;
  for(i=0; i<n; i++){
    x1 += *x; x2+= (*x)*(*x);
    y1 += *y; y2+= (*y)*(*y);
    xy += (*x)*(*y); x++; y++;
  }
  float r=(n*xy-y1*x1)/sqrt((n*x2-x1*x1)*(n*y2-y1*y1));
  return(r);
}

double Mean_hydro(float *P, float *hydro, int n){
  double ave=0; int a;
  for(a=0; a<n; a++)ave+=hydro[a]*P[a];
  return(ave);
}
