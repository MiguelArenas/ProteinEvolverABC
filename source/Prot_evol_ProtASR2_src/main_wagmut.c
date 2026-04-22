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

int CpG=1;

int Read_files(char ***files, int argc, char **argv);
int Read_matrix(double **matrix, double *freq_aa, char *infile);
char *Read_float(float *tmp, char *s);
extern void Fit_nuc_freq(float *mut_par, int *num_aa,
			 char **codon, char *coded_aa, int CpG);
void Change_AA_order(int *iaa, char *AA);
void Print_matrix(float **exchange, double *freq_aa, int *iwag, char *name);
static void Compute_exchange_mut(float **exchange,
				 float *mut_par, float tt_ratio,
				 float TWONUC, float k_CpG, int CpG,
				 char **codon, char *coded_aa);
static void Normalize_flux(float **exchange, float *f, int n);

static void Update_F(float *F_mut, float *S,
		     char *cod, int c, int a,
		     int j, int nt, int n, char *cod2,
		     float *w, float tt_ratio, float k_CpG,
		     char **codon, char *coded_aa);
float Compare_flux(float **m1, float **m2, float *f, int n);

main(int argc, char **argv){
  char **files;
  int n=Read_files(&files, argc, argv), i;

  double *freq_aa=malloc(20*sizeof(double));
  for(i=0; i<20; i++)freq_aa[i]=0;
  double **matrix=Allocate_mat2_d(20, 20);
  double norm=0;
  float **tmp=Allocate_mat2_f(20, 20);
  for(i=0; i<n; i++){
    norm+=Read_matrix(matrix, freq_aa, files[i]); 
  }

  // Compute matrix
  int a,b;
  for(a=0; a<20; a++){
    freq_aa[a]/=norm;
    for(b=0; b<a; b++){
      WAG[b][a]=WAG[a][b];
      if(matrix[a][b]==0)matrix[a][b]=matrix[b][a];
      if(matrix[b][a]==0)matrix[b][a]=matrix[a][b];
    }
  }
  int *iwag=malloc(20*sizeof(int));
  Change_AA_order(iwag, AA_WAG);
  float *fw=malloc(20*sizeof(float));
  float **exchange=Allocate_mat2_f(20, 20);
  for(a=0; a<20; a++){
    int ia=iwag[a];
    fw[ia]=fWAG[a];
    for(b=0; b<20; b++){
      int ib=iwag[b];
      exchange[ia][ib]=WAG[a][b]*(fWAG[a]*fWAG[b]*norm)/matrix[ia][ib];
      if(b<a)printf("%.3f ", exchange[ia][ib]);
    }
    printf("\n");
  }
  printf("%.0f aa read\n", norm);

  // Print matrix
  Normalize_flux(exchange, fw, 20);
  Print_matrix(exchange, freq_aa, iwag, "wagmut");

  // Fit matrix
  double err_min=100000;
  float **exchange_mut=Allocate_mat2_f(20,20);
  float mut_par[9], f1, f2, f3, f4, tt, k_CpG, s=1.02;
  float mut_par_opt[9], tt_opt=0, R, R_opt=0, k_CpG_opt=1;
  float k_CpG_min=1, k_CpG_max=1;

  // Fitting nucleotides
  int num_aa[20], L=2000;
  for(a=0; a<20; a++)num_aa[a]=L*fw[a]; // Fit WAG
  // for(a=0; a<20; a++)num_aa[a]=L*freq_aa[a]; // Fit PDB
  Fit_nuc_freq(mut_par, num_aa, codon, coded_aa, CpG);
  printf("Fitted nucleotide frequencies:  ");
  for(a=0; a<4; a++)printf("%c %.3f ", Nuc_code(a), mut_par[a]);
  if(CpG)printf("f_CpG= %.3f", mut_par[4]); printf("\n");

  int NPAR=4; 
  if(CpG){k_CpG_min=0.7; k_CpG_max=4; NPAR=5;}
  for(R=0; R<=0.4; R+=0.01){
    printf("R= %.3f\n", R);
    /*for(f1=0.12; f1<=0.4; f1+=0.005){
      for(f2=0.12; f2<=0.4; f2+=0.005){
      for(f3=0.12; f3<=0.4; f3+=0.005){
      f4=1-f1-f2-f3; if(f4<0)continue;
      mut_par[0]=f1; mut_par[1]=f2; 
      mut_par[2]=f3; mut_par[3]=f4; */
    for(tt=0.7; tt<3.0; tt+=0.1){
      for(k_CpG=1; k_CpG<=3; k_CpG+=0.1){
	Compute_exchange_mut(exchange_mut, mut_par, tt, R, k_CpG, CpG,
			     codon,coded_aa);
	Normalize_flux(exchange_mut, fw, 20);
	float err=Compare_flux(exchange_mut,exchange,fw,20);
	if(err<err_min){
	  err_min=err; tt_opt=tt; R_opt=R; k_CpG_opt=k_CpG;
	  for(a=0; a<NPAR; a++)mut_par_opt[a]=mut_par[a];
	  printf("%.6f  %.3f %.2f %.2f %.3f %.3f %.3f %.3f",
		 sqrt(err), R, k_CpG, tt_opt,
		 mut_par_opt[0], freq_opt[1], freq_opt[2], freq_opt[3]);
	  if(CpG)printf("  %.2f", freq_opt[4]); printf("\n"); 
	}
      }
    }
  }
  /*    }
	}
	}*/

  printf("Optimal mutation model:\n");
  printf("err=%.6f  TWONUC= %.3f k_CPG= %.3f tt= %.3f ",
	 sqrt(err_min), R_opt, k_CpG_opt, tt_opt);
  printf("%c %.3f %c %.3f %c %.3f %c %.3f",
	 Nuc_code(0), freq_opt[0], Nuc_code(1), freq_opt[1],
	 Nuc_code(2), freq_opt[2], Nuc_code(3), freq_opt[3]); 
  if(CpG)printf("  f_CPG= %.2f", freq_opt[4]); printf("\n"); 
  Compute_exchange_mut(exchange_mut, freq_opt, tt_opt, R_opt, k_CpG_opt,
		       CpG, codon,coded_aa);
  Normalize_flux(exchange_mut, fw, 20);

  float freq_mut[20];
  Mutational_distribution(freq_mut, mut_par, codon, coded_aa, CpG);
  printf("WAG distribution:\n");
  for(a=0; a<20; a++)printf("%.3f\t", fw[a]); printf("\n");
  printf("Fitted mutational distribution:\n");
  for(a=0; a<20; a++)printf("%.3f\t", freq_mut[a]); printf("\n");
  printf("PDB distribution:\n");
  for(a=0; a<20; a++)printf("%.3f\t", freq_aa[a]); printf("\n");
  for(a=0; a<20; a++)printf("%c\t", Amin_code(a)); printf("\n");


  for(a=0; a<20; a++)freq_aa[a]=freq_mut[a];
  Print_matrix(exchange_mut, freq_aa, iwag, "mutmut");

  return(0);

}

void Print_matrix(float **exchange, double *freq_aa, int *iwag, char *what)
{
  int a,b;
  char name[100];
  sprintf(name, "%s.txt", what);
  FILE *file_mat=fopen(name, "w");
  printf("Writing %s\n", name);
  for(a=1; a<20; a++){
    for(b=0; b<a; b++)exchange[b][a]=exchange[a][b];
  }
  for(a=1; a<20; a++){
    float *e=exchange[iwag[a]];
    for(b=0; b<a; b++)fprintf(file_mat, "%.4f\t", e[iwag[b]]);
    fprintf(file_mat, "\n");
  }
  fprintf(file_mat, "\n");
  for(a=0; a<20; a++)fprintf(file_mat, "%.4f\t", freq_aa[iwag[a]]);
  fprintf(file_mat, "\n");
  for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_WAG[a]);
  fprintf(file_mat, "\n");
  for(a=0; a<20; a++)fprintf(file_mat, "%d\t", iwag[a]);
  fprintf(file_mat, "\n");
  fclose(file_mat);

  sprintf(name, "%s.h", what);
  file_mat=fopen(name, "w");
  printf("Writing %s\n", name);
  fprintf(file_mat, "float %s[20][20]={\n", what);
  for(a=0; a<20; a++){
    float *e=exchange[iwag[a]];
    fprintf(file_mat, "{");
    for(b=0; b<a; b++)fprintf(file_mat, "%.4f,", e[iwag[b]]);
    for(b=a; b<19; b++)fprintf(file_mat, "0,");
    if(a<19){fprintf(file_mat, "0},\n");}
    else{fprintf(file_mat, "0}\n");}
  }
  fprintf(file_mat, "};\n");

  fprintf(file_mat, "float fPDB[20]={");
  for(a=0; a<19; a++){
    fprintf(file_mat, "%.4f,", freq_aa[iwag[a]]);
  }
  fprintf(file_mat, "%.4f};\n", freq_aa[iwag[19]]);
  fclose(file_mat);
}

int Read_files(char ***files, int argc, char **argv){
  if(argc<2){
    printf("ERROR, file name must be provided\n"); exit(8);
  }
  FILE *file_in=fopen(argv[1], "r");
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", argv[1]); exit(8);
  }
  int n=0, i; char string[1000];
  while(fgets(string, sizeof(string), file_in)!=NULL)n++;
  fclose(file_in);
  if(n==0){
    printf("ERROR, no file names read in %s\n", argv[1]); exit(8);
  }


  file_in=fopen(argv[1], "r");
  (*files)=malloc(n*sizeof(char *)); n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    (*files)[n]=malloc(100*sizeof(char));
    sscanf(string, "%s", (*files)[n]); n++;
  }
  fclose(file_in);

  printf("%d files to read\n", n);
  return(n);
}

int Read_matrix(double **matrix, double *freq_aa, char *infile)
{
  FILE *file_in=fopen(infile, "r");
  if(file_in==NULL){
    printf("WARNING, file %s nor found\n", infile); return(0);
  }

  char string[1000], aa[20][3]; int read=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='A'){
      sscanf(string, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
	     aa[0], aa[1], aa[2], aa[3], aa[4], 
	     aa[5], aa[6], aa[7], aa[8], aa[9],
	     aa[10],aa[11],aa[12],aa[13],aa[14], 
	     aa[15],aa[16],aa[17],aa[18],aa[19]);
      read=1; break;
    }
  }
  fclose(file_in);

  if(read==0){
    printf("WARNING, string with amino acid codes not found in %s\n",infile);
    return(0);
  }

  int *iaa=malloc(20*sizeof(int)), a, b;
  for(a=0; a<20; a++)iaa[a]=Code_AA(aa[a][0]);
  //for(a=0; a<20; a++)printf("%c\t", aa[a][0]); printf("\n");
  //for(a=0; a<20; a++)printf("%d\t", iaa[a]); printf("\n");

  int L=0; float tmp;
  file_in=fopen(infile, "r"); a=1;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='L'){
      sscanf(string+2, "%d", &L);
    }else if(a<20){
      double *m=matrix[iaa[a]]; char *s=string;
      for(b=0; b<a; b++){
	s=Read_float(&tmp, s); m[iaa[b]]+= tmp;
      }
      a++;
    }else if(string[0]!='\n'){
      char *s=string;
      for(b=0; b<20; b++){
	s=Read_float(&tmp, s); freq_aa[iaa[b]]+= L*tmp;
      }
      break;
    }
  }
  fclose(file_in);
  printf("%3d residues and %d aa read in %s\n", L, a, infile);
  return(L);
}

char *Read_float(float *tmp, char *s){
  char *t=s;
  sscanf(s, "%f", tmp);
  while((*t!='\t')&&(*t!=' ')&&(*t!='\n'))t++;
  while((*t=='\t')||(*t==' '))t++;
  return(t);
}

void Change_AA_order(int *iaa, char *AA){
  int a; for(a=0; a<20; a++)iaa[a]=Code_AA(AA[a]);
}

void Compute_exchange_mut(float **S_mut,
			  float *mut_par, float tt_ratio,
			  float TWONUC, float k_CpG, int CpG,
			  char **codon, char *coded_aa)
{
  // S(a,b)=sum_{c\in C(a) c'\in C(b)} w_c*w_c'*S(c,c')/(W_a*W_b)
  int a, b, c;
  for(a=0; a<20; a++)for(b=0; b<20; b++)S_mut[a][b]=0;
  double *norm=malloc(20*sizeof(double));
  for(a=0; a<20; a++)norm[a]=0;
  float w[64];
  if(CpG){
    for(c=0; c<64; c++)w[c]=Weight_codon_CpG(codon[c], mut_par);
  }else{
    for(c=0; c<64; c++)w[c]=Weight_codon(codon[c], mut_par);
  }
  float r=1; if(TWONUC>0)r=TWONUC; // rate

  for(c=0; c<64; c++){
    if(coded_aa[c]=='*')continue; // Stop codon
    a=Code_AA(coded_aa[c]); norm[a]+=w[c];
    //printf("%2d a=%c b= ", c, Amin_code(a));
    float *Fa=S_mut[a], rc=r*w[c];
    int j, n; char cod2[3], *cod=codon[c];
    for(j=0; j<3; j++)cod2[j]=cod[j];
    for(j=0; j<3; j++){
      char nj=cod[j]; int nt=Transition(nj);
      for(n=0; n<4; n++){
	if(Nuc_code(n)==nj)continue;
	// cod2[j]=n
	float S=r;
	Update_F(Fa, &S, cod, c, a, j, nt, n, cod2, w,
		 tt_ratio, k_CpG, codon, coded_aa);
	if(TWONUC>0){
	  int j2, n2, j3, n3;
	  for(j2=j+1; j2<3; j2++){
	    int nt2=Transition(cod[j2]);
	    for(n2=0; n2<4; n2++){
	      if(Nuc_code(n2)==cod[j2])continue;
	      float S2=S*r;
	      Update_F(Fa, &S2, cod, c, a, j2, nt2, n2, cod2, w,
		       tt_ratio, k_CpG, codon, coded_aa);
	      for(j3=j2+1; j3<3; j3++){
		int nt3=Transition(cod[j3]);
		for(n3=0; n3<4; n3++){
		  if(Nuc_code(n3)==cod[j3])continue;
		  float S3=S2*3;
		  Update_F(Fa, &S3, cod, c, a, j3, nt3, n3, cod2, w,
			   tt_ratio, k_CpG, codon, coded_aa);
		  cod2[j3]=cod[j3];
		}
	      }
	      cod2[j2]=cod[j2];
	    }
	  }
	}
	cod2[j]=cod[j];
      } // End n loop
    } // End j loop
  } // End c loop
  
  // Impose symmetry
  double sum=0;
  for(a=0; a<20; a++)sum+=norm[a];
  for(a=0; a<20; a++)norm[a]/=sum;
  for(a=0; a<20; a++){
    for(b=0; b<a; b++){
      S_mut[a][b]=(S_mut[a][b]+S_mut[b][a])/(2*norm[a]*norm[b]);
      S_mut[b][a]=S_mut[a][b];
    }
  }

  // Impose normalization Q_aa=-sum_b Q_ab
  /*for(a=0; a<20; a++){
    double sum=0;
    for(b=0; b<20; b++){
      if(a!=b)sum+=S_mut[a][b]*norm[b];
    }
    S_mut[a][a]=-sum/norm[a];
    }*/
}

void Update_F(float *F_mut, float *S,
	      char *cod, int c, int a,
	      int j, int nt, int n, char *cod2,
	      float *w, float tt_ratio, float k_CpG,
	      char **codon, char *coded_aa)
{
  cod2[j]=Nuc_code(n);
  int c2=Code_codon(cod2, codon);
  if(coded_aa[c2]=='*')return; // Stop codon
  int b=Code_AA(coded_aa[c2]); if(b==a)return;
  //printf("%c ", Amin_code(b));
  if(n==nt){
    (*S)*=tt_ratio;
    if((strncmp(cod, "CG", 2)==0)||(strncmp(cod+1, "CG", 2)==0)){
      (*S)*=k_CpG;
    }
  }
  F_mut[b]+=w[c]*w[c2]*(*S);
}

void Normalize_flux(float **exchange, float *f, int n)
{
  double sum=0;
  int i, j;
  for(i=1; i<n; i++){
    for(j=0; j<i; j++){
      sum+=exchange[i][j]*f[i]*f[j];
    }
  }
  sum *=2;
  for(i=1; i<n; i++){
    for(j=0; j<i; j++)exchange[i][j]/=sum;
  }
}

float Compare_flux(float **m1, float **m2, float *f, int n){
  double sum=0, e; int i, j;
  for(i=1; i<n; i++){
    for(j=0; j<i; j++){
      e=m1[i][j]-m2[i][j]; sum+=(e*e)*f[i]*f[j];
    }
  }
  return(2*sum);
}

