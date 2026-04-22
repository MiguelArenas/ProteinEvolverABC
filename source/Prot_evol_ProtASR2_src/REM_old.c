#include "REM.h"
#include "energy_BKV.h"  // Contact energy parameters
#include "protein3.h"
#include "allocate.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "output.h"
#include "Sec_str_comp.h"

// Parameters for selecting compact structures
#define NAA 21 // Number of amino acid types
#define C_ASYMPT 4.00  // Asymptotic value of Nc/N
#define C_SLOPE 8.07   // Nc/N ~ C_ASYMPT - C_SLOPE N^(-1/3)
#define C_THR 1.0      // Maximum deviation of Nc/N from typical value
                       // Very important parameter!!!
int E_36=0; // WARNING!!! If E_36=1 and C_THR is large, there will be problems
int E343_2=0;
float T0=-1;
int LEN=0;
int LEN_MAX=300;

#define L_MAX 100000
#define NCMAX 30    // Maximum number of contacts per residue

// Common parameters defined in .h files
//float **interactions;
struct protein prot[N_PROT_MAX], target; int N_PROT;

double *Cont_freq=NULL, NC1=0, NC2=0, NC3=0;
double *C221_ij=NULL, *C232_i=NULL, C242=0;
double *C321_ij=NULL, **C332_ij=NULL, *C342_ij=NULL;
double **C343_2_ij=NULL, *C343_3=NULL, *C353_i=NULL, C363=0;
double *nc1i=NULL, *nc2i=NULL, *nc3i=NULL, *CNc=NULL;
double *NC1nc1i=NULL, *NC2nc1i=NULL, *NC1nc2i=NULL;
double **Cnn=NULL, **ncnc=NULL, **C1n;
double nc2_sum, nc12_sum;

void Initialize_REM(int *LEN, int L, char *file_str);
void Contact_range(int *NCmin, int *NCmax, int L);
void Compute_E343_2(struct REM *E, short *aa_seq);
int Read_structures(char *file_prot, struct protein *prot);
void Threading_energies(double *EC1, double *EC2, double *EC3, double *U2,
			short *aa_seq, int len_nat);
float G_misfold_threading(struct REM *E);

/*****************************************************************/

void Initialize_E_REM(struct REM *E, int L, int REM,
		      float T, float S_C, float S_U,
		      char *file_str, int ini)
{
  // Temperature
  if(Econt_T==NULL)Econt_T=Allocate_mat2_f(NAA, NAA);
  if(T!=T0){
    printf("Dividing contact energies by T= %.2f\n", T);
    T0=T; int i, j;
    for(i=0; i<NAA; i++){
      for(j=0; j<NAA; j++)Econt_T[i][j]=Econt[i][j]/T;
    }
    Initialize_E_loc(T, SEC_STR);
  }
  // Statistics of misfolded structures
  if(L != LEN && LEN!=LEN_MAX){
    Initialize_REM(&LEN, L, file_str);
  }
  if(ini){free(E->c1U1); free(E->c332U2);}
  E->L=L;
  E->c1U1=malloc(L*sizeof(double));
  E->c332U2=malloc(L*sizeof(double));
  E->REM=REM;
  E->T=T;
  E->S_C=S_C;
  E->S_U=S_U;
  //E->S_C=sC0+sC1*L;
  //E->S_U=sU1*L;
  Set_to_zero(E);
}

void Set_to_zero(struct REM *E)
{
  E->E_nat=0;
  E->E_loc=0;
  E->E1=0;
  E->E22=0;
  E->E23=0;
  E->E321=0; E->E332=0;
  E->E342=0; E->E343_2=0; E->E343_3=0;
  E->E353=0; E->E363=0;
  E->e342=0;
  for(int i=0; i<E->L; i++){
    E->c1U1[i]=0; E->c332U2[i]=0;
  }
}

void Initialize_REM(int *LEN, int L_pdb, char *file_str)
{
  if(*LEN==L_pdb)return;

  // Empty_memory
  if((*LEN)){
    if(Cont_freq)free(Cont_freq);
    if(C221_ij)free(C221_ij);
    if(C232_i)free(C232_i);
    if(nc1i)free(nc1i);
    if(nc2i)free(nc2i);
    if(nc3i)free(nc3i);
    if(NC1nc1i)free(NC1nc1i);
    if(NC2nc1i)free(NC2nc1i);
    if(NC1nc2i)free(NC1nc2i);
    if(CNc)free(CNc);
    if(C321_ij)free(C321_ij);
    if(C332_ij)Empty_matrix_d(C332_ij, *LEN);
    if(C342_ij)free(C342_ij);
    if(C343_2_ij)Empty_matrix_d(C343_2_ij, *LEN);
    if(C343_3)free(C343_3);
    if(Cnn)Empty_matrix_d(Cnn, *LEN);
    if(C1n)Empty_matrix_d(C1n, *LEN);
    if(ncnc)Empty_matrix_d(ncnc, *LEN);
  }
  int L=L_pdb; if(L>LEN_MAX)L=LEN_MAX;
  *LEN=L;

  // Allocate
  int i;
  NC1=0; NC2=0; NC3=0;
  Cont_freq =malloc(L_pdb*sizeof(double));
  C221_ij   =malloc(L_pdb*sizeof(double));
  C232_i    =malloc(L_pdb*sizeof(double));
  nc1i =malloc(L_pdb*sizeof(double));
  nc2i =malloc(L_pdb*sizeof(double));
  nc3i =malloc(L_pdb*sizeof(double));
  NC1nc1i =malloc(L_pdb*sizeof(double));
  NC2nc1i =malloc(L_pdb*sizeof(double));
  NC1nc2i =malloc(L_pdb*sizeof(double));
  CNc  =malloc(L_pdb*sizeof(double));

  C321_ij = malloc(L_pdb*sizeof(double));
  C332_ij = Allocate_mat2_d(L_pdb,L_pdb);
  C342_ij = malloc(L_pdb*sizeof(double));
  C343_2_ij=Allocate_mat2_d(L_pdb,L_pdb);
  C343_3 = malloc(L_pdb*sizeof(double));
  C353_i = malloc(L_pdb*sizeof(double));

  Cnn  =Allocate_mat2_d(L_pdb,L_pdb);
  ncnc =Allocate_mat2_d(L_pdb,L_pdb);
  C1n  =Allocate_mat2_d(L_pdb,L_pdb);

  for(i=0; i<L_pdb; i++){
    Cont_freq[i]=0; CNc[i]=0;
    nc1i[i]=0; nc2i[i]=0; nc3i[i]=0;
    NC1nc1i[i]=0; NC1nc2i[i]=0; NC2nc1i[i]=0;
  }

  N_PROT=Read_structures(file_str, prot);
  int NCmin, NCmax;
  Contact_range(&NCmin, &NCmax, L);
  printf("### Contact statistics\n");
  printf("%d contact matrices found in file %s\n", N_PROT, file_str);
  printf("Selecting fragments of %d residues with %d < NC < %d\n",
	 L, NCmin, NCmax);
  printf("(NC/L ~ %.2f - %.2f^(-1/3) +/- %.2f)\n", C_ASYMPT,C_SLOPE,C_THR);

  /* Generate alternative structures by threading */
  int *nci=malloc(L*sizeof(int));
  short j_prot=0, first, res1, *res2, j;
  long n_str=0, Noncompact=0;
  for(j_prot=0; j_prot< N_PROT; j_prot++){

    short len2=prot[j_prot].length;
    if(len2 < L)continue;
    short **contact=prot[j_prot].contact;

    for(first=0; first<=len2-L; first++){
      
      /* New fragment */
      int last=first+L, Nc=0;
      for(i=0; i<L; i++)nci[i]=0;
      for(res1=first; res1<last; res1++){
	res2=contact[res1]; i=res1-first;
	int res2_min=res1+IJ_MIN;
	while((*res2>=0)&&(*res2 < last)){
	  if(*res2 >= res2_min){
	    nci[i]++; nci[*res2-first]++;
	  }
	  res2++;
	}
	Nc+=nci[i]; 
      }
      Nc/=2;
      // Exclude non compact conformations and average
      if((Nc<NCmin)||(Nc>NCmax)){Noncompact++; continue;}
      n_str++;
      double Nc2=Nc*Nc;
      NC1+=Nc; NC2+=Nc2; NC3+=(double)(Nc2*Nc);
      for(i=0; i<L; i++){
	long n=nci[i];nc1i[i]+=n;
	long n2=n*n;  nc2i[i]+=n2;
	long n3=n2*n; nc3i[i]+=n3;
	NC1nc1i[i]+=n*Nc; NC1nc2i[i]+=n2*Nc;
	NC2nc1i[i]+=Nc2*n;
	if(ncnc)for(j=i+1; j<L; j++)ncnc[i][j]+=n*nci[j];
      } 
      for(res1=first; res1<last; res1++){
	res2=contact[res1]; i=res1-first;
	//if(i>=L)break;
	while((*res2>=0)&&(*res2 < last)){
	  if(*res2 >= res1+IJ_MIN){
	    j=*res2-first; // j>i
	    Cont_freq[j-i]++;
	    CNc[j-i]+=Nc;
	    if(Cnn){
	      Cnn[i][j]+=nci[i]*nci[j];
	      C1n[i][j]+=nci[i];
	      C1n[j][i]+=nci[j];
	    }
	  }
	  res2++;
	}
      }
    }
  }
  float discard=Noncompact/(float)(Noncompact+n_str);
  printf("%ld fragments used for statistics\n", n_str);
  printf("Fraction of discarded contact matrices: %.3f\n", discard);
  printf("Minimum contact distance IJ_MIN= %d\n", IJ_MIN);
  free(nci);

  /* End of statistics, compute averages */
  NC1/=n_str; NC2/=n_str; NC3/=n_str;
  nc2_sum=0; nc12_sum=0;
  for(i=0; i<L; i++){
    nc1i[i]/=n_str;
    nc2i[i]/=n_str;
    nc3i[i]/=n_str;
    NC1nc1i[i]/=n_str;
    NC1nc2i[i]/=n_str;
    NC2nc1i[i]/=n_str;
    nc2_sum  += nc1i[i]*nc1i[i];
    nc12_sum += nc2i[i];
    if(ncnc){
      for(j=i+1; j<L; j++){
	ncnc[i][j]/=n_str;
	ncnc[j][i]=ncnc[i][j];
	Cnn[i][j]/=n_str;
	Cnn[j][i]=Cnn[i][j];
	C1n[i][j]/=n_str;
	C1n[j][i]/=n_str;
      }
    }
  }
  int Nl=1, l;
  double sum_c2=0;
  for(l=L-1; l>=1; l--){
    Cont_freq[l]/=(n_str*Nl);
    CNc[l]/=(n_str*Nl);
    sum_c2 += Cont_freq[l]*Cont_freq[l]*Nl;
    Nl++;
  }

  double sum_c2_i[L];
  for(i=0; i<L; i++)sum_c2_i[i]=0;
  for(i=0; i<L; i++){
    for(int j=0; j<L; j++){
      int l=i-j; if(l<IJ_MIN)continue;
      float c2=Cont_freq[l]*Cont_freq[l];
      sum_c2_i[i]+=c2;
      sum_c2_i[j]+=c2;
    }
  }

  double NC1_2= NC1*NC1, NC1_3=NC1_2*NC1;

  // Compute coefficients of the free energy, 2nd mom
  C242 = ((NC2-nc2_sum)/(NC1_2-nc12_sum)-1)/2;
  double C232_ave=0;
  for(i=0; i<L; i++){
    // Site i in contact with two other sites 
    C232_i[i] = ((nc2i[i]-nc1i[i])/(nc1i[i]*nc1i[i]-sum_c2_i[i])-1)/2;
    C232_ave += C232_i[i];
    C232_i[i]-= C242; // C242/4
  }
  C232_ave/=L;
  for(l=IJ_MIN; l<L; l++){
    double c1=Cont_freq[l], c2 = c1*c1;
    C221_ij[l] = (c1-c2)/2 - 2*c2*C232_ave;
  }

  // Compute coefficients of the free energy, 3rd mom
  // 363           //(ij)(kl)(mn)  mult=1 -> NULL ->353
  // 353(i);       //(ij)(ik)(lm)  mult=3 -> 363  -> 34X
  // 343_3(i);     //(ij)(ik)(il)  mult=1 -> 353  -> 332 
  // 343_2(i,j);   //(ij)(ik)(jl), (ij)(ik)(jk) mult=6 -> 353 -> 332
  // 342_ij[L];    //(ij)^2(kl)  mult=3 -> 353 -> 332
  // 332_ij(i,j);  //(ij)^2(ik)  mult=3 -> 34X -> 321
  // 321_ij[(i-j)  //(ij)^3 mult=1 -> 332 -> NULL */

  // i-j: A321 and A342
  double A342_ij[L]; 
  double A332_sum[L], norm332[L];
  for(l=0; l<L; l++){
    if(l>=IJ_MIN){
      double c1=Cont_freq[l], c2 = c1*c1, c3=c2*c1;
      C321_ij[l] = (c1-3*c2+2*c3);            // B321
      A342_ij[l] = (CNc[l]-NC1*c1)*(1-2*c1);  // A342
    }else{
      C321_ij[l] = 0;
      A342_ij[l] = 0;
    }
    A332_sum[l]=0;
    norm332[l]=0;
  }

  // Pairwise: A332 (not symm) and A343_2 (symm)
  double A332_ij[L][L];
  double A343_2[L][L], W343_2[L][L];
  for(i=0; i<L; i++){
    for(j=0; j<L; j++){
      A332_ij[i][j]=0;
      A343_2[i][j]=0;
      W343_2[i][j]=0;
    }
  }
  for(i=0; i<L; i++){
    float nc1=nc1i[i];
    for(j=i+IJ_MIN; j<L; j++){
      int l=j-i;
      double c1=Cont_freq[l];
      A332_ij[i][j]=(C1n[i][j]-c1*nc1)*(1-2*c1);
      A332_ij[j][i]=(C1n[j][i]-c1*nc1i[j])*(1-2*c1);
      
      A332_sum[l]+=(A332_ij[i][j]+A332_ij[j][i]);
      norm332[l]+=2;

      A343_2[i][j]=Cnn[i][j]-c1*ncnc[i][j]-
       nc1*C1n[i][j]-nc1i[j]*C1n[j][i]+2*c1*nc1*nc1i[j];
      W343_2[i][j]=nc1*nc1i[j];
      A343_2[j][i]=A343_2[i][j];
      W343_2[j][i]=W343_2[i][j];
    }
  }
  for(l=IJ_MIN; l<L; l++)A332_sum[l]/=norm332[l];

  double C332_sum[L];
  for(l=IJ_MIN; l<L; l++){C332_sum[l]=0; norm332[l]=0;}
  double W342=NC1;
  double A353_sum=0, W353_sum=0;
  double B353_ave=0;
  for(i=0; i<L; i++){
    float nc1=nc1i[i], nc2=nc1*nc1, nc3=nc2*nc1;
    float c_ave=sum_c2_i[i]/nc1;
    double A343_2_sum=0, W343_2_sum=0;
    double A342_sum=0, W342_sum=0;
    for(j=0; j<L; j++){
      int l=abs(j-i);
      if(l<IJ_MIN){
	C332_ij[i][j]=0;
	C343_2_ij[i][j]=0;
      }else{
	double c1=Cont_freq[l];
	C332_ij[i][j]=(A332_ij[i][j]-C321_ij[l])/(nc1-c1); // B332
	C332_ij[i][j]/=2;
	C332_sum[l]+=C332_ij[i][j]; norm332[l]++;

	A342_sum+=A342_ij[l]; W342_sum+=W342;
	
	A343_2_sum += A343_2[i][j];
	W343_2_sum += W343_2[i][j];
	C343_2_ij[i][j]=(A343_2[i][j]-A332_ij[i][j]-A332_ij[j][i])/
	  (W343_2[i][j]-nc1i[i]-nc1i[j]);                // B343_2
      }
    }
    double A332_sum_i=(nc2i[i]-nc2)*(1-2*c_ave);
    double A343_3=nc3i[i]-3*nc2i[i]*nc1+2*nc3, W343_3=nc3;
    C343_3[i]=(A343_3-A332_sum_i)/              // B343_3
      (W343_3-nc1*sum_c2_i[i]);
    C343_3[i]/=6;

    double A353_i=NC1nc2i[i]-NC1*nc2i[i]-2*nc1*(NC1nc1i[i]-NC1*nc1);
    double W353_i=NC1*nc2;
    A353_sum+=A353_i; W353_sum+=W353_i;
    C353_i[i]=(A353_i-A342_sum-2*A343_2_sum-A343_3)/ 
      (W353_i-W342_sum-2*W343_2_sum-W343_3); // B353
    C353_i[i]/=6;
    B353_ave+=C353_i[i];
  }
  B353_ave/=L;

  double A363=(NC3-3*NC2*NC1+2*NC1_3), W363=NC1_3;
  C363 = (A363-3*A353_sum)/(W363-3*W353_sum); 
  C363/=6;  // C363

  for(l=IJ_MIN; l<L; l++){
    C321_ij[l] /= 6;          // C321
    //C332_sum[l]/=norm332[l];
    //C321_ij[l] -= C332_sum[l];
    C342_ij[l]=(A342_ij[l]-A332_sum[l])/(W342);
    C342_ij[l] /= 2;
    //C342_ij[l]-=2*B353_ave; // C342
  }
  /*float C36=3*C363;
  for(i=0; i<L; i++){
    for(j=0; j<L; j++){
      C332_ij[i][j]-=  // C332
	(C342_ij[abs(i-j)]+C343_2_ij[i][j]+C343_3[i]+C343_3[j]);
      C343_2_ij[i][j]-=2*C353_i[i];  // C343_2
    }
    C343_3[i]-=C353_i[i];  // C343_3
    C353_i[i]-=C36;          // C353   
    }*/

  double Z2= NC1*NC1, Z3=NC1*NC1*NC1;
  NC3=NC3-3*NC1*NC2+2*Z3;
  NC2=NC2-Z2;
  printf("<NC>= %.0f M2(NC)/<NC>^2= %.3g M3(NC)/<NC>^3= %.3g\n",
	 NC1, NC2/Z2, NC3/Z3);

  char nameout[400];
  sprintf(nameout, "Contact_statistics_%d.dat", L);
  printf("Writing %s\n", nameout);
  FILE *file_out=fopen(nameout, "w");
  fprintf(file_out, "# Contact statistics for protein of length %d\n", L);
  fprintf(file_out, "# Number of fragments used in statistics: %ld\n", n_str);
  fprintf(file_out, "# Fraction of discarded contact matrices: %.3f\n",discard);
  fprintf(file_out, "# (NC/L ~ %.2f - %.2f^(-1/3) +/- %.2f)\n",
	  C_ASYMPT,C_SLOPE,C_THR);
  fprintf(file_out, "# <NC>= %.2f /L = %.3g\n", NC1, NC1/L);
  fprintf(file_out, "# <NC^2>-<NC>^2= %.2f  /<NC>= %.3f\n", NC2, NC2/NC1);
  fprintf(file_out, "# <(NC-<NC>)^3>= %.2f  /<NC>= %.3f\n", NC3, NC3/NC1);
  NC3/=(6*Z3);
  NC2/=(2*Z2);
  fprintf(file_out, "# [<NC^2>/<NC>^2-1]/2= %.2g  [<(NC/<NC>-1)^3>]/6= %.2g\n",
	  NC2, NC3);

  fprintf(file_out, "#l Cont_freq(l) C221(l) C232(i)");
  fprintf(file_out, " C321(l) C332(i,4) C342(l) C343_2(l) C343_3(l) C353(i)\n");

  for(l=IJ_MIN; l<L; l++){
    fprintf(file_out, "%3d %.5f %.5f %.5f ",
	    l, Cont_freq[l], C221_ij[l], C232_i[l]);
    fprintf(file_out,  "%.5f %.5f %.5f %.5f %.5f %.5f\n",
	    C321_ij[l], C332_ij[l-IJ_MIN][l], C342_ij[l],
	    C343_2_ij[l-IJ_MIN][l], C343_3[l], C353_i[l]);   

  }
  fprintf(file_out, "#l Cont_freq(l) C221(l) C232(i)\n");

  fclose(file_out);
  return;
}

float Compute_DG_overT_threading(struct REM *E,
				short *aa_seq, int **C_nat, int *i_sec)
{

  Set_to_zero(E);

  // Native energy
  E->E_loc=Compute_E_loc(i_sec, aa_seq, E->L);
  int i, j;
  for(i=0; i<E->L; i++){
    float *U=Econt_T[aa_seq[i]];
    for(j=0; j<E->L; j++){
      if((j>i)&&(C_nat[i][j]))E->E_nat+=U[aa_seq[j]];
    }
  }

  // Misfolding energy
  double U2;
  Threading_energies(&E->E1, &E->E2, &E->E3, &U2, aa_seq, E->L);
  E->G_misf=G_misfold_threading(E);
  float DeltaG_overT=DeltaG((float)(E->E_nat), E->G_misf, E->S_U);
  return(DeltaG_overT);
}

float G_misfold_threading(struct REM *E)
{
  if(E->REM==0){
    E->E1=0; E->E2=0; E->E3=0;
  }else if(E->REM==1){
    E->E2=0; E->E3=0;
  }else if(E->REM==2){
    E->E3=0;
  }
  E->Tf=Compute_Tfreeze(E);
  float T=1; if(E->Tf>T){T=E->Tf;}
  double H=E->E1-E->E2/T;
  if(E->REM==3)H+=E->E3/(T*T);
  return(H-E->S_C);
}

float Compute_DG_overT_contfreq(struct REM *E,
				short *aa_seq, int **C_nat, int *i_sec)
{
  Set_to_zero(E);
  E->E_loc=Compute_E_loc(i_sec, aa_seq, E->L);
  int i, j;
  for(i=0; i<E->L; i++){
    if(aa_seq[i]<0)continue;
    float *U=Econt_T[aa_seq[i]];
    for(j=0; j<E->L; j++){
      if(aa_seq[j]<0)continue;
      if((j>i)&&(C_nat[i][j]))E->E_nat+=U[aa_seq[j]];
      if(E->REM==0)continue;
      int l=abs(i-j); if(l<IJ_MIN)continue;
      if(l>=LEN_MAX)continue; //l=LEN_MAX-1;
      float c=Cont_freq[l];
      float Uij=U[aa_seq[j]], cU=c*Uij;
      E->c1U1[i]+=cU;
      if(E->REM >= 2){
	float Uij2=Uij*Uij;
	if(j>i){
	  E->E22  += C221_ij[l]*Uij2;
	  if(E->REM==3){
	    E->E321 += C321_ij[l]*Uij2*Uij;
	    E->e342 += C342_ij[l]*Uij2;
	  }
	}
	if(E->REM==3)
	  E->c332U2[i]+= C332_ij[i][j]*Uij2;
      }
    }
    E->E1  += E->c1U1[i];
  }
  E->E1/=2;


  if(E343_2) Compute_E343_2(E, aa_seq);
  
  E->G_misf=G_misfold(E);
  float DeltaG_overT=DeltaG((float)(E->E_nat), E->G_misf, E->S_U);
  return(DeltaG_overT);
}

float Mutate_DG_overT_contfreq(struct REM *E, short *aa_seq,
			       int **C_nat, int *is,
			       int res_mut, short aa_new)
{
  int i_sec=is[res_mut];
  int *Cnat_i=C_nat[res_mut], j, l;
  float *inter_mut=Econt_T[aa_new];
  float *inter_wt=Econt_T[aa_seq[res_mut]];
  double c1U1=0, c332U2=0;
  float dEl=Delta_E_loc(i_sec, aa_seq[res_mut], aa_new);
  for(j=0; j<E->L; j++){
    l=abs(j-res_mut);
    if(l<IJ_MIN)continue;
    if(l>=LEN_MAX)continue; //l=LEN_MAX-1;
    float U_mut=inter_mut[aa_seq[j]];
    float U_wt=inter_wt[aa_seq[j]];
    float d=U_mut-U_wt;
    if(Cnat_i[j])E->E_nat+=d;
    if(E->REM==0)continue;
    float c=Cont_freq[l], cd=c*d;
    c1U1+=cd;
    E->c1U1[j]+=cd;
    if(E->REM==1)continue;
    float U2_mut=U_mut*U_mut, U2_wt=U_wt*U_wt;
    float d2=U2_mut-U2_wt;
    float cd2=c*d2;
    E->E22 += C221_ij[l]*d2;
    if(E->REM==3){
      E->E321  += C321_ij[l]*(U2_mut*U_mut-U2_wt*U_wt);
      c332U2+= C332_ij[res_mut][j]*d2;
      E->c332U2[j]+=C332_ij[res_mut][j]*cd2;
      E->e342  += C342_ij[l]*d2;
    }
  }
  E->E_nat+=dEl;
  if(E->REM){
    E->E1+=c1U1;
    E->c1U1[res_mut]+=c1U1;
    if(E->REM==3)E->c332U2[res_mut]+=c332U2;
    if(E343_2) Compute_E343_2(E, aa_seq);
  }
  E->G_misf=G_misfold(E);
  float DeltaG_overT=DeltaG((float)(E->E_nat), E->G_misf, E->S_U);
  return(DeltaG_overT);
}

float G_misfold(struct REM *E)
{
  if(E->REM==0){
    E->E1=0; E->E2=0; E->E3=0;
  }else if(E->REM==1){
    E->E2=0; E->E3=0;
  }else if(E->REM>=2){
    E->E23=0;
    E->E332=0; E->E343_3=0; E->E353=0;
    for(int i=0; i<E->L; i++){
      if(i>=LEN_MAX)continue;
      float u12 = E->c1U1[i]*E->c1U1[i];
      E->E23  += C232_i[i]*u12;
      if(E->REM==3){
	E->E332  += E->c1U1[i]*E->c332U2[i];
	E->E343_3+= E->c1U1[i]*u12*C343_3[i];
	E->E353  += C353_i[i]*u12;
      }
    }
    E->E24=(C242)*E->E1*E->E1;
    E->E2 = E->E22+E->E23+E->E24;
    if(E->REM==3){
      E->E342   = E->e342*E->E1;
      E->E353 *= E->E1;
      E->E363=(C363)*E->E1*E->E1*E->E1;
      E->E3=  E->E321 + E->E332 + E->E342 + E->E343_3 + E->E353 + E->E363;
      if(E343_2) E->E3 += E->E343_2;
    }
  }
  E->Tf=Compute_Tfreeze(E);
  float T=1; if(E->Tf>T){T=E->Tf;}
  double H=E->E1-E->E2/T;
  if(E->REM==3)H+=E->E3/(T*T);
  return(H-E->S_C);
}


float DeltaG(float E_nat, float G_misf, float S_U)
{
  // DeltaG/T =log(exp((G_N-G_U)/T)+exp((G_N-G_M)/T))
  //       = G_N/T +log(exp(-G_U/T)+exp(-G_M/T))
  // Note that G_U/T =-S_U
  double fU=exp(E_nat+S_U);   // 
  double fm=exp(E_nat-G_misf);
  float DeltaG_overT=log(fU+fm);
  return(DeltaG_overT);
}

float DeltaG_old(float E_nat, float G_misf, float S_U)
{
  // DeltaG/T = G_N/T+log(exp(-G_U/T)+exp(-G_M/T))
  //          = (G_N-G_M)/T + log(1+exp((G_M-G_U)/T))
  // Note that G_U/T =-S_U
  double H=G_misf+S_U, delta;
  if(H>30){delta=H;}
  else if(H<-30){delta=0;}
  else{delta=log(1+exp(H));}
  float DeltaG_overT=E_nat-G_misf+delta;
  return(DeltaG_overT);
}

float Compute_Tfreeze(struct REM *E){
  // Freezing temperature: S_C-E2/T^2+2*E3/T^3=0
  float Tf=sqrt(E->E2/E->S_C);
  if(E->REM<3)return(Tf);
  if((E->E3>0)&&(E->E3 >(E->E2*Tf/sqrt(27))))return(3*E->E3/E->E2);
  float K3=2*E->E3, K2=E->E2, S3=3*E->S_C;
  for(int it=0; it<20; it++){
    float T2=Tf*Tf;
    float f=E->S_C*T2*Tf-K2*Tf+K3;
    float f1=S3*T2-K2;
    if(fabs(f)<0.001)return(Tf);
    Tf-=(f/f1);
  }
  return(Tf); 
}

void Compute_E343_2(struct REM *E, short *aa_seq){
  double sum=0;
  for(int i=0; i<E->L; i++){
    if(aa_seq[i]<0)continue; if(i>=LEN_MAX)continue;
    float *U=Econt_T[aa_seq[i]]; double sum_i=0;
    for(int j=i+IJ_MIN; j<E->L; j++){
      if(aa_seq[j]<0)continue;
      if(j>=LEN_MAX)continue;
      sum_i += C343_2_ij[i][j]*U[aa_seq[j]]*E->c1U1[j];
    }
    sum += sum_i*E->c1U1[i];
  }
  E->E343_2=2*sum;
}

float Compute_hydro(short *aa_seq, int L)
{
  double h=0; int n=0;
  for(int i=0; i<L; i++){
    int a=aa_seq[i];
    if((a>=0)&&(a<20)){h+=hydro[a]; n++;}
  }
  if(n)h/=n;
  return(h);
}



float Print_DG_contfreq(struct REM *E, char *name)
{
  FILE *file_out=fopen(name, "w");
  printf("Writing %s\n", name);

  E->E24=(C242)*E->E1*E->E1;
  G_misfold(E); // Compute T freeze!
  float G_misf_freeze=E->E1-E->E2/E->Tf-E->Tf*E->S_C;
  if(E->REM==3)G_misf_freeze+=E->E3/(E->Tf*E->Tf);
  fprintf(file_out, "# Third moment of energy used: ");
  if(E->REM==3){fprintf(file_out,"YES\n");}
  else{fprintf(file_out,"NO\n");}
  fprintf(file_out,"# Tfreeze= %.2f\n", E->Tf);
  fprintf(file_out,
	  "# E_nat/L= %.4f Gfreeze/TL= %.4f E1/L= %.3g E2/L= %.2g E3/L= %.2g\n",
	  E->E_nat/E->L, G_misf_freeze/E->L,
	  E->E1/E->L, E->E2/E->L, E->E3/E->L);
  double nc1=0;
  for(int l=IJ_MIN; l<E->L; l++)nc1+=Cont_freq[l]*(E->L-l);
  fprintf(file_out, "# NC1/L= %.3f NC2/L= %.2g NC3/L= %.2g nc1/L= %.3f\n",
	  NC1/E->L, NC2/E->L, NC3/E->L, nc1/E->L);
  fprintf(file_out, "# L= %d S_U/L= %.3f S_C/L= %.3f\n",
	  E->L, E->S_U/E->L, E->S_C/E->L);
  fprintf(file_out, "# Temp. DeltaG  G_misf-G_unf\n");
  float T_INI=0.2, T_END=1.2, T_STEP=0.01;
  float T0=-1, DG0=100000;
  float DeltaG, G_misf;
  for(float T=T_INI; T<T_END; T+=T_STEP){
    if(T<E->Tf){
      G_misf=G_misf_freeze; // Freezing
    }else{
      G_misf=E->E1-E->E2/T-T*E->S_C;
      if(E->REM==3)G_misf+=E->E3/(T*T);
    }
    DeltaG=E->E_nat-G_misf+T*log(1+exp(G_misf/T+E->S_U));
    if(DeltaG>-100)fprintf(file_out, "%.2f\t%.2f\t%.2f\n",
			   T*E->T, DeltaG, G_misf+T*E->S_U);
    if(DeltaG < DG0){DG0=DeltaG; T0=T*E->T;}
  }
  fclose(file_out);
  return(T0);
}

int **Get_target(char *file_str, char *name_tar, int *len_tar)
{
  FILE *file_in=fopen(file_str, "r");
  char string[200], name[20], dumm[4];
  int N_prot=0, length, n_cont;
  int i, res1, res2, nc[L_MAX];
  struct contact *cont_list; short **contact=NULL;

  if(file_in==NULL){
    printf("ERROR, file %s not found (structures)\n", file_str); exit(8);
  }

  while(fgets(string, sizeof(string), file_in)!=NULL){

    if(strncmp(string, "#  All", 6)==0)continue;
    sscanf(string,"%s%d%d%s", dumm,&length,&n_cont,name);
    if(string[0]!='#'){
      printf("No protein start symbol at %d. Read:\n%s",N_prot+1,string);
      exit(8);
    }
    N_prot++;

    /* Check if target */
    if(strncmp(name, name_tar, 4)==0){
      *len_tar=length;
      strcpy(target.name,name);
      target.n_cont=n_cont;
      if(target.contact){
	for(i=0; i<target.length; i++)free(target.contact[i]);
	free(target.contact);
      }
      contact=malloc((length)*sizeof(short *));
      target.contact=contact;
      for(i=0; i<length; i++){
	contact[i]=malloc(NCMAX*sizeof(short)); nc[i]=0;
      }
      if(target.cont_list)free(target.cont_list);
      cont_list=malloc((n_cont+1)*sizeof(struct contact));
      target.cont_list=cont_list;
      target.length=length;

      /* Read contact map */
      for(i=0; i<n_cont; i++){
	fgets(string, sizeof(string), file_in);
	sscanf(string,"%d%d",&res1,&res2);
	cont_list->res1=res1; cont_list->res2=res2; cont_list++;
	contact[res1][nc[res1]]=res2; nc[res1]++;
	//contact[res2][nc[res2]]=res1; nc[res2]++;
      }
      for(i=0; i<length; i++)contact[i][nc[i]]=-1;
      cont_list->res1=-1;
      break;
    }else{
      // Other structure; skipped
      for(i=0; i<n_cont; i++)fgets(string, sizeof(string), file_in);
    }
  }
  fclose(file_in);

  if(contact==NULL){
    printf("ERROR in Get_target, pdb %s not found in file %s, %d proteins\n",
	   name_tar, file_str, N_prot); exit(8);
  }
  int **C_nat=Fill_C_nat(length, contact);
  return(C_nat);
}

int **Fill_C_nat(int length, short **contact)
{
  int i, j, Nc=0;
  // Allocate native contact matrix
  int **C_nat=malloc(length*sizeof(int *));
  for(i=0; i<length; i++){
    C_nat[i]=malloc(length*sizeof(int));
    for(j=0; j<length; j++)C_nat[i][j]=0;
  }
  for(i=0; i<length; i++){
    short *Ci=contact[i];
    while(*Ci>=0){
      C_nat[i][*Ci]=1; C_nat[*Ci][i]=1; Ci++; Nc++;
    }
  }
  printf("Allocating native contact matrix, L=%d Nc=%d\n", length, Nc);
  return(C_nat);
}

void Test_contfreq(struct REM *E, short *aa_seq, char *name)
{
  // WARNING, you have to make sure that Initialize_REM has been run before
  double EC1, EC2, EC3, U2; int L=E->L;
  Threading_energies(&EC1, &EC2, &EC3, &U2, aa_seq, L);
  char s[5000];
  sprintf(s, "T= %.2f L= %d S_U=%.2f S_C=%.2f\n", E->T, E->L, E->S_U, E->S_C);
  sprintf(s, "%sE_nat/L = %.4f\n", s, E->E_nat/L);
  sprintf(s, "%sTHREADING results REM= %d\n", s, E->REM);
  sprintf(s, "%s<U^2>= %.3g\n", s, U2);
  sprintf(s, "%s<EC>= %.3g           /L= %8.3g\n", s, EC1, EC1/L);
  sprintf(s, "%s<(EC-<EC>)^2>/2= %6.2f  /L= %8.3g\n", s, EC2/2, EC2/(2*L));
  sprintf(s, "%s<(EC-<EC>)^3>/6= %6.2f  /L= %8.3g\n", s, EC3/6, EC3/(6*L));
  sprintf(s, "%sDG_misfold/L = %.4g\n", s,
	  (EC1-EC2/(2*E->T))/L); //+EC3/(6*E->T*E->T)

  E->E24=(C242)*E->E1*E->E1;
  G_misfold(E);

  U2=0; int i, j; double nc=0;
  for(i=0; i<L; i++){
    float *U=Econt_T[aa_seq[i]];
    for(j=i+IJ_MIN; j<L; j++){
      float Uij=U[aa_seq[j]];
      int l=j-i; if(l>=LEN_MAX)continue; //l=LEN_MAX-1;
      U2+=Cont_freq[l]*Uij*Uij; nc+=Cont_freq[l];
    }
  }

  sprintf(s, "%sCONTFREQ results:\n", s);
  sprintf(s, "%s<U^2>= %.3g\n", s, U2/nc);
  sprintf(s, "%sE1= %.3g\n", s, E->E1);
  sprintf(s, "%s<(EC-<EC>)^2>/2= %6.2f  E22= %.2g E23= %.2g E24= %.2g\n",
	 s, E->E2, E->E22, E->E23, E->E24);
  sprintf(s,"%s<(EC-<EC>)^3>/6= %6.2f  ", s, E->E3);
  sprintf(s,"%sE321= %.2g E332= %.2g\n", s, E->E321, E->E332);
  sprintf(s,"%s\tE342= %.2g E343_2= %.2g E343_3= %.2g",
	  s, E->E342, E->E343_2, E->E343_3);
  sprintf(s,"%s\tE353=%.2g E363= %.2g\n", s, E->E353, E->E363);
  sprintf(s,"%sE36 used: ", s);
  if(E_36){sprintf(s,"%sYES\n", s);}else{sprintf(s,"%sNO\n", s);}
  sprintf(s, "%sDG_misfold/L = %.4g\n", s,
	  (E->E1-E->E2/(2*E->T))/L); //+E->E3/(6*E->T*E->T));
  printf("%s", s);

  FILE *file_out=Output_file(name, "Threading", "dat");
  fprintf(file_out, "%s", s);
  fclose(file_out);

}

void Threading_energies(double *EC1, double *EC2, double *EC3, double *U2,
			short *aa_seq, int len_pdb)
{
  int len_dom, ndom;
  if(len_pdb>LEN_MAX){
    ndom=len_pdb/LEN_MAX;
    if(len_pdb > ndom*LEN_MAX)ndom++;
    len_dom=len_pdb/ndom;
  }else{
    len_dom=len_pdb; ndom=1;
  }

  //int N_PROT=Read_structures(file_str, prot);
  short j_prot=0, first;
  short res1, *res2;
  long n_str=0, Noncompact=0;

  int NCmin, NCmax;
  Contact_range(&NCmin, &NCmax, len_dom);
  *EC1=0; *EC2=0; *EC3=0; *U2=0;
 
  /* Generate alternative structures by threading */
  printf("### Contact statistics\n");
  printf("%d contact matrices found in file %s\n", N_PROT, FILE_STR);
  printf("Selecting fragments of %d residues with %d < NC < %d\n",
	 len_dom, NCmin, NCmax);
  printf("(NC/L ~ %.2f - %.2f^(-1/3) +/- %.2f)\n", C_ASYMPT,C_SLOPE,C_THR);
  double Nc_sum=0;

  for(j_prot=0; j_prot< N_PROT; j_prot++){

    short len2=prot[j_prot].length;
    if(len2 < len_dom)continue;
    short **contact=prot[j_prot].contact;

    for(first=0; first<=len2-len_dom; first++){
      
      /* New fragment */
      int last=first+len_dom, Nc=0;
      for(res1=first; res1<last; res1++){
	res2=contact[res1];
	while((*res2>=0)&&(*res2 < last)){
	  if(*res2 >= res1+IJ_MIN)Nc++;
	  res2++;
	}
      }
      // Exclude non compact conformations
      if((Nc<NCmin)||(Nc>NCmax)){Noncompact++; continue;}
      n_str++;

      // Align structure with each of the domains
      int ini=0; double e1=0, e2=0;
      for(int idom=0; idom<ndom; idom++){
	int shift=ini-first;
	for(res1=first; res1<last; res1++){
	  res2=contact[res1];
	  float *U=Econt_T[aa_seq[res1+shift]], Uij;
	  while((*res2>=0)&&(*res2 < last)){
	    if(*res2 >= res1+IJ_MIN){
	      Uij=U[aa_seq[*res2+shift]];
	      e1+=Uij; e2+=Uij*Uij;
	    }
	    res2++;
	  }
	}
	Nc_sum+=Nc; 
	ini+=len_dom;
      } // end sequence domain
      (*U2)+=e2;
      (*EC1)+=e1; e2=e1*e1; (*EC2)+=e2; (*EC3)+=e2*e1;
    } // end fragment 
  } // end protein

  /* End of statistics, compute averages */
  (*EC1)/=n_str; (*EC2)/=n_str; (*EC3)/=n_str;
  (*EC3)=(*EC3)-3*(*EC1)*(*EC2)+2*(*EC1)*(*EC1)*(*EC1);
  (*EC2)=(*EC2)-(*EC1)*(*EC1);
  (*U2)/=Nc_sum;
  // Normalize by the length of the sampled fragments
  float norm=(float)len_pdb/(ndom*len_dom);
  (*EC1)*=norm; (*EC2)*=norm; (*EC3)*=norm;

  return;
}

int Read_structures(char *file_prot, struct protein *prot)
{
  FILE *file_in=fopen(file_prot, "r");
  if(file_in==NULL){
    printf("ERROR, I do not find the file with contact matrices named %s\n",
	   file_prot); exit(8);
  }

  char string[200], name[20], dumm[4];
  int N_prot=0, length_old=10000, length, n_cont;
  int i, res1, res2;
  struct contact *cont_list; short **contact, nc[L_MAX];
  
  N_prot=0;

  fgets(string, sizeof(string), file_in);
  while(fgets(string, sizeof(string), file_in)!=NULL){

    sscanf(string,"%s%d%d%s", dumm,&length,&n_cont,name);
    if(string[0]!='#'){
      printf("No protein start symbol at %d\n",N_prot+1);
      printf("%s\n",string);
      exit(8);
    }
    if(length>length_old){
      printf("Bad ordered proteins: N_prot=%d, length=%d, previous=%d\n",
	     N_prot+1, length, length_old);
      exit(8);
    }

    length_old=length;
    cont_list=malloc((n_cont+1)*sizeof(struct contact));
    contact=malloc((length)*sizeof(short *));
    for(i=0; i<length; i++){
      contact[i]=malloc(NCMAX*sizeof(short)); nc[i]=0;
    }

    // Empty memory, if not empty
    if(prot[N_prot].cont_list)free(prot[N_prot].cont_list);
    if(prot[N_prot].contact){
      for(i=0; i<prot[N_prot].length; i++)free(prot[N_prot].contact[i]);
      free(prot[N_prot].contact);
    }

    prot[N_prot].length=length; prot[N_prot].n_cont=n_cont;
    prot[N_prot].contact=contact;
    prot[N_prot].cont_list=cont_list;
    strcpy(prot[N_prot].name,name);
    N_prot++;

    /* Read contact map */
    for(i=0; i<n_cont; i++){
      fgets(string, sizeof(string), file_in);
      sscanf(string,"%d%d",&res1,&res2);
      cont_list->res1=res1; cont_list->res2=res2; cont_list++;
      contact[res1][nc[res1]]=res2; nc[res1]++;
      contact[res2][nc[res2]]=res1; nc[res2]++;
    }
    for(i=0; i<length; i++)contact[i][nc[i]]=-1;
    cont_list->res1=-1;
    if(N_prot==N_PROT_MAX-1)break;
  }
  fclose(file_in);
  prot[N_prot].length=0;

  return(N_prot);
}

 void Store(struct state *ptr, struct protein *protein, int first,
	    float energy, int n_cont){
   ptr->prot_ptr=protein;
  ptr->first=first;
  ptr->energy=energy;
  ptr->n_cont=n_cont;
 }

void Contact_range(int *NCmin, int *NCmax, int L){
  float c=C_ASYMPT - C_SLOPE*pow((float)L, -0.3333333);
  int NC=c*L, thr=C_THR*L;
  (*NCmin)=NC-thr;
  (*NCmax)=NC+thr;
}

void Copy_E_REM(struct REM *E2, struct REM *E1)
{
  E2->E_nat=E1->E_nat;
  E2->E1=E1->E1;
  E2->E22=E1->E22;
  if(E1->REM==3){
    E2->E321=E1->E321;
    E2->E332=E1->E332;
    E2->E342=E1->E342;
  }
  for(int i=0; i<E1->L; i++){
    E2->c1U1[i]=E1->c1U1[i];
    E2->c332U2[i]=E1->c332U2[i];
  }
}

void Empty_E_REM(struct REM *E){
  free(E->c1U1); free(E->c332U2);
} 
