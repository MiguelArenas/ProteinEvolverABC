/*
  REM computation of the free energy of the misfolded ensemble
  E1 = sum_ij <Cij>Uij   <U>=E1/NC (weighted mean energy)
  E2 = sum_ijkl (<CijCkl>-<Cij><Ckl>)UijUkl
  c1U1(i)= sum_j <Cij>Uij

  E2=sum_ijkl (<CijCkl>/<Cij><Ckl>-1)<Cij>Uij <Ckl>Ukl
    =
     sum_ij (<Cij>-<Cij>^2) Uij^2 (ij=kl)
    +sum_ij sum_kl (<CijCkl>/<Ckl>-<Cij>)/L^2 *Uij *E1
    -sum_ij (1-<Cij>) Uij *E1
    +sum_ij sum_kl (<CikCjl>/<Cik><Cjl>-1)/L^2 *c1U1(i)*c1U1(j)
    =
     sum_ij (<Cij>-<Cij>^2) Uij^2
    -sum_ij (1-<Cij>) Uij *E1/L^2 (negligible)
    +[sum_ij (<Cij*Nc>/<Nc>-<Cij>) *Uij]*E1
    + sum_ij (<ni*nj>/<ni><nj>-1) *c1U1(i)*c1U1(j)

*/

static int LEN_MAX=250; // Max. l for REM calculations
int dj_max=0;    // Max value for loops d,d+l

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
#include "externals.h"

// Parameters for selecting compact structures
#define NAA 21 // Number of amino acid types
#define C_ASYMPT 4.00  // Asymptotic value of Nc/N
#define C_SLOPE 8.07   // Nc/N ~ C_ASYMPT - C_SLOPE N^(-1/3)
#define C_THR 1.0      // Maximum deviation of Nc/N from typical value
                       // Very important parameter!!!
int E_36=0; // WARNING!!! If E_36=1 and C_THR is large, there will be problems
int E343_2=0;
float T0=-1;
int Ini_VBS=1;

// Global variables
//int LEN=0;
//int DTAIL_MAX=30;
//int LEN_MAX=300;

#define L_MAX 100000
#define NCMAX 30    // Maximum number of contacts per residue

// Common parameters defined in .h files
//float **interactions;
struct protein prots[N_PROT_MAX], target; int N_PROT;

// NEW: Concact statistics
int DM1; // DTAIL_MAX+1;
// <Nc>, <Nc^2>, <Nc^3> as a function of x(L)=L^(-1/3)
double Slope_x; // <Nc/L>=Nc1L0*(1+Slope_x*x(L))
double Y_ave;  // Nc/L=y*(1+Slope_x*x(L)) Y_ave=<y>
// <(y-<y>)^2>~exp(Y2L0+Y2L1*log(L)), corr=r2
// <(y-<y>)^3>~exp(Y3L0+Y3L1*log(L)), corr=r3
double z_expo;  // z_expo=(r2*Y2L1/2+r3*Y3L1/3)/(r2+r3)
double z2_ave, z3_ave, zs2, zs3;
//  <(y-<y>)^n> ~ L^(n*z_expo) => (y-<y>) ~ z*L^z_expo
//  z=(y-<y>)*L^(-z_expo) 
// z2_ave= <z^2> z3_ave= <z^3>

// Length-independent variables
//float *Cnc1=NULL;  // <C_ij(nc_i+nc_j)/2>
float *Cont_f=NULL;  // <C_ij>=y*Cont_f(|i-j|)
float *Cz1=NULL;     // <Cz>
float *Cz2=NULL;     // <Cz^2>
float *nc_ave=NULL;  // <nc[i]>/(1+Slope_x*x(L))
float *nc2_ave=NULL;
float *nc_z=NULL;    // <nc[i]z>
float *nc_z2=NULL;   // <nc[i]z*z>
float **nc_nc=NULL;  // <nc[i]nc[i+l]>/(1+Slope_x*x(L))^2
float **nc_nc_z=NULL; // <nc[i]nc[i+l]Nc>/(1+Slope_x*x(L))^3

// fit parameters
float Cf_0=0, Cf_1=0; 
float Cz1_0=0, Cz1_1=0;
float Cz2_0=0, Cz2_1=0;
float *nc_nc_0=NULL, *nc_nc_1=NULL;
float *nc_nc_z_0=NULL, *nc_nc_z_1=NULL;
float nc2L_0, nc2L_1;
float nc1L_0, nc1L_1;
double Slope_x2;
float *Corr_nc_nc=NULL;
float *Corr_nc_nc_z=NULL;

//float Cnc1_0=0, Cnc1_1=0;
// sum_i <(nci)^2>/L~(nc2L0+nc2L1*L^(-1/3))
// nc2~ sum_i (nc[i]/(1+Slope_x*x(L))^2 ~ nc2L0+nc2L1*x[L]
//
// Length rescaled variables
// Nc/L=(<y>+zL^(-z_expo))(1+Slope_x*L^(-1/3)) <z>=1
float NC1=0, NC2=0, NC3=0; // moments of the number of contacts
float *Cont_L=NULL;
float *C2_L=NULL;
float *C3_L=NULL;
float *nc_L=NULL;
float *nc2_L=NULL;
float *nc_Nc_L=NULL;
float *nc_Nc2_L=NULL;
float **nc_nc_L=NULL;
float **nc_nc_Nc_L=NULL;
float *Cont_1L=NULL;
float *Cont_2L=NULL;
float *Cont_3L=NULL;
//float *Cnc1_L=NULL;

// Obsolete
/*double *C221_ij=NULL, *C232_i=NULL, C242=0;
double *Cont_freq=NULL;
double *C321_ij=NULL, **C332_ij=NULL, *C342_ij=NULL;
double **C343_2_ij=NULL, *C343_3=NULL, *C353_i=NULL, C363=0;
double *nc1i=NULL, *nc2i=NULL, *nc3i=NULL, *CNc=NULL;
double *NC1nc1i=NULL, *NC2nc1i=NULL, *NC1nc2i=NULL;
double **Cnn=NULL, **ncnc=NULL, **C1n;
double nc2_sum, nc12_sum; */

int init_REM=1;
int Initialize_REM(char *file_str, char *nameout);
void Rescale_REM(int L_pdb, char *nameout);
static float Corr_coeff(float *off, float *slope, float *xx, float *yy, int n);
void Contact_range(int *NCmin, int *NCmax, int L);
void Compute_E343_2(struct REM *E, short *aa_seq);
int Read_structures(char *file_prot, struct protein *prot);
void Threading_energies(double *E_ave,
			double *EC1, double *EC2, double *EC3, double *U2,
			double *Nc1, double *Nc2, double *Nc3, 
			short *aa_seq, float **C_nat, int len_nat,
			float *T, int nT);
float G_misfold_threading(struct REM *E);
int tail(int i, int L);
void Count_sec(int *N_helix, int *N_strand, char *c_sec, int L);

/*****************************************************************/

void Initialize_E_REM(struct REM *E, int L, int REM,
		      float T, float S_C, float S_U,
		      char *file_str)
{
  // Temperature
  if(Econt_T==NULL){
    Econt_T=Allocate_mat2_f(NAA, NAA);
    printf("Dividing contact energies by T= %.2f\n", T);
    T0=T; int i, j;
    for(i=0; i<NAA; i++){
      AA_code[i]=AA_BKV[i];
      hydro[i]=hspec[i];
      for(j=0; j<NAA; j++)Econt_T[i][j]=EBKV[i][j]/T;
    }
    Initialize_E_loc(T, SEC_STR);
  }

  // Statistics of misfolded structures
  char nameout[100]="REM.txt";
  if(init_REM){
    N_PROT=Initialize_REM(file_str, nameout); init_REM=0;
  }
  if(L != LEN){
    printf("Rescaling REM\n");
    Rescale_REM(L, nameout); LEN=L; 
  }

  E->L=L;
  //if(E->c1U1){free(E->c1U1); E->c1U1=NULL;}
  E->c1U1=malloc(L*sizeof(double));
  E->REM=REM;
  E->T=T;
  E->S_C=S_C;
  E->S_U=S_U;
  //E->c332U2=malloc(L*sizeof(double));
  //E->S_C=sC0+sC1*L;
  //E->S_U=sU1*L;
  Set_to_zero(E);
}

void Set_to_zero(struct REM *E)
{
  E->E_nat=0;
  E->E_loc=0;
  E->E1=0;
  E->E2=0;
  E->E3=0;
  E->E2cont1=0;
  E->E2cont2=0;
  E->E2site1=0;
  E->E2site2=0;
  E->E3cont1=0;
  E->E3cont2=0;
  E->E3cont3=0;
  E->E3site1=0;
  E->E3site2=0;
  for(int i=0; i<E->L; i++)E->c1U1[i]=0;
  /*E->E321=0; E->E332=0;
  E->E342=0; E->E343_2=0; E->E343_3=0;
  E->E353=0; E->E363=0;
  E->e342=0;
  for(int i=0; i<E->L; i++)E->c332U2[i]=0;*/
}


int Initialize_REM(char *file_str, char *nameout)
{
  // Allocate
  int L_max=LEN_MAX, i;
  DM1=DTAIL_MAX+1;
  Cont_f =malloc(L_max*sizeof(float));
  Cz1  =malloc(L_max*sizeof(float));
  Cz2  =malloc(L_max*sizeof(float));
  //Cnc1  =malloc(L_max*sizeof(float));
  for(i=0; i<L_max; i++){Cont_f[i]=0; Cz1[i]=0; Cz2[i]=0;}
  nc_z=malloc(L_max*sizeof(float));
  nc_z2=malloc(L_max*sizeof(float));
  nc_ave=malloc(L_max*sizeof(float));
  nc2_ave=malloc(L_max*sizeof(float));
  nc_nc=malloc(DM1*sizeof(float *));
  nc_nc_z=malloc(DM1*sizeof(float *));
  for(i=0; i<=DTAIL_MAX; i++){
    nc_nc[i] =malloc(L_max*sizeof(float));
    nc_nc_z[i] =malloc(L_max*sizeof(float));
  }
  nc_nc_0=malloc(DM1*sizeof(float));
  nc_nc_1=malloc(DM1*sizeof(float));
  Corr_nc_nc=malloc(DM1*sizeof(float));
  nc_nc_z_0=malloc(DM1*sizeof(float));
  nc_nc_z_1=malloc(DM1*sizeof(float));
  Corr_nc_nc_z=malloc(DM1*sizeof(float));


  // Read proteins and prepare output
  int N_prot=Read_structures(file_str, prots);

  FILE *file_out=fopen(nameout, "w");
  printf("Writing contact statistics in file %s\n", nameout);

  char out[1000];
  sprintf(out,
	  "### Contact statistics\n"
	  "#%d contact matrices found in file %s\n", N_prot, file_str);
  printf("%s", out); fprintf(file_out, "%s", out);

  // Scaling of number of contacts
  printf("Fit of number of contacts:\n");
  float Slope_tmp=C_SLOPE/C_ASYMPT; 
  float Ncp[N_prot], Xp[N_prot], logLp[N_prot];
  int np=0, select[N_prot];
  for(i=0; i< N_prot; i++){
    struct protein *prot=prots+i;
    float L=prot->length;
    float x=pow(L, -0.333333), c=(prot->n_cont/L);
    if(fabs(c-C_ASYMPT*(1-Slope_tmp*x))>C_THR){
      select[i]=0;
    }else{
      select[i]=1;
      logLp[np]=log(L); Xp[np]=x; Ncp[np]=c;
      np++;
    }
    //printf("Prot %d: %.0f %.3f %.3f %d\n", i, L, c, x, select[i]);
  }
  printf("%d proteins stored\n", np);

  // Fits of Nc(L), Nc^2(L) etc.
  sprintf(out,
	  "# Fitting number of contacts versus surface to volume scaling.\n"
	  "# %d proteins kept and %d discarded due to large deviations "
	  "|Nc/L-%.2f(1-%.2fL^(-1/3)|>%.2f\n",
	  np, N_prot-np, C_ASYMPT,Slope_tmp,C_THR);
  printf("%s", out); fprintf(file_out, "%s", out);

  // First moment
  float Nc1L0, Nc1L1;    // <Nc>~L(NC1L0+NC1L1*L^(-1/3))
  float r=Corr_coeff(&Nc1L0, &Nc1L1, Xp, Ncp, np);
  Slope_x=Nc1L1/Nc1L0; 
  sprintf(out,

	  "# <Nc/L> ~ %.2f + %.1fL^(-1/3)  r= %.3f\n", Nc1L0, Nc1L1, r);
  printf("%s", out); fprintf(file_out, "%s", out);

  // Normalize Nc/L by fitted value
  sprintf(out, "# y = (NC/L)/(1 + %.2fL^(-1/3))\n", Slope_x);
  printf("%s", out); fprintf(file_out, "%s", out);
  float Yp[np]; Y_ave=0;
  for(i=0; i< np; i++){
    Yp[i]=Ncp[i]/(1+Slope_x*Xp[i]);
    Y_ave+=Yp[i];
  } 
  Y_ave/=np;
  float off, slope;
  r=Corr_coeff(&off, &slope, Xp, Yp, np);
  sprintf(out,"# Slope= %.3f  <y>= %.3f\n"
	  "# Fit: y ~ %.3f + %.1f*L^(-1/3)  r= %.3f\n",
	  Slope_x, Y_ave,off,slope,r);
  printf("%s", out); fprintf(file_out, "%s", out);

  // Second and third moment with power law
  float logYp[np];
  for(i=0; i< np; i++)logYp[i]=log(fabs(Yp[i]-Y_ave));
  float YL0, YL1, r2;  // <(y-<y>)^2>~exp(Y2L0+Y2L1*log(L))
  r2=Corr_coeff(&YL0, &YL1, logLp, logYp, np);
  z_expo=YL1; if(fabs(r2)<0.5)z_expo=0;
  // z_expo=-0.40 instead of -0.48 gives better <NC^2>
  sprintf(out,
	  "# <(y-<y>)^2> ~ (%.3f*L^(%.3f))^2  r= %.3f\n"
	  "# Scaling exponent: %.3f\n",
	  exp(YL0), YL1, r2, z_expo);
  printf("%s", out); fprintf(file_out, "%s", out);

  // Scale-less variable Z 
  float Zp[np]; double z1_ave=0; z2_ave=0; z3_ave=0; 
  for(i=0; i< np; i++){
    float z=(Yp[i]-Y_ave)*exp(-z_expo*logLp[i]);
    Zp[i]=z; z1_ave+=z; z2_ave+=z*z; z3_ave+=z*z*z;
  }
  z1_ave/=np; z2_ave/=np; z3_ave/=np;
  r=Corr_coeff(&off, &slope, Xp, Zp, np);
  sprintf(out,
	  "# z=(y-<y>)*L^%.2f ~ %.2f + %.2fL^(-1/3)  r= %.3f\n"
	  "# <z>= %.2g <z^2>^(1/2)= %.3g <z^3>^(1/3)= %.3g <z^3>= %.3g\n",
	  -z_expo,off,slope,r,z1_ave,sqrt(z2_ave),
	  pow(fabs(z3_ave),1/3)*(z3_ave/fabs(z3_ave)), z3_ave);
  printf("%s", out); fprintf(file_out, "%s", out);

  /**************** Statistics of contacts ****************/
  sprintf(out,"# Statistics of contacts, %d proteins L<%d\n",np, L_max);
  printf("%s", out); fprintf(file_out, "%s", out);

  long Cont_num[L_max], Cont_norm[L_max],   //ncCont_num[L_max], 
    z_Cont_num[L_max], z2_Cont_num[L_max];
  double num_nc_nc[DM1][L_max], nc_nc_sum[DM1][L_max], 
    nc_nc_z_sum[DM1][L_max];
  double num_nc[L_max], nc_sum[L_max], nc2_sum[L_max],
    nc_z_sum[L_max], nc_z2_sum[L_max];
  for(i=0; i<L_max; i++){
    Cont_num[i]=0; Cont_norm[i]=0; //ncCont_num[i]=0; 
    z_Cont_num[i]=0; z2_Cont_num[i]=0;
    num_nc[i]=0;
    nc_sum[i]=0;
    nc2_sum[i]=0;
    nc_z_sum[i]=0;
    nc_z2_sum[i]=0;
  }
  for(int j=0; j<DM1; j++){
    for(i=0; i<L_max; i++){
      num_nc_nc[j][i]=0;
      nc_nc_sum[j][i]=0;
      nc_nc_z_sum[j][i]=0;
    }
  }

  float nc2[N_prot]; int kp=0;
  for(int ip=0; ip< N_prot; ip++){
    if(select[ip]==0)continue;
    struct protein *prot=prots+ip;
    int L=prot->length;
    int nc[L], ncont_l[L];
    for(i=0; i<L; i++){nc[i]=0; ncont_l[i]=0;}
    struct contact *cont=prot->cont_list;
    for(int ic=0; ic<prot->n_cont; ic++){
      int l=cont->res2-cont->res1;
      if(l>=IJ_MIN){
	nc[cont->res1]++; nc[cont->res2]++;
	if(l<L_max)ncont_l[l]++;
      }
      cont++;
    }
    /*// <c_ij ni>
    for(int ic=0; ic<prot->n_cont; ic++){
      struct contact *cont=prot->cont_list+ic;
      int l=cont->res2-cont->res1; if(l<IJ_MIN)continue;
      ncCont_num[l]+=(nc[cont->res1]+nc[cont->res2]);
      }*/
    // Sum protein to global counters
    float z=Zp[kp];
    for(i=0; i<L; i++){
      if(i==L_max)break;
      Cont_norm[i]+=(L-i);
      Cont_num[i]+=ncont_l[i];
      z_Cont_num[i]+=ncont_l[i]*z;  // Average contacts per residue
      z2_Cont_num[i]+=ncont_l[i]*z*z;
    }
    
    float core=(1+Slope_x*Xp[kp]); // Ycore=Yp[kp]*core;
    float ci[L];
    for(i=0; i<L; i++){
      ci[i]=nc[i]/core;
      // distance from tail
      int di=tail(i, L);
      num_nc[di]++;
      nc_sum[di]+=ci[i];
      nc2_sum[di]+=ci[i]*ci[i];
      float cz=ci[i]*z;
      nc_z_sum[di]+=cz;
      nc_z2_sum[di]+=cz*z;
      for(int j=i; j>=0; j--){
	int l=i-j;
	if(l>=L_max)break;
	int dj=tail(j, L), d;
	if(di<=dj){d=di;}else{d=dj;}
	if(d>DTAIL_MAX){d=DTAIL_MAX;}
	nc_nc_sum[d][l]+=ci[i]*ci[j];
	nc_nc_z_sum[d][l]+=cz*ci[j];
	num_nc_nc[d][l]++;
      }
    } // End pairs
    // sum_i (ni/L)^2
    float n2=0; for(i=0; i<L; i++)n2+=nc[i]*nc[i];
    nc2[kp]=n2/L; kp++;
  } // end proteins
  if(kp!=np){
    printf("ERROR, wrong number of proteins %d instead of %d\n",
	   kp, np); exit(8);
  }
  
  // sum_i <(ni/L)^2>~(nc2L_0+nc2L_1*L^(-1/3))
  r=Corr_coeff(&nc2L_0, &nc2L_1, Xp, nc2, np);
  sprintf(out, "# sum_i nci^2/L ~ %.3g + %.3gL^(-1/3)  r= %.3f\n",
	  nc2L_0, nc2L_1, r);
  printf("%s", out); fprintf(file_out, "%s", out);
  for(i=0; i<np; i++)nc2[i]=sqrt(nc2[i]);
  r=Corr_coeff(&nc1L_0, &nc1L_1, Xp, nc2, np);
  Slope_x2=nc1L_1/nc1L_0;
  sprintf(out,"# sqrt(sum_i nci^2/L) ~ %.3g + %.3gL^(-1/3)  r= %.3f\n"
	  "# Slope2= %.3g/%.3g= %.3g\n",
	  nc1L_0, nc1L_1, r, nc1L_1, nc1L_0, Slope_x2);
  printf("%s", out); fprintf(file_out, "%s", out);

 
  for(i=0; i<L_max; i++){
    float norm=Cont_norm[i];
    //Cnc1[i]=(float)ncCont_num[i]/(2*norm);
    Cont_f[i]=(float)Cont_num[i]/norm;
    Cz1[i]=(float)z_Cont_num[i]/norm;
    Cz2[i]=(float)z2_Cont_num[i]/norm;
    //if((i>=4)&&(i<30))printf("%.4f %.3f\n", Cont_f[i], CNc1[i]);
    if((i>=4)&&(Cont_f[i]<0)){
      printf("ERROR, negative Cont_freq:\n");
      printf("%d cf= %.4f %ld %ld Ncf= %.4f Nc2f= %.4f\n",
	     i, Cont_f[i], Cont_num[i], Cont_norm[i], Cz1[i], Cz2[i]);
      exit(8);
    }
  }

  dj_max=L_max-1;
  for(int d=0; d<L_max; d++){
    if(num_nc[d]<100){dj_max=d-1; break;}
    nc_ave[d]=nc_sum[d]/num_nc[d];
    nc2_ave[d]=nc2_sum[d]/num_nc[d];
    nc_z[d]=nc_z_sum[d]/num_nc[d]; 
    nc_z2[d]= nc_z2_sum[d]/num_nc[d];
  }
  for(int d=0; d<=DTAIL_MAX; d++){
    for(int l=0; l<L_max; l++){
      int dj=d+l;
      if(dj<=dj_max && num_nc_nc[d][l]>=100){
	nc_nc[d][l]=nc_nc_sum[d][l]/num_nc_nc[d][l];
	nc_nc_z[d][l]=nc_nc_z_sum[d][l]/num_nc_nc[d][l];
	// Center the correlation
	nc_nc[d][l]-=nc_ave[d]*nc_ave[dj];
	nc_nc_z[d][l]-=(nc_ave[d]*nc_z[dj]+nc_ave[dj]*nc_z[d])/2;
      }else{
	if(dj<dj_max)dj_max=dj;
	nc_nc[d][l]=0;
	nc_nc_z[d][l]=0;
      }
    }
  }
  printf("Max. value of d+l= %d\n", dj_max);

  // Fit for larger distance
  int L2=L_max/2, n=0;
  for(i=L2; i<L_max; i++){
    if(Cont_norm[i]<100){break;} n++;
  }
  int m=n;

  float xx[L_max], x[L_max], yy[L_max]; n=0;
  for(i=L2; i<(L2+m); i++){
    xx[n]=log(i); yy[n]=log(Cont_f[i]); n++;
  }
  r=Corr_coeff(&off, &slope, xx, yy, n);
  Cf_0=off; Cf_1=slope;
  sprintf(out,"# <C(l)> ~ %.2g*l^(%.2f)  r= %.3f\n",exp(off),slope,r);
  printf("%s", out); fprintf(file_out, "%s", out);
  

  //n=0; for(i=L2; i<(L2+m); i++){yy[n]=log(Cnc1[i]); n++;}
  //r=Corr_coeff(&off, &slope, xx, yy, n);
  //Cnc1_0=off; Cnc1_1=slope;
  //sprintf(out,"# <nci*C(i+l)> ~ %.2f*l^(%.2f) r= %.3f\n",exp(off),slope,r);
  //printf("%s", out); fprintf(file_out, "%s", out);

  n=0; for(i=L2; i<(L2+m); i++){yy[n]=Cz1[i]; n++;}
  r=Corr_coeff(&off, &slope, xx, yy, n);
  if(fabs(r)>0.25){Cz1_0=off; Cz1_1=slope;}
  else{Cz1_1=0; Cz1_0=Cz1[L2+n-1];}
  sprintf(out,"# <C(l)*z(Nc/L)> ~ %.3g+%.3g*log(l)  r= %.3f"
	  " Fit para: %.2g %.2g\n", off, slope,r, Cz1_0, Cz1_1);
  printf("%s", out); fprintf(file_out, "%s", out);

  n=0; for(i=L2; i<(L2+m); i++){yy[n]=Cz2[i]; n++;}
  r=Corr_coeff(&off, &slope, xx, yy, n);
  if(fabs(r)>0.25){Cz2_0=off; Cz2_1=slope;}
  else{Cz2_1=0; Cz2_0=Cz2[L2+n-1];}
  sprintf(out,"# <C(l)*z(Nc/L)^2> ~ %.3g+%.3g*log(l)  r= %.3f"
	  " Fit para: %.2g %.2g\n", off, slope, r, Cz2_0, Cz2_1);
  printf("%s", out); fprintf(file_out, "%s", out);

  for(int d=0; d<=DTAIL_MAX; d++){
    n=0; int nan=0;
    for(i=10; i<dj_max; i++){
      if(isnan(nc_nc[d][i])){nan=1;}
      else if(nc_nc[d][i]){
	yy[n]=(nc_nc[d][i]/(nc_ave[d]*nc_ave[i])); x[n]=i;
	//if(yy[n]>0){yy[n]=log(yy[n]);}else{continue;} x[n]=log(i);
	n++;
      }
    }
    r=Corr_coeff(&off, &slope, x, yy, n);
    if(isnan(r) || nan){
      printf("ERROR, fit of nc_nc[%d] is nan\n", d); 
      for(i=0; i<dj_max; i++){
	if(isnan(nc_nc[d][i]))printf("%d %d %.2g\n", d, i, nc_nc[d][i]);
      }
    }
    Corr_nc_nc[d]=r;
    if(fabs(r)>0.25){nc_nc_0[d]=off; nc_nc_1[d]=slope;}
    else{nc_nc_0[d]=nc_nc[d][L2+n-1]; nc_nc_1[d]=0;}
    if(d<15 || d==DTAIL_MAX){
      sprintf(out,"# <nc[%d]nc[%d+l]>/<nc[%d]><nc[%d+l]>-1 "
	      //"~ (%.3g+%.3g*l)  r= %.3f n=%d\n",
	      "~ exp(%.3g+%.3g*log(l))  r= %.3f n=%d\n",
	      d, d, d, d, off, slope, r, n);
      printf("%s", out); fprintf(file_out, "%s", out);
    }
  }

  //
  for(int d=0; d<=DTAIL_MAX; d++){
    n=0; int nan=0;
    for(i=10; i<dj_max; i++){
      if(isnan(nc_nc_z[d][i])){nan=1;}
      else if(nc_nc_z[d][i]){
	yy[n]=(nc_nc_z[d][i]); x[n]=i; n++;
      }
    }
    r=Corr_coeff(&off, &slope, x, yy, n);
    if(isnan(r)|| nan){
      printf("ERROR, fit of nc_nc_z[%d] is nan\n", d); 
      for(i=0; i<dj_max; i++){
	if(isnan(nc_nc_z[d][i]))printf("%d %d %.2g\n", d, i, nc_nc_z[d][i]);
      }
    }
    Corr_nc_nc_z[d]=r;
    if(fabs(r)>0.25){nc_nc_z_0[d]=off; nc_nc_z_1[d]=slope;}
    else{nc_nc_z_1[d]=0; nc_nc_0[d]=nc_nc_z[d][L2+n-1];}
    if(d<15 || d==DTAIL_MAX){
      sprintf(out,"# <nc[%d]nc[%d+l]z>-<nc><nc*z> "
	      "~ (%.3g+%.3g*l)  r= %.3f n=%d\n",
	      d, d, off,slope,r,n);
      printf("%s", out); fprintf(file_out, "%s", out);
    }
  }

  fprintf(file_out,
	  "#l=|i-j| <C> <Cz> <Cz^2> <n0*nl>/<n0><nl> <ni*nj>/<ni><nl> sum\n");
  for(i=0; i<50; i++){
    int k=i; if(k>dj_max)k=dj_max;
    int k2=i+DTAIL_MAX; if(k2>dj_max)k2=dj_max;
    fprintf(file_out, "%d\t%.4g\t%.3g\t%.3g\t%.3g\t%.3g\t%ld\n",
	    i, Cont_f[i], Cz1[i], Cz2[i],
	    nc_nc[0][i]/(nc_ave[0]*nc_ave[k]),
	    nc_nc[DTAIL_MAX][i]/(nc_ave[DTAIL_MAX]*nc_ave[k2]),
	    Cont_norm[i]);
  }
  fclose(file_out);

  printf("Contact statistics written in %s\n", nameout);
  return(N_prot);
}

void Rescale_REM(int L, char *nameout)
/*void Normalize_cont_freq(float *Cont_L, float *Cont_f,
			 float *Cnc1_L, float *Cnc1,
			 float *CNc1_L, float *Cz1,
			 float *CNc2_L, float *Cz2,
			 float **nc_nc_L, float *nc_nc,
			 int L)*/
{
  printf("Rescaling contacts for L= %d\n", L);

  double x=pow(L,-0.33333333); // x(L)
  float core=(1+Slope_x*x);
  float Lcore=L*core;
  float Lc2=Lcore*Lcore;

  float scale=pow(L, z_expo); // scale=1
  float sc2=scale*scale;
  zs2=z2_ave*sc2; zs3=z3_ave*sc2*scale;

  double cc1=Y_ave;     // <(Nc/Lcore)>
  double cc12=cc1*cc1;   // <(Nc/Lcore)>^2
  NC1=Y_ave*Lcore;     // <Nc>
  NC2=zs2*Lc2;       // <(Nc-<Nc>)^2>
  NC3=zs3*Lc2*Lcore; // <(Nc-<Nc>)^3>

  int l, d;
  //if(Cnc1_L)free(Cnc1_L);
  if(Cont_L)free(Cont_L);
  if(C2_L)free(C2_L);
  if(C3_L)free(C3_L);
  if(Cont_1L)free(Cont_1L);
  if(Cont_2L)free(Cont_2L);
  if(Cont_3L)free(Cont_3L);
  if(nc_L)free(nc_L);
  if(nc_Nc_L)free(nc_Nc_L);
  if(nc_Nc2_L)free(nc_Nc2_L);
  if(nc_nc_L){
    for(d=0; d<=DTAIL_MAX; d++)free(nc_nc_L[d]);
    free(nc_nc_L);
  }
  if(nc_nc_Nc_L){
    for(d=0; d<=DTAIL_MAX; d++)free(nc_nc_Nc_L[d]);
    free(nc_nc_Nc_L);
  }

  // Allocate
  //Cnc1_L=malloc(L*sizeof(float));
  Cont_L=malloc(L*sizeof(float));
  Cont_1L=malloc(L*sizeof(float));
  Cont_2L=malloc(L*sizeof(float));
  Cont_3L=malloc(L*sizeof(float));
  C2_L=malloc(L*sizeof(float));
  C3_L=malloc(L*sizeof(float)); 
  nc_L=malloc(DM1*sizeof(float));
  nc2_L=malloc(DM1*sizeof(float));
  nc_Nc_L=malloc(DM1*sizeof(float));
  nc_Nc2_L=malloc(DM1*sizeof(float));
  nc_nc_L =(float **)malloc(DM1*sizeof(float *));
  nc_nc_Nc_L =(float **)malloc(DM1*sizeof(float *));
  for(d=0; d<=DTAIL_MAX; d++){
    nc_nc_L[d] =malloc(L*sizeof(float));
    nc_nc_Nc_L[d] =malloc(L*sizeof(float));
  }

  // Store length rescaled variables
  for(l=0; l<LEN_MAX; l++){
    if(l>=L)break;
    Cont_L[l]=Cont_f[l];
    C2_L[l]=Cz1[l];
    C3_L[l]=Cz2[l];
    for(d=0; d<=DTAIL_MAX; d++){
      nc_nc_L[d][l]=nc_nc[d][l];
      nc_nc_Nc_L[d][l]=nc_nc_z[d][l];
    }
  }
  if(L>LEN_MAX){
    printf("Performing fits for l>=%d\n", LEN_MAX);
    double log_l[L]; for(l=LEN_MAX; l<L; l++)log_l[l]=log(l);
    for(l=LEN_MAX; l<L; l++){
      Cont_L[l]= exp(Cf_0+Cf_1*log_l[l]);
      C2_L[l]=Cz1_0+Cz1_1*log_l[l];
      C3_L[l]=Cz2_0+Cz2_1*log_l[l];
      for(d=0; d<=DTAIL_MAX; d++){
	//nc_nc_L[d][l]=exp(nc_nc_0[d]+nc_nc_1[d]*log_l[l]);
	nc_nc_L[d][l]=nc_ave[d]*nc_ave[LEN_MAX-1]*(nc_nc_0[d]+nc_nc_1[d]*l);
	//nc_nc_Nc_L[d][l]=exp(nc_nc_z_0[d]+nc_nc_z_1[d]*log_l[l]);
	nc_nc_L[d][l]=nc_nc_z_0[d]+nc_nc_z_1[d]*l;
      }
    }
  }

  // Normalization
  double sum_c=0, sum_C2=0, sum_C3=0, sum_nl=0; //sum_c_nc=0, 
  double sum_nc_nc[DTAIL_MAX+1], sum_nc_nc_Nc[DTAIL_MAX+1];
  for(d=0; d<=DTAIL_MAX; d++){
    sum_nc_nc[d]=0; sum_nc_nc_Nc[d]=0;
  }
  for(l=0; l<L; l++){
    int nl=L-l;
    sum_nl+=nl;
    sum_c+=Cont_L[l]*nl;
    sum_C2 +=C2_L[l]*nl;
    sum_C3 +=C3_L[l]*nl;
    //sum_c_nc+=Cnc1_L[l]*nl;
    for(d=0; d<=DTAIL_MAX; d++){
      int nld;
      if(d==DTAIL_MAX){nld=nl-2*DTAIL_MAX; if(nld<0)continue;}
      else if(d>l){nld=2;}
      else{nld=1;}
      sum_nc_nc[d]+=nc_nc_L[d][l]*nld;
      sum_nc_nc_Nc[d]+=nc_nc_Nc_L[d][l]*nld;
    }
  }

  float norm_c=   NC1/sum_c; // sum_l C(l)=Nc1
  float norm_C2= zs2*Lcore/(cc1*sum_C2);
  float norm_C3= zs3*Lcore/(cc12*sum_C3);
  float norm_nc_nc[DTAIL_MAX+1];
  float norm_nc_nc_Nc[DTAIL_MAX+1];
  float c_scale=cc1*scale;
  for(d=0; d<=DTAIL_MAX; d++){
    nc_L[d]=nc_ave[d]*core;
    nc2_L[d]=nc2_ave[d]*core;
    nc_Nc_L[d]=(nc_ave[d]*cc1+nc_z[d]*scale)*core*Lcore;
    norm_nc_nc[d]=2*nc_z[d]*scale*L/(sum_nc_nc[d]*nc_ave[d]); // 
    nc_Nc2_L[d]=
      (nc_ave[d]*cc12+2*c_scale*nc_z[d]+sc2*nc_z2[d])*core*Lc2;
    norm_nc_nc_Nc[d]=2*nc_z2[d]*sc2*Lcore/
      (sum_nc_nc_Nc[d]*nc_ave[d]*cc1);  // *cc1 
  }
  for(l=0; l<L; l++){
    Cont_L[l]*=norm_c;
    if(Cont_L[l]){
      C2_L[l]=(C2_L[l]*norm_C2)/Cont_L[l];
      C3_L[l]=(C3_L[l]*norm_C3)/Cont_L[l]; //-2*C2_L[l];
      Cont_1L[l]=1-Cont_L[l];
      Cont_2L[l]=C2_L[l]*(1-2*Cont_L[l]);
      Cont_3L[l]=1-3*Cont_L[l]+2*Cont_L[l]*Cont_L[l];
    }else{
      C2_L[l]=0; C3_L[l]=0;
      Cont_1L[l]=0; Cont_2L[l]=0; Cont_3L[l]=0;
    }
    for(d=0; d<=DTAIL_MAX; d++){
      int j=d+l; if(j>DTAIL_MAX)j=DTAIL_MAX;
      nc_nc_L[d][l]*=(norm_nc_nc[d]/nc_ave[j]); 
      nc_nc_Nc_L[d][l]*=(norm_nc_nc_Nc[d]/nc_ave[j]);
    }
  }

  printf("Writing %s\n", nameout);
  FILE *file_out=fopen(nameout, "a");

  double NC1_2= NC1*NC1, NC2_all=NC2+NC1_2;
  fprintf(file_out, "# Contact statistics for proteins of length %d\n", L);
  fprintf(file_out, "# <NC>= %.3g /L = %.3g\n", NC1, NC1/L);
  fprintf(file_out, "# <NC^2>-<NC>^2= %.3g  /<NC>= %.3g\n", NC2, NC2/NC1);
  fprintf(file_out, "# <(NC-<NC>)^3>= %.3g  /<NC>= %.3g\n", NC3, NC3/NC1);
  fprintf(file_out, "# [<NC^2>/<NC>^2-1]= %.2g  [<(NC/<NC>-1)^3>]= %.2g\n",
	  NC2/(NC1_2), NC3/(NC1_2*NC1));

  fprintf(file_out, "#i\t<nc_i>\t<nc_i*Nc>/<Nc>\t<nc_i*Nc^2>/<Nc^2>"
	  "\tr(nc_i*nc_i+l,log(l))\toffs\tslope"
	  "\tr(nc_i*nc_i+l*z,log(l))\toffs\tslope\n");
  for(d=0; d<=DTAIL_MAX; d++){
    fprintf(file_out, "%d\t%.3g\t%.3g\t%.3g"
	    "\t%.3g\t%.3g\t%.3g" "\t%.3g\t%.3g\t%.3g\n",
	    d, nc_L[d], nc_Nc_L[d]/NC1, nc_Nc2_L[d]/(NC2_all),
	    Corr_nc_nc[d], nc_nc_0[d], nc_nc_1[d],
	    Corr_nc_nc_z[d], nc_nc_z_0[d], nc_nc_z_1[d]);
  }
  fprintf(file_out, "#l\t<C>\t<C(Nc-<Nc>)>\t<(C-<C>)(Nc-<Nc>)^2>"
	  "\t<nc_i*nc_j>/<nc_i><Nc_j>-1(i=0,j=l)\tsame(i=%d,j=%d+l)"
	  "\t<nc_i*nc_j*Nc>/<nc_i><nc_j><Nc>-1(i=0,j=l)\tsame(i=%d,j=%d+l)\n",
	  DTAIL_MAX, DTAIL_MAX, DTAIL_MAX, DTAIL_MAX);
  for(l=0; l<100; l++){
    if(l>=L)break;
    fprintf(file_out, "%3d\t%.4f\t%.2g\t%.2g\t%.3g\t%.3g\t%.3g\t%.3g\n",
	    l, Cont_L[l], C2_L[l], C3_L[l],
	    nc_nc_L[0][l], nc_nc_L[DTAIL_MAX][l],
	    nc_nc_Nc_L[0][l], nc_nc_Nc_L[DTAIL_MAX][l]); 
  }
  fclose(file_out);
  printf("End rescaling\n");

}

/*// Empty_memory
  if(Cont_freq)free(Cont_freq);
  if(nc1i)free(nc1i);
  if(nc2i)free(nc2i);
  if(nc3i)free(nc3i);
  if(NC1nc1i)free(NC1nc1i);
  if(NC2nc1i)free(NC2nc1i);
  if(NC1nc2i)free(NC1nc2i);
  if(CNc)free(CNc);
  if(C321_ij)free(C321_ij);
  if(C332_ij)Empty_matrix_d(C332_ij, LEN);
  if(C342_ij)free(C342_ij);
  if(C343_2_ij)Empty_matrix_d(C343_2_ij, LEN);
  if(C343_3)free(C343_3);
  if(Cnn)Empty_matrix_d(Cnn, LEN);
  if(C1n)Empty_matrix_d(C1n, LEN);
  if(ncnc)Empty_matrix_d(ncnc, LEN);

  // Allocate
  int i;
  NC1=0; NC2=0; NC3=0;
  Cont_freq =malloc(L*sizeof(double));
  nc1i =malloc(L*sizeof(double));
  nc2i =malloc(L*sizeof(double));
  nc3i =malloc(L*sizeof(double));
  NC1nc1i =malloc(L*sizeof(double));
  NC2nc1i =malloc(L*sizeof(double));
  NC1nc2i =malloc(L*sizeof(double));
  CNc  =malloc(L*sizeof(double));

  C321_ij = malloc(L*sizeof(double));
  C332_ij = Allocate_mat2_d(L,L);
  C342_ij = malloc(L*sizeof(double));
  C343_2_ij=Allocate_mat2_d(L,L);
  C343_3 = malloc(L*sizeof(double));
  C353_i = malloc(L*sizeof(double));

  Cnn  =Allocate_mat2_d(L,L);
  ncnc =Allocate_mat2_d(L,L);
  C1n  =Allocate_mat2_d(L,L);

  for(i=0; i<L; i++){
    Cont_freq[i]=0; CNc[i]=0;
    C232_i[i]=0; C221_ij[i]=0;
    nc1i[i]=0; nc2i[i]=0; nc3i[i]=0;
    NC1nc1i[i]=0; NC1nc2i[i]=0; NC2nc1i[i]=0;
  }

  // End of statistics, compute averages
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
      for(int j=i+1; j<L; j++){
	ncnc[i][j]/=n_str;
	ncnc[j][i]=ncnc[i][j];
	Cnn[i][j]/=n_str;
	Cnn[j][i]=Cnn[i][j];
	C1n[i][j]/=n_str;
	C1n[j][i]/=n_str;
      }
    }
  } 
  int Nl=1;
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
  // 321_ij[(i-j)  //(ij)^3 mult=1 -> 332 -> NULL 

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
  double A332_ij[L][L]; int j;
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
  float C36=3*C363;
  for(i=0; i<L; i++){
    for(j=0; j<L; j++){
      C332_ij[i][j]-=  // C332
	(C342_ij[abs(i-j)]+C343_2_ij[i][j]+C343_3[i]+C343_3[j]);
      C343_2_ij[i][j]-=2*C353_i[i];  // C343_2
    }
    C343_3[i]-=C353_i[i];  // C343_3
    C353_i[i]-=C36;          // C353   
    }
  return;
} */


float Compute_DG_overT_contfreq(struct REM *E,
				short *aa_seq, float **C_nat,
				int *i_sec, char *c_sec, int VBS)
{  
  Set_to_zero(E); int L=E->L;
  Compute_E_nat(E, aa_seq, C_nat, i_sec, c_sec);
  if(E->REM==0)goto compute_dg;

  int i, j;
  for(i=0; i<L; i++){
    if(aa_seq[i]<0)continue;
    float *U=Econt_T[aa_seq[i]];
    for(j=i+IJ_MIN; j<L; j++){
      if(aa_seq[j]<0)continue;
      float Uij=U[aa_seq[j]];
      int l=j-i;
      float cU=Cont_L[l]*Uij;
      E->c1U1[i]+=cU;
      E->c1U1[j]+=cU;
      if(E->REM >= 2){
	E->E2cont1+=cU*Uij*Cont_1L[l]; // (<Cij>-<Cij>^2)*Uij ok
	E->E2cont2+=cU*C2_L[l];
	if(E->REM==3){
	  float cU2=cU*Uij;
	  E->E3cont1+=cU2*Uij*Cont_3L[l];
	  E->E3cont2+=cU2*Cont_2L[l];
	  E->E3cont3+=cU*C3_L[l];
	}
      }
    }
    E->E1  += E->c1U1[i];
  }
  E->E1/=2;
  //printf("DBG: E_nat, E1, E2 computed in DG\n");

  //if(E343_2) Compute_E343_2(E, aa_seq);
 compute_dg:
  E->G_misf=G_misfold(E);
  //printf("DBG: G_misfold computed in DG\n");

  float DeltaG_overT=DeltaG((float)(E->E_nat), E->G_misf, E->S_U);
  //printf("DBG: DeltaG computed in DG\n");

  if(Ini_VBS){
    printf("Computing native energy and misfold energy with Contfreq\n");
    printf("L= %d E_nat= %.3g E1= %.3g E2= %.3g G_misf= %.3g REM=%d\n",
	   E->L, E->E_nat, E->E1, E->E2, E->G_misf, E->REM);
    Ini_VBS=0;
  }
  return(DeltaG_overT);
}


void Mutate_Enat(struct REM *E, short *aa_seq,
		  float **C_nat, int *is,
		  int res_mut, short aa_new)
{
  int aa_old=aa_seq[res_mut]; if(aa_old<0)return;
  float dEl=Delta_E_loc(is[res_mut], aa_old, aa_new);
  E->E_loc_misf+=Mutate_E_loc_misf(aa_old, aa_new);
  E->E_nat+=dEl;

  float *Cnat_i=C_nat[res_mut];
  float *inter_mut=Econt_T[aa_new];
  float *inter_wt=Econt_T[aa_old];
  for(int j=0; j<E->L; j++){
    int l=abs(j-res_mut);
    if(l<IJ_MIN || aa_seq[j]<0 || Cnat_i[j]==0)continue;
    E->E_nat+=Cnat_i[j]*(inter_mut[aa_seq[j]]-inter_wt[aa_seq[j]]);
  }
  // Disulfide bonds
  if(AA_BKV[aa_old]=='C'){
    for(int j=0; j<E->L; j++){
      if(Cnat_i[j] && AA_BKV[aa_seq[j]]=='C'){ // Dis destroyed
	N_disulf_seq--; int l=abs(res_mut-j);
	Len_disulf_seq-=l;
	Log_len_disulf_seq-=log(l);
      }
    }
  }
  if(AA_BKV[aa_new]=='C'){
    for(int j=0; j<E->L; j++){
      if(Cnat_i[j] && AA_BKV[aa_seq[j]]=='C'){ // Dis created
	N_disulf_seq++; int l=abs(res_mut-j);
	Len_disulf_seq+=l;
	Log_len_disulf_seq+=log(l);
      }
    }
  }


}

float Mutate_DG_overT_contfreq(struct REM *E, short *aa_seq,
			       float **C_nat, int *is, char *c_sec,
			       int res_mut, short aa_new)
{
  if(aa_new>=20){ //aa_new<0 || 
    printf("WARNING, mutated amino acid is %d\n", aa_new);
    return(100);
  }

  int aa_old=aa_seq[res_mut]; if(aa_old<0)return(0);
  Mutate_Enat(E, aa_seq, C_nat, is, res_mut, aa_new);
  if(E->REM==0)goto mutate_DG;

  int j, l;
  float *inter_mut=Econt_T[aa_new];
  float *inter_wt=Econt_T[aa_old];
  double c1U1=0; //c332U2=0;
  for(j=0; j<E->L; j++){
    l=abs(j-res_mut);
    if(l<IJ_MIN)continue;
    if(aa_seq[j]<0)continue;
    float U_mut=inter_mut[aa_seq[j]];
    float U_wt=inter_wt[aa_seq[j]];
    float dU=U_mut-U_wt;
    float cd=Cont_L[l]*dU;
    c1U1+=cd;
    E->c1U1[j]+=cd;
    if(E->REM>=2){
      float cdU2=cd*(U_mut+U_wt);
      E->E2cont1+= cdU2*Cont_1L[l];
      E->E2cont2+= cd*C2_L[l];
      if(E->REM==3){
	E->E3cont1+=cd*(U_mut*U_mut*U_mut-U_wt*U_wt*U_wt)*Cont_3L[l];
	E->E3cont2+=cdU2*Cont_2L[l];
	E->E3cont3+=cd*C3_L[l];
      }
    }
    /*float U2_mut=U_mut*U_mut, U2_wt=U_wt*U_wt;
    float dU2=U2_mut-U2_wt;
    float cd2=c*dU2;
    E->E21 += C221_ij[l]*d2;
    if(E->REM==3){
      E->E321  += C321_ij[l]*(U2_mut*U_mut-U2_wt*U_wt);
      c332U2+= C332_ij[res_mut][j]*d2;
      E->c332U2[j]+=C332_ij[res_mut][j]*cd2;
      E->e342  += C342_ij[l]*d2;
      }*/
  }
  E->E1+=c1U1;
  E->c1U1[res_mut]+=c1U1;
  //if(E->REM==3)E->c332U2[res_mut]+=c332U2;
  //if(E343_2) Compute_E343_2(E, aa_seq);
 mutate_DG:
  E->G_misf=G_misfold(E);
  float DeltaG_overT=DeltaG((float)(E->E_nat), E->G_misf, E->S_U);
  return(DeltaG_overT);
}

float G_misfold(struct REM *E)
{
  if(E->REM>=2){
    E->E2site1=0; E->E3site1=0;
    E->E2site2=0; E->E3site2=0;
    int L=E->L;
    for(int i=0; i<L; i++){
      int di=tail(i, L);
      if(di>DTAIL_MAX){di=DTAIL_MAX;}
      double E22i=0, E32i=0;
      for(int j=i+1; j<L; j++){
	int l=j-i;
	int dj=tail(j, L), d;
	if(di<=dj){d=di;}else{d=dj;}
	if(d>DTAIL_MAX){d=DTAIL_MAX;}
	E22i+=nc_nc_L[d][l]*E->c1U1[j];
	if(E->REM==3)E32i+=nc_nc_Nc_L[d][l]*E->c1U1[j];
      }
      E->E2site1+=nc_nc_L[di][0]*E->c1U1[i]*E->c1U1[i];
      E->E2site2+=2*E22i*E->c1U1[i];
      if(E->REM==3){
	E->E3site1+=nc_nc_Nc_L[di][0]*E->c1U1[i]*E->c1U1[i];
	E->E3site2+=2*E32i*E->c1U1[i];
      }
    }
    //E->E22/=2;
    //E->E32/=2;
    E->E2 = E->E2cont1+(E->E2cont2*E->E1+E->E2site1+E->E2site2);
    E->E3 = E->E3cont1+(E->E3cont3*E->E1+E->E3site1+E->E3site2+E->E3cont2)
      *E->E1;

    /*if(E->REM==3){
      E->E342  = E->e342*E->E1;
      E->E353 *= E->E1;
      E->E363=(C363)*E->E1*E->E1*E->E1;
      E->E3 = E->E321 + E->E332 + E->E342 + E->E343_3 + E->E353 + E->E363;
      if(E343_2) E->E3 += E->E343_2;
      }*/
  }
  E->Tf=Compute_Tfreeze(E);
  float T=1; if(T<E->Tf){T=E->Tf;} //E->T=T;
  double H=E->E1-E->E2/(2*T);
  if(E->REM==3)H+=E->E3/(6*T*T);
  return(H-E->S_C);
}

float Compute_Tfreeze(struct REM *E){
  // Freezing temperature: S_C-E2/(2*T^2)+E3/(3*T^3)=0
  float Tf=sqrt(E->E2/(2*E->S_C));
  if(E->REM<3)return(Tf);
  if((E->E3>0)&&(E->E3 >(E->E2*Tf/sqrt(3))))//return(E->E3/E->E2);
    return(0.01);
  float K3=E->E3/3, K2=E->E2/2, S3=3*E->S_C;
  for(int it=0; it<20; it++){
    float T2=Tf*Tf;
    float f=E->S_C*T2*Tf-K2*Tf+K3;
    float f1=S3*T2-K2;
    if(fabs(f)<0.001)return(Tf);
    Tf-=(f/f1);
  }
  return(Tf); 
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

void Compute_E_nat(struct REM *E, short *aa_seq, float **C_nat,
		   int *i_sec, char *c_sec)
{
  // Native energy
  E->E_loc=Compute_E_loc(i_sec, aa_seq, E->L);
  E->E_loc_misf=Compute_E_loc_misf(aa_seq, E->L);
  Count_sec(&N_helix_seq, &N_strand_seq, c_sec, E->L);

  E->E_nat=0;
  for(int i=0; i<E->L; i++){
    if(aa_seq[i]<0)continue;
    float *U=Econt_T[aa_seq[i]];
    for(int j=i+IJ_MIN; j<E->L; j++){
      if(aa_seq[j]<0 || C_nat[i][j]==0)continue;
      E->E_nat+=C_nat[i][j]*U[aa_seq[j]];
    }
  }

  // Disulfides
  N_disulf_seq=0; Len_disulf_seq=0; Log_len_disulf_seq=0;
  for(int i=0; i<N_disulf; i++){
    int i1=Disulf_res1[i], i2=Disulf_res2[i];
    if(AA_BKV[aa_seq[i1]]=='C' &&
       AA_BKV[aa_seq[i2]]=='C'){
      N_disulf_seq++; int l=abs(i2-i1);
      Len_disulf_seq+=l;
      Log_len_disulf_seq+=log(l);
    }
  }
  //printf("N_disulf_seq= %d\n", N_disulf_seq);
  return;
}

void Count_sec(int *N_helix, int *N_strand, char *c_sec, int L)
{
  // Interactions among secondary strand elements
  // For each consecutive string of H or E, sum n-1
  *N_helix=0; *N_strand=0;
  for(int i=0; i<L; i++){
    int j=-1;
    if(c_sec[i]=='H'){
      j=i+1; while(c_sec[j]=='H'){(*N_helix)++; j++;}
    }else if(c_sec[i]=='E'){
      j=i+1; while(c_sec[j]=='E'){(*N_strand)++; j++;}
    }
    if(j>0)i=j-1;
  }
}

float Compute_DG_overT_threading(struct REM *E,
				 short *aa_seq, float **C_nat,
				 int *i_sec, char *c_sec)
{

  Set_to_zero(E);
  Compute_E_nat(E, aa_seq, C_nat, i_sec, c_sec);

  // Misfolding energy
  double U2, Nc1, Nc2, Nc3; 
  Threading_energies(NULL, &E->E1, &E->E2, &E->E3, &U2, &Nc1, &Nc2, &Nc3,
		     aa_seq, C_nat, E->L, NULL, 0);
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
  float T=1; if(E->Tf>T){T=E->Tf;} //E->T=T;
  double H=E->E1-E->E2/(2*T);
  if(E->REM==3)H+=E->E3/(6*T*T);
  return(H-E->S_C);
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

void Compute_E343_2(struct REM *E, short *aa_seq){
  double sum=0;
  for(int i=0; i<E->L; i++){
    if(aa_seq[i]<0)continue;
    //float *U=Econt_T[aa_seq[i]];
    double sum_i=0;
    for(int j=i+IJ_MIN; j<E->L; j++){
      if(aa_seq[j]<0)continue;
      if(j>=LEN_MAX)continue;
      //sum_i += C343_2_ij[i][j]*U[aa_seq[j]]*E->c1U1[j];
    }
    sum += sum_i*E->c1U1[i];
  }
  //E->E343_2=2*sum;
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

  //E->E24=(C242)*E->E1*E->E1;
  G_misfold(E); // Compute T freeze!
  fprintf(file_out, "# Third moment of energy used: ");
  if(E->REM==3){fprintf(file_out,"YES\n");}
  else{fprintf(file_out,"NO\n");}
  fprintf(file_out, "# Temperature factor: %.2f\n", E->T);
  fprintf(file_out, "# E2/L= %.3g E3/L= %.3g s_C= %.3g\n",
	  E->E2/E->L, E->E3/E->L, E->S_C/E->L);
  fprintf(file_out,"# Tfreeze= %.2f\n", E->T*E->Tf);
  float G_misf_freeze=E->E1-E->E2/(2*E->Tf)-E->Tf*E->S_C;
  if(E->REM==3)G_misf_freeze+=E->E3/(6*E->Tf*E->Tf);
  fprintf(file_out,
	  "# E_nat/L= %.4f Gfreeze/TL= %.4f E1/L= %.3g\n",
	  E->E_nat/E->L, G_misf_freeze/E->L, E->E1/E->L);
  double nc1=0;
  for(int l=IJ_MIN; l<E->L; l++)nc1+=Cont_L[l]*(E->L-l);
  fprintf(file_out, "# NC1/L= %.3f NC2/L= %.2g NC3/L= %.2g nc1/L= %.3f\n",
	  NC1/E->L, NC2/E->L, NC3/E->L, nc1/E->L);
  fprintf(file_out, "# L= %d S_U/L= %.3f S_C/L= %.3f\n",
	  E->L, E->S_U/E->L, E->S_C/E->L);
  fprintf(file_out, "# Temp. DeltaG(nat-nonat)  G_misf-G_unf\n");
  float T_STEP=0.01, T=E->Tf-2*T_STEP, T0=-1;
  float DeltaG, G_misf, DG0=100000;
  int k_nonat=0, k_max=200;
  //for(T=T_INI; T<T_END; T+=T_STEP){
  for(int k=0; k<k_max; k++){
    if(T<E->Tf){
      G_misf=G_misf_freeze; // Freezing
    }else{
      G_misf=E->E1-E->E2/(2*T)-T*E->S_C;
      if(E->REM==3)G_misf+=E->E3/(6*T*T);
    }
    DeltaG=E->E_nat-G_misf+T*log(1+exp(G_misf/T+E->S_U));
    float G_nonat=G_misf+T*E->S_U;
    if(DeltaG>-100)fprintf(file_out, "%.3f\t%.2f\t%.2f\n",
			   T*E->T, DeltaG, G_nonat);
    if(DeltaG < DG0){DG0=DeltaG; T0=T*E->T;}
    if(G_nonat>=0){k_nonat++; if(k_nonat==4)break;}
    k++; T+=T_STEP;
  }
  fclose(file_out);
  return(T0);
}

float **Get_target(char *file_str, char *name_tar, int *len_tar)
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
  float **C_nat=Fill_C_nat(length, contact, 1);
  return(C_nat);
}

float **Fill_C_nat(int length, short **contact, int num_chain)
{
  // Allocate native contact matrix
  float **C_nat=malloc(length*sizeof(float *));
  int i, j, Nc=0;
  for(i=0; i<length; i++){
    C_nat[i]=malloc(length*sizeof(float));
    for(j=0; j<length; j++)C_nat[i][j]=0;
  }
  for(i=0; i<length; i++){
    short *Ci=contact[i];
    while(*Ci>=0){
      C_nat[i][*Ci]++; Ci++; Nc++;
    }
  }
  Nc/=num_chain;
  for(i=0; i<length; i++){
    for(j=0; j<i; j++){
      float C=(C_nat[i][j]+C_nat[j][i]);
      if(num_chain != 1)C/=num_chain;
      C_nat[i][j]=C; C_nat[j][i]=C;
    }
  }

  printf("Allocating native contact matrix, L=%d Nc=%d\n", length, Nc);
  return(C_nat);
}

void Test_contfreq(struct REM *E, short *aa_seq, float **C_nat,
		   int *i_sec, char *c_sec, char *name)
{
  // WARNING, you have to make sure that Initialize_REM has been run before
  FILE *file_out=Output_file(name, "Threading", "dat");
  char s[5000];

  int nT=16; float T_ini=0.75, T_fact=1.15, Ts[nT];
  Ts[0]=1; for(int j=1; j<nT; j++){Ts[j]=T_ini; T_ini*=T_fact;}
  double E_ave[nT];
  double EC1, EC2, EC3, U2, Nc1, Nc2, Nc3; int L=E->L;
  Threading_energies(E_ave, &EC1, &EC2, &EC3, &U2, &Nc1, &Nc2, &Nc3,
		     aa_seq, C_nat, L, Ts, nT);

  sprintf(s,
	  "T= %.2f L= %d S_U=%.2f S_C=%.2f\n" 
	  "THREADING results REM= %d\n"
	  "<Nc/L>= %.2g <(Nc-<Nc>)^2>/<Nc>^2 = %.2g "
	  "<(Nc-<Nc>)^3>/<Nc>^3 = %.2g\n"
	  "<U^2>= %.3g\n<EC>= %.3g           /L= %8.3g\n"
	  "<(EC-<EC>)^2>/2= %6.2f  /L= %8.3g\n"
	  "<(EC-<EC>)^3>/6= %6.2f  /L= %8.3g\n"
	  "DG_misfold/L = %.4g\n"
	  "Boltzmann averages: %.4g T= %.2g\n#T <E>\n",
	  E->T, E->L, E->S_U, E->S_C,
	  E->REM, 
	  Nc1/L, (Nc2/(Nc1*Nc1)-1),
	  (Nc3/(Nc1*Nc1*Nc1)-3*Nc2/(Nc1*Nc1)+2),
	  U2, EC1, EC1/L, EC2/2, EC2/(2*L),
	  EC3/6, EC3/(6*L), (EC1-EC2/(2*E->T))/L, E_ave[0]/L, Ts[0]);

  printf("%s", s);
  fprintf(file_out, "%s", s);
  for(int j=1; j<nT; j++){
    fprintf(file_out, "%.3g\t%.4g\n", Ts[j], E_ave[j]/L);
  }

  int i, j; double nc=0;
  for(i=0; i<L; i++){
    if(aa_seq[i]<0)continue;
    float *U=Econt_T[aa_seq[i]];
    for(j=i+IJ_MIN; j<L; j++){
      if(aa_seq[j]<0)continue;
      float Uij=U[aa_seq[j]];
      int l=j-i; if(l>=LEN_MAX)continue; //l=LEN_MAX-1;
      U2+=Cont_L[l]*Uij*Uij; nc+=Cont_L[l];
    }
  }

  int REM=E->REM; E->REM=2;
  Compute_DG_overT_contfreq(E, aa_seq, C_nat, i_sec, c_sec, 1);
  float T2=1; if(T2<E->Tf){T2=E->Tf;}
  float G2=E->G_misf;

  E->REM=3;
  Compute_DG_overT_contfreq(E, aa_seq, C_nat, i_sec, c_sec, 1);
  float T3=1; if(T3<E->Tf){T3=E->Tf;}
  float G3=E->G_misf;
  E->REM=REM;

  float H2=(E->E1-E->E2/T2)/L;
  float H3=(E->E1-E->E2/T3+E->E3/(2*T3*T3))/L;

  sprintf(s,
	  "E_nat/L = %.4f G_unfolded/L= %.4f\n"
	  "CONTFREQ results:\n"
	  "<Nc/L>= %.2g <(Nc-<Nc>)^2>/<Nc>^2 = %.2g "
	  "<(Nc-<Nc>)^3>/<Nc>^3 = %.2g\n"
	  "<U^2>= %.3g\n<EC>= %.3g\n"
	  "<(EC-<EC>)^2>/2= %.4g\t"
	  "E2cont1= %.3g E2cont2*E1= %.3g E2site1= %.3g E2site2= %.3g\n"
	  "<(EC-<EC>)^3>/6= %.4g\t"
	  "E3cont1= %.3g E3cont2*E1= %.3g E3cont3*E1*E1= %.2g "
	  "E3site1*E1= %.2g E3site2*E1= %.2g\n"
	  "Tfreeze= %.2g (REM2) %.2g (REM3)\n"
	  "DG_misfold/L = %.4g (REM2) %.4g (REM3)\n"
	  "Boltzmann average E/L: %.4g (REM2) %.4g (REM3)\n",
	  //"<(EC-<EC>)^3>/6= %6.2f  E321= %.2g E332= %.2g\n"
	  //"\tE342= %.2g E343_2= %.2g E343_3= %.2g\tE353=%.2g E363= %.2g\n"
	  //"E36 used: %d\n"
	  E->E_nat/L, E->S_U/L,
	  NC1/L, zs2/(Y_ave*Y_ave), zs3/(Y_ave*Y_ave*Y_ave),
	  //NC2/(NC1*NC1), NC3/(NC1*NC1*NC1),
	  //Y_ave, zs2/(Y_ave*Y_ave), zs3/(Y_ave*Y_ave*Y_ave),
	  U2/nc, E->E1,
	  E->E2/2, E->E2cont1, E->E2cont2*E->E1, E->E2site1, E->E2site2, 
	  E->E3/6, E->E3cont1, E->E3cont2*E->E1, E->E3cont3*E->E1*E->E1,
	  E->E3site1*E->E1, E->E3site2*E->E1,
	  //E->E21, E->E332, E->E342, E->E343_2, E->E343_3,
	  //E->E353, E->E363, E_36,
	  T2, T3, G2/L, G3/L, H2, H3);
  printf("%s", s);
  fprintf(file_out, "%s", s);
  
  fclose(file_out);
}

void Threading_energies(double *E_ave,
			double *EC1, double *EC2, double *EC3, double *U2,
			double *Nc1, double *Nc2, double *Nc3,
			short *aa_seq, float **C_nat, int len_pdb,
			float *T, int nT)
{
  float q_thr=0.5; // Max overlap
  int LEN_THR=100;

  //int N_PROT=Read_structures(file_str, prots);
  long n_str=0, Noncompact=0;

  float Nc_nat=0;
  if(C_nat){
    for(int i=0; i<len_pdb; i++)
      for(int j=i+IJ_MIN; j<len_pdb; j++)
	if(C_nat[i][j])Nc_nat+=C_nat[i][j];
  }

  double E_norm[nT];
  if(E_ave)for(int j=0; j<nT; j++){E_norm[j]=0; E_ave[j]=0;}
  *EC1=0; *EC2=0; *EC3=0; *U2=0;
  *Nc1=0; *Nc2=0; *Nc3=0;
 
  /* Generate alternative structures by threading */
  //printf("### Contact statistics\n");
  //printf("%d contact matrices found in file %s\n", N_PROT, FILE_STR);
  printf("Threading seq on fragments >= %d residues and ", LEN_THR);
  printf("NC/L ~ %.2f - %.2fL^(-1/3) +/- %.2f\n", C_ASYMPT,C_SLOPE,C_THR);

  for(short j_prot=0; j_prot< N_PROT; j_prot++){

    short len2=prots[j_prot].length;
    if(len2 < LEN_THR)continue;
    int ndom;
    if(len_pdb>len2){
      ndom=len_pdb/len2;
      if(len_pdb > ndom*len2)ndom++;
    }else{
      ndom=1;
    }
    int len_dom=len_pdb/ndom;
    int NCmin, NCmax;
    Contact_range(&NCmin, &NCmax, len_dom);

    short **contact=prots[j_prot].contact;

    short res1, *res2;
    for(int first=0; first<=len2-len_dom; first++){
      
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
      if((Nc<NCmin)){Noncompact++; continue;} //||(Nc>NCmax)
      n_str++;

      // Align structure with each of the ndom domains of the threaded sequence
      int ini=0; // Initial site of the threaded sequence 
      double e1=0, e2=0; int overlap=0; Nc=0;
      for(int idom=0; idom<ndom; idom++){
	int shift=ini-first;
	for(res1=first; res1<last; res1++){
	  int i=res1+shift; if(aa_seq[i]<0)continue;
	  float *U=Econt_T[aa_seq[i]], Uij;
	  res2=contact[res1];
	  while((*res2>=0)&&(*res2 < last)){
	    if(*res2 >= res1+IJ_MIN){
	      int j=*res2+shift;
	      if(aa_seq[j]>=0){
		if(C_nat && C_nat[i][j])overlap++;
		Uij=U[aa_seq[j]];
		e1+=Uij; e2+=Uij*Uij; Nc++;
	      }
	    }
	    res2++;
	  }
	}
	ini+=len_dom;
      } // end sequence domain
      if(overlap > q_thr*sqrt(Nc*Nc_nat))continue;
      (*Nc1)+=Nc; float c2=Nc*Nc;
      (*Nc2)+=c2; (*Nc3)+=c2*Nc;
      (*U2)+=e2;
      (*EC1)+=e1; e2=e1*e1; (*EC2)+=e2; (*EC3)+=e2*e1;
      if(E_ave){
	for(int j=0; j<nT; j++){
	  double w=exp(-e1/T[j]); E_ave[j]+=e1*w; E_norm[j]+=w;
	}
      }
    } // end fragment 
  } // end protein

  /* End of statistics, compute averages */
  (*EC1)/=n_str; (*EC2)/=n_str; (*EC3)/=n_str;
  (*EC3)=(*EC3)-3*(*EC1)*(*EC2)+2*(*EC1)*(*EC1)*(*EC1);
  (*EC2)=(*EC2)-(*EC1)*(*EC1);
  (*U2)/=(*Nc1);
  (*Nc1)/=n_str; (*Nc2)/=n_str; (*Nc3)/=n_str;
  
  if(E_ave){
    for(int j=0; j<nT; j++){E_ave[j]/=E_norm[j];}
  }
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

    prot[N_prot].length=length;
    prot[N_prot].n_cont=n_cont;
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
  E2->DeltaG=E1->DeltaG;
  E2->E_nat=E1->E_nat;
  E2->E1=E1->E1;
  E2->E2=E1->E2;
  E2->E2cont1=E1->E2cont1;
  E2->E2cont2=E1->E2cont2;
  E2->E2site1=E1->E2site1;
  E2->E2site2=E1->E2site2;
  E2->E3=E1->E3;
  E2->E3cont1=E1->E3cont1;
  E2->E3cont2=E1->E3cont2;
  E2->E3cont3=E1->E3cont3;
  E2->E3site1=E1->E3site1;
  E2->E3site2=E1->E3site2;
  for(int i=0; i<E1->L; i++){
    E2->c1U1[i]=E1->c1U1[i];
    //E2->c332U2[i]=E1->c332U2[i];
  }
}

void Empty_E_REM(struct REM *E){
  if(E->c1U1){free(E->c1U1); E->c1U1=NULL;}
  //free(E->c332U2);
} 

float Corr_coeff(float *offset, float *slope, float *xx, float *yy, int n){
  double x1=0, x2=0, y1=0, y2=0, xy=0;
  int i; float *x=xx, *y=yy;
  for(i=0; i<n; i++){
    x1 += *x; x2+= (*x)*(*x);
    y1 += *y; y2+= (*y)*(*y);
    xy += (*x)*(*y); x++; y++;
  }
  float nx2=n*x2-x1*x1, ny2=n*y2-y1*y1;
  *slope=(n*xy-y1*x1); if(nx2)(*slope)/=nx2;
  *offset=(y1-(*slope)*x1)/n;
  float r=(*slope);
  if((nx2)&&(ny2))r*=sqrt(nx2/ny2);
  return(r);
}

int tail(int i, int L){
  int d=L-1-i;
  if(i<d){return(i);}else{return(d);}
}
