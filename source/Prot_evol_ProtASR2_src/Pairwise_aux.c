#include "pairwise.h"
#include "allocate.h"
#include "protein4.h"
#include "energy_BKV.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "Pairwise_aux.h"
#include "optimization.h"
#include "REM_new.h"

// Parameters for selecting compact structures
#define C_OFFSET 4.0  // Asymptotic value of Nc/N
#define C_SLOPE 8.0   // Nc/N ~ C_OFFSET - C_SLOPE N^(-1/3)
//#define C_SLOPE 8.2   // Nc/N ~ C_OFFSET - C_SLOPE N^(-1/3)
#define C_THR 1.0      // Maximum deviation of Nc/N from typical value
                       // Very important parameter!!!
int IJ_MIN=4, L_MAX;
float ERR_THR=0.04; // Mean accepted error for Normalize_Q ERR_THR/40=0.001

// Energy parameters
int CP_Error=0;


// Concact statistics
float *Cont_f, *Cz1, *Cz2, *Cnc1;
float *nc_nc;
// <Nc>, <Nc^2>, <Nc^3>
//float *Nc1L, *Nc2L, *Nc3L;
float *nc2L;
float *Nc2L_center, *Nc3L_center;
double y1_ave, z_expo, Slope_x;

void Normalize_cont_freq(float *Cont_norm, float *Cont_f,
			 float *Cnc1_norm, float *Cnc1,
			 float *nc_nc_norm, float *nc_nc,
			 float *CNc1_norm, float *Cz1,
			 float *CNc2_norm, float *Cz2,
			 float *Cont_indir_norm, float *Cont_f_indir,
			 double Nc1L, double Nc2L, double Nc3L,
			 double NcL_indir, int L);

/*float Compute_cont_freq(int ij, int L);
float Compute_cont_freq_Nc(int ij, int L);
float Compute_cont_freq_Nc2(int ij, int L);
float Compute_cont_freq_indir(int ij, int L);*/



// Pairwise probabilities Q, predicted
float qmin=0.000000001, lmin; // Minimum value of q
void Compute_Qprime(float *Qprime, struct pair_class  *pair, int cont,
		    float Lambda);
float Normalize_Qprime(float *Q_pred, struct pair_class  *pair);
void Update_q(float *qnew, float *Q, float *zi, float *zj);
float Update_z_naive(float *zi, float *zj, float *Q, struct pair_class *pair);
float Update_z_Newton(float *zi, float *zj, float *yi, float *yj,
		      float *Q, struct pair_class *pair);
float Error_norm(float *Q, struct pair_class *pair);
void Print_error_Q(struct pair_class *pair, float *Q,
		   float * zi, float * zj, float * ziopt, float * zjopt);
int printconv=0;


// Q observed
double Marginalize_N2(float *P1_i, float *P1_j, double *N2_obs);
void Pairwise_Q(float *Q_obs, float *P1_i, float *P1_j,
		  double *N2_obs, double n);

// Scores
void Compute_logarithm(float *log_Q, float *Q, float qmin, float lmin, int Na2);
float Compute_specific_mut_inf(struct pair_class *pair);
float dKL_2(float *log_Qpred, struct pair_class *pair);
float RME_2(float *log_Q_pred, struct pair_class *pair);
float Likelihood_2(float *Qpred, struct pair_class *pair);
static float Corr_coeff(float *offset, float *slope, float *x, float *y, int n);
static float Corr_coeff_w(float *offset, float *slope,
			  float *x, float *y, double *w, int n);

// Indirect contacts
float *Cont_freq_indir;
float NcL0_indir, NcL1_indir;    // <Nc_indirect>~L(NCL0+NCL1*L^(-1/3))
float *NcL_indir;
float Score_indir(struct pair_class *pair_p, int Nij_class,
		  int LS, int comp, float f_indir);
float Indirect_Q_PDB(float *Q, 
		     struct pair_class *pair_ij,
		     struct pair_class *pairs,
		     int Nij_class, int contact, float f);
float Indirect_Q_ali(float *Q, int i_ij, struct pair_class *pairs, float f);
void Indirect_Q_k(double *P_ind, int Naa,
		  struct pair_class *pair_ij,
		  struct pair_class *pair_ik,
		  struct pair_class *pair_jk);
float *Most_lik_Q(struct pair_class *pair);
float Compute_indirect(float *Q, struct pair_class *pair, double *P, float f);
static int Direct_contact(int *ic, int i, int j,
			  struct contact *cont_list, int Nc);

/*void Indirect_Q_previous(float *Qprime, int Naa,
			 struct pair_class *pair_ij,
			 struct pair_class *pair_ik,
			 struct pair_class *pair_jk);
void Indirect_Q_internal(float *Qprime, int Naa,
			 struct pair_class *pair_ij,
			 struct pair_class *pair_ik,
			 struct pair_class *pair_jk);
void Indirect_Q_null(float *Qprime, int Naa,
		     struct pair_class *pair_ij,
		     struct pair_class *pair_ik,
		     struct pair_class *pair_jk);
*/


// Auxiliary
static void Copy_vec(float *v1, float *v2, int N);
int Check_AA(short i_aa, int Naa);

/************************ Contact frequencies ***************************/

// Contact statistics

void Initialize_REM_2(int L_max, int ij_min,
		      struct protein *prot_ptr, int N_prot)
{
  IJ_MIN=ij_min;
  L_MAX=L_max;

  // Allocate
  Cont_f =malloc(L_max*sizeof(float));
  Cnc1  =malloc(L_max*sizeof(float));
  Cz1  =malloc(L_max*sizeof(float));
  Cz2  =malloc(L_max*sizeof(float));
  nc_nc =malloc(L_max*sizeof(float));

  printf("### Contact statistics\n%d contact matrices read\n", N_prot);

  /* Scaling of contacts */
  float Slope_tmp=C_SLOPE/C_OFFSET; 
  float Ncp[N_prot], Xp[N_prot], logLp[N_prot];
  int np=0, select[N_prot], i;
  for(i=0; i< N_prot; i++){
    struct protein *prot=prot_ptr+i;
    float L=prot->nres;
    float x=pow(L, -0.333333), c=(prot->Nc/L);
    if(fabs(c-C_OFFSET*(1-Slope_tmp*x))>C_THR){
      select[i]=0;
    }else{
      select[i]= 1;
      logLp[np]=log(L); Xp[np]=x; Ncp[np]=c;
      np++;
    }
  }
  // Fits of Nc(L), Nc^2(L) etc.
  printf("Fitting number of contacts versus surface to volume scaling.\n");
  printf("%d proteins kept and %d omitted due to large deviations\n",
	 np, N_prot-np);
  // First moment
  float Nc1L0, Nc1L1;    // <Nc>~L(NC1L0+NC1L1*L^(-1/3))
  float r=Corr_coeff(&Nc1L0, &Nc1L1, Xp, Ncp, np);
  printf("<Nc/L> ~ %.2f + %.1fL^(-1/3)  r= %.3f\n", Nc1L0, Nc1L1, r);
  Slope_x=Nc1L1/Nc1L0;

  // Normalize Nc/L by fitted value
  float Yp[np]; y1_ave=0;
  printf("y = NC/L(1 + %.2fL^(-1/3))\n", Slope_x);
  for(i=0; i< np; i++){
    Yp[i]=Ncp[i]/(1+Slope_x*Xp[i]);
    y1_ave+=Yp[i];
  } 
  y1_ave/=np;
  float coeff, slope;
  r=Corr_coeff(&coeff, &slope, Xp, Yp, np);
  printf("y ~ %.3f + %.1f*L^(-1/3)  r= %.3f\n", coeff, slope, r);
  printf("<y>= %.3f\n", y1_ave);

  // Second and third moment with power law
  float logY2p[np], logY3p[np];
  for(i=0; i< np; i++){
    float ly=log(fabs(Yp[i]-y1_ave));
    logY2p[i]=2*ly; logY3p[i]=3*ly;
  }   
  float Y2L0, Y2L1, r2;  // <(y-<y>)^2>~exp(Y2L0+Y2L1*log(L))
  r2=Corr_coeff(&Y2L0, &Y2L1, logLp, logY2p, np);
  printf("<(y-<y>)^2> ~ exp(%.2f + %.1flog(L))  r= %.3f\n", Y2L0,Y2L1,r2);
  float Y3L0, Y3L1, r3;  // <(y-<y>)^3>~exp(Y3L0+Y3L1*log(L))
  r3=Corr_coeff(&Y3L0, &Y3L1, logLp, logY3p, np);
  printf("<(y-<y>)^3> ~ exp(%.2f + %.1flog(L))  r= %.3f\n", Y3L0,Y3L1,r3);
  r2*=r2; r3*=r3;
  z_expo=(r2*Y2L1/2+r3*Y3L1/3)/(r2+r3);
  printf("Average scaling exponent: %.3f\n", z_expo);

  // Compute Z score
  float Zp[np]; double z2_ave=0, z3_ave=0; 
  for(i=0; i< np; i++){
    float z=(Yp[i]-y1_ave)*exp(-z_expo*logLp[i]);
    Zp[i]=z; z2_ave+=z*z; z3_ave+=z*z*z;
  }
  z2_ave/=np; z3_ave/=np;
  r=Corr_coeff(&coeff, &slope, Xp, Zp, np);
  printf("z=(y-<y>)*L^%.2f ~ %.2f + %.2fL^(-1/3)  r= %.3f\n",
	 -z_expo,coeff,slope,r);
  printf("<z^2>= %.3g\n", z2_ave);

  /**************** Statistics of contacts ****************/
  long Cont_num[L_max], Cont_norm[L_max], nc_nc_num[L_max],
    ncCont_num[L_max], z_Cont_num[L_max], z2_Cont_num[L_max];
  for(i=0; i<L_max; i++){
    Cont_num[i]=0; Cont_norm[i]=0;
    ncCont_num[i]=0; nc_nc_num[i]=0;
    z_Cont_num[i]=0; z2_Cont_num[i]=0;
  }

  float nc2[N_prot]; int kp=0;
  for(int ip=0; ip< N_prot; ip++){
    if(select[ip]==0)continue;
    struct protein *prot=prot_ptr+ip;
    int L=prot->nres;
    int nc[L], ncont_l[L];
    for(i=0; i<L; i++){nc[i]=0; ncont_l[i]=0;}
    for(int ic=0; ic<prot->Nc; ic++){
      struct contact *cont=prot->cont_list+ic;
      ncont_l[cont->res2-cont->res1]++;
      nc[cont->res1]++; nc[cont->res2]++;
    }
    // <c_ij ni>
    for(int ic=0; ic<prot->Nc; ic++){
      struct contact *cont=prot->cont_list+ic;
      ncCont_num[cont->res2-cont->res1]+=(nc[cont->res1]+nc[cont->res2]);
    }
    // Sum protein to global counters
    float z=Zp[kp];
    float scale2=(1+Slope_x*Xp[kp]); scale2*=scale2;
    //float z=Yp[kp];
    for(i=0; i<L; i++){
      if(i>=L_max)break;
      Cont_norm[i]+=L-i;
      Cont_num[i]+=ncont_l[i];
      z_Cont_num[i]+=ncont_l[i]*z;  // Average contacts per residue
      z2_Cont_num[i]+=ncont_l[i]*z*z;
      for(int j=0; j<(L-i); j++){
	if(j >= L_max)break;
	nc_nc_num[i]+=nc[j]*nc[j+i]/scale2;
      }
    } // End pairs
    // sum_i (ni/L)^2
    float n2=0; for(i=0; i<L; i++)n2+=nc[i]*nc[i]/scale2;
    nc2[kp]=n2/L; kp++;
  } // end proteins
  if(kp!=np){
    printf("ERROR, wrong number of proteins %d instead of %d\n",
	   kp, np); exit(8);
  }

  float nc2L0, nc2L1;  // sum_i <(ni/L)^2>~(nc2L0+nc2L1*L^(-1/3))
  r=Corr_coeff(&nc2L0, &nc2L1, Xp, nc2, np);
  printf("sum_i ni^2/L ~ %.2f + %.1fL^(-1/3)  r= %.3f\n", nc2L0, nc2L1, r);
  nc2L0=0; for(i=0; i<np; i++)nc2L0+=nc2[i]; nc2L0/=np;
  printf("Average sum_i ni^2/L: %.2f\n", nc2L0);

  //printf("Cont_freq <NcCij>\n");
  for(i=0; i<L_max; i++){
    float norm=Cont_norm[i];
    Cont_f[i]=(float)Cont_num[i]/norm;
    Cnc1[i]=(float)ncCont_num[i]/(2*norm);
    Cz1[i]=(float)z_Cont_num[i]/norm;
    Cz2[i]=(float)z2_Cont_num[i]/norm;
    nc_nc[i]=nc_nc_num[i]/norm;
    //if((i>=4)&&(i<30))printf("%.4f %.3f\n", Cont_f[i], CNc1[i]);
    if((i>=4)&&(Cont_f[i]<0)){
      printf("ERROR, negative Cont_freq:\n");
      printf("%d cf= %.4f %ld %ld Ncf= %.4f Nc2f= %.4f\n",
	     i, Cont_f[i], Cont_num[i], Cont_norm[i], Cz1[i], Cz2[i]);
      exit(8);
    }
  }

  nc2L=malloc(L_max*sizeof(float));
  Nc1L=malloc(L_max*sizeof(float));
  Nc2L=malloc(L_max*sizeof(float));
  Nc3L=malloc(L_max*sizeof(float));
  Nc2L_center=malloc(L_max*sizeof(float));
  Nc3L_center=malloc(L_max*sizeof(float));
  //printf("#L N<Nc> <Nc^2> <Nc^2>-<Nc>^2\n");
  for(i=4; i<L_max; i++){
    double x=pow(i,-0.33333333);
    float core=(1+Slope_x*x);
    //nc2L[i]=(nc2L0+nc2L1*x)*i;
    nc2L[i]=(nc2L0*core*core)*i;
    double scale=pow(i, z_expo);
    float Lcore=core*i;
    Nc1L[i]=y1_ave*Lcore;
    float y2=y1_ave*y1_ave, z2=z2_ave*scale*scale;
    Nc2L[i]=(y2+z2)*Lcore*Lcore;
    float Lscale=Lcore*scale;
    Nc2L_center[i]=z2_ave*Lscale*Lscale;
    Nc3L_center[i]=z3_ave*Lscale*Lscale*Lscale;
    Nc3L[i]=y1_ave*(y2+3*z2)*Lcore*Lcore*Lcore+Nc3L_center[i];
    //if((i>90)&&(i<120))
    //printf("%d %.3g %.3g %.3g\n", i, Nc1L[i], Nc2L[i], Nc2L_center[i]);
  }
  //exit(8);

}

void Indirect_cont_stat(int L_max, int ij_min,
			struct protein *prot_ptr, int N_prot)
{
  // Allocate
  int i;
  Cont_freq_indir =malloc(L_max*sizeof(float));
  long *Cont_num=malloc(L_max*sizeof(long));
  long *Cont_norm=malloc(L_max*sizeof(long));
  for(i=0; i<L_max; i++){
    Cont_num[i]=0; Cont_norm[i]=0;
  }

  /* Statistics of contacts */
  float *Nc1=malloc(N_prot*sizeof(float));
  float *Nc2=malloc(N_prot*sizeof(float));
  float *Xp=malloc(N_prot*sizeof(float));
  for(short ip=0; ip< N_prot; ip++){
    struct protein *prot=prot_ptr+ip;
    float c=prot->Nc_indirect/(float)prot->nres;
    Nc1[ip]=c; Nc2[ip]=c*c;
    Xp[ip]=pow(prot->nres, -0.333333);
    if(prot->indirect_contacts==NULL)continue;
    struct contact *cont=prot->indirect_contacts;
    for(i=0; i<prot->Nc; i++){
      int l=cont->res2-cont->res1; cont++;
      if(l<L_max)Cont_num[l]++; 
    }
    for(i=0; i<L_max; i++){
      int nl=prot->nres-i;
      if(nl<=0)break;
      Cont_norm[i]+=nl;
    }
  }
  printf("Cont_freq_indirect:\n");
  for(i=0; i<L_max; i++){
    Cont_freq_indir[i]=(float)Cont_num[i]/(float)Cont_norm[i];
    if((i>=4)&&(i<40))printf("%d %.4f\n", i, Cont_freq_indir[i]);
    if((i>=4)&&(Cont_freq_indir[i]<0)){
      printf("ERROR, negative Cont_freq:\n");
      printf("%d cf= %.4f %ld %ld \n",
	     i, Cont_freq_indir[i], Cont_num[i], Cont_norm[i]);
      exit(8);
    }
  }
  float r=Corr_coeff(&NcL0_indir, &NcL1_indir, Xp, Nc1, N_prot);
  printf("NC_indirect/L ~ %.2f + %.2f^(-1/3)  r= %.3f\n",
	 NcL0_indir, NcL1_indir, r);

  NcL_indir=malloc(L_max*sizeof(float));
  for(i=10; i<L_max; i++){
    float x=pow(i,-0.33333);
    NcL_indir[i]=(NcL0_indir+NcL1_indir*x)*i;
  }
}


/*
float Compute_cont_freq(int ij, int L){
  double norm=0;
  for(int l=IJ_MIN; l<L; l++)norm+=Cont_f[l]*(L-l);
  return(Cont_f[ij]*Nc1L[L]/norm);
}

float Compute_cont_freq_Nc(int ij, int L){
  double norm=0;
  for(int l=IJ_MIN; l<L; l++)norm+=CNc1[l]*(L-l);
  return(CNc1[ij]*Nc2L[L]/norm);
}

float Compute_cont_freq_Nc2(int ij, int L){
  double norm=0;
  for(int l=IJ_MIN; l<L; l++)norm+=CNc2[l]*(L-l);
  return(CNc2[ij]*Nc3L[L]/norm);
}

float Compute_cont_freq_indir(int ij, int L){
  double norm=0;
  for(int l=IJ_MIN; l<L; l++)norm+=Cont_freq_indir[l]*(L-l);
  return(Cont_freq_indir[ij]*NcL_indir[L]/norm);
}
*/

void Normalize_cont_freq(float *Cont_norm, float *Cont_f,
			 float *Cnc1_norm, float *Cnc1,
			 float *nc_nc_norm, float *nc_nc,
			 float *CNc1_norm, float *Cz1,
			 float *CNc2_norm, float *Cz2,
			 float *Cont_indir_norm, float *Cont_f_indir,
			 double Nc1, double Nc2, double Nc3,
			 double Nc_indir, int L)
{
  double sum_c=0, sum_c_nc=0, sum_c_z=0,
    sum_c_z2=0, sum_nc_nc=0, sum_nl=0, sum_c_indir=0;
  int l;
  for(l=IJ_MIN; l<L; l++){
    int nl=L-l;
    sum_nl+=nl;
    sum_c+=Cont_f[l]*nl;
    sum_c_nc+=Cnc1[l]*nl;
    sum_c_z+=Cz1[l]*nl;
    sum_c_z2+=Cz2[l]*nl;
    sum_nc_nc+=nc_nc[l]*nl;
    if(Cont_f_indir)sum_c_indir+=Cont_f_indir[l]*nl;
  }
  float Nc12=Nc1*Nc1;
  float norm_c_z= (Nc2-Nc12)/sum_c_z;
  float norm_c_z2=(Nc3-Nc12*Nc1)/sum_c_z2;
  float norm_c=Nc1/sum_c;
  float norm_c_nc=nc2L[L]/sum_c_nc;
  float norm_nc_nc=2*Nc2/sum_nc_nc; // sum_i<j n_i*nj=4Nc/2!!
  float norm_c_indir=0;
  if(Cont_f_indir)norm_c_indir=Nc_indir/sum_c_indir;
  float nc=2*Nc1/L, nc2=2*Nc1*Nc1/sum_nl;

  // float cc=Nc1L*Nc1L/norm_nl;
  for(l=IJ_MIN; l<L; l++){
    Cont_norm[l]=Cont_f[l]*norm_c;
    nc_nc_norm[l]=nc_nc[l]*norm_nc_nc-nc2;
    CNc1_norm[l]=Cz1[l]*norm_c_z;
    CNc2_norm[l]=Cz2[l]*norm_c_z2-2*Nc1*CNc1_norm[l];
    Cnc1_norm[l]=Cnc1[l]*norm_c_nc-Cont_norm[l]*nc;
    if(Cont_f_indir)Cont_indir_norm[l]=Cont_f_indir[l]*norm_c_indir;
  }
}

/***************************** Classes *******************************/

struct pair_class
**Allocate_pairs_PDB(struct prot_class **prot_class,
		     int *Np_class, int *Nij_class,
		     int N_ij, int  N_U, int *N_L,
		     int **L_bin, int L_min, int L_max, int L_STEP)
{
  //int Naa=20; // No gaps for PDB structures
  *N_L=(L_max-L_min)/L_STEP;
  if(*N_L<=0){
    printf("ERROR, negative number of bins of length %d\n", *N_L);
    exit(8);
  }
  *L_bin=malloc((*N_L)*sizeof(int));
  int L=L_min;
  for(int i_L=0; i_L<(*N_L); i_L++){
    L+=L_STEP; (*L_bin)[i_L]=L;
  }

  //Na2=Naa*Naa;
  int Np=N_U*(*N_L), Npair=N_ij*3;
  *Np_class=Np; *Nij_class=Npair;
  *prot_class=malloc(Np*sizeof(struct prot_class));
  struct pair_class **pair_class=malloc(Np*sizeof(struct pair_class *));
  for(int i_L=0; i_L<(*N_L); i_L++){
    for(int i_U=0; i_U<N_U; i_U++){
      int ip=i_L*N_U+i_U;
      struct prot_class *prot=*prot_class+ip;
      prot->i_L=i_L;
      prot->i_U=i_U;
      Initialize_prot_class(prot);

      struct pair_class **pair_p=pair_class+ip;
      *pair_p=malloc(Npair*sizeof(struct pair_class));
      struct pair_class *pair=*pair_p;
      for(int ij=0; ij<N_ij; ij++){
	for(int C_nat=0; C_nat<3; C_nat++){
	  Allocate_pair(pair, prot);
	  pair->C_nat=C_nat;
	  pair++;
	} // close C_nat
      } // close ij
    } // close i_U
  } // close i_L
  return(pair_class);
}

void  Initialize_prot_class(struct prot_class *prot){
  prot->Np=0;
  prot->L=0;
  prot->Nc=0;
  prot->U_ave=0;
  for(int a=0; a<Naa; a++)prot->U1[a]=0;
  prot->U_norm=0;
  prot->Q_global=malloc(Na2*sizeof(float));
  prot->log_Q_glob=malloc(Na2*sizeof(float));
  prot->P1_i=malloc(Naa*sizeof(float));
  prot->P1_j=malloc(Naa*sizeof(float));
  for(int c=0; c<3; c++)prot->Npair[c]=0;
}

void Allocate_pair(struct pair_class *pair, struct prot_class *prot)
{
  pair->prot_class=prot;
  // Set to zero
  pair->ij=0;
  pair->Cont_freq=0;
  pair->nc_nc=0;
  pair->Cont_freq_nc=0;
  pair->Cont_freq_Nc=0;
  pair->Cont_freq_Nc2=0;
  pair->Cont_freq_indir=0;
  pair->E_cont=0;
  pair->min_dist=-1;
  pair->n=0;

  pair->N2_obs=malloc(Na2*sizeof(double));
  for(int a=0; a<Na2; a++)pair->N2_obs[a]=0;
  pair->P1i_obs=malloc(Naa*sizeof(float));
  pair->P1j_obs=malloc(Naa*sizeof(float));
  pair->Q_obs=malloc(Na2*sizeof(float));
  pair->log_Qobs=malloc(Na2*sizeof(float));
  pair->Q_pred_0=malloc(Na2*sizeof(float));
  pair->Q_pred_1=malloc(Na2*sizeof(float));
  pair->Q_pred_all=malloc(Na2*sizeof(float));
  //pair->P2_null=malloc(Na2*sizeof(float));
}

int Sum_pairs(struct prot_class *prot_c, struct pair_class *pair_c,
	      struct protein *protein, int L_seq,
	      int **C_mat, int **C2, int N_ij, int *ij_bin, int ij_min, 
	      int **label_ij)
{
  int Npairs=0, L=protein->nres;
  prot_c->L+=L; 
  prot_c->Nc+=Nc1L[L];
  float wp=protein->Nc;
  prot_c->U_ave+=wp*protein->U_ave;
  for(int a=0; a<Naa; a++)prot_c->U1[a]+=wp*protein->U1[a];
  prot_c->U_norm+=wp;
  prot_c->Np++;

  // Normalization of contacts
  float Cont_norm[L], Cnc1_norm[L], CNc1_norm[L], CNc2_norm[L],
    Cont_indir_norm[L], nc_nc_norm[L];
  Normalize_cont_freq(Cont_norm, Cont_f,
		      Cnc1_norm, Cnc1,
		      nc_nc_norm, nc_nc,
		      CNc1_norm, Cz1,
		      CNc2_norm, Cz2,
		      Cont_indir_norm, Cont_freq_indir,
		      Nc1L[L], Nc2L[L], Nc3L[L], NcL_indir[L], L);

  int k=0, na=-1; short *i_aa=protein->i_aa;
  struct pair_class *pair;
  for(int i=0; i<L_seq; i++){
    if(i_aa){
      if(Check_AA(i_aa[i], Naa)<0){continue;}else{na=i_aa[i]*Naa;}
    }
    for(int j=i+ij_min; j<L_seq; j++){
      if(i_aa && (Check_AA(i_aa[j], Naa)<0))continue;
      int ij=j-i; if(ij>=L)ij=L-1;
      if(C_mat){
	// PDB structures
	k=Find_ij_class(ij, C_mat[i][j], C2[i][j], N_ij, ij_bin);
	pair=pair_c+k;
	pair->i=0; pair->j=1;
      }else if(label_ij){
	// MSA
	pair=pair_c+k;
	pair->i=i; pair->j=j;
	// Remind to eliminate columns not present in target
	label_ij[i][j]=k;
	label_ij[j][i]=k; k++;
      }else{
	printf("ERROR in set_pair\n"); exit(8);
      }
      pair->n++;
      pair->ij+=ij;
      pair->Cont_freq+=Cont_norm[ij];
      pair->Cont_freq_nc+=Cnc1_norm[ij];
      pair->Cont_freq_Nc+=CNc1_norm[ij];
      pair->Cont_freq_Nc2+=CNc2_norm[ij];
      pair->nc_nc+=nc_nc_norm[ij];
      pair->Cont_freq_indir+=Cont_indir_norm[ij];
      if(i_aa){
	pair->N2_obs[na+i_aa[j]]++;
	pair->E_cont+=E_cont_T[i_aa[i]][i_aa[j]];
      }
    }
  }
  return(Npairs);
}

void Protein_stat(struct prot_class *prot, struct pair_class *pairs,
		  int Nij_class, int ij_min, int L_max, int I_WEIGHT)
{ 
  if(prot->Np==0)return;
  int a;
  if(prot->Np>1){
    prot->L/=prot->Np;
    prot->Nc/=prot->Np;
  }
  prot->U_ave/=prot->U_norm;
  for(a=0; a<Naa; a++)prot->U1[a]/=prot->U_norm;
  if(Naa==21)prot->U1[20]=E_cont_gap[20][20];

  double norm=0; int j;
  struct pair_class *pair=pairs;
  double norm_c=0, norm_nc_nc=0, norm_c_nc=0,
    norm_c_Nc=0, norm_c_Nc2=0, norm_nl=0;
  for(j=0; j<Nij_class; j++){
    prot->Npair[pair->C_nat]+=pair->n;
    if(pair->C_nat==1){pair->min_dist=4.4;}
    else if(pair->C_nat==2){pair->min_dist=6.0;}
    else{pair->min_dist=15;}
    if(I_WEIGHT==0){pair->w=1;}
    else if(I_WEIGHT==1){pair->w=pair->n;}
    else if(I_WEIGHT==2){pair->w=log(pair->n);}
    if(pair->C_nat==2)pair->w=0; // Omit indirect contacts!
    norm+=pair->w;
    //if(pair->n==1){pair++; continue;}

    pair->E_cont/=pair->n;
    pair->ij/=pair->n;
    if(pair->ij < ij_min){
      printf("ERROR, too small <|i-j|>= %.1f n=%.0f\n",pair->ij, pair->n);
      exit(8);
    }else if(pair->ij >= L_max){
      printf("Protein stat |i-j|= %.1f -> %d\n", pair->ij, L_max-1);
      pair->ij=L_max-1;
    }
    float nl= pair->n/prot->Np; // Average number of pairs per protein
    norm_nl+=nl;
    pair->Cont_freq /= pair->n,
    norm_c+=pair->Cont_freq*nl;
    pair->nc_nc /= pair->n;
    norm_nc_nc+=pair->nc_nc*nl;
    pair->Cont_freq_nc /= pair->n;
    norm_c_nc+=pair->Cont_freq_nc*nl;
    pair->Cont_freq_Nc /= pair->n;
    norm_c_Nc+=pair->Cont_freq_Nc*nl;
    pair->Cont_freq_Nc2 /= pair->n;
    norm_c_Nc2+=pair->Cont_freq_Nc2;
    pair->Cont_freq_indir/=pair->n;
    pair++;
  }
  for(j=0; j<Nij_class; j++)pairs[j].w/=norm;

  pair=pairs;
  int L=prot->L;
  norm_c=Nc1L[L]/norm_c;
  norm_nc_nc=2*Nc2L_center[L]/norm_nc_nc;
  norm_c_nc=0.5*(nc2L[L]-4*Nc1L[L]*Nc1L[L]/L)/norm_c_nc;
  norm_c_Nc=Nc2L_center[L]/norm_c_Nc;
  norm_c_Nc2=Nc3L_center[L]/norm_c_Nc2;
  for(j=0; j<Nij_class; j++){
    pair->Cont_freq *= norm_c;
    pair->nc_nc *= norm_nc_nc;
    pair->Cont_freq_nc *=norm_c_nc;
    pair->Cont_freq_Nc *=norm_c_Nc;
    pair->Cont_freq_Nc2 *=norm_c_Nc2;
    Compute_Phi(pair, prot);
    pair++;
  }
}

/****************** Predicted Q *****************************/

void Compute_Phi(struct pair_class *pair, struct prot_class *prot)
{
  pair->Phi0=0; pair->Phi11=0; pair->Phi01=0; pair->Phi2=0;
  if(REM){
    //float Cij_ave=0.10;
    pair->Phi0=-pair->Cont_freq;
    //if(prot->U_ave>=0)return;
    if(REM>=2){
      //if(pair->Cont_freq_Nc > pair->Cont_freq_nc)
      pair->Phi0+=(pair->Cont_freq_Nc-pair->Cont_freq_nc)*prot->U_ave; //
      pair->Phi11=0.5*pair->nc_nc/TEMP; 
      float Cij_2=pair->Cont_freq*(1-pair->Cont_freq);
      pair->Phi01=0.5*(pair->Cont_freq_nc-Cij_2)/TEMP;
      pair->Phi2=0.5*Cij_2/TEMP;
      if(REM>=3){
	//pair->Phi0+=pair->Cont_freq_Nc2*prot->U_ave*prot->U_ave;
	pair->Phi11+=0.5*pair->nc_nc*Nc1L[(int)prot->L]*prot->U_ave/(TEMP*TEMP);
      }
    }
    if(pair->Phi0 > 10){
      int L=prot->L;
      printf("WARNING, very large Phi= %.2f\n", pair->Phi0);
      printf("L= %d <U>= %.2g ij= %.2g ", L, prot->U_ave, pair->ij);
      printf("cij= %.2f\n", pair->Cont_freq);
      printf("<Nc>= %.0f <Cij*Nc>-<Cij><Nc>= %.2f\n",
	     Nc1L[L], pair->Cont_freq_Nc);
      if(REM>=3)
	printf("<Nc^2>= %.0f Cij_Nc2= %.2g\n",
	       Nc2L[L], pair->Cont_freq_Nc2);
    }
  }
}

void Compute_Qprime(float *Qprime, struct pair_class  *pair, int cont,
		    float Lambda)
{
  /*printf("Computing beta, Phi1= %.3g Phi2= %.3g Lambda= %.3f\n",
    Phi1,Phi2,Lambda);*/

  float LPhi0=Lambda*(cont+pair->Phi0), LPhi11=Lambda*pair->Phi11;
  float LPhi01=Lambda*pair->Phi01, LPhi2=Lambda*pair->Phi2;
  float *U1=pair->prot_class->U1;

  for(int a=0; a<Naa; a++){
    int ab=a*Naa; float *q=Qprime+ab; 
    for(int b=0; b<=a; b++){
      double beta=
	E_cont_gap[a][b]*LPhi0;
      if(REM>=2)
	beta+=E_cont_gap[a][b]*(LPhi2*E_cont_gap[a][b]+ LPhi01*(U1[a]+U1[b]))
	  +LPhi11*U1[a]*U1[b];
      *q=exp(-beta); if(b<a){int ba=b*Naa+a; Qprime[ba]=*q;} q++;
    }
  }
}

float Normalize_Qprime(float *Q, struct pair_class  *pair)
{
  if(pair->n==0){
    printf("ERROR, nothing to normalize\n"); exit(8);
  }
  int VBS=0; // Verbose

  // Optimization parameters
  int IT_MAX=400, IT1=IT_MAX/2, iter, nini=0, a, b;
  float EPS=0.0004; // Mean increment = EPS/40= 0.00001
  float error, e_thr=ERR_THR; // Mean error = e_thr/30 = 0.0005
  float e_min=1000, e_min_0, e_ini=e_min;

  // Dynamic variables
  float zi[Naa], zj[Naa], ziopt[Naa], zjopt[Naa], yi[Naa], yj[Naa];
  // Initialization: Normalize Q'
  float *q =Q, *p; double sum=0, z;
  for(a=0; a<Naa; a++){
    z=0; p=pair->P1j_obs; for(b=0; b<Naa; b++){z+=(*q)*(*p); q++; p++;}
    yi[a]=z; sum+=pair->P1i_obs[a]*z;
  }
  for(a=0; a<Naa; a++){
    z=0; q=Q+a; p=pair->P1i_obs;
    for(b=0; b<Naa; b++){z+=(*q)*(*p); q+=Naa; p++;}
    yj[a]=z;
  }
  q=Q; for(a=0; a<Na2; a++){(*q)/=sum; q++;}
  // Initialize zi=1, zj=1
  for(a=0; a<Naa; a++){zi[a]=yi[a]/sum; zj[a]=yj[a]/sum;}

  float d=0;
  for(iter=0; iter<IT_MAX; iter++){

    if(iter<IT1){
      d=Update_z_Newton(zi, zj, yi, yj, Q, pair);
    }else if(iter==IT1){
      e_min_0=e_min;
      for(a=0; a<Naa; a++){zi[a]=1; zj[a]=1;}
      d=Update_z_naive(zi, zj, Q, pair);
    }else{
      d=Update_z_naive(zi, zj, Q, pair);
    }

    // Compute error
    float qnew[Na2];
    Update_q(qnew, Q, zi, zj);
    error=Error_norm(qnew, pair);
    if(VBS)printf("e= %.4g\n", error);   
    if(error<e_min){
      e_min=error; for(a=0; a<Naa; a++){ziopt[a]=zi[a]; zjopt[a]=zj[a];}
    }
    if(error<e_thr)break;
   
    // Change zk if no convergence
    if(d<EPS){
      int it=iter; if(it>=IT1)it-=IT1;
      float x=(float)it/IT1, f=20/x; int k=2*Naa*x;
      if(VBS)
	printf("z converges but Q does not. Changing z %d factor %.2g\n", k,f);
      if(nini==0){e_ini=e_min;} nini++;
      for(a=0; a<Naa; a++){zi[a]=1; zj[a]=1;}
      if(k<Naa){zj[k]*=f;}else{zi[k-Naa]*=f;}
    }

  }
  // Use optimal z always
  for(a=0; a<Naa; a++){zi[a]=ziopt[a]; zj[a]=zjopt[a];}
  
  // Final prediction of Q
  // xi = zi/Pi
  Update_q(Q, Q, zi, zj);

  // Test normalization condition
  error=Error_norm(Q, pair);
  if(error > e_thr){
    printf("ERROR, Q_pred is not well normalized err= %.4g\n", error);
    printf("%d initializations e_ini= %.4g err_Newton= %.4g\n",
	   nini, e_ini, e_min_0);
    printf("prot: L=%.0f U=%.4f pair: ij= %.1f <Cij>= %.2g C_nat= %d  ",
	   pair->prot_class->L, pair->prot_class->U_ave,
	   pair->ij, pair->Cont_freq, pair->C_nat);
    printf("i= %d j= %d\n", pair->i, pair->j);
    float qmax=Q[0], qmin=qmax; int c;
    for(c=1; c<Na2; c++){
      if(Q[c]>qmax){qmax=Q[c];}else if(Q[c]<qmin){qmin=Q[c];}
    }
    printf("Q: Max= %.2g Min= %.2g\n", qmax, qmin);
    if(VBS && (printconv==0)){
      Print_error_Q(pair, Q, zi,zj, ziopt, zjopt);
    }
    printconv++;
  } // end of test 
  return(error);
}

float Update_z_naive(float *zi, float *zj, float *Q, struct pair_class *pair)
{

  /* Q_ij(a,b)= Q'_ij(a,b)xi(a)xj(b)
   sum_b Q'_ij(a,b)xi(a)xj(b)P_j(b)=1 forall i same for j =>
   1/xi(a)=sum_b Q'_ij(a,b)xj(b)Pj(b)
   We fix the scale imposing that sum_a xi(a)/20 = 1

   This normalization is necessary because of the invariance
   xi(a)=k*xi(a), xj(a)=(1/k)*xj(a). We choose the gauge so that <xi>=1
   Normalizing both xi and xj makes dKL and lik much worse
   Not normalizing any of them makes lik and dKL a bit better
   but it makes Lambda less robust.
   The best option seems to be to normalize only xi */

  float dsumi=0, dsumj=0; double z;
  float sumi=0, sumj=0, zinew[Naa], zjnew[Naa], pz[Naa], *p, *q; 
  int a, b;

  // Update zj
  for(a=0; a<Naa; a++)pz[a]=zi[a]*pair->P1i_obs[a];
  for(b=0; b<Naa; b++){
    z=0; q=Q+b; p=pz; for(a=0; a<Naa; a++){z+=(*q)*(*p); q+=Naa; p++;}
    z=1/z; zjnew[b]=z; sumj+=z*pair->P1j_obs[b];
  }
  if(sumj==0){printf("ERROR j in normalize\n");}
  else{
    for(a=0; a<Naa; a++){
      //zjnew[a]/=sumj;
      float d=fabs(zj[a]-zjnew[a]); dsumj+=d; zj[a]=zjnew[a];
    }
  }

  // Update and normalize zi
  for(a=0; a<Naa; a++)pz[a]=zj[a]*pair->P1j_obs[a];
  q=Q;
  for(a=0; a<Naa; a++){
    z=0; p=pz; for(b=0; b<Naa; b++){z+=(*q)*(*p); q++; p++;}
    z=1/z; zinew[a]=z; sumi+=z*pair->P1i_obs[a]; 
  }
  if(sumi==0){printf("ERROR i in normalize\n");}
  else{
    dsumi=0;
    for(a=0; a<Naa; a++){
      zinew[a]/=sumi;
      float d=fabs(zi[a]-zinew[a]); dsumi+=d; zi[a]=zinew[a];
    }
  }
  return(dsumi+dsumj);
}

float Update_z_Newton(float *zi, float *zj, float *yi, float *yj,
		      float *Q, struct pair_class *pair)
{
  /* Q_ij(a,b)= Q'_ij(a,b)xi(a)xj(b)
   F_ia == X_ia sum_b Q'_ia,jb P_jb X_jb  - 1  == X_ia*Y_ia -1
   dF_ia/dX_ic = delta_ac Y_ia
   dF_ia/dX_jb = X_ia Q'_ia,jb P_jb
   0=F_ia(X0)~ F_ia(X)+ dF_ia/dX_ia*(X0_ia-X_ia)+sum_b dF_ia/dX_jb*(X0_jb-X_jb)
   = (X_ia*Y_ia-1)+Y_ia*(X0_ia-X_ia)+X_ia*sum_b Q'_ia,jb P_jb*(X0_jb-X_jb)
   = X0_ia*Y_ia -  X_ia*Y_ia + X_ia*Y_ia(X0_j) -1 =>
   X0_ia = X_ia - (X_ia*Y_ia(X0_j)-1)/Y_ia(X_j) 
  */

  float dsumi=0, dsumj=0; double z; int a, b;
  float sumi=0, sumj=0, zinew[Naa], zjnew[Naa], yinew[Naa], yjnew[Naa];
  float pz[Naa], *p, *q, d; 
  
  // Update zj
  for(a=0; a<Naa; a++)pz[a]=zi[a]*pair->P1i_obs[a];
  for(b=0; b<Naa; b++){
    z=0; q=Q+b; p=pz; for(a=0; a<Naa; a++){z+=(*q)*(*p); q+=Naa; p++;}
    yjnew[b]=z; zjnew[b]=zj[b]-(zj[b]*yjnew[b]-1)/yj[b];
    sumj+=zjnew[b]*pair->P1j_obs[b];
  }
  if(sumj==0){printf("ERROR j in normalize\n");}
  else{
    for(a=0; a<Naa; a++){
      //zjnew[a]/=sumj;
      d=fabs(zj[a]-zjnew[a]); dsumj+=d; zj[a]=zjnew[a]; yj[a]=yjnew[a];
    }
  }
  
  // Update and normalize zi
  for(a=0; a<Naa; a++)pz[a]=zj[a]*pair->P1j_obs[a];
  q=Q;
  for(a=0; a<Naa; a++){
    z=0; p=pz; for(b=0; b<Naa; b++){z+=(*q)*(*p); q++; p++;}
    yinew[a]=z; zinew[a]=zi[a]-(zi[a]*yinew[a]-1)/yi[a];
    sumi+=zinew[a]*pair->P1i_obs[a];
  }
  if(sumi==0){printf("ERROR i in normalize\n");}
  else{
    for(a=0; a<Naa; a++){
      //zinew[a]/=sumi;
      d=fabs(zi[a]-zinew[a]); dsumi+=d; zi[a]=zinew[a]; yi[a]=yinew[a];
    }
  }
  return(dsumi+dsumj);
}

void Update_q(float *qnew, float *Q, float *zi, float *zj)
{
  int a, b; float *q=qnew, *qold=Q; 
  for(a=0; a<Naa; a++){
    for(b=0; b<Naa; b++){
      (*q)=(*qold)*zi[a]*zj[b]; q++; qold++;
    }
  }
}

float Error_norm(float *Q, struct pair_class *pair)
{
  float err_i=0, err_j=0, *q; int a, b;
  for(a=0; a<Naa; a++){
    double P1=0; q=Q+a*Naa; float *p=pair->P1i_obs;
    for(b=0; b<Naa; b++){P1+=(*q)*(*p); q++; p++;}
    P1-=1; err_i+=fabs(P1);
    P1=0; q=Q+a; p=pair->P1j_obs;
    for(b=0; b<Naa; b++){P1+=(*q)*(*p); q+=Naa; p++;}
    P1-=1; err_j+=fabs(P1);
  }
  return(err_i+err_j);
}

void Print_error_Q(struct pair_class *pair, float *Q,
		   float * zi, float * zj, float * ziopt, float * zjopt)
{
  int a;
  printf("Phi0= %.2g Phi01= %.2g Phi11= %.2g Phi2= %.2g\n",
	 pair->Phi0, pair->Phi01, pair->Phi11, pair->Phi2);
  printf("Q_prime= ");
  for(a=0; a<Na2; a++)printf(" %.2g", Q[a]); printf("\n");
  printf("P_i= ");
  for(a=0; a<Naa; a++)printf(" %.3f", pair->P1i_obs[a]); printf("\n");
  printf("P_j= ");
  for(a=0; a<Naa; a++)printf(" %.3f", pair->P1j_obs[a]); printf("\n");
  printf("z_i= ");
  for(a=0; a<Naa; a++)printf(" %.2g", zi[a]); printf("\n");
  printf("z_j= ");
  for(a=0; a<Naa; a++)printf(" %.2g", zj[a]); printf("\n");
  printf("z_i^opt= ");
  for(a=0; a<Naa; a++)printf(" %.2g", ziopt[a]); printf("\n");
  printf("z_j^opt= ");
  for(a=0; a<Naa; a++)printf(" %.2g", zjopt[a]); printf("\n");
  for(a=0; a<Na2; a++)if(isnan(Q[a]))exit(8);
}

/******************** Auxiliary routines *********************/

float **Energy_over_T(float **E_cont_gap, int Naa, float TEMP)
{
  float **E_cont_T=malloc(Naa*sizeof(float *));
  for(int a=0; a<Naa; a++){
    E_cont_T[a]=malloc(Naa*sizeof(float));
    for(int b=0; b<Naa; b++)E_cont_T[a][b]=E_cont_gap[a][b]/TEMP;
  }
  return(E_cont_T);
}

float Average_energy_over_T(float *U1, short *i_aa, int L, float TEMP)
{
  // Computes the average contact energy for a given sequence

  // Normalize contact frequencies: not necessary

  float norm[Naa]; int a;
  for(a=0; a<Naa; a++){U1[a]=0; norm[a]=0;}
  double E_ave=0, Zc=0;
  for(int l=IJ_MIN; l<L; l++){
    float cij=Cont_f[l]; //cij=1;
    for(int i=0; i<(L-l); i++){
      int a=i_aa[i], b=i_aa[i+l];
      float cU=cij*E_cont_gap[a][b];
      E_ave+=cU; Zc+=cij;
      U1[a]+=cU; norm[a]+=cij;
      U1[b]+=cU; norm[b]+=cij;
    }
  }
  //float sT=sqrt(TEMP);
  for(a=0; a<Naa; a++)if(norm[a])U1[a]/=(norm[a]); // *sT
  return(E_ave/(TEMP*Zc));
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


float Corr_coeff_w(float *offset, float *slope,
		   float *xx, float *yy, double *ww, int n){
  double x1=0, x2=0, y1=0, y2=0, xy=0, w1=0;
  int i; float *x=xx, *y=yy; double *w=ww;
  for(i=0; i<n; i++){
    float wx=(*w)*(*x), wy=(*w)*(*y);
    x1 += wx; x2+= wx*(*x);
    y1 += wy; y2+= wy*(*y);
    xy += wx*(*y); w1+=*w;
    x++; y++; w++;
  }

  double nx2=w1*x2-x1*x1, ny2=w1*y2-y1*y1;
  *slope=(w1*xy-y1*x1); if(nx2>0)(*slope)/=nx2;
  *offset=(y1-(*slope)*x1)/w1;
  float r=(*slope); 
  if((nx2>=0)&&(ny2>0)){
    r*=sqrt(nx2/ny2);
  }else{
      printf("ERROR in Corr_coeff, negative variance: ");
      printf("X2= %.2g Y2= %.2g W1= %.2g\n", nx2/w1, ny2/w1, w1);
      return(-100);
  }
  return(r);
}

void Copy_vec(float *v1, float *v2, int N){
  for(int a=0; a<N; a++)v1[a]=v2[a];
}

void Fill_C_mat(int **C, int nres, struct contact *cont_list, int Nc){
  int i, j; struct contact *cl=cont_list;
  for(i=0; i<nres; i++)for(j=0; j<nres; j++)C[i][j]=0;
  for(i=0; i<Nc; i++){
    C[cl->res1][cl->res2]=1; C[cl->res2][cl->res1]=1; cl++; 
  }
}

int Check_AA(short i_aa, int Naa){
  if((i_aa<0)||(i_aa >=Naa)){
    printf("WARNING, wrong amino acid code %d\n", i_aa);
    return(-1);
  }
  return(0);
}

int Find_p_class(int nres, float U_ave,
		 int N_L, int *L_bin, int N_U, float *U_bin)
{
  // Place L and then U
  int i=Find_bin_i(nres, N_L, L_bin);
  int j=Find_bin_f(U_ave,N_U, U_bin);
  return(i*N_U+j);
}

int Find_ij_class(int ij, int C_mat, int C2, int N_ij, int *ij_bin)
{
  // Place ij and then C_nat
  int i=Find_bin_i(ij, N_ij, ij_bin),j;
  if(C_mat){j=1;}
  else if(C2){j=2;}
  else{j=0;}
  return(i*3+j);
}

int Find_bin_i(int L, int N_L, int *L_bin){
  int i; for(i=0; i<N_L; i++)if(L<=L_bin[i])return(i);
  return(N_L-1);
}

int Find_bin_f(float L, int N_L, float *L_bin){
  int i; for(i=0; i<N_L; i++)if(L<=L_bin[i])return(i);
  return(N_L-1);
}

/*********************** Observed Q *****************************/

double Marginalize_N2(float *P1_i, float *P1_j, double *N2_obs)
{  
  double n=0; int a, b; 
  for(a=0; a<Naa; a++){
    double P1=0, *N2=N2_obs+a*Naa;
    for(b=0; b<Naa; b++){P1+=*N2; N2++;}
    P1_i[a]=P1; n+=P1;
    P1=0; N2=N2_obs+a;
    for(b=0; b<Naa; b++){P1+=*N2; N2+=Naa;}
    P1_j[a]=P1;
  }
  if(n==0)return(0);
  float Psumi=0, Psumj=0;
  for(a=0; a<Naa; a++){
    P1_i[a]/=n;
    P1_j[a]/=n;
    Psumi+=P1_i[a];
    Psumj+=P1_j[a];
  }
  if((fabs(Psumi-1.0)>0.001)||(fabs(Psumj-1.0)>0.001)){
    printf("ERROR, wrongly normalized P1: %.4f (i) %.4f (j)\n",
	   Psumi, Psumj); exit(8);
  }
  return(n);
}

void Pairwise_Q(float *Q_obs, float *P1_i, float *P1_j,
		  double *N2_obs, double n)
{
  // n=sum_c N2[c]
  int a, b, ab=0;
  for(a=0; a<Naa; a++){
    float na=P1_i[a]*n;
    for(b=0; b<Naa; b++){
      if((P1_i[a]>0)&&(P1_j[b]>0)){
	Q_obs[ab]=N2_obs[ab]/(na*P1_j[b]);
      }else{
	if(N2_obs[ab]!=0.0){
	  printf("ERROR: Wrong probabilities P1i= %.4f P1j= %.4f",
		 P1_i[a], P1_j[b]);
	  printf(" N2= %.0f\n", N2_obs[ab]); exit(8);
	}
	Q_obs[ab]=1.0;
      }
      if(isnan(Q_obs[ab])){
	printf("ERROR in Pairwise Q, Q[%d]= %.2g\n", ab, Q_obs[ab]);
	exit(8);
      }
      ab++;
    }
  }
}

/************************** Main computation **************************/

float Score(struct pair_class *pair_p, int Nij_class, float Lambda,
	    int LS, int comp, float *f_indir)
{
  float *Q, e, e_thr=ERR_THR;
  //struct prot_class *prot=pair_p->prot_class;

  // Predict Q
  double norm=0; int noconv=0; printconv=0;
  for(int j=0; j<Nij_class; j++){
    struct pair_class *pair=pair_p+j;
    if(pair->n==0)continue;
    
    for(int contact=0; contact<=1; contact++){

      if(contact==1){Q=pair->Q_pred_1;} // Direct contact
      else{Q=pair->Q_pred_0;}

      // Q'(a,b)=Phi1*U(a,b)+Phi2*U(a,b)^2
      Compute_Qprime(Q, pair, contact, Lambda);
      // Multiply times non-specific Q
      for(int c=0; c<Na2; c++)Q[c]*=pair->prot_class->Q_global[c];      
      // Normalization: sum_b P_j(b)Q_ij(a,b)=1. Test normalization
      e=Normalize_Qprime(Q, pair);
      if(e>e_thr){printf("attempted contact= %d\n", contact); noconv++;}
      pair->score.log_lik[contact]=Likelihood_2(Q, pair);
      if((contact==0)||(pair->score.log_lik[1]> pair->score.log_lik[0])){
	pair->error=e;
      }

    } // end contact
    norm+=pair->w;
  } // end ij

  if(noconv)
    printf("WARNING, %d pairs over %d were not properly normalized\n",
	   noconv, 2*Nij_class);
  
  float f, Score;
  if(INDIRECT==0)*f_indir=0;
  if(comp){
    f=*f_indir; Score=Score_indir(pair_p, Nij_class, LS, comp, f)/norm;
    return(Score);
  }

  f=0; Score=Score_indir(pair_p, Nij_class, LS, comp, f)/norm;
  if(INDIRECT==0)return(Score);

  // Optimize f_indir
  int ITER_MAX=20, ITER_MIN=1, iter; 
  float f0=f, S0=Score, S_opt=S0, f_opt=f;
  printf("Optimizing f_indir=%.3g Score=%.5g\n", f, Score);

  f=*f_indir/2; Score=Score_indir(pair_p, Nij_class, LS, comp, f)/norm;
  float f1=f, S1=Score; if(Score>S_opt){S_opt=Score; f_opt=f;}
  printf("Optimizing f_indir=%.2g Score=%.5g\n", f, Score);

  f=*f_indir; Score=Score_indir(pair_p, Nij_class, LS, comp, f)/norm;
  float f2=f, S2=Score; if(Score>S_opt){S_opt=Score; f_opt=f;}
  printf("Optimizing f_indir=%.2g Score=%.5g\n", f, Score);

  float f_high=f*100;
  for(iter=0; iter<ITER_MAX; iter++){
    f=Find_max_quad(f0, f1, f2, S0, S1, S2, 0, f_high);
    Score=Score_indir(pair_p, Nij_class, LS, comp, f)/norm;
    printf("Optimizing f_indir=%.2g Score=%.5g\n", f, Score);
    if(Score > S_opt){f_opt=f; S_opt=Score;}
    else if(iter >= ITER_MIN){break;}
    
    if(f < f0){
      f2=f1; S2=S1; f1=f0; S1=S0; f0=f; S0=Score;
    }else if(f < f1){
      f2=f1; S2=S1; f1=f; S1=Score;
    }else if(f < f2){
      f0=f1; S0=S1; f1=f; S1=Score;
    }else{
      f0=f1; S0=S1; f1=f2; S1=S2; f2=f; S2=Score;
    }
  }
  if(iter==ITER_MAX)printf("WARNING, f_indir optimization did not converge\n");
  //S_opt=Score_indir(pair_p, Naa, Nij_class, comp, f)/norm;
  printf("Optimal f_indir= %.2g Lambda= %.3f Score= %.4g\n",
	 f_opt, Lambda, S_opt);
  *f_indir=f_opt;
  return(S_opt);
}

float Score_indir(struct pair_class *pair_p, int Nij_class,
		  int LS, int comp, float f_indir)
{
  if(f_indir)printf("Computing indirect correlations, f=%.3g\n", f_indir);
  double Score=0; float Qtmp[Na2], *Q, e;
  for(int j=0; j<Nij_class; j++){
    struct pair_class *pair=pair_p+j;
    if(pair->n==0)continue;

    // Indirect contacts
    for(int contact=0; contact<=1; contact++){      
      if(contact){Q=pair->Q_pred_1;}
      else{Q=pair->Q_pred_0;}
      if(f_indir){  //&& (pair->C_nat>=2)
	Copy_vec(Qtmp, Q, Na2); // indirect contact
	if(PDB){
	  e=Indirect_Q_PDB(Qtmp, pair, pair_p, Nij_class, contact, f_indir);
	}else{
	  e=Indirect_Q_ali(Qtmp, j, pair_p, f_indir);
	}
	// Update likelihood
	pair->score.log_lik[contact]=Likelihood_2(Qtmp, pair);
	Q=Qtmp;
      } 
      // Store predicted distributions
      if((contact==0)||
	 (pair->score.log_lik[contact] > pair->score.log_lik_opt)){
	if(f_indir)pair->error=e;
	pair->score.c_pred=contact;
	pair->score.log_lik_opt=pair->score.log_lik[contact];
	Copy_vec(pair->Q_pred_all, Q, Na2);
      }
    } // end contact
    
    // Compute scores
    float log_Qpred[Na2];
    if(comp || LS){
      Compute_logarithm(log_Qpred, pair->Q_pred_all, qmin, lmin, Na2);
      // Kullback-Leibler divergence
      if(comp || (LS==1))pair->score.d_KL=dKL_2(log_Qpred, pair);
      // Relative Mean Error
      if(comp || (LS==2))pair->score.RME=RME_2(log_Qpred, pair);
      // Contact probability
      //pair->score.lik_all=Likelihood_2(Q_all, pair);
      //if(Ncontact>2)c2ij=prot->Cont_freq_indir[(int)pair->ij];
    }
    
    // Score to be optimized:
    if(LS==0){
      Score+=pair->w*pair->score.log_lik_opt;
    }else if(LS==1){
      Score-=pair->w*pair->score.d_KL;
    }else{
      Score-=pair->w*pair->score.RME;
    }
   
  } // End loop on ij

  //if(LS==0)Score/=Lambda;
  return(Score);
}

float **Compute_E_cont_gap(int Naa){
  if((Naa!=20)&&(Naa!=21)){
    printf("ERROR, Naa=%d but instead of 20 a.a. or 21 (20+gap)\n", Naa);
    exit(8);
  }
  // Opening a gap is scored as the minimum energy of an a.a.
  float **E=Allocate_mat2_f(Naa, Naa), E_gap[Naa];
  int a,b;
  for(a=0; a<20; a++){
    E_gap[a]=0;
    for(b=0; b<20; b++){
      E[a][b]=Econt[a][b];
      if(E[a][b]<E_gap[a])E_gap[a]=E[a][b];
    }
  }
  if(Naa==21){
    float E_gg=100;
    for(a=0; a<20; a++){
      E[20][a]=-E_gap[a]; E[a][20]=-E_gap[a];
      if((a==0)||(E_gap[a]>E_gg))E_gg=E_gap[a];
    }
    E[20][20]=-E_gg;
  }
  return(E);
}

/************************* Indirect contacts *******************/

float Indirect_Q_ali(float *Q, int i_ij, struct pair_class *pairs,
		    float f_indir)
{
  if(label_ij==NULL){
    printf("ERROR, no label available for pairs of aligned residues\n");
    exit(8);
  }
  struct pair_class *pair_ij=pairs+i_ij;
  int i=pair_ij->i, j=pair_ij->j;
  double P_ind[Na2], Pk[Na2]; for(int c=0; c<Na2; c++)P_ind[c]=0;
  for(int k=0; k<pair_ij->prot_class->L; k++){
    if((k==i)||(k==j))continue;
    int ik, jk;
    if(k>i){ik=label_ij[i][k];}else{ik=label_ij[k][i];}
    if(k>j){jk=label_ij[j][k];}else{jk=label_ij[k][j];}
    if((ik<0)||(jk<0))continue;
    Indirect_Q_k(Pk, Naa, pair_ij, pairs+ik, pairs+jk);
    for(int c=0; c<Na2; c++)P_ind[c]+=Pk[c];
  }
  float e=Compute_indirect(Q, pair_ij, P_ind, f_indir);
  return(e);
}

float Indirect_Q_PDB(float *Q,
		     struct pair_class *pair_ij,
		     struct pair_class *pairs,
		     int Nij_class, int contact,
		     float f_indir)
{
  // Find most likely intermediate contact
  int N_ij=Nij_class/3, ik, c;

  double P_ind[Na2], Pk[Na2], norm=0;
  for(c=0; c<Na2; c++)P_ind[c]=0;
  for(ik=0; ik<Nij_class; ik++){
    struct pair_class *pair_ik=pairs+ik;
    if(pair_ik->C_nat>1)continue;
    int l_ik, jk; float w;

    // i<=k<=j
    l_ik=pair_ij->ij+pair_ik->ij;
    if(l_ik<pair_ij->prot_class->L){
      jk=Find_ij_class(l_ik, 1, 0, N_ij, ij_bin); // contact
      Indirect_Q_k(Pk, Naa, pair_ij, pairs+ik, pairs+jk);
      w=(pairs[ik].n*pairs[jk].n); norm+=w; //w*=L;
      for(c=0; c<Na2; c++)P_ind[c]+=w*Pk[c];

      jk=Find_ij_class(l_ik, 0, 0, N_ij, ij_bin); // no contact
      Indirect_Q_k(Pk, Naa, pair_ij, pairs+ik, pairs+jk);
      w=(pairs[ik].n*pairs[jk].n); norm+=w; //w*=L;
      for(c=0; c<Na2; c++)P_ind[c]+=w*Pk[c];
      }

    l_ik=pair_ij->ij-pair_ik->ij;
    if(l_ik<IJ_MIN)continue;

    jk=Find_ij_class(l_ik, 1, 0, N_ij, ij_bin); // contact
    Indirect_Q_k(Pk, Naa, pair_ij, pairs+ik, pairs+jk);
    w=(pairs[ik].n*pairs[jk].n); norm+=w; //w*=L;
    for(c=0; c<Na2; c++)P_ind[c]+=w*Pk[c];

    jk=Find_ij_class(l_ik, 0, 0, N_ij, ij_bin); // no contact
    Indirect_Q_k(Pk, Naa, pair_ij, pairs+ik, pairs+jk);
    w=(pairs[ik].n*pairs[jk].n); norm+=w; //w*=L;
    for(c=0; c<Na2; c++)P_ind[c]+=w*Pk[c];
    
  }
  for(c=0; c<Na2; c++)P_ind[c]/=norm;
  float e=Compute_indirect(Q, pair_ij, P_ind, f_indir);
  return(e);
}

float Compute_indirect(float *Q, struct pair_class *pair,
		      double *P, float f_indir)
{
  for(int c=0; c<Na2; c++){
    if(isnan(P[c])){
      printf("ERROR in Compute_indirect! c=%d\n", c); exit(8);
    }
    Q[c]*=(1+P[c]*f_indir);
  }
  float e=Normalize_Qprime(Q, pair);
  return(e);
}

void Indirect_Q_k(double *P_ind, int Naa,
		  struct pair_class *pair_ij,
		  struct pair_class *pair_ik,
		  struct pair_class *pair_jk)
{
  double *q=P_ind;
  float *Qik=Most_lik_Q(pair_ik), *Qjk=Most_lik_Q(pair_jk);
  int a, b, c;

  if((pair_ij->ij>=pair_ik->ij)&&(pair_ij->ij>=pair_jk->ij)){
    // i < k < j
    float *P1k=pair_ik->P1j_obs;
    for(a=0; a<Naa; a++){
      int na=a*Naa;
      for(b=0; b<Naa; b++){
	double Q=0; int ac=na, cb=b;  
	for(c=0; c<Naa; c++){
	  Q+=P1k[c]*(Qik[ac]-1)*(Qjk[cb]-1); ac++; cb+=Naa;
	}
	(*q)=Q; q++;
      }
    }
  }else if((pair_ik->ij<pair_jk->ij)){ // k < i < j
    float *P1k=pair_ik->P1i_obs;
    for(a=0; a<Naa; a++){
      for(b=0; b<Naa; b++){
	double Q=0; int ca=a, cb=b;	  
	for(c=0; c<Naa; c++){
	  Q+=P1k[c]*(Qik[ca]-1)*(Qjk[cb]-1); ca+=Naa; cb+=Naa;
	}
	(*q)=Q; q++;
      }
    }
  }else{ // i < j < k
    float *P1k=pair_ik->P1j_obs;
    for(a=0; a<Naa; a++){
      int na=Naa*a;
      for(b=0; b<Naa; b++){
	double Q=0; int ac=na, bc=b*Naa;	  
	for(c=0; c<Naa; c++){
	  Q+=P1k[c]*(Qik[ac]-1)*(Qjk[bc]-1); ac++; bc++;
	}
	(*q)=Q; q++;
      }
    }
  }

}

float *Most_lik_Q(struct pair_class *pair)
{
  if(pair->score.log_lik[1]>pair->score.log_lik[0]){
    return(pair->Q_pred_1);
  }else{
    return(pair->Q_pred_0); 
  }
  return(NULL);
}

void Indirect_contacts(int **C2, int **k_indirect, int **l_indirect,
		       struct contact **indirect_cont, int *Nc_indirect,
		       int nres, int **C1)
{
  int ic=0;
  int i, j, k;
  for(i=0; i<nres; i++){
    int *Ci=C1[i];
    for(j=i+1; j<nres; j++){
      int Csum=0, *Cj=C1[j], lmax=0, kmax=-1;
      for(k=0; k<nres; k++){
	if((Ci[k])&&(Cj[k])){
	  Csum++;
	  int l1=abs(i-k), l2=abs(j-k);
	  if(l1<l2){if(l1>lmax){lmax=l1; kmax=k;}}
	  else{if(l2>lmax){lmax=l2; kmax=k;}}
	}
      }
      C2[i][j]=Csum; C2[j][i]=Csum;
      if(Csum){
	k_indirect[i][j]=kmax; l_indirect[i][j]=lmax;
	k_indirect[j][i]=kmax; l_indirect[j][i]=lmax;
	if(Ci[j]==0)ic++;
      }
    }
  }
  int Nic=ic; ic=0;
  (*indirect_cont)=malloc(Nic*sizeof(struct contact));
  for(i=0; i<nres; i++){
    for(j=i+1; j<nres; j++){
      if((C2[i][j])&&(C1[i][j]==0)){
	if(ic>=Nic){
	  printf("ERROR, too many indirect contacts\n");
	  printf("ic= %d i= %d j= %d nres= %d\n", ic, i, j, nres);
	  exit(8);
	}
	(*indirect_cont)[ic].res1=i;
	(*indirect_cont)[ic].res2=j; ic++;
      }
    }
  } 
  *Nc_indirect=ic;
}

struct contact *Indirect_contact_list(int *Nc_indirect,
				      struct contact *cont_list,
				      int Nc, int nres)
{
  int N_ind_max=250, i, j;
  short **ind_list=malloc(nres*sizeof(short *));
  short *n_ind=malloc(nres*sizeof(short));
  for(i=0; i<nres; i++){
    n_ind[i]=0; ind_list[i]=malloc(N_ind_max*sizeof(short));
  }
  int dir=0;
  struct contact *c1=cont_list;
  for(int ic=0; ic<Nc; ic++){
    struct contact *c2=c1+1;
    for(int jc=ic+1; jc<Nc; jc++){
      if(c1->res1==c2->res1){
	i=c1->res2; j=c2->res2;
      }else if(c1->res2==c2->res1){
	i=c1->res1; j=c2->res2;
      }else if(c1->res2==c2->res2){
	i=c1->res1; j=c2->res1;
      }else{
	goto next;
      }
      if(Direct_contact(&dir, i, j, cont_list, Nc))goto next;
      ind_list[i][n_ind[i]]=j; n_ind[i]++;
      if(n_ind[i]>N_ind_max){
	printf("ERROR, too many indirect contacts > %d per res\n",n_ind[i]);
	printf("Nc= %d ic= %d jc= %d i=%d\n", Nc, ic, jc, i); exit(8);
      }
    next:
      c2++;
    }
    c1++;
  }
  *Nc_indirect=0; int ic=0;
  for(i=0; i<nres; i++)(*Nc_indirect)+=n_ind[i];
  struct contact *indirect=malloc((*Nc_indirect)*sizeof(struct contact));
  struct contact *ind=indirect;
  for(i=0; i<nres; i++){
    for(j=0; j<n_ind[i]; j++){
      ind->res1=i; ind->res2=ind_list[i][j]; ind++; ic++;
    }
  }
  *Nc_indirect=ic;
  for(i=0; i<nres; i++)free(ind_list[i]);
  free(ind_list); free(n_ind);
  return(indirect);
}

int Direct_contact(int *ic, int i, int j, struct contact *cont_list, int Nc)
{
  struct contact *cont=cont_list+*ic;
  for(int c=*ic; c<Nc; c++){
    if(cont->res1 > i){
      return(0);
    }else if(cont->res1 == i){
      if(cont->res2 > j){return(0);}
      if(cont->res2 == j){*ic=c; return(1);}
    }
    cont++;
  }
  return(0);
}

/************************ Scores *******************************/

float RME_2(float *log_Qpred, struct pair_class *pair)
{
  double d2=0, dd=0; //d1=0, norm=0; 
  for(int c=0; c<Na2; c++){
    if(pair->N2_obs[c]){
      float w=pair->N2_obs[c];
      dd+=w*fabs(pair->log_Qobs[c]-log_Qpred[c]);
      d2+=w*fabs(pair->log_Qobs[c]);
    }
  }
  return(dd/d2);
}

void Compute_logarithm(float *log_Q, float *Q, float qmin, float lmin, int Na2)
{
  float *q=Q, *lq=log_Q;
  for(int c=0; c<Na2; c++){
    if(*q>qmin){*lq=log(*q);}else{*lq=lmin;} q++; lq++;
  }
}

float Likelihood_2(float *Q_pred, struct pair_class *pair)
{
  double logP=0; float l;
  for(int c=0; c<Na2; c++){
    if(Q_pred[c]<qmin){l=lmin;}
    else{l=log(Q_pred[c]);}
    logP+=pair->N2_obs[c]*l;
  }
  return(logP);
}

float dKL_2(float *log_Qpred, struct pair_class *pair)
{
  double d=0;
  for(int c=0; c<Na2; c++){
    if(pair->N2_obs[c]){
      d+=pair->N2_obs[c]*(pair->log_Qobs[c]-log_Qpred[c]);
    }
  }
  return(d/pair->n);
}

float Compute_specific_mut_inf(struct pair_class *pair){
  double M2=0; float *log_Q_glob=pair->prot_class->log_Q_glob;
  for(int c=0; c<Na2; c++){
    if(pair->N2_obs[c])
      M2+=pair->N2_obs[c]*(pair->log_Qobs[c]-log_Q_glob[c]);
  }
  return(M2/pair->n);
}

float Global_mut_inf(struct pair_class *pair, int Na2){
  double M2=0;
  float *log_Q_glob=pair->prot_class->log_Q_glob;
  for(int c=0; c<Na2; c++){
    if(pair->N2_obs[c])M2+=pair->N2_obs[c]*log_Q_glob[c];
  }
  return(M2/pair->n);
}

void Contact_probabilities(struct pair_class *pairs, int Nij_class)
{
  /* float C_EXP=0.5; // P(C_ij=1)=<C_ij>^C_EXP, 0<=C_EXP<=1
  // int APC=1; // Perform Average Pair Correction on d_lik?
  // float SD=S_DEPEND; // Number of dependent sequences */

  // Likelihood difference for all pairs and Cij
  int i; struct pair_class *pair=pairs;
  double P_cont_pred=0, P_cont_obs=0;
  for(i=0; i<Nij_class; i++){
    pair=pairs+i;
    float cij=pair->Cont_freq;
    pair->score.d_log_lik=(pair->score.log_lik[1]-pair->score.log_lik[0]);
    if(pair->score.d_log_lik < 0){
      float fc=cij*exp(pair->score.d_log_lik);
      pair->score.P_C_nat=fc/(fc+(1-cij));
    }else{
      float f=exp(-pair->score.d_log_lik);
      pair->score.P_C_nat=cij/(cij+f*(1-cij));
    }
    P_cont_pred+=pair->n*pair->score.P_C_nat;
    P_cont_obs+=pair->n*pair->Cont_freq;
  }
  printf("Computing contact probabilities.\n");

  // Normalize so that predicted and observed contacts are equal?
  double normP=P_cont_obs/P_cont_pred;
  if((P_cont_pred > 0) && (normP>0)){
    printf("Correction factor: %.2g\n",normP);
    for(i=0; i<Nij_class; i++)pairs[i].score.P_C_nat*=normP;
  }
  // Important: reduce prob if E_cont >0 or Mut_inf<0

  double Mut_inf_1=0, Mut_inf_2=0;
  double E_1=0, E_2=0;
  for(i=0; i<Nij_class; i++){
    pair=pairs+i;
    Mut_inf_1+=pair->score.mut_inf;
    Mut_inf_2+=pair->score.mut_inf*pair->score.mut_inf;
    E_1+=pair->E_cont;
    E_2+=pair->E_cont*pair->E_cont;
  }
  E_1/= Nij_class; E_2-=Nij_class*E_1*E_1;
  E_2=sqrt(E_2/(Nij_class-1));
  Mut_inf_1/= Nij_class; Mut_inf_2-= Nij_class*Mut_inf_1*Mut_inf_1;
  Mut_inf_2=sqrt(Mut_inf_2/(Nij_class-1));
  float E_thr=E_1+E_2, M_thr=Mut_inf_1+Mut_inf_2;

  pair=pairs; 
  for(i=0; i<Nij_class; i++){
    if(pair->E_cont > E_thr)
      pair->score.P_C_nat *= exp(-(pair->E_cont-E_thr)/E_2);
    if(pair->score.mut_inf<M_thr)
      pair->score.P_C_nat *= exp((pair->score.mut_inf-M_thr)/Mut_inf_2);
    if(pair->score.P_C_nat>0.5){pair->score.c_pred=1;}
    else{pair->score.c_pred=0;}
    pair++;
  }
}

void Prot_results(struct prot_class *prot_class, struct pair_class *pairs,
		  int Nij_class, float Lambda, float f_indir, int L_tar)
{
  prot_class->Lambda=Lambda;
  prot_class->f_indir=f_indir;
  //int Na2=Naa*Naa;

  for(int j=0; j<Nij_class; j++){
    struct pair_class *pair=pairs+j;
    pair->score.r=
      Corr_coeff_w(&(pair->score.offset), &(pair->score.slope),
		   pair->Q_pred_all, pair->Q_obs, pair->N2_obs, Na2);
    pair->score.mut_inf=Compute_specific_mut_inf(pair);
  }
  APC_pairs(pairs, Nij_class, L_tar, "mut_inf");
  if(pairs->C_nat<0)return;
  double norm[3];
  for(int c=0; c<3; c++){
    prot_class->r[c]=0;
    prot_class->lik[c]=0;
    prot_class->d_KL[c]=0;
    prot_class->RME[c]=0;
    norm[c]=0;
  }
  for(int j=0; j<Nij_class; j++){
    struct pair_class *pair=pairs+j;
    prot_class->r[pair->C_nat]+=pair->w*pair->score.r;
    prot_class->lik[pair->C_nat]+=pair->w*pair->score.log_lik_opt;
    prot_class->d_KL[pair->C_nat]+=pair->w*pair->score.d_KL;
    prot_class->RME[pair->C_nat]+=pair->w*pair->score.RME;
    norm[pair->C_nat]+=pair->w;
  }
  for(int c=0; c<3; c++){
    if(norm[c]==0)continue;
    prot_class->r[c]/=norm[c];
    prot_class->d_KL[c]/=norm[c];
    prot_class->lik[c]/=norm[c];
    prot_class->RME[c]/=norm[c];
  }
}

/******************* Pairwise statistics ***************************/

void Pairwise_statistics(struct pair_class *pairs, int Npairs)
{
  lmin=log(qmin);

  // Observed Q
  int empty=0;
  printf("Computing observed and specific mutual information\n");
  for(int i=0; i<Npairs; i++){
    struct pair_class *pair=pairs+i;
    if(1)Symmetrize_mat(pair->N2_obs);
    pair->n=Marginalize_N2(pair->P1i_obs, pair->P1j_obs, pair->N2_obs);
    if(pair->n==0){empty++; continue;}
    // Q_ij(a,b)=P_ij(a,b)/P_i(a)P_j(b)
    Pairwise_Q(pair->Q_obs,pair->P1i_obs,pair->P1j_obs,pair->N2_obs,pair->n);
    //Normalize_Qprime(pair->Q_obs, pair); // Not needed
    Compute_logarithm(pair->log_Qobs, pair->Q_obs, qmin, lmin, Na2);
  }
  printf("%d classes with data and %d empty\n",Npairs-empty, empty);
  if(empty==Npairs){printf("No classes left, exiting\n"); exit(8);}
}

void Compute_Q_global(struct prot_class *prot,
		      struct pair_class *pairs, int Npairs)
{
  int a;

  // global Q
  printf("Computing unspecific mutual information Q_global\n");
  double N2_global[Na2];
  for(a=0; a<Na2; a++)N2_global[a]=0;
  struct pair_class *pair;
  for(int i=0; i<Npairs; i++){
    pair=pairs+i;
    // Exclude short range loops: local optimum at 33 only slightly better 
    // if(pair->ij < 33)continue; 
    for(a=0; a<Na2; a++)N2_global[a]+=pair->N2_obs[a];
  }
  // Symmetrize N2_global
  if(1)Symmetrize_mat(N2_global);

  double n=Marginalize_N2(prot->P1_i, prot->P1_j, N2_global);
  Pairwise_Q(prot->Q_global, prot->P1_i, prot->P1_j, N2_global, n);
  // Q_global(a,b)= N2(a,b)/P1(a)P1(b)

}

void Compute_Q_global_all(struct prot_class *prots, int Np_class,
			  struct pair_class **pairs, int Npairs,
			  int iprot, int all_prot)
{
  if((all_prot)&&(iprot))return;

  // global Q
  printf("Computing unspecific mutual information Q_global for all prots\n");
  double N2_global[Na2]; int a;
  for(a=0; a<Na2; a++)N2_global[a]=0;

  for(int ip=0; ip<Np_class; ip++){
    if((all_prot==0)&&(ip!=0))continue;
    struct pair_class *pair=pairs[ip];
      for(int i=0; i<Npairs; i++){
	// Exclude short range loops: local optimum at 33 only slightly better 
	// if(pair->ij < 33)continue; 
	for(a=0; a<Na2; a++)N2_global[a]+=pair->N2_obs[a];
	pair++;
      }
  }
  // Symmetrize N2_global
  if(1)Symmetrize_mat(N2_global);

  struct prot_class *prot_res=prots;
  if(all_prot)prot_res=prots+iprot;
  double n=Marginalize_N2(prot_res->P1_i, prot_res->P1_j, N2_global);
  Pairwise_Q(prot_res->Q_global, prot_res->P1_i, prot_res->P1_j, N2_global, n);
  // Q_global(a,b)= N2(a,b)/P1(a)P1(b)
  if(all_prot){
    for(int ip=1; ip<Np_class; ip++){
      Copy_vec((prots+ip)->Q_global, prot_res->Q_global, Na2);
      Copy_vec((prots+ip)->P1_i, prot_res->P1_i, Naa);
      Copy_vec((prots+ip)->P1_j, prot_res->P1_j, Naa);
    }
  }
}

void Symmetrize_mat(double *X)
{
  for(int a=0; a<Naa; a++){
    int ab=a*Naa, ba=a;
    for(int b=0; b<a; b++){
      float xx=0.5*(X[ab]+X[ba]);
      X[ab]=xx; X[ba]=xx; ab++; ba+=Naa;
    }
  }
}

void Update_Q_global_all(struct prot_class *prots, int Np_class,
			 struct pair_class **pairs, int Npairs,
			 int iprot, int all_prot)
{
  if(all_prot && iprot)return;
  int a;

  // global Q
  printf("Updating unspecific mutual information Q_global\n");
  double Qnew_global[Na2], norm=0;
  for(a=0; a<Na2; a++)Qnew_global[a]=0;

  for(int ip=0; ip<Np_class; ip++){
    if((all_prot==0)&&(ip!=0))continue;
    struct pair_class *pair=pairs[ip];
    for(int i=0; i<Npairs; i++){
      for(a=0; a<Na2; a++)
	if(pair->Q_pred_all[a])
	  Qnew_global[a]+=pair->w*(pair->Q_obs[a]/pair->Q_pred_all[a]);
      norm+=pair->w;
      pair++;
    }
  }

  // Symmetrize N2_global
  Symmetrize_mat(Qnew_global);
  struct prot_class *prot_res=prots;
  if(all_prot)prot_res=prots+iprot;

  for(a=0; a<Na2; a++)prot_res->Q_global[a]*=(Qnew_global[a]/norm);
  struct pair_class pair_prot;
  pair_prot.n=1;
  pair_prot.P1i_obs=prot_res->P1_i;
  pair_prot.P1j_obs=prot_res->P1_j;
  Normalize_Qprime(prot_res->Q_global, &pair_prot);
  if(all_prot){
    for(int ip=1; ip<Np_class; ip++){
      Copy_vec((prots+ip)->Q_global, prot_res->Q_global, Na2);
    }
  }

}

void Pairwise_statistics_old(struct pair_class *pairs, int Npairs,
			     struct prot_class *prot, int L_tar, int PDB)
{
  int i, a, b, c;
  lmin=log(qmin);

 // global Q
  printf("Computing unspecific mutual information Q_global\n");
  double N2_global[Na2];
  for(a=0; a<Na2; a++)N2_global[a]=0;
  struct pair_class *pair;
  for(i=0; i<Npairs; i++){
    pair=pairs+i;
    for(a=0; a<Na2; a++)N2_global[a]+=pair->N2_obs[a];
  }
  // Symmetrize N2_global
  if(1)Symmetrize_mat(N2_global);

  double n=Marginalize_N2(prot->P1_i, prot->P1_j, N2_global);
  Pairwise_Q(prot->Q_global, prot->P1_i, prot->P1_j, N2_global, n);
  // Q_global(a,b)= N2(a,b)/P1(a)P1(b)

  // Observed Q
  int empty=0; float Q_global[Na2], P2_null[Na2];
  printf("Computing observed and specific mutual information\n");
  for(i=0; i<Npairs; i++){
    pair=pairs+i;
    pair->n=Marginalize_N2(pair->P1i_obs, pair->P1j_obs, pair->N2_obs);
    if(pair->n==0){empty++; continue;}
    Copy_vec(Q_global, prot->Q_global, Na2);
    Normalize_Qprime(Q_global, pair);
    c=0; double norm=0;
    for(a=0; a<Naa; a++){
      for(b=0; b<Naa; b++){
	P2_null[c]=Q_global[c]*pair->P1i_obs[a]*pair->P1j_obs[b];
	norm+=P2_null[c]; c++;
      }
    }
    if(norm<=0){
      printf("ERROR norm=%.2g, ij= %.0f\n",norm,pair->ij/pair->n); exit(8);
    }
    // Q specific = P2/P2_null
    for(c=0; c<Na2; c++){
      P2_null[c]/=norm;
      if(P2_null[c] >0){
	pair->Q_obs[c]=pair->N2_obs[c]/(P2_null[c]*pair->n);
      }else{
	pair->Q_obs[c]=1.0;
      }
    }
    Compute_logarithm(pair->log_Qobs, pair->Q_obs, qmin, lmin, Na2);
    pair->score.mut_inf=Compute_specific_mut_inf(pair);
  }
  printf("%d classes with data and %d empty\n",Npairs-empty, empty);

  if(L_tar<=0)return;
  APC_pairs(pairs, Npairs, L_tar, "mut_inf");
}

void APC_pairs(struct pair_class *pairs, int Npairs, int L, char *what)
{
  // Compute Average Product Correction (APC) to mutual information or P_cont
  printf("Computing Average Product Correction to %s\n", what);
  int i_what, i;
  if(strncmp(what, "mut_inf", 7)==0){i_what=0;}
  else if(strncmp(what, "d_lik", 5)==0){i_what=1;}
  else{printf("ERROR in APC, undefined %s\n", what); exit(8);}

  double average[L], norm[L];
  for(i=0; i<L; i++){
    average[i]=0;
    norm[i]=0;
  } 
  struct pair_class *pair=pairs; float wy;
  for(i=0; i<Npairs; i++){
    if(i_what==0){
      wy=pair->w*pair->score.mut_inf;
    }else{
      wy=pair->w*exp(pair->score.d_log_lik);
    }
    average[pair->i]+=wy; norm[pair->i]+=pair->w;
    average[pair->j]+=wy; norm[pair->j]+=pair->w;
    pair++;
  }

  double ave_all=0, norm_all=0;
  for(i=0; i<L; i++){
    if(norm[i]){
      ave_all+=average[i];
      norm_all+=norm[i];
      average[i]/=norm[i];
      if(i_what)average[i]=log(average[i]);
    }
  }
  if(norm_all)ave_all/=norm_all;
  if(i_what)ave_all=log(ave_all);
  for(i=0; i<Npairs; i++){
    pair=pairs+i;
    float APC=average[pair->i]*average[pair->j];
    if(ave_all)APC/=ave_all;
    if(i_what==0){
      pair->score.mut_inf-=APC;
    }else{
      pair->score.d_log_lik-=APC;
    }
  }
}

void Print_Cont_stat(int L, int ij_min)
{
  char name_out[100]; sprintf(name_out, "Cont_stat_%d.dat", L);
  FILE *file_out=fopen(name_out, "w");
  printf("Writing contact statistics in file %s\n", name_out);

  float Cont_norm[L], nc_nc_norm[L], CNc1_norm[L], CNc2_norm[L], Cnc1_norm[L];
  Normalize_cont_freq(Cont_norm, Cont_f,
		      Cnc1_norm, Cnc1,
		      nc_nc_norm, nc_nc,
		      CNc1_norm, Cz1,
		      CNc2_norm, Cz2,
		      NULL, NULL,
		      Nc1L[L], Nc2L[L], Nc3L[L], 0, L);

  fprintf(file_out, "#|i-j| 2=<Cij> 3=<Cij*Nc>-<Cij><Nc> 4=<ni*nj>-<ni*nj>");
  fprintf(file_out, " 5=<Cij*Nc^2>-<Cij><Nc^2>-2(<Cij*Nc>-<Cij><Nc>)");
  fprintf(file_out, " 6=<Cij*ni>-<Cij><ni>\n");
  for(int l=ij_min; l<L; l++){
    fprintf(file_out, "%d\t%.4f", l, Cont_norm[l]);
    fprintf(file_out, "\t%.4g", CNc1_norm[l]);
    fprintf(file_out, "\t%.4g", nc_nc_norm[l]);
    fprintf(file_out, "\t%.4g", CNc2_norm[l]);
    fprintf(file_out, "\t%.4g", Cnc1_norm[l]);
    fprintf(file_out, "\n");
  }
  fclose(file_out);
}


