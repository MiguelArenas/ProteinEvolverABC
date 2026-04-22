#include "REM.h"
//#include "energy_BKV.h"  // Contact energy parameters
#include "protein3.h"
#include "allocate.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "meanfield.h"
extern float Entropy(float *P, int n);

#define IT_MAX 25
#define EPS 0.0000001
double norm_lik=0;

float COEFF=1;
float **Econt_T=NULL, T_ratio=1; // Not used anymore


int INI_MF=0;

static void Average_ene(float **u1_ja, float **u2_ja, float **u3_ja,
			float **P_MF_ia, int L, int Naa);
void  Sum_ene(float **psi1_ja, float **psi2_ja,
	      float **Upsi1_ja, float **Upsi2_ja,
	      float **U2psi1_ja, float **u1U_ja,
	      float **u1_ja, float **u2_ja,
	      float **P_MF_ia, double *Cont_freq, 
	      int L, int Naa);
static double Compute_DG_overT(struct REM *E, double *T_eff, double *factor,
			       float **u1_ja, float **u2_ja, float **u3_ja,
			       float **psi1_ja,float **psi2_ja,
			       float **P_MF_ia, int L,
			       int **C_nat, int *i_sec, double *Cont_freq);
static void Update_P(float **P_new_ia, float *P_mut_a,
		     float **u1_ja, float **u2_ja, float **u3_ja,
		     float **psi1_ja, float **psi2_ja,
		     float **Upsi1_ja, float **Upsi2_ja,
		     float **U2psi1_ja, float **u1U_ja,
		     struct REM *E, int L, int Naa,
		     int **C_nat, int *i_sec, double *Cont_freq,
		     double T_eff, double factor, float Lambda);
static void Update_P_nat(float **P_new_ia, float *P_mut_a, float **u1_ja,
			 int L, int Naa, int **C_nat, int *i_sec, float Lambda);

static float Test_convergence(float **P1_ia, float **P2_ia, int L, int Naa);
void Boltzmann_distr(float *P, float beta, float *hydro, int N);
void Initialize_P(float **P_MF_ia, int **C_nat, int *i_sec, double *Cont_freq,
		  float *hydro, int L, float Lambda);
static void Allocate_mean_field(float ***P_new_ia, float ***P_opt_ia,
				float ***u1_ja, float ***u2_ja, float ***u3_ja,
				float ***psi1_ja, float ***psi2_ja,
				float ***Upsi1_ja, float ***Upsi2_ja,
				float ***U2psi1_ja, float ***u1U_ja, int L);
static void Empty_mean_field(float **P_new_ia, float **P_opt_ia,
			     float **u1_ja, float **u2_ja, float **u3_ja,
			     float **psi1_ja, float **psi2_ja,
			     float **Upsi1_ja, float **Upsi2_ja,
			     float **U2psi1_ja, float **u1U_ja, int L);

/*****************************************************************************/

int meanfield(float **P_MF_ia, struct MF_results *res,
	      //double *DG_ave, float *Tfreeze, float *Kul_Leib
	      float *P_mut_a, float **P_ini_ia,
	      int **C_nat, int *i_sec, int L,
	      float Lambda, struct REM E_wt, int ALL) //Npop
{
  // WARNING: E_wt must be initialized before
  // Initialization
  // Initialize_REM(&S_C, &S_U, &LEN, L, file_str, sC0, sC1, sU1);

  printf("Mean field distribution, Lambda=%.3f\n", Lambda);

  // Allocate
  float **P_new_ia, **P_opt_ia, **u1_ja, **u2_ja, **u3_ja,
    **psi1_ja, **psi2_ja, **Upsi1_ja, **Upsi2_ja, **U2psi1_ja, **u1U_ja;
  Allocate_mean_field(&P_new_ia, &P_opt_ia, &u1_ja, &u2_ja, &u3_ja,
		      &psi1_ja, &psi2_ja, &Upsi1_ja, &Upsi2_ja,
		      &U2psi1_ja, &u1U_ja, L);

  // Initialize distribution
  if(P_ini_ia){
    Copy_P(P_MF_ia, P_ini_ia, L, 20);
  }else{
    Initialize_P(P_MF_ia, C_nat, i_sec, Cont_freq, hydro, L, Lambda);
  }

  struct REM E;
  Initialize_E_REM(&E, E_wt.L, E_wt.REM, E_wt.T, E_wt.S_C, E_wt.S_U,
		   FILE_STR, 0);
  double DG, T_eff, DG_ini=0;


  // Iteration
  int convergence=0, iter, i; double factor=1;
  float conv, free_e, score, score_min=10000000, DGopt=0, KLD=0;
  float KLD_opt, Tf_opt; //T_old
  printf("#iteration DeltaG Tfreeze factor Kul-Leib ");
  if(ALL){printf("KL+Lambda*DG");}else{printf("KL+Lambda*E_nat");}
  printf(" convergence (Lambda=%.2f T=%.2f SC=%.2f SU=%.2f)\n",
	 Lambda, E_wt.T, E_wt.S_C, E_wt.S_U);
  for(iter=0; iter<IT_MAX; iter++){
    Average_ene(u1_ja, u2_ja, u3_ja, P_MF_ia, L, 20);
    Sum_ene(psi1_ja, psi2_ja, Upsi1_ja, Upsi2_ja, U2psi1_ja, u1U_ja,
	    u1_ja, u2_ja, P_MF_ia, Cont_freq, L, 20);
    DG=Compute_DG_overT(&E, &T_eff, &factor, u1_ja, u2_ja, u3_ja,
			psi1_ja, psi2_ja, P_MF_ia, L, C_nat, i_sec,
			Cont_freq);
    //if(iter==0)T_old=T_eff;
    //T_new=(T_old+T_eff)/2; T_old=T_eff;
    float T_new=T_eff;

    if(ALL){
      if(isnan(factor)){printf("WARNING, factor is nan\n"); break;}
      Update_P(P_new_ia, P_mut_a, u1_ja, u2_ja, u3_ja, psi1_ja, psi2_ja,
	       Upsi1_ja, Upsi2_ja, U2psi1_ja, u1U_ja, &E, 
	       L, 20, C_nat, i_sec, Cont_freq, T_new, factor, Lambda);
    }else{
      Update_P_nat(P_new_ia, P_mut_a, u1_ja, L, 20, C_nat, i_sec, Lambda);
    }
    conv=Test_convergence(P_new_ia, P_MF_ia, L, 20);
    KLD=0; for(i=0; i<L; i++)KLD+=KL(P_MF_ia[i], P_mut_a, 20);
    if(ALL){free_e = KLD+COEFF*Lambda*DG;}
    else{free_e = KLD+Lambda*E.E_nat;}
    score=conv; 
    //score= free_e;
    if(score < score_min){
      score_min = score;
      DGopt=DG; KLD_opt=KLD; Tf_opt=E.Tf;
      Copy_P(P_opt_ia, P_MF_ia, L, 20);
    }

    Copy_P(P_MF_ia, P_new_ia, L, 20);
    printf("%2d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2g\n",
	   iter, DG, E.Tf, factor, KLD, free_e, sqrt(conv));
    if((conv<EPS)&&(fabs(DG-DG_ini)<0.01)){convergence=1; break;}
    DG_ini=DG;
  }
  if(convergence==0){
    Copy_P(P_MF_ia, P_opt_ia, L, 20);
    DG=DGopt; KLD=KLD_opt; E.Tf=Tf_opt;
    printf("WARNING, mean field computation did not converge in %d", iter+1);
    printf(" steps.\nUsing distribution with minimum ");
    //printf("Kullback-Leiber div. + Lambda*DG= %.3g\n", score_min);
    printf("convergence error %.3g\n", score_min);
    printf("DG= %.3f\n", DG);
  }
  res->DG=DG;
  res->KL_mut=KLD/L;
  res->Tf=E.Tf;
  res->Lambda[0]=Lambda;
  res->h=Hydro_ave(P_MF_ia, hydro, L);

  Empty_mean_field(P_new_ia, P_opt_ia, u1_ja, u2_ja, u3_ja, psi1_ja, psi2_ja,
		   Upsi1_ja, Upsi2_ja, U2psi1_ja, u1U_ja, L);
  Empty_E_REM(&E); 
  return(convergence);
}



void Average_ene(float **u1_ja, float **u2_ja, float **u3_ja,
		 float **P_MF_ia, int L, int Naa)
{
  int i, a, b;
  for(i=0; i<L; i++){
    float *p=P_MF_ia[i];
    float *u1j=u1_ja[i], *u2j=u2_ja[i], *u3j=u3_ja[i];
    for(a=0; a<Naa; a++){
      float *U=Econt_T[a], x;
      u1j[a]=0; u2j[a]=0; u3j[a]=0; 
      for(b=0; b<Naa; b++){
	x=U[b]*p[b]; u1j[a]+=x;
	x*=U[b]; u2j[a]+=x;
	x*=U[b]; u3j[a]+=x;
      }
    }
  }
}

void  Sum_ene(float **psi1_ja, float **psi2_ja,
	      float **Upsi1_ja, float **Upsi2_ja,
	      float **U2psi1_ja, float **u1U_ja, 
	      float **u1_ja, float **u2_ja,
	      float **P_MF_ia, double *Cont_freq,
	      int L, int Naa)
{
  int i, j, a, b;
  for(i=0; i<L; i++){
    float *psi1=psi1_ja[i], *psi2=psi2_ja[i];
    for(a=0; a<Naa; a++){psi1[a]=0; psi2[a]=0;}
    for(j=0; j<L; j++){
      int l=abs(i-j); if(l<IJ_MIN)continue;
      float c=Cont_freq[l];
      for(a=0; a<Naa; a++){
	psi1[a]+= c*u1_ja[j][a];
	psi2[a]+= c*u2_ja[j][a];
      }
    }
  }

  for(i=0; i<L; i++){
    float *Upsi1=Upsi1_ja[i];
    float *Upsi2=Upsi2_ja[i], *U2psi1=U2psi1_ja[i];
    float *psi1=psi1_ja[i], *psi2=psi2_ja[i], *u1U=u1U_ja[i];
    float *p=P_MF_ia[i];
    for(a=0; a<Naa; a++){
      Upsi1[a]=0;  Upsi2[a]=0; U2psi1[a]=0; u1U[a]=0;
      float *U=Econt_T[a];
      for(b=0; b<Naa; b++){
	float pU=p[b]*U[b];
	Upsi1[a] += pU*psi1[b];
	Upsi2[a] += pU*psi2[b];
	U2psi1[a]+= pU*psi1[b]*U[b];
	u1U[a] +=  pU*u1_ja[i][b];
      }
    }
  }
}

double Compute_DG_overT(struct REM *E, double *T_eff, double *factor,
			float **u1_ja, float **u2_ja, float **u3_ja,
			float **psi1_ja, float **psi2_ja,
			float **P_MF_ia,
			int L, int **C_nat, int *i_sec,
			double *Cont_freq)
{
  Set_to_zero(E);
  double E12_sum=0; 
  double u12_c2[20];
  int i, j, a;
  for(i=0; i<L; i++){
    float *p=P_MF_ia[i];
    if(SEC_STR){
      float *El=E_loc_over_T[i_sec[i]];
      for(a=0; a<20; a++)E->E_loc+=p[a]*El[a];
    }
    double E1ii=0;
    for(a=0; a<20; a++)u12_c2[a]=0;
    for(j=0; j<L; j++){
      int l=abs(i-j); if(l<IJ_MIN)continue;
      double U1ij=0; float *u1j=u1_ja[j];
      for(a=0; a<20; a++)U1ij+=p[a]*u1j[a];
      if((j>i)&&(C_nat[i][j]))E->E_nat+=U1ij;
      if(E->REM==0)continue;
      float c=Cont_freq[l];
      E1ii +=c*U1ij;
      if(E->REM > 1){
	float c2=c*c;
	for(a=0; a<20; a++)u12_c2[a]+=c2*u1j[a]*u1j[a];
	if(j>i){    
	  double U2ij=0; float *u2j=u2_ja[j];
	  for(a=0; a<20; a++)U2ij+=p[a]*u2j[a];
	  E->E22 += C221_ij[l]*U2ij;
	  if(E->REM==3){
	    double U3ij=0; float *u3j=u3_ja[j];
	    for(a=0; a<20; a++)U3ij+=p[a]*u3j[a];
	    E->E321+= C321_ij[l]*U3ij;
	    E->e342+= C342_ij[l]*U2ij;
	    E->c332U2[i]+= C332_ij[i][j]*U2ij;
	    E->c332U2[j]+= C332_ij[j][i]*U2ij;
	  }
	}
      }
    }
    if(E->REM==0)continue;
    E1ii/=2;
    E->E1 += E1ii;
    E->c1U1[i]=E1ii;
    if(E->REM > 1){
      E12_sum += E1ii*E1ii;
      //E12_sum += E1ii*E->c1U1[i];
      double E11=0; float *psi1=psi1_ja[i];
      for(a=0; a<20; a++)E11+=p[a]*(psi1[a]*psi1[a]-u12_c2[a]);
      E->E23  += C232_i[i]*E11;
      /*if(E->REM==3){
	double E12=0; float *psi2=psi2_ja[i];
	for(a=0; a<20; a++)E12+=p[a]*psi1[a]*psi2[a];
	E332 += M332_i[i]*E12;
	}*/
    }    
  }

  /*if(E->REM > 1){E->E24=(C242)*E->E1*E->E1-E12_sum;}
    else{E->E24=0;}*/
  E->G_misf=G_misfold(E);
  if(isnan(E->G_misf)){
    printf("WARNING, G_misf= %.2g Tf= %.2g E1= %.2g E2= %.2g ",
	   E->G_misf, E->Tf, E->E1, E->E2);
    printf("E32= %.2g E342= %.2g E332= %.2g E363= %.2g\n",
	   E->E321, E->E342, E->E332, E->E363);
    for(i=0; i<10; i++){
      printf("i= %d P= ", i);
      for(a=0; a<20; a++)printf(" %.2f", P_MF_ia[i][a]);
      printf("\n");
    }
  }
  float DeltaG_overT=DeltaG((float)E->E_nat, E->G_misf, E->S_U);
  if(E->REM==1)E->Tf=0;
  if(E->Tf>1){*T_eff=E->Tf;}else{*T_eff=1;}
  E->G_misf+=E->S_U;
  if(E->G_misf < -30){*factor=1;}
  else if(E->G_misf > 30){*factor=0;}
  else{*factor=1./(1+exp(E->G_misf));}
  return(DeltaG_overT);
}

void Update_P_nat(float **P_new_ia, float *P_mut_a, float **u1_ja,
		  int L, int Naa, int **C_nat, int *i_sec, float Lambda)
{
  int i, j, a;
  //float x=exp(DG), Lambda=(x/(1+x))*Npop;
  //printf("lambda= %.2g\n", Lambda);
  for(i=0; i<L; i++){
    float *p=P_new_ia[i]; double norm=0;
    for(a=0; a<Naa; a++){
      double h=0;
      if(SEC_STR)h=E_loc_over_T[i_sec[i]][a];
      for(j=0; j<L; j++){
	if(C_nat[i][j])h+=u1_ja[j][a];
      }
      p[a]=P_mut_a[a]*exp(-Lambda*h);
      norm+=p[a];
    }
    for(a=0; a<Naa; a++)p[a]/=norm;
  }
}

void Update_P(float **P_new_ia, float *P_mut_a,
	      float **u1_ja, float **u2_ja, float **u3_ja,
	      float **psi1_ja, float **psi2_ja,
	      float **Upsi1_ja, float **Upsi2_ja,
	      float **U2psi1_ja, float **u1U_ja,
	      struct REM *E, int L, int Naa, int **C_nat, int *i_sec,
	      double *Cont_freq, double T_ratio, double factor,
	      float Lambda)
{
  int i, j, a;
  //float x=exp(DG), Lambda=(x/(1+x))*Npop;
  //printf("lambda= %.2g\n", Lambda);
  float T2=T_ratio*T_ratio;
  float dE3X3=3*C363*E->E1*E->E1;
  float h_a[20];

  for(i=0; i<L; i++){
    float h_min=1000;
    for(a=0; a<Naa; a++){
      double hnat=0;
      if(SEC_STR)hnat=E_loc_over_T[i_sec[i]][a];
      float h22=0, h23=0, h24_j=0;
      double h32=0, h34=0, u12_c2=0;
      for(j=0; j<L; j++){
	if(C_nat[i][j])hnat+=u1_ja[j][a];
	if(E->REM==0)continue;
	int l=i-j; if(l<0)l=-l; if(l<IJ_MIN)continue;
	float c=Cont_freq[l];
	if(E->REM > 1){
	  h22 += C221_ij[l]*u2_ja[j][a];
	  if(j<LEN_MAX)
	    h23 += C232_i[j]*4*c*(Upsi1_ja[j][a]-c*u1U_ja[j][a]);
	  float x=c*u1_ja[j][a]; u12_c2 += x*x;
	  h24_j += x*E->c1U1[j];
	  if(E->REM==3){
	    h32 += C321_ij[l]*u3_ja[j][a];
	    // h33 += C332_i[j]*c*(Upsi2_ja[j][a]+U2psi1_ja[j][a]);
	    h34 += C342_ij[l]*u2_ja[j][a];
	  }
	}
      }
      double h = hnat;
      if(E->REM){
	float psi=psi1_ja[i][a];
	double h1=-psi;
	if(E->REM > 1){
	  double h24=C242*((E->E1-E->c1U1[i])*psi-h24_j);
	  h1+=(h22+h23+C232_i[i]*(psi*psi-u12_c2)+h24)/T_ratio; // / or *?
	  if(E->REM==3){
	    h1 -= (h32+h34*E->E1+(E->E342+dE3X3)*psi)/T2;  // / or *?
	    //+h33+C332_i[i]*psi1_ja[i][a]*psi2_ja[i][a]
	  }
	}
	h += factor*h1;
      }
      h_a[a]=h;
      if((a==0)||(h< h_min))h_min=h;
    }
    double norm=0;
    float *p=P_new_ia[i];
    for(a=0; a<Naa; a++){
      double h=h_a[a]-h_min;
      p[a]=P_mut_a[a]*exp(-Lambda*h);
      norm+=p[a];
      if(isnan(p[a])){
	printf("ERROR, P_MF vanishes!\n");
	printf("SITE %d out of %d, amm %d local field %.2g ", i,L,a,h);
	printf(" norm= %.2g f=%.2g\n", norm, factor);
	exit(8);
      }
    }
    if(norm <= 0){
      printf("ERROR, P_MF vanishes!\n");
      printf(" norm= %.2g f=%.2g\n", norm, factor);
      exit(8);
    }
    for(a=0; a<Naa; a++)p[a]/=norm;
  }
}


float Test_convergence(float **P1_ia, float **P2_ia, int L, int Naa)
{
  double delta=0; int i, a;
  for(i=0; i<L; i++){
    float *p1=P1_ia[i], *p2=P2_ia[i], d;
    for(a=0; a<Naa; a++){
      d=p1[a]-p2[a]; delta+=d*d;
    }
  }
  return(delta/(Naa*L));
}

void Copy_P(float **P1_ia, float **P2_ia, int L, int Naa)
{
  int i, a;
  for(i=0; i<L; i++){
    float *p1=P1_ia[i], *p2=P2_ia[i];
    for(a=0; a<Naa; a++)p1[a]=p2[a];
  }
}
 
void Divide_by_Pmut(float **P_MF_ia, float *P_mut_a, int L, int Naa){
  int i, a;
  for(i=0; i<L; i++){
    for(a=0; a<Naa; a++)P_MF_ia[i][a]/=P_mut_a[a];
  }
}

void Initialize_P(float **P_MF_ia, int **C_nat, int *i_sec,
		  double *Cont_freq, float *hydro, int L, float Lambda)
{
  printf("Initializing amino acid distribution as exp(L*(c_i-cave)/s*h[a])\n");
  // Compute average native contact matrix
  float csum[L];
  double ave=0; int i, j;
  for(i=0; i<L; i++){
    float c=0;
    for(j=0; j<L; j++){
      int l=abs(i-j); if(l<IJ_MIN)continue;
      if(C_nat[i][j])c++;
      c-=Cont_freq[l];
    }
    csum[i]=c; ave+=c;
  }
  ave/=L;

  // Compute BETA
  for(i=0; i<L; i++){
    float beta=Lambda*(csum[i]-ave);
    Boltzmann_distr(P_MF_ia[i], beta, hydro, 20);
    if(SEC_STR){
      float *El=E_loc_over_T[i_sec[i]], *p=P_MF_ia[i];
      double sum=0; int a;
      for(a=0; a<20; a++){
	p[a]*=exp(-Lambda*El[a]); sum+=p[a];
      }
      for(a=0; a<20; a++)p[a]/=sum;
    }
  }

}

void Boltzmann_distr(float *P, float beta, float *hydro, int N)
{
  double sum=0; int a;
  for(a=0; a<N; a++){
    P[a]=exp(beta*hydro[a]);
    sum+=P[a];
  }
  for(a=0; a<N; a++)P[a]/=sum;
}

void Test_distr(//float *DG, float *lik, float *Tf,
		struct MF_results *MF_res, float **P_MF_ia,
		float **f_reg_ia, float **f_msa_ia, float *wi, int L,
		int **C_nat, int *i_sec, struct REM E_wt)
{
  double T_ratio, factor;
  float **P_new_ia, **P_opt_ia, **u1_ja, **u2_ja, **u3_ja,
    **psi1_ja, **psi2_ja, **Upsi1_ja, **Upsi2_ja, **U2psi1_ja, **u1U_ja;
  Allocate_mean_field(&P_new_ia, &P_opt_ia, &u1_ja, &u2_ja, &u3_ja,
		      &psi1_ja, &psi2_ja, &Upsi1_ja,
		      &Upsi2_ja, &U2psi1_ja, &u1U_ja, L);
  Average_ene(u1_ja, u2_ja, u3_ja, P_MF_ia, L, 20);
  Sum_ene(psi1_ja, psi2_ja, Upsi1_ja, Upsi2_ja, U2psi1_ja, u1U_ja,
	  u1_ja, u2_ja, P_MF_ia, Cont_freq, L, 20);
  struct REM E;
  Initialize_E_REM(&E, E_wt.L, E_wt.REM, E_wt.T, E_wt.S_C, E_wt.S_U,
		   FILE_STR, 0);
  MF_res->DG=Compute_DG_overT(&E, &T_ratio, &factor, u1_ja, u2_ja, u3_ja,
			     psi1_ja, psi2_ja, P_MF_ia, L, C_nat, i_sec,
			     Cont_freq);
  MF_res->Tf=E.Tf;
 
  Empty_mean_field(P_new_ia, P_opt_ia, u1_ja, u2_ja, u3_ja, psi1_ja, psi2_ja,
		   Upsi1_ja, Upsi2_ja, U2psi1_ja, u1U_ja, L);
  Empty_E_REM(&E);
  MF_res->h=Hydro_ave(P_MF_ia, hydro, L);
  
  //MF_res->entropy=Average_entropy(P_MF_ia, L);
  //MF_res->lik=Compute_lik(P_MF_ia, L, f_reg_ia);
  Compute_score(MF_res, P_MF_ia, L, f_reg_ia, f_msa_ia, wi);

}

void Compute_score(struct MF_results *res,
		   float **P_MF_ia, int L,
		   float **f_reg_ia, float **f_msa_ia, float *wi)
{
  if(norm_lik==0){
    for(int i=0; i<L; i++)norm_lik+=wi[i];
  }

  Entr_reg=0;
  double log_lik=0, KL_reg=0, entropy=0;
  for(int i=0; i<L; i++){
    if(wi[i]==0)continue;
    Entr_reg+=wi[i]*Entropy(f_reg_ia[i], 20);
    float ll_i=0, ll_mod_i=0, ent_i=0;
    for(int a=0; a<20; a++){
      float l_mod=log(P_MF_ia[i][a]);
      ent_i+=P_MF_ia[i][a]*l_mod;  // Model entropy
      ll_i+=f_reg_ia[i][a]*l_mod;  // Likelihood of data
      ll_mod_i+=P_MF_ia[i][a]*log(f_reg_ia[i][a]); // Likelihood of model
    }
    entropy-=wi[i]*ent_i;
    KL_reg+=wi[i]*(ent_i - ll_mod_i);
    log_lik += wi[i]*ll_i;
  }
  Entr_reg/=norm_lik;
  res->entropy=entropy/norm_lik;
  res->lik_reg=log_lik/norm_lik;
  res->KL_mod=-Entr_reg-res->lik_reg;
  res->KL_reg=KL_reg/norm_lik;
  res->lik_MSA=Compute_lik(P_MF_ia, L, f_msa_ia);
  //res->score=-res->KL_reg;
  //res->score=res->lik_reg;
  res->score=-res->KL_mod-res->KL_reg;
}

float Compute_lik(float **P_MF_ia, int L, float **freq_ia)
{
  if(norm_lik==0){
    for(int i=0; i<L; i++){
      for(int a=0; a<20; a++)norm_lik+=freq_ia[i][a];
    }
  }

  double log_lik=0; int i; float eps=0.00000001;
  for(i=0; i<L; i++){
    for(int a=0; a<20; a++){
      if(freq_ia[i][a]==0)continue;
      float p=P_MF_ia[i][a];
      if(p<eps)p=eps;
      log_lik+=freq_ia[i][a]*log(p);
    }
  }
  return(log_lik/norm_lik);
}

void Allocate_mean_field(float ***P_new_ia, float ***P_opt_ia,
			 float ***u1_ja, float ***u2_ja, float ***u3_ja,
			 float ***psi1_ja, float ***psi2_ja,
			 float ***Upsi1_ja, float ***Upsi2_ja,
			 float ***U2psi1_ja, float ***u1U_ja, int L)
{
  *P_new_ia=Allocate_mat2_f(L, 20);
  *P_opt_ia=Allocate_mat2_f(L, 20);
  *u1_ja=Allocate_mat2_f(L, 20);
  *u2_ja=Allocate_mat2_f(L, 20);
  *u3_ja=Allocate_mat2_f(L, 20);
  *psi1_ja=Allocate_mat2_f(L, 20);
  *psi2_ja=Allocate_mat2_f(L, 20);
  *Upsi1_ja=Allocate_mat2_f(L, 20);
  *Upsi2_ja=Allocate_mat2_f(L, 20);
  *U2psi1_ja=Allocate_mat2_f(L, 20);
  *u1U_ja=Allocate_mat2_f(L, 20);
}

void Empty_mean_field(float **P_new_ia, float **P_opt_ia,
		      float **u1_ja, float **u2_ja, float **u3_ja,
		      float **psi1_ja, float **psi2_ja,
		      float **Upsi1_ja, float **Upsi2_ja,
		      float **U2psi1_ja, float **u1U_ja, int L)
{
  Empty_matrix_f(P_new_ia, L);
  Empty_matrix_f(u1_ja, L);
  Empty_matrix_f(u2_ja, L);
  Empty_matrix_f(u3_ja, L);
  Empty_matrix_f(psi1_ja, L);
  Empty_matrix_f(psi2_ja, L);
  Empty_matrix_f(Upsi1_ja, L);
  Empty_matrix_f(Upsi2_ja, L);
  Empty_matrix_f(U2psi1_ja, L);
  Empty_matrix_f(u1U_ja, L);
}

float KL(float *P, float *Q, int n){ // Kullback-Leibler divergence
  int i; double KL=0;
  for(i=0; i<n; i++)if(P[i])KL+=P[i]*log(P[i]/Q[i]);
  return(KL);
}
