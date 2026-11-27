#include "thread_contfreq.h"
#include "energy_BKV.h"  // Contact energy parameters
#include "protein3.h"
#include "allocate.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// Parameters for calculating alpha_REM
#define AA 0.1
#define BB 4.0
#define Q0 0.1
#define L_MAX 10000
#define NCMAX 24    // Maximum number of contacts per residue

// Common parameters defined in .h files
//float **interactions;
struct protein prot[N_PROT_MAX], target;
float TEMP;
float sC1, sC0, SC, s0, Conf_entropy, K2Thr, TEMP2;
int LEN;

//float conf_entropy;
int INI_STAT=0;
float **C_nat=NULL;
double **Cont_freq=NULL, **C2=NULL, **C23=NULL, **C3=NULL;
double A2, A3;

void Contact_statistics(double ***Cont_freq, double ***C2, double ***C23,
			double ***C3, double *A2, double *A3,
			int len_nat, char *file_str);
float G_misf_2_TEMP(double E1, double E2, double A2);
float G_misf_3_TEMP(double E1, double E2, double A2,
		    double E23, double E3, double A3);
float Compute_Tfreeze_3(float K2, float K3, float SC);


float Compute_DG_REM3(short *aa_seq, int L, double *E_nat,
		      double *E1_misf, double *E2_misf,
		      double *E23_misf, double *E3_misf,
		      char *file_str)
{
  if(INI_STAT==0){
    Contact_statistics(&Cont_freq, &C2, &C23, &C3, &A2, &A3, L, file_str);
    SC=sC0+sC1*L; TEMP2=TEMP*TEMP; K2Thr=SC*TEMP2;
    Conf_entropy=s0*L; LEN=L;
    INI_STAT=1;
  }

  double Enat=0, E1=0, E2=0, E23=0, E3=0;
  float U, U2, U3, Gmisf, DeltaG;
  int i, j;

  for(i=0; i<L; i++){
    if(aa_seq[i]<0)continue;
    float *U_inter=Econt[aa_seq[i]];
    for(j=i+IJ_MIN; j<L; j++){
      if(aa_seq[j]<0)continue;
      U=U_inter[aa_seq[j]];
      E1+=Cont_freq[i][j]*U; U2=U*U;
      E2+=C2[i][j]*U2;
      E23+=C23[i][j]*U2; U3=U2*U;
      E3+=C3[i][j]*U3;
      if(C_nat[i][j])Enat+=C_nat[i][j]*U;
    }
  }
  (*E1_misf)=E1;
  (*E2_misf)=E2;
  (*E3_misf)=E3;
  (*E23_misf)=E23;
  (*E_nat)=Enat+Conf_entropy*TEMP;
  Gmisf=G_misf_3_TEMP(E1, E2, A2, E23, E3, A3);
  DeltaG=(Enat-Gmisf)/TEMP;
  return(DeltaG);
}

float Mutate_DG_REM3(short *aa_seq, int L, int res_mut, short aa_new,
		     double *E_nat, double *E1, double *E2,
		     double *E23, double *E3)
{
  if(aa_new<0 || aa_new>=20){
    printf("WARNING in Mutate_DG_REM3, aa_new= %d\n", aa_new);
    return(100);
  }
  if(aa_seq[res_mut]<0)return(0);
  float *inter_mut=Econt[aa_new], U_mut;
  float *inter_wt=Econt[aa_seq[res_mut]], U_wt;
  float *Cnat_i=C_nat[res_mut]; int j;
  double *C1_i=Cont_freq[res_mut], *C2_i=C2[res_mut];
  double *C23_i=C23[res_mut], *C3_i=C3[res_mut];
  for(j=0; j<L; j++){
    if(aa_seq[j]<0)continue;
    if(abs(j-res_mut)<IJ_MIN)continue;
    U_mut=inter_mut[aa_seq[j]];
    U_wt=inter_wt[aa_seq[j]];
    (*E1)+=C1_i[j]*(U_mut-U_wt);
    float U2_mut=U_mut*U_mut, U2_wt=U_wt*U_wt;
    (*E2)+=C2_i[j]*(U2_mut-U2_wt);
    (*E23)+=C23_i[j]*(U2_mut-U2_wt);
    float U3_mut=U2_mut*U_mut, U3_wt=U2_wt*U_wt;
    (*E3)+=C3_i[j]*(U3_mut-U3_wt);
    if(Cnat_i[j])(*E_nat)+=(U_mut-U_wt);
  }
  float Gmisf=G_misf_3_TEMP(*E1, *E2, A2, *E23, *E3, A3);
  float DeltaG=((*E_nat)-Gmisf)/TEMP;
  return(DeltaG);
}

float G_misf_3_TEMP(double E1, double E2, double A2,
		    double E23, double E3, double A3)
{
  float K2=A2*E1*E1+E2, S;
  float K3=A3*E1*E1*E1+E3+E1*E23;
  float Tfreeze=Compute_Tfreeze_3(K2, K3, SC);
  if(TEMP<Tfreeze){K2/=Tfreeze; K3/=(Tfreeze*Tfreeze); S=SC*Tfreeze;}
  else{K2/=TEMP; K3/=TEMP2; S=SC*TEMP;}
  float Gmisf=E1-K2-S+K3;
  return(Gmisf);
}

float Compute_DG_REM2(short *aa_seq, int L, double *E_nat,
		      double *E1_misf, double *E2_misf,
		      char *file_str)
{
  if(INI_STAT==0){
    Contact_statistics(&Cont_freq, &C2, &C23, &C3, &A2, &A3, L, file_str);
    SC=sC0+sC1*L; TEMP2=TEMP*TEMP; K2Thr=SC*TEMP2;
    Conf_entropy=s0*L; LEN=L;
    printf("Conf_entropy= %.2f\n", Conf_entropy);
    INI_STAT=1;
  }

  double Enat=0, E1=0, E2=0;
  float U, U2, Gmisf, DeltaG;
  int i, j;

  for(i=0; i<L; i++){
    if(aa_seq[i]<0)continue;
    float *U_inter=Econt[aa_seq[i]];
    for(j=i+IJ_MIN; j<L; j++){
      if(aa_seq[j]<0)continue;
      U=U_inter[aa_seq[j]];
      E1+=Cont_freq[i][j]*U; U2=U*U;
      E2+=C2[i][j]*U2;
      if(C_nat[i][j])Enat+=C_nat[i][j]*U;
    }
  }
  (*E1_misf)=E1;
  (*E2_misf)=E2;
  (*E_nat)=Enat+Conf_entropy*TEMP;
  Gmisf=G_misf_2_TEMP(E1, E2, A2);
  DeltaG=((*E_nat)-Gmisf)/TEMP;
  return(DeltaG);
}

float G_misf_2_TEMP(double E1, double E2, double A2){
  float K2=A2*E1*E1+E2, T;
  if(K2>K2Thr){T=sqrt(K2/SC);} // Freezing
  else{T=TEMP;}
  float Gmisf=E1-K2/T-SC*T;
  return(Gmisf);
}

float Mutate_DG_REM2(short *aa_seq, int L, int res_mut, short aa_new,
		     double *E_nat, double *E1, double *E2)
{
  if(aa_new<0 || aa_new>=20){
    printf("WARNING in Mutate_DG_REM2, aa_new= %d\n", aa_new);
    return(100);
  }
  int aa_old=aa_seq[res_mut];
  if(aa_old<0 || aa_old>=20){
    printf("WARNING in Mutate_DG_REM2, aa_wt= %d at res %d\n",
	   aa_old, res_mut); return(0);
  }
  float *inter_mut=Econt[aa_new], U_mut;
  float *inter_wt=Econt[aa_old], U_wt;
  float *Cnat_i=C_nat[res_mut];
  double *C1_i=Cont_freq[res_mut], *C2_i=C2[res_mut];
  for(int j=0; j<L; j++){
    if(abs(j-res_mut)<IJ_MIN)continue;
    if(aa_seq[j]<0)continue;
    U_mut=inter_mut[aa_seq[j]];
    U_wt=inter_wt[aa_seq[j]];
    (*E1)+=C1_i[j]*(U_mut-U_wt);
    (*E2)+=C2_i[j]*(U_mut*U_mut-U_wt*U_wt);
    if(Cnat_i[j])(*E_nat)+=(U_mut-U_wt);
  }
  float Gmisf=G_misf_2_TEMP(*E1, *E2, A2);
  float DeltaG=((*E_nat)-Gmisf)/TEMP;
  return(DeltaG);
}


float Print_DG_REM2(double E_nat, double E1, double E2, char *name)
{
  FILE *file_out=fopen(name, "w");
  printf("Writing %s\n", name);
  float T, T_INI=0.1, T_END=5, T_STEP=0.1;
  float T0, DG0=100000;
  float K2=A2*E1*E1+E2, K2S;
  float Tfreeze=sqrt(K2/SC);
  float E_NAT=E_nat-Conf_entropy*TEMP;
  fprintf(file_out,
	  "# E_nat/L= %.3f K2/L= %.4f SC/L= %.3f Tfreeze= %.3f\n",
	  E_NAT/LEN, K2/LEN, SC/LEN, Tfreeze);
  fprintf(file_out, "# Temp. DeltaG/L (REM2)\n");
  for(T=T_INI; T<T_END; T+=T_STEP){
    if(T<Tfreeze){K2S=K2/Tfreeze+SC*Tfreeze;}  // Freezing
    else{K2S=K2/T+SC*T;}
    float DeltaG=(E_NAT-E1+K2S+Conf_entropy*T);
    fprintf(file_out, "%.1f %.5f\n", T, DeltaG/LEN);
    if(fabs(DeltaG)<fabs(DG0)){
      if((DeltaG<0)||(DG0>0)){DG0=DeltaG; T0=T;}
    }
  }
  fclose(file_out);
  return(T0);
}

float Print_DG_REM3(double E_nat, double E1, double E2,
		     double E23, double E3, char *name){
  FILE *file_out=fopen(name, "w");
  printf("Writing %s\n", name);
  float T, T_INI=0.1, T_END=5, T_STEP=0.1;
  float K2=A2*E1*E1+E2, K3S;
  float K3=A3*E1*E1*E1+E3+E1*E23;
  float Tfreeze=Compute_Tfreeze_3(K2, K3, SC);
  float Tfreeze2=Tfreeze*Tfreeze;
  float T0, DG0=100000;
  float E_NAT=E_nat-Conf_entropy*TEMP;
  fprintf(file_out,
	  "# Enat/L= %.3f K2/L= %.4f SC/L= %.3f K3/L= %.3f Tfreeze= %.3f\n",
	  E_NAT/LEN, K2/LEN, SC/LEN, K3/LEN, Tfreeze);
  fprintf(file_out, "# Temp. DeltaG/L (REM3)\n");
  for(T=T_INI; T<T_END; T+=T_STEP){
    if(T<Tfreeze){K3S=K2/Tfreeze+SC*Tfreeze-K3/Tfreeze2;}  // Freezing     
    else{K3S=K2/T+SC*T-K3/(T*T);}
    float DeltaG=(E_NAT-E1+K3S+Conf_entropy*T);
    fprintf(file_out, "%.1f %.5f\n", T, DeltaG/LEN);
    if(fabs(DeltaG)<fabs(DG0)){
      if((DeltaG<0)||(DG0>0)){DG0=DeltaG; T0=T;}
    }
  }
  fclose(file_out);
  return(T0);
}

float Compute_Tfreeze_3(float K2, float K3, float SC){
  if(K3>0)return(0);
  float T=sqrt(K2/SC), T_new, f, f1; int it;
  for(it=0; it<20; it++){
    f=SC*T*T*T-K2*T+2*K3;
    if(fabs(f)<0.0001)break;
    f1=3*SC*T*T-K2;
    T_new=T-f/f1;
  }
  return(T); 
}


int Get_target(char *file_str, char *name_tar, int *len_tar)
{
  FILE *file_in=fopen(file_str, "r");
  char string[200], name[20], dumm[4];
  int N_prot=0, length, n_cont, i_res, j;
  int i, res1, res2, nc[L_MAX], found=0;
  struct contact *cont_list; short **contact;

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
      target.length=length;
      cont_list=malloc((n_cont+1)*sizeof(struct contact));
      contact=malloc((length)*sizeof(short *));
      target.contact=contact;
      target.cont_list=cont_list;
      for(i=0; i<length; i++){
	contact[i]=malloc(NCMAX*sizeof(short)); nc[i]=0;
      }

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
      found=1;
      break;
    }else{
      // Other structure; skipped
      for(i=0; i<n_cont; i++)fgets(string, sizeof(string), file_in);
    }
  }
  fclose(file_in);

  if(found==0){
    printf("ERROR in Get_target, pdb %s not found in file %s, %d proteins\n",
	   name_tar, file_str, N_prot); exit(8);
  }
  Fill_C_nat(length, contact, 1);
  return(1);
}


void Contact_statistics(double ***Cont_freq, double ***C2, double ***C23,
			double ***C3, double *A2, double *A3,
			int len_nat, char *file_str)
{
  int N_PROT=Read_structures(file_str, prot);
  short j_prot=0, first, i, j;
  short res1, *res2, s1, s2;
  long n_str=0;

  // Allocate
  struct contact *Cont_list=malloc((40*len_nat)*sizeof(struct contact));
  double **Cont_Nc=Allocate_mat2_d(len_nat, len_nat, "Cont_Nc");
  (*Cont_freq)=Allocate_mat2_d(len_nat, len_nat, "Cont_freq");
  (*C2)=Allocate_mat2_d(len_nat, len_nat, "C2");
  (*C3)=Allocate_mat2_d(len_nat, len_nat, "C3");
  (*C23)=Allocate_mat2_d(len_nat, len_nat, "C23");
  double NC1=0, NC2=0, NC3=0, Ncc;

  /* Generate alternative structures by threading */
  for(j_prot=0; j_prot< N_PROT; j_prot++){

    short len2=prot[j_prot].length;
    if(len2 < len_nat)continue;
    short **contact=prot[j_prot].contact;

    for(first=0; first<=len2-len_nat; first++){
      
      /* New structure */
      short last=first+len_nat;
      short Nc=0; struct contact *cont=Cont_list;
      for(res1=first; res1<last; res1++){
	res2=contact[res1]; s1=res1-first;
	while((*res2>=0)&&(*res2 < last)){
	  if(*res2 >= res1+IJ_MIN){
	    cont->res1=s1;
	    cont->res2=*res2-first;
	    Nc++; cont++;
	  }
	  res2++;
	}
      }
      // Statistics
      n_str++;
      NC1+=Nc; Ncc=Nc*Nc; NC2+=Ncc; Ncc*=Nc; NC3+=Ncc;
      cont=Cont_list;
      for(i=0; i<Nc; i++){
	(*Cont_freq)[cont->res1][cont->res2]++;
	Cont_Nc[cont->res1][cont->res2]+=Nc;
	cont++;
      }
    }
  }
  free(Cont_list);

  /* End of statistics, compute averages */
  double eta2=0, eta3=0, NC_C2=0, c, c2, c3;
  for(i=0; i<len_nat; i++){
    for(j=i+IJ_MIN; j<len_nat; j++){
      (*Cont_freq)[i][j]/=n_str;
      (*Cont_freq)[j][i] = (*Cont_freq)[i][j];
      Cont_Nc[i][j]/=n_str;
      c=(*Cont_freq)[i][j]; c2=c*c;
      eta2 += c2; eta3+=c2*c;
      NC_C2 += Cont_Nc[i][j]*c;
    }
  }
  NC1/=n_str; NC2/=n_str; NC3/=n_str; 
  Empty_matrix_d(Cont_Nc, len_nat);
  
  double NC1_2=NC1*NC1;
  double K1= NC1-eta2;
  double K2= NC2-NC1_2;
  double K3= NC3-3*NC2*NC1+2*NC1_2*NC1;
  double B2= (K2-K1)/(1.-eta2/NC1_2);
  double B33= K3-3*K2+6*NC_C2-6*eta2*(NC1+1)+2*NC1+4*eta3;
  double B32= K2-2*NC_C2+2*NC1*eta2-NC1+3*eta2-2*eta3;
  
  double C2_1=0.5;
  double C2_2=0.5*(1+(K2-K1)/NC1_2);
  double C23_1=3*B32/(NC1-eta2/NC1);
  double C23_2=3*B33/NC1_2;
  double C3_2=3*(1+B32/NC1_2);
  double C3_3=2*(1+B33/(NC1_2*NC1));
  
  for(i=0; i<len_nat; i++){
    for(j=0; j<len_nat; j++){
      c=(*Cont_freq)[i][j]; c2=c*c; c3=c2*c;
      (*C2)[i][j]=C2_1*c-C2_2*c2;
      (*C23)[i][j]=C23_1*c-C23_2*c2;
      (*C3)[i][j]=c-C3_2*c2+C3_3*c3;
    }
  }
  *A2=B2/(2*NC1_2);
  *A3=B33/(6*(NC1_2*NC1-3*eta2*NC1));
  
  char nameout[400];
  sprintf(nameout, "Contact_statistics_%d.dat", len_nat);
  printf("Writing %s\n", nameout);
  FILE *file_out=fopen(nameout, "w");
  fprintf(file_out, "Contact statistics for protein of length %d\n",len_nat);
  fprintf(file_out, "<NC>= %.2f\n", NC1);
  fprintf(file_out, "<NC^2>-<NC>^2= %.2f  /<NC>= %.3f\n", K2, K2/NC1);
  fprintf(file_out, "<(NC-<NC>)^3>= %.2f  /<NC>= %.3f\n", K3, K3/NC1);
  fprintf(file_out, "sum_ij <C_ij>^2= %.2f\n", eta2); 
  fprintf(file_out, "sum_ij <C_ij>^3= %.2f\n", eta3);
  fprintf(file_out, "B2= %.3f B32= %.3f B33= %.3f\n", B2, B32, B33);
  fprintf(file_out,
	  "C2_2= %.3f C23_1= %.3f C23_2= %.3f C3_2= %.3f C3_3= %.3f\n",
	  C2_2, C23_1, C23_2, C3_2, C3_3);
  fprintf(file_out, "A2= %.6f A3= %.6f\n", *A2, *A3);
  fclose(file_out);
  Empty_matrix_d(Cont_Nc, len_nat);

  return;
}



Read_structures(char *file_prot, struct protein *prot)
{
  FILE *file_in=fopen(file_prot, "r");
  char string[200], name[20], dumm[4];
  int N_prot=0, length_old=10000, length, n_cont, i_res;
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
