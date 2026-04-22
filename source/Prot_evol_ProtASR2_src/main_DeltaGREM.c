/* Program DeltaGREM
   Computes contact free energy (native-unfolding-misfolding)
   of a set of sequences with respect to the best of a list of
   target structures. Sequences can also be given as mutations.
*/

float LMIN=0.5; // Minimum length for selecting alignments
int PRINT_MUT=1; // Print details of mutations?
char FILE_STR[200];

#define FILE_CODE_DEF "gen_code_ATGC.in"

#include "REM.h"
#include "coord.h"
#include "mut_del.h"
#include "alignments.h"
#include "gen_code.h"
#include "allocate.h"
#include "protein3.h"
#include "read_pdb.h"
#include "random3.h"           /* Generating random numbers */
#include "mutation.h"
#include "codes.h"
#include "input.h"
#include "get_para_DeltaGREM.h"
#include "Sec_str_all.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

// External variables
int DTAIL_MAX=20;
float cont_thr=4.5;
int LEN; // For REM calculations
char AA_code[21];
float hydro[21];
int Verbose=0;
//int L_diso=0;

int ini_print;

/* 
   Give as input a PDB file and either an alignment file or a file with
   mutations. Output is the DeltaG of wild type and mutant sequences.
*/

#define N_CHAR 300          // Max. length of file names

struct result{
  float E_nat_1;      // Mean of Gnat over native models
  float E_nat_2;      // Mean square of Gnat over native models
  float E_nat_min;    // Minimum value of Gnat over native models
  int pdbbest;      // Best model where minimum is attained
  float seqid_opt;  // Largest seqid over possible models
  float seqid_mod;  // Seqid for best energy model
  float hydro;      // Average hydrophobicity
  int n, L;         // Number of models and length of optimal model
  int nmut;         // Number of mutations
  int N_disulf;     // Number of native disulfide bonds
  int Len_disulf;   // Sum of length of native disulfide bonds
  float Log_len_disulf; // Sum of log length of native disulfide bonds
  int N_helix;      // Helical interactions in native structure
  int N_strand;     // Strand interactions in native structure

  //float DeltaG_1;   // Mean of DeltaG over native models
  //float DeltaG_2;   // Mean square of DeltaG over native models
  //float DeltaG_min; // Minimum of DeltaG over native models

};

int N_disulf_seq=0;
int Len_disulf_seq=0;
float Log_len_disulf_seq=0;
int N_helix_seq=0;
int N_strand_seq=0;

// Input codes
int Remove_gaps(char **MSA, int n_seq, int L_ali, char **name,
		char *node, int c, char gap);

// Output
struct result *Allocate_res(int n);
void Initialize_res(struct result *);
void Record(struct result *r, struct REM E, int i_pdb, float seq_id, int L);
char *Get_name(char *file_name, char *chain);
void Print_result(struct result *r, double *DDG_1, double *DDG_2, int *n,
		  char *name, char **name_pdb, float DG_wt, float G_nat_wt,
		  float hydro_wt, float fit_wt, struct REM E, FILE *file_out);
void Print_mut(short *aa_seq, short *aa_seq0, float **C_nat, int len_amm,
	       FILE *file_mut_out);
void Print_DG(struct REM E, struct result *r, short *aa_seq, int *i_sec,
	      char *c_sec, char *name_pdb, int ini, char *name);
void Empty_prot(struct protein *prot);
static void Modify_name_del_interface(struct deletion *del, struct protein prot,
				      struct residue *res, int L);
static void Modify_name_mut_interface(struct mutation *mut, struct protein prot,
				      struct residue *res, int L);

static void Print_interface(struct protein prot, short *aa_seq,
			    char *name_pdb, char *chain);
extern int 
Find_site(char amm, int rmut, char chain, struct residue *res, int N_res);

// Input parameters
// A: Input files defining the protein
static char dir_out[100];
// B: Thermodynamic parameters
float sC1=0.065, sC0=0, sU1=0.140; // Configuration entropy
int REM=2;   // Use 1st (1), 2nd (2) and 3rd (3) moment of misfolded energy 
float S_C=0, S_U=0;
float TEMP=0.5; // default
// C2: Mean-field model
//int MEANFIELD=1;  // Compute site-specific mean-field distributions
//int OPT_LAMBDA=1; // 1=Optimize Lambda by maximum likelihood
//float LAMBDA=0.9; // Lambda parameter for meanfield if OPT_LAMBDA=0
//char MODEL[N_CHAR]="ALL"; // Type of mean-field model (ALL, NAT, DG)
//float DG_OPT=-1;  // DG target of the optimization if MODEL="DG"
// D1: Mutation model, P_mut
//#define MUTPAR 10
//float mut_par[MUTPAR], tt_ratio=1, kCpG=1;
// D1: Mutation model, exchangeability
//char MATRIX[40]="WAG"; // Empirical exchangeability matrix
//char EXCHANGE='F'; // exchangeability model. M=MUT F=FLUX Q=RATE E=EXCH
//float TWONUC=0; // Rate of two to one nucleotide substitutions  if MUT
// E: Output
//int FORMAT=1;   // PAML format

// Derived data
unsigned long iran;
static int len_amm;
float DG_thr;

float Econt[21][21];
float **Econt_T=NULL;
float SEC_STR=1; // Coefficient for local interactions

int main(int argc, char **argv){


  /***********************
          INPUT
  ************************/
  // Input files
  //char Input_dir[N_CHAR];
  char **file_pdb, **chain_pdb;
  char *file_ali=NULL, *file_mut=NULL;

  /***********************
         SEQUENCES
  ************************/
  //int nuc_mut, res_mut=0, aa_new;

  /***********************
         Wild Type
  ************************/
  //float fitness_wt;

  /***********************
         MUTANT
  ************************/
  //float fitness_mut;

  /***********************
          DUMMIES
  ************************/
  int i, j;

  /***********************
          OUTPUT
  ************************/

  /******************** Input operations   ************************/
  int N_pdb=Get_para(argc, argv, &file_pdb, &chain_pdb,
		     FILE_STR, &TEMP, &sU1, &sC0, &sC1, &REM, &SEC_STR,
		     &file_ali, &file_mut, dir_out);
  printf("%d PDB files to read\n", N_pdb);
  printf("Temperature= %.3f\n", TEMP);

  // alignment
  int n_seq=0, L_ali=0, *selected=NULL;
  char **name_seq;
  char **MSA=Read_MSA(&n_seq, &L_ali, &name_seq, &selected, file_ali, LMIN);
  if(MSA==NULL){
    printf("No multiple alignment given, checking mutations\n");
    file_ali=NULL;
  }else{
    Remove_gaps(MSA, n_seq, L_ali, name_seq, "node", 4, 'L');
  }
 
  // Mutations-deletions
  int Nmut=0; struct mutation *mutations=NULL;
  int Ndel=0; struct deletion *deletions=NULL;
  FILE *file_mut_out=NULL;
  if(file_mut){
    Nmut=Read_mut(&mutations, file_mut, &deletions, &Ndel);
    printf("%d mutations, %d deletions\n", Nmut, Ndel);
    if(PRINT_MUT && Nmut){
      char name_mut[900];
      sprintf(name_mut, "%s_mutations.dat", file_mut);
      file_mut_out=fopen(name_mut, "w");
      printf("Writing %s\n", name_mut);
    }
  }
  if( n_seq &&((Nmut)||(Ndel))){
      printf("ERROR, it is not allowed to input MSA (%d seq.) and ", n_seq);
      printf("mutations (%d a.a. changes and %d indels)\n", Nmut, Ndel);
      exit(8);
  }

  /**************************  Allocate  ****************************/
  char **name_pdb=malloc(N_pdb*sizeof(char *));
  float *seq_id=NULL; int *seq_L=NULL, i_pdb;
  if(n_seq){
    seq_id=malloc(n_seq*sizeof(float));
    seq_L=malloc(n_seq*sizeof(int));
  }
  struct result *result_seq=Allocate_res(n_seq);
  struct result *result_mut=Allocate_res(Nmut);
  struct result *result_del=Allocate_res(Ndel);
  struct result *result_wt=Allocate_res(1);
  struct REM E_ref;
  float DG_wt, E_nat_wt, hydro_wt, fit_wt;
  
  int n_pdb=0;
  for(i_pdb=0; i_pdb<N_pdb; i_pdb++){
    int N_mod=Count_models_PDB(file_pdb[i_pdb]), imod;
    if(N_mod<=0){
      printf("WARNING, no PDB file found: %s\n", file_pdb[i_pdb]);
      continue;
    }
    printf("PDB file %s, %d models\n", file_pdb[i_pdb], N_mod);
    struct protein *models=malloc(N_mod*sizeof(struct protein));
    struct residue *res=NULL; struct res_short *res_short; int *res_index;
    len_amm=Get_pdb(models, &res, &res_short, file_pdb[i_pdb],
		    chain_pdb[i_pdb], &res_index, N_mod);
    name_pdb[i_pdb]=Get_name(file_pdb[i_pdb], chain_pdb[i_pdb]);

    if(len_amm<=0){
      printf("WARNING, no chain found: %s chain %s\n",
	     file_pdb[i_pdb], chain_pdb[i_pdb]); continue;
    }

    target=models[0];
    float len_str=target.L_PDB/target.num_chain;
    //L_diso=len_amm-len_str;
    // len_amm=target.length are all residues, including disordered ones
    // target.L_PDB is the number of structured residues, one per each chain
    //S_C=sC0+len_amm*sC1; S_U=len_amm*sU1;
    S_C=sC0+len_str*sC1; S_U=len_str*sU1;

    short *aa_seq0=target.aa_seq;
    short aa_seq[len_amm];
    int ali_seq[len_amm];
    int k_pdb; // MSA sequence that corresponds to wild-type

    if(target.chain_all > 1){
      for(i=0; i<Nmut; i++){
	Modify_name_mut_interface(mutations+i, target, res, target.L_PDB);
      }
      for(i=0; i<Ndel; i++){
	Modify_name_del_interface(deletions+i, target, res, target.L_PDB);
      }
    }

    if(MSA==NULL){
      // If not given, construct MSA with wt and mutants
      k_pdb=0;
      n_seq=1+Nmut+Ndel;
      L_ali=len_amm;
      MSA=malloc(n_seq*sizeof(char *));
      name_seq=malloc(n_seq*sizeof(char *));
      seq_L=malloc(n_seq*sizeof(int));
      seq_id=malloc(n_seq*sizeof(float));
      file_ali=malloc(40*sizeof(char));
      strcpy(file_ali, name_pdb[i_pdb]);
      int n, i, j;
      for(n=0; n<n_seq; n++){
	MSA[n]=malloc(L_ali*sizeof(char));
	name_seq[n]=malloc(100*sizeof(char));
	seq_L[n]=len_amm;
      }
      n=0; seq_id[0]=1.0000; strcpy(name_seq[n], name_pdb[i_pdb]);
      for(j=0; j<L_ali; j++)MSA[0][j]=Amin_code(aa_seq0[j]); // wt
      for(j=0; j<L_ali; j++)ali_seq[j]=j; // Align wild type with MSA

      int discard=0;
      for(i=0; i<Nmut; i++){
	printf("Mutation %s\n", mutations[i].name);
	n=i+1;
	strcpy(name_seq[n], mutations[i].name);
	for(j=0; j<len_amm; j++)aa_seq[j]=aa_seq0[j];
	if(Construct_mut(aa_seq, mutations+i, res, len_amm)<0){
	  discard++; continue;
	}
	for(j=0; j<L_ali; j++)MSA[n][j]=Amin_code(aa_seq[j]); // mut
	seq_id[n]=1-mutations[i].nmut/len_amm;
      }
      int discard2=0; int start=1+Nmut-discard;
      for(i=0; i<Ndel; i++){
	printf("Deletion %s %d\n", deletions[i].name, deletions[i].Ldel);
	n=i+start;
	strcpy(name_seq[n], deletions[i].name);
	seq_L[n]-=deletions[i].Ldel;
	seq_id[n]-=deletions[i].Ldel;
	seq_id[n]/=len_amm;
	for(j=0; j<len_amm; j++)aa_seq[j]=aa_seq0[j];
	if(Construct_del(aa_seq, deletions+i, res, len_amm)<0){
	  discard2++; continue;
	}
	for(j=0; j<L_ali; j++)MSA[n][j]=Amin_code(aa_seq[j]); // del
      }
      if(discard || discard2){
	printf("WARNING, %d mutations over %d and %d deletions over %d "
	       "have been discarded because of sequence discrepancies\n",
	       discard, Nmut, discard2, Ndel);
	n_seq-=(discard+discard2);
      }
      result_seq=Allocate_res(n_seq);
      Nmut=0; Ndel=0;
      // end build MSA
    }else{
      // Align PDB with MSA
      char seq_PDB[len_amm];
      for(i=0; i<len_amm; i++)seq_PDB[i]=Amin_code(aa_seq0[i]);
      k_pdb=Find_seq(ali_seq,seq_L,seq_id, seq_PDB,len_amm,
		     MSA,n_seq,L_ali);
      if(k_pdb<0){
	printf("WARNING, PDB sequence %s chain %s not found in MSA %s\n",
	       file_pdb[i_pdb], chain_pdb[i_pdb], file_ali);
      }
    }


    // Contact energy computations
    struct REM E_wt; E_wt.L_str=len_str;
    Initialize_E_REM(&E_wt, len_amm, REM, TEMP, S_C, S_U, FILE_STR, 0);
    struct REM E_mut; E_mut.L_str=len_str;
    if(Nmut || Ndel || n_seq){
      Initialize_E_REM(&E_mut, len_amm, REM, TEMP, S_C, S_U, FILE_STR, 0);
    }
    if(target.chain_all > 1){
      Print_interface(target,target.aa_seq,name_pdb[i_pdb],chain_pdb[i_pdb]);
    }

    Initialize_res(result_wt);

    for(imod=0; imod<N_mod; imod++){
 
      /******************** Folding stability **********************/
      int *i_sec=models[imod].i_sec;
      char *c_sec=models[imod].sec_str;
      float **C_nat=Fill_C_nat(len_amm, models[imod].contact,
			       models[0].num_chain);
      //float **Dist_nat=Dist_nat(len_amm, models[imod].dist);

      if(imod==0){
	n_pdb++;
	Test_contfreq(&E_wt, aa_seq0, C_nat, i_sec, c_sec, name_pdb[i_pdb]);
	E_wt.DeltaG=
	  //Compute_DG_overT_contfreq(&E_wt, aa_seq0, C_nat, i_sec, c_sec);
	  Compute_DG_overT_threading(&E_wt, aa_seq0, C_nat, i_sec, c_sec);
	result_wt[0].hydro=Compute_hydro(aa_seq0, len_amm);
      }else{
	Compute_E_nat(&E_wt, aa_seq0, C_nat, i_sec, c_sec);
      }
      // Wild type
      Record(result_wt, E_wt, i_pdb, 1.000, len_amm);
 
      if(MSA && k_pdb>=0){
	for(i=0; i<n_seq; i++){
	  for(j=0; j<len_amm; j++){
	    if(ali_seq[j]>=0){aa_seq[j]=Code_AA(MSA[i][ali_seq[j]]);}
	    else{aa_seq[j]=-1;} // gap, previously it was 20
	  }
	  if(imod==0){
	    E_mut.DeltaG=
	      Compute_DG_overT_threading(&E_mut,aa_seq,C_nat,i_sec,c_sec);
	    E_mut.DeltaG+=(sU1)*(seq_L[i]-len_str); //+sC1
	    result_seq[i].hydro=Compute_hydro(aa_seq, len_amm);
	  }else{
	    Compute_E_nat(&E_mut, aa_seq, C_nat, i_sec, c_sec);
	  }
	  Record(result_seq+i, E_mut, i_pdb, seq_id[i], seq_L[i]);
	}
      }

      // Mutation
      for(i=0; i<Nmut; i++){
	printf("Mutation %s\n", mutations[i].name);
	for(j=0; j<len_amm; j++)aa_seq[j]=aa_seq0[j];
	if(Construct_mut(aa_seq, mutations+i, res, len_amm)<0)continue;
	if(PRINT_MUT)
	  Print_mut(aa_seq, aa_seq0, C_nat, len_amm, file_mut_out);
	if(imod==0){
	  E_mut.DeltaG=
	    //Compute_DG_overT_contfreq(&E_mut, aa_seq, C_nat, i_sec, c_sec);
	    Compute_DG_overT_threading(&E_mut, aa_seq, C_nat, i_sec, c_sec);
	  result_mut[i].hydro=Compute_hydro(aa_seq, len_amm);
	}else{
	  Compute_E_nat(&E_mut, aa_seq, C_nat, i_sec, c_sec);
	}
	// Statistics
	float id=1.-mutations[i].nmut/(float)len_amm;
	Record(result_mut+i, E_mut, i_pdb, id, len_amm);
	result_mut[i].nmut=mutations[i].nmut;
      }
      
      // Deletion
      for(i=0; i<Ndel; i++){
	int del=deletions[i].Ldel;
	printf("Deletion %s %d\n", deletions[i].name, del);
	for(j=0; j<len_amm; j++)aa_seq[j]=aa_seq0[j];
	if(Construct_del(aa_seq, deletions+i, res, len_amm)<0)continue;
	E_mut.DeltaG=
	  //Compute_DG_overT_contfreq(&E_mut, aa_seq, C_nat, i_sec, c_sec);
	  Compute_DG_overT_threading(&E_mut, aa_seq, C_nat, i_sec, c_sec);
	E_mut.DeltaG-=(sU1)*del;
	if(imod==0)result_del[i].hydro=Compute_hydro(aa_seq, len_amm);
	// Statistics
	float id=1.-del/(float)len_amm;
	Record(result_del+i, E_mut, i_pdb, id, len_amm-del);
      }
      for(i=0; i<len_amm; i++){free(C_nat[i]);} free(C_nat);
    } // end models
    Print_DG(E_wt, result_wt, aa_seq0, models[0].i_sec, models[0].sec_str,
	     name_pdb[i_pdb], 1, "WT");

    if(n_pdb==1){
      E_ref=E_wt;
      E_nat_wt=result_wt[0].E_nat_min;
      DG_wt=E_wt.DeltaG+E_nat_wt-E_wt.E_nat;;
      hydro_wt=result_wt[0].hydro;
      fit_wt=-log(1+exp(DG_wt));
    }

    for(imod=0; imod<N_mod; imod++)Empty_prot(models+imod);
    free(models);
    free(res);
  }
  if(n_pdb==0){
    printf("ERROR, no pdb file found\n"); exit(8);
  }
  
  /************************* Output files ***************************/
  char nameout[300], para[80]="";
  sprintf(para, "T%.2f_SU%.3f_SC%.3f", TEMP, sU1, sC1);
  if(file_ali){ // MSA
    sprintf(nameout, "%s%s", dir_out, file_ali);
  }else if((Nmut) || (Ndel)){
    sprintf(nameout, "%s%s", dir_out, file_mut);
  }else{
    sprintf(nameout, "%sProts", dir_out);
  }
  strcat(nameout, "_DeltaG.dat");
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);

  char header[2000];
  sprintf(header,
	  "# T= %.2f sU1= %.3f sC1= %.3f sC0= %.3f REM= %d\n"
	  "# %d PDB files used\n", TEMP, sU1, sC1, sC0, REM, n_pdb);
  if((Nmut==0)&&(Ndel==0)){
    strcat(header, "#L\tDelta_G/L(best)\tsigma(Delta_G)/L\tG_nat/N");
    if(n_pdb>1)strcat(header, "\tseq.id(best)");
    strcat(header, "\tseq.id(model)\tN_disulf\tLen_disulf\tLog_len_disulf"
	   "\tN_helix\tN_strand\thydro\tlog(f)\ttemplate\tsequence");
  }else{
    strcat(header, "#L\tDDG\tsigma(DDG)\tDG_nat");
    if(n_pdb>1)strcat(header,"\tseq.id(best)");
    strcat(header,"\tseq.id(model)\tN_disulf\tLen_disulf\tLog_len_disulf"
	   "\tN_helix\tN_strand\thydro\tlog(f)\tDeltaG/nmut\ttemplate\tmut");
  }

  fprintf(file_out, "%s\n", header);

  double Dh_1=0, Dh_2=0, h, DDG_1=0, DDG_2=0; int num=0;
  Print_result(result_wt, &DDG_1, &DDG_2, &num, "WT",
	       name_pdb, 0, 0, 0, 0, E_ref, file_out);
  for(i=0; i<n_seq; i++){
    Print_result(result_seq+i, &DDG_1, &DDG_2, &num, name_seq[i],
		 name_pdb, DG_wt, E_nat_wt, hydro_wt, fit_wt, E_ref,
		 file_out);
    h=result_seq[i].hydro; Dh_1+=h; Dh_2+=h*h;
  }
  for(i=0; i<Nmut; i++){
    Print_result(result_mut+i, &DDG_1, &DDG_2, &num, mutations[i].name,
		 name_pdb, DG_wt, E_nat_wt, hydro_wt, fit_wt, E_ref,
		 file_out);
    h=result_mut[i].hydro; Dh_1+=h; Dh_2+=h*h;
  }
  for(i=0; i<Ndel; i++){
    result_del[i].nmut=1;
    Print_result(result_del+i, &DDG_1, &DDG_2, &num, deletions[i].name,
		 name_pdb, DG_wt, E_nat_wt, hydro_wt, fit_wt, E_ref, file_out);
    h=result_del[i].hydro; Dh_1+=h; Dh_2+=h*h;
  }

  if(num>2){
    DDG_1/=num;
    DDG_2=sqrt((DDG_2/num-DDG_1*DDG_1)/(num-1));
    DDG_1-=DG_wt;
    fprintf(file_out, "# Mean (D)DG= %.3f S.E.M.= %.3f\n", DDG_1, DDG_2);
    printf("%s Mean (D)DG= %.3f Error= %.3f\n", file_mut, DDG_1, DDG_2);
    Dh_1/=num;
    Dh_2=sqrt((Dh_2/num-Dh_1*Dh_1)/(num-1));
    Dh_1-=hydro_wt;
    fprintf(file_out, "# Mean (D)hydro= %.3g S.E.M.= %.2g\n", Dh_1, Dh_2);
  }
  fclose(file_out);

  return(0);
}
 
void Empty_prot(struct protein *prot){
  if(prot->length==0)return;
  if(prot->contact){
    for(int i=0; i<prot->length; i++){
      free(prot->contact[i]);
    }
    free(prot->contact);
    prot->contact=NULL;
  }
  if(prot->cont_list){
    free(prot->cont_list); prot->cont_list=NULL;
  }
  if(prot->aa_seq){
    free(prot->aa_seq); prot->aa_seq=NULL;
  }
  if(prot->sec_str){
    free(prot->sec_str); prot->sec_str=NULL;
  }
  if(prot->i_sec){
    free(prot->i_sec); prot->i_sec=NULL;
  }
}


char *Get_name(char *file_name, char *chain){
  char *name=malloc(100*sizeof(char)), *n=name;
  char *s=file_name, *s1=s;
  while(*s!='\0'){if(*s=='/')s1=s; s++;} s1++;
  while(*s1!='\0'){if(*s1=='.')break; *n=*s1; n++; s1++;}
  *n='\0';
  strcat(name, chain);
  return(name);
}

int Remove_gaps(char **MSA, int n_seq, int L_ali, char **name,
		char *node, int c, char cgap){
  int iseq, i, Lgap=0, inode=0, ini=1, undecided=0, skip=0, mseq=0;
  for(i=0; i<L_ali; i++){
    int n_gap=0, n_c=0;
    for(iseq=0; iseq<n_seq; iseq++){
      if(strncmp(name[iseq], node, c)==0)continue;
      if(MSA[iseq][i]=='-'){n_gap++;}
      else if(MSA[iseq][i]==cgap){n_c++;}
      if(i==0)mseq++;
    }
    if(n_gap==0)continue;    // No gap present in real sequence
    if(n_c)undecided++;      // Not clear decision
    if(n_c > n_gap){skip++; continue;} // More characters than gaps in column
    for(iseq=0; iseq<n_seq; iseq++){
      if(strncmp(name[iseq], node, c)!=0)continue;
      if(MSA[iseq][i]==cgap){MSA[iseq][i]='-'; Lgap++;}
      if(ini)inode++;
    }
    if(ini)ini=0;
  }
  printf("%d undecided columns over %d, %d of them treated as gap\n",
	 undecided, L_ali, undecided-skip);
  printf("%d gap characters introduced in %d internal nodes\n",
	 Lgap, inode);
  //exit(8);

  char nameout[50]; strcpy(nameout, "newalignment.fasta");
  FILE *file_out=fopen(nameout, "w");
  for(iseq=0; iseq<n_seq; iseq++){
    fprintf(file_out, ">%s\n", name[iseq]);
    for(i=0; i<L_ali; i++)fprintf(file_out, "%c", MSA[iseq][i]);
    fprintf(file_out, "\n");
  }
  fclose(file_out);
  printf("Writing %s\n", nameout);
  return(0);
}

void Initialize_res(struct result *r){
  //r->DeltaG_1=0;
  //r->DeltaG_2=0;
  //r->DeltaG_min=0;
  r->E_nat_1=0;
  r->E_nat_2=0;
  r->E_nat_min=0;
  r->pdbbest=-1;
  r->seqid_opt=0;
  r->seqid_mod=0;
  r->n=0;
  r->nmut=0;
  r->N_disulf=0;
  r->Len_disulf=0;
  r->Log_len_disulf=0;
  r->N_helix=0;
  r->N_strand=0;
}

struct result *Allocate_res(int n)
{
  if(n==0)return(NULL);
  struct result *res=malloc(n*sizeof(struct result)), *r=res;
  int i;
  for(i=0; i<n; i++){
    Initialize_res(r);
    r++;
  }
  return(res);
}

void Record(struct result *r, struct REM E, int i_pdb, float seq_id, int L)
{
  r->E_nat_1+=E.E_nat;
  r->E_nat_2+=E.E_nat*E.E_nat;
  if((r->n==0)||(E.E_nat<r->E_nat_min)){
    r->E_nat_min=E.E_nat;
    r->pdbbest=i_pdb;
    r->seqid_mod=seq_id; r->L=L;
  }
  if(N_disulf_seq>r->N_disulf){
    r->N_disulf=N_disulf_seq;
    r->Len_disulf=Len_disulf_seq;
    r->Log_len_disulf=Log_len_disulf_seq;
  }
  if(N_helix_seq>r->N_helix){
    r->N_helix=N_helix_seq;
  }
  if(N_strand_seq>r->N_strand){
    r->N_strand=N_strand_seq;
  }

  /*r->DeltaG_1+=E.DeltaG;
  r->DeltaG_2+=E.DeltaG*E.DeltaG;
  if((r->n==0)||(E.DeltaG<r->DeltaG_min)){
    r->DeltaG_min=E.DeltaG; r->G_nat=E.E_nat;
    r->pdbbest=i_pdb; r->seqid_mod=seq_id; r->L=L;
    }*/

  if((r->n==0)||(seq_id > r->seqid_opt))r->seqid_opt=seq_id;
  (r->n)++;
}

void Print_result(struct result *r, double *DDG_1, double *DDG_2, int *n,
		  char *name, char **name_pdb, float DG_wt, float E_nat_wt,
		  float hydro_wt, float fit_wt, struct REM E,
		  FILE *file_out)
{
  if(r->n==0)return;

  //int L=E.L_str; // Only structured length
  //float S_C=sC1*L+sC0;

  float E_nat=r->E_nat_min;
  float DG=E.DeltaG+E_nat-E.E_nat;
  float fit=-log(1+exp(DG));
  if(DDG_1){
    (*n)++;
    *DDG_1+=DG;
    *DDG_2+=DG*DG;
  }

  r->E_nat_1/=r->n;
  float s=r->E_nat_2-r->n*r->E_nat_1*r->E_nat_1;
  if(r->n>1)s=sqrt(s/(r->n-1));

  fprintf(file_out, "%3d", r->L);
  if(DG_wt==0){ // wild type
    fprintf(file_out, "\t%.4g", DG/r->L);
    fprintf(file_out, "\t%.2g", s/r->L);
    fprintf(file_out, "\t%.4g", r->E_nat_min/r->L);
    if(r->n>1)fprintf(file_out, "\t%.3f", r->seqid_opt);
    fprintf(file_out, "\t%.3f", r->seqid_mod);
    fprintf(file_out, "\t%d", r->N_disulf);
    fprintf(file_out, "\t%d", r->Len_disulf);
    fprintf(file_out, "\t%.4g", r->Log_len_disulf);
    fprintf(file_out, "\t%d", r->N_helix);
    fprintf(file_out, "\t%d", r->N_strand);
    fprintf(file_out, "\t%.4f", r->hydro);
    fprintf(file_out, "\t%.3g", fit);
    fprintf(file_out, "\t%.3f", 0.00);
  }else{
    float DDG=DG-DG_wt;
    fprintf(file_out, "\t%.4g", DDG);
    fprintf(file_out, "\t%.2g", s);
    fprintf(file_out, "\t%.4g", E_nat-E_nat_wt);
    if(r->n>1)fprintf(file_out, "\t%.3f", r->seqid_opt);
    fprintf(file_out, "\t%.3f", r->seqid_mod);
    fprintf(file_out, "\t%d", r->N_disulf);
    fprintf(file_out, "\t%d", r->Len_disulf);
    fprintf(file_out, "\t%.4g", r->Log_len_disulf);
    fprintf(file_out, "\t%d", r->N_helix);
    fprintf(file_out, "\t%d", r->N_strand);
    fprintf(file_out, "\t%.4f", r->hydro);
    fprintf(file_out, "\t%.3g", fit);
    if(r->nmut)DDG/=r->nmut;
    fprintf(file_out, "\t%.3f", DDG);
  }

  fprintf(file_out, "\t%s", name_pdb[r->pdbbest]);
  fprintf(file_out, "\t%s", name);
  fprintf(file_out, "\n");
}

void Print_DG(struct REM E, struct result *r, short *aa_seq, int *i_sec,
		char *sec_str, char *name_pdb, int ini, char *name)
{
  char nameout[200];
  sprintf(nameout, "%s_DG.dat", name_pdb);
  char mode[4];
  if(ini){strcpy(mode, "w");}else{strcpy(mode, "a");}
  FILE *file_out=fopen(nameout, mode);
  if(ini){
    printf("Writing %s\n", nameout);
    if(SEC_STR){
      fprintf(file_out, "# seq:     ");
      for(int i=0; i<E.L; i++)fprintf(file_out, "%c", Amin_code(aa_seq[i]));
      fprintf(file_out, "\n");
      fprintf(file_out, "# sec_str: ");
      for(int i=0; i<E.L; i++)fprintf(file_out, "%c", sec_str[i]);
      fprintf(file_out, "\n");
      fprintf(file_out, "# sec_str: ");
      for(int i=0; i<E.L; i++)fprintf(file_out, "%c", SEC_EL[i_sec[i]]);
      fprintf(file_out, "\n");

    }
    fprintf(file_out,"# Extensive terms are divided times L\n");
    fprintf(file_out,"# T= %.2f S_C= %.3g*L+%.3g\n", TEMP, sC1, sC0);
    fprintf(file_out,
	    "#E_nat\t2=E_loc\t3=E_loc_misf"
	    "\t4=N_disulf\t5=Len_disulf\t6=Log_len_disulf"
	    "\t7=N_helix\t8=N_strand"
	    "\t9=E1\t10=E2/2\t11=E3/6\t12=Tf2\t13=Tf3"
	    "\t14=L\t15=DG(T,S_C)\tmut\n");
  }

  //float L=E.L_str; // Only structured length
  float L=E.L; // Only structured length
  float S_C=sC1*L+sC0;
  float E_nat=r->E_nat_min;
  float DG=E.DeltaG+E_nat-E.E_nat;

  fprintf(file_out,
	  "%.4g\t%.4g\t%.4g"
	  "\t%.4g\t%.4g\t%.4g"
	  "\t%.4g\t%.4g"
	  "\t%.4g\t%.4g\t%.4g"
	  "\t%.2g\t%.2g\t"
	  "%.0f\t%.4g\t%s\n",
	  E_nat/L, E.E_loc/L, E.E_loc_misf/L,
	  N_disulf_seq/L, Len_disulf_seq/L, Log_len_disulf_seq/L,
	  N_helix_seq/L, N_strand_seq/L,
	  E.E1/L, E.E2/(2*L), E.E3/(6*L),
	  sqrt(E.E2/(2*S_C)), E.Tf,
	  L, DG/L, name);

  fclose(file_out);
}

void Print_mut(short *aa_seq, short *aa_seq0, float **C_nat, int len_amm,
	       FILE *file_out)
{
  if(file_out==NULL)return;
  int i, j;
  for(i=0; i<len_amm; i++){
    if(aa_seq[i]==aa_seq0[i])continue;
    float *U=Econt[aa_seq[i]];
    float *U_old=Econt[aa_seq0[i]];
    for(j=0; j<len_amm; j++){
      if(C_nat[i][j]==0)continue;
      fprintf(file_out, "%c%c-%c%c-%d ",
	      Amin_code(aa_seq0[i]), Amin_code(aa_seq0[j]),
	      Amin_code(aa_seq[i]), Amin_code(aa_seq[j]), j);
      fprintf(file_out, " %.3f ", U[aa_seq[j]]-U_old[aa_seq0[j]]);
    }
  }
  fprintf(file_out, "\n");
}

void Modify_name_mut_interface(struct mutation *mut, struct protein prot,
			       struct residue *res, int L)
{
  if(mut->interface)return;
  for(int k=0; k<mut->nmut; k++){
    int pos=Find_site(mut->amm1_mut[k],mut->rmut[k],mut->chain[k], res, L);
    if(pos<1)continue;
    struct inter_chain *ic=prot.inter_chain;
    for(int j=0; j<prot.N_inter_chain; j++){
      if(ic->res1==pos || ic->res2==pos){
	char tmp[20];
	sprintf(tmp, "interface%c%c", ic->chain1, ic->chain2);
	strcat(mut->name, tmp);
	mut->interface=1;
	return;
      }
      ic++;
    }
  }
  return;
}

void Modify_name_del_interface(struct deletion *del, struct protein prot,
			       struct residue *res, int L)
{
  if(del->interface)return;
  for(int k=0; k<del->ndel; k++){
    int pos1=Find_site(del->aa1[k], del->res1[k], del->chain[k], res, L);
    int pos2=Find_site(del->aa2[k], del->res2[k], del->chain[k], res, L);
    if((pos1<0)||(pos2<0))continue;
    struct inter_chain *ic=prot.inter_chain;
    for(int j=0; j<prot.N_inter_chain; j++){
      if((ic->res1>=pos1 && ic->res1<= pos2) ||
	 (ic->res2>=pos1 && ic->res2<= pos2)){
	char tmp[20];
	sprintf(tmp, "interface%c%c", ic->chain1, ic->chain2);
	strcat(del->name, tmp);
	del->interface=1;
	return;
      }
      ic++;
    }
  }
  return;
}

void Print_interface(struct protein prot, short *aa_seq,
		     char *name_pdb, char *chain)
{
  if(prot.chain_all <= 1)return;
  float E_THR=-1;
  /* The threshold used in CAFFE is -1.2/0.485 Rigid body entropy in Kcal/mol
     divided by the conversion factor from contact energies to Kcal/mol */
  int nc=prot.chain_all, i, j;
  float E_inter[nc][nc]; int ncont[nc][nc];
  for(i=0; i<nc; i++)for(j=i+1; j<nc; j++){E_inter[i][j]=0; ncont[i][j]=0;}

  struct inter_chain *ic=prot.inter_chain;
  for(int k=0; k<prot.N_inter_chain; k++){
    i=ic->ichain1; j=ic->ichain2;
    if(j<i){int tmp=j; j=i; i=tmp;}
    if(i<0 || j>=nc){
      printf("ERROR in interface contact %d out of %d, "
	     "chains %d-%c and %d-%c of %d\n",
	     k, prot.N_inter_chain,
	     ic->ichain1, ic->chain1, ic->ichain2, ic->chain2, nc);
      exit(8);
    }
    if(ic->res1 <0 || ic->res1 >= prot.length ||
       ic->res2 <0 || ic->res2 >= prot.length){
      printf("ERROR in interface contact %d, res1= %d res2=%d of %d\n",
	     j, ic->res1, ic->res2, prot.length);
      exit(8);
    }
    E_inter[i][j]+=Econt_T[aa_seq[ic->res1]][aa_seq[ic->res2]];
    ncont[i][j]++;
    ic++;
  }

  char name_out[200];
  sprintf(name_out, "%s_complex.dat", name_pdb);
  FILE *file_out=fopen(name_out, "w");
  printf("Writing %s\n", name_out);
  fprintf(file_out, "Protein %s %d chains %d different chains\n",
	  name_pdb, nc, prot.nchain);
  int nclust=nc, iclust[nc], numclust[nc];
  for(i=0; i<nc; i++){iclust[i]=i; numclust[i]=1;}
  for(i=0; i<nc; i++){
    for(j=i+1; j<nc; j++){
      fprintf(file_out, "chains %c %c %2d contacts E_cont= %.3g ",
	      chain[i], chain[j], ncont[i][j], E_inter[i][j]);
      if(E_inter[i][j] < E_THR){
	fprintf(file_out, "Form complex\n");
	int cmin=iclust[i], cmax=iclust[j];
	if(cmin==cmax)continue;
	if(cmax<cmin){int tmp=cmin; cmin=cmax; cmax=tmp;}
	for(int k=0; k<nc; k++)if(iclust[k]==cmax)iclust[k]=cmin;
	numclust[cmin]+=numclust[cmax]; numclust[cmax]=0;
	nclust--;
      }else{
	fprintf(file_out, "Do not form complex\n");
      }
    }
  }
  fprintf(file_out, "%d complexes found\n", nclust);
  for(i=0; i<nc; i++){
    if(numclust[i]==0)continue;    
    fprintf(file_out, "Complex %d %d chains: ",i, numclust[i]);
    for(j=0; j<nc; j++)if(iclust[j]==i)fprintf(file_out, "%c", chain[j]);
    fprintf(file_out, "\n");
  }
  fclose(file_out);

}
