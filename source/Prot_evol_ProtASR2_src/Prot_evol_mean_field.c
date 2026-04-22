/* Program dna_evol.
   Needs a target structure and a starting DNA sequence.
   Attempts DNA mutations with given mutational bias.
   alpha, Z, energy are computed for the mutated sequence.
   The mutation is accepted or not according to selection criteria.
   The base composition is evaluated at every time step.
*/

int REMUT=0;  // Repeat the fit of mutation parameters?
int PRINT_TN=1; // Print Tajima-Nei divergence?

#define FILE_CODE_DEF "gen_code_ATGC.in"
// #define FILE_ENE_DEF  "energy.in"
#define VBT 0 // Verbatim

#include "REM.h"
#include "coord.h"
#include "gen_code.h"
#include "allocate.h"
#include "protein3.h"
#include "read_pdb.h"
#include "random3.h"           /* Generating random numbers */
#include "mutation.h"
#include "codes.h"
#include "input.h"
#include "get_para_pop_dyn.h"
#include "subst_exhaustive.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include "Print_meanfield.h"
#include "fits.h"
#include "meanfield.h"  //NEW! 

int ini_print;

/* New version 17/03/2014:
   Now energy parameters and genetic codes are given in header files
   energy_BKV.h gen_code.h
   It is possible to give as input either a PDB file (extension .pdb is
   needed) or a file with a nucleotide sequence and the name of the target
   pdb, which must be contained in the defeult input file structures.in
*/

//#define FILE_IN "Prot_evol.in"
#define N_CHAR 300          // Max. length of file names
#define EXT_AVE "_ave.dat"  // Extension for output file with averages
#define EXT_DNA "_dna.dat"  // Extension for output file with dna statistics
#define EXT_OUT "_stab.dat" // Extension for output file with folding
                            // thermodynamics and fitness


// Input codes
void Read_ene_par(char *, float **interactions);
int Read_ene_new(char *, float **interactions);
char *Read_sequence(int *len_dna, char *inseq);
unsigned long randomgenerator(void);
float *Get_counts(short *seq, int L, int NAA);

// Output
void Output_name(char *file_name, char *dir_out, char *prot_name,
		 float TEMP, float sU1, float sC1,  int MEANFIELD,
		 char *MODEL, float LAMBDA, int OPT_LAMBDA, //NEW
		 int NEUTRAL, int N_pop, float *mut_par);
FILE *open_file(char *, char *, short *, int, char *fit_def);
int Print_dna(char *, FILE *, int);
void Print_ave(FILE *file_ave,
	       long it_sum, long t_indip,
	       int N_pop,
	       double f_sum, double f_dev,
	       double E_sum, double E_dev,
	       double DG_sum, double DG_dev,
	       long num_syn_subst, long num_aa_subst,
	       long num_syn_mut, long num_aa_mut,
	       float seq_entr,
	       struct load mut_load,
	       struct load trans_load);
void Print_TN_div(short *aa_seq, short *aa_seq0, int L,
		  int num_aa_subst, FILE *file_out);

void Print_matrix(struct protein target);
FILE *Open_summary(char *name_file, struct REM E_wt,
		   float *mut_par, short *aa_seq);
void Get_mean(double *ave, double *err, float sum, float dev,
	      long it_sum, float t_indep);
void Record_evo(float *fit_evo, double fitness,
		float *DG_evo,  double DG,
		float *Enat_evo, double E_nat, int *nsam);
void Print_final(char *name_file, long it_sum,
		 float TEMP, float sU1, float sC1, float sC0,
		 int MEANFIELD, float Lambda, int N_pop,
		 float *mut_par, float tt_ratio,
		 //
		 float *fit_evo, float *DG_evo, 
		 float *Enat_evo, int nsam, int step,
		 double *f_ave, double f_dev, 
		 double *E_ave, double E_dev,
		 double *DG_ave, double DG_dev,
		 float seq_entr, double seq_entr_dev,
		 struct load *mut_load, struct load *trans_load,
		 long num_syn_subst, long num_aa_subst,
		 long num_syn_mut, long num_aa_mut,
		 float *dN_dS, float *accept,
		 long **nuc_evo, int len_dna, int sample);

// Calculations
void Initialize_load(struct load *load);
void Sum_loads(struct load *load, struct load *load1);
void WildType_Mutants(float *Lambda_opt, double *DG_ave, 
		      float **P_WT_ia, float *P_mut_a,
		      int **C_nat, int *i_sec, short *aa_seq,
		      struct REM *E_wt, char *name_file, FILE *file_summ);
void Compute_P_WT(float **P_WT_ia, float Lambda,
		  float **DG_mut, float *P_mut_a, int L);
int Selection(float fitness, float fitness_old, int N_pop);
int Detailed_balance(float *p, int xend, int xini);
float Sequence_entropy(double **aa_distr, int L);
float Sequence_entropy_mut(float *mut_par, char **codon, char *coded_aa);
void Compute_freq_codons(float *mut_par, float *freq_aa,
			 char **codon, char *coded_aa);
void Compute_load(double *Tload_sum, double *Tload_dev,
		  double *Mload_sum, double *Mload_dev, int *Nload,
		  struct REM *E_wt, int **C_nat, int *i_sec, short *aa_seq,
		  float fitness_wt, char *dna_seq, int len_dna, 
		  char **codon, char *coded_aa);

// Input parameters
// A: Input files defining the protein
static char file[N_CHAR], seq_name[N_CHAR], dir_out[N_CHAR];
static char seq_name[N_CHAR];
// B: Thermodynamic parameters
float sC1=0.065, sC0=0, sU1=0.140; // Configuration entropy
int REM=2;   // Use 1st (1), 2nd (2) and 3rd (3) moment of misfolded energy 
float TEMP, S_C, S_U;
float SEC_STR=1; // Coefficient for local interactions
// C: Selection model
static long IT_MAX=0; // Number of iterations
static int NEUTRAL=1; // Neutral (1) versus continous (0) fitness
float DG_THR_COEFF=0.95; // Threshold in DeltaG for neutral selection
static int N_pop=100;  // Population size for continuous selection model
// C2: Mean-field model
int MEANFIELD=1;  // Compute site-specific mean-field distributions
int OPT_LAMBDA=1; // 1=Optimize Lambda by maximum likelihood
float LAMBDA=0.9; // Lambda parameter for meanfield if OPT_LAMBDA=0
char MODEL[N_CHAR]="ALL"; // Type of mean-field model (ALL, NAT, DG)
float DG_OPT=-1;  // DG target of the optimization if MODEL="DG"
// D1: Mutation model, P_mut
int CpG=1;       // Enhance mutation rate at CpG dinucleotides
int GET_FREQ=2;  // 0= Get nucleotide frequencies from input
                 // 1= Fit nucleotide frequencies from AA sequence
                 // 2= Get P_mut[a] from fit plus AA sequence
                 // 3= Get P_mut[a] from AA sequence
float mut_par[MUTPAR], tt_ratio=4, kCpG=4;
// D1: Mutation model, exchangeability
char MATRIX[40]="WAG"; // Empirical exchangeability matrix
char EXCHANGE='F'; // exchangeability model. M=MUT F=FLUX Q=RATE E=EXCH
float TWONUC=0; // Rate of two to one nucleotide substitutions  if MUT
// E: Output
int FORMAT=1;   // PAML format
int PRINT_E=0; // Print exchangeability matrix for all sites?

// Derived data
unsigned long iran;
static int len_amm, len_dna;
static int count[4];
static float rate[4];
float DG_thr;

int main(int argc, char **argv){


  /***********************
          INPUT
  ************************/
  // Input files
  //char Input_dir[N_CHAR];
  char file_pdb[N_CHAR], file_seq[N_CHAR]="";
  char chain[]="\0\0\0";
  char name_file[N_CHAR];
  char FILE_CODE[N_CHAR];
  char *file_ali=NULL;

  // Genetic code
  // char *codon[64], coded_aa[64], name_code[200];

  /***********************
         SEQUENCES
  ************************/
  char *dna_seq=NULL;
  short *nuc_seq=NULL;

  /***********************
         Wild Type
  ************************/
  double fitness_wt;

  /***********************
          DUMMY
  ************************/
  int i, j, a;

  /******************** Input operations   ************************/
  TEMP=0.5; // default
  sprintf(FILE_CODE, FILE_CODE_DEF);
  Get_para(argc, argv, file_pdb, chain, file_seq,
	   FILE_STR, &TEMP, &sU1, &sC0, &sC1, &REM, &SEC_STR, &REMUT,
	   &IT_MAX, &NEUTRAL, &MEANFIELD, &OPT_LAMBDA, &LAMBDA, &DG_OPT,
	   &N_pop, &GET_FREQ, mut_par, &tt_ratio, &kCpG, &TWONUC,
	   &EXCHANGE, MATRIX, &FORMAT, &PRINT_E, &file_ali, MODEL, dir_out);

  // Random numbers
  iran=randomgenerator();
  InitRandom( (RANDOMTYPE)iran);

  // Compute contact matrix from PDB file
  struct residue *res; int imod=-1;
  len_amm=Get_pdb(&target, &res, file_pdb, chain, imod);
  int **C_nat, *i_sec; short *aa_seq=NULL;
  if(len_amm > 0){
    aa_seq=target.aa_seq;
    i_sec=target.i_sec;
    C_nat=Fill_C_nat(len_amm, target.contact);
    Print_matrix(target);
  }else{
    // or read contact matrix from precomputed matrices
    C_nat=Get_target(FILE_STR, file_pdb, &len_amm);
    if(len_amm>0){
      printf("contact matrix %s found in %s, %d residues\n",
	     file_pdb, FILE_STR, len_amm);
    }else{
      printf("contact matrix %s not found in %s\n",file_pdb, FILE_STR);
      exit(8);
    }
  }
    
    
  // DNA sequence
  if(file_seq[0]!='\0'){
    // Read from file
    char inputseq[N_CHAR]; strcpy(inputseq, file_seq);
    if(Check_file(inputseq)==0){
      printf("WARNING, DNA sequence file %s not found\n", file_seq);
    }else{
      dna_seq=Read_sequence(&len_dna, inputseq);
      if((dna_seq==NULL)||(len_dna==0)||
	 ((aa_seq!=NULL)&&(len_amm>0)&&
	  (Compare_amm_dna(dna_seq,len_dna,aa_seq,len_amm,codon,coded_aa)==0))){
	len_dna=0;
      }
    }
  }

  // If AA sequence not given, get it from DNA sequence
  if(aa_seq==NULL){
    if(len_dna){
      len_amm=len_dna/3;
      aa_seq=malloc(len_amm*sizeof(short));
      Translate_new(dna_seq, aa_seq, len_amm, codon, coded_aa);
    }else{
      printf("ERROR, no sequence given\n"); exit(8);
    }
  }

  // If DNA sequence not given, get it from AA sequence
  int idna=1;
  if(len_dna<=0){
    idna=0;
    printf("Extracting DNA sequence from AA sequence\n");
    dna_seq=Extract_dna(&len_dna, len_amm, aa_seq, codon, coded_aa);
  }
  nuc_seq=malloc(len_dna*sizeof(short));
  for(i=0; i<len_dna; i++)nuc_seq[i]=Code_nuc(dna_seq[i]);
  float *num_dna=NULL;
  if(idna)num_dna=Get_counts(nuc_seq, len_dna, 4);
    

    
    
    
    // Miguel
    FILE *fpmy;
    fpmy = fopen ( "PDBseq.txt", "w" );
    if (fpmy==NULL) {fputs ("File error of PDBseq.txt",stderr); exit (1);}
    
    printf("\n\nNumberOfSites=%d\n", len_amm);
    printf("Sequence="); // Miguel
    for(i=0; i<len_amm; i++){printf("%c", AMIN_CODE[aa_seq[i]]);}
    printf("\n\n\n");
    
    fprintf(fpmy, "NumberOfSites=%d\n", len_amm);
    fprintf(fpmy, "Sequence="); // Miguel
    for(i=0; i<len_amm; i++){fprintf(fpmy, "%c", AMIN_CODE[aa_seq[i]]);}
    
    fclose ( fpmy );
    // Miguel
    
    
    

  /**************************  End INPUT  ****************************/

  /*********** Compute background amino acid frequencies P_mut *******/

  /* Obtain mutational model:
     if (GET_FREQ==0), from input data
     if (idna) Obtain nucleotide frequencies from input DNA sequence
     else get mutation parameters fitting amino acid frequencies
  */
  float P_mut_a[20], *num_aa=Get_counts(aa_seq, len_amm, 20);
  if(file_ali)Sum_aa(num_aa, file_ali);
  if(kCpG<tt_ratio)kCpG=tt_ratio;
  mut_par[4]=kCpG; mut_par[5]=tt_ratio; mut_par[6]=TWONUC;
  float **Q_cod=Allocate_mat2_f(64, 64), P_cod[64];
  Get_mut_par(mut_par, P_mut_a, P_cod, Q_cod, GET_FREQ,
	      num_aa, len_amm, num_dna, target.name, 0);
  if(GET_FREQ){
    kCpG=mut_par[4]; tt_ratio=mut_par[5]; TWONUC=mut_par[6];
  }

  /************************* Output files ***************************/
  Output_name(name_file, dir_out, target.name, TEMP, sU1, sC1, MEANFIELD,
	      MODEL, LAMBDA, OPT_LAMBDA, NEUTRAL, N_pop, mut_par);

  /******************** Thermodynamics of wild type **********************/
  S_C=sC0+len_amm*sC1; S_U=len_amm*sU1;
  struct REM E_wt;
  Initialize_E_REM(&E_wt, len_amm, REM, TEMP, S_C, S_U, FILE_STR, 0);

  E_wt.DeltaG=Compute_DG_overT_contfreq(&E_wt, aa_seq, C_nat, i_sec);
  printf("T= %.1f sU1= %.2f sC1= %.2f sC0= %.2f L=%d\n",
	 TEMP, sU1, sC1, sC0, len_amm);
  printf("DeltaG/T= %.2f\n", E_wt.DeltaG);

  char nameout[200];
  sprintf(nameout, "%s%s_DeltaG.dat", dir_out, target.name);
  double T0=Print_DG_contfreq(&E_wt, nameout);
  printf("Temperature with smallest DG: T=%.2f\n", T0);
  Test_contfreq(&E_wt, aa_seq, name_file);
  if(E_wt.DeltaG > 0){
    printf("WARNING, unstable target structure!\n");
    if(NEUTRAL){
      printf("This condition is illegal with NEUTRAL evolution, ");
      printf("please change the input file!\n");
      exit (8);
    }
  }

  // Mean field computation.
  float **P_MF_ia=NULL;
  if(MEANFIELD){
    FILE *file_summ=Open_summary(target.name, E_wt, mut_par, aa_seq);

    P_MF_ia=Allocate_mat2_f(len_amm, 20);
    float **P_WT_ia=Allocate_mat2_f(len_amm, 20), Lambda;
    double DG_MF=0, DG_WT=0;
    struct MF_results mut_res;

    int NMUT=1; if(REMUT && GET_FREQ)NMUT=2; // Change mutation model    
    for(int iter=0; iter<NMUT; iter++){
 
      // Mutation model
      for(i=0; i<len_amm; i++)for(a=0; a<20; a++)P_MF_ia[i][a]=P_mut_a[a];
      mut_res.h=0; for(a=0; a<20; a++)mut_res.h+=P_mut_a[a]*hydro[a];
      Test_distr(&mut_res, P_MF_ia, aa_seq, len_amm, C_nat, i_sec, E_wt);
      Print_results(mut_res, "mut", file_summ);

      // Pairwise model
      WildType_Mutants(&Lambda, &DG_WT, P_WT_ia, P_mut_a,
		       C_nat, i_sec, aa_seq, &E_wt, name_file, file_summ);


      // Mean-field model
      if(OPT_LAMBDA==0){ 
	Fixed_Lambda(&DG_MF, P_MF_ia, P_mut_a, C_nat, i_sec, aa_seq,
		     E_wt, LAMBDA, name_file);
      }else{
	Optimize_Lambda(&LAMBDA, &DG_MF, P_MF_ia, P_mut_a, C_nat, i_sec,aa_seq,
			E_wt, DG_OPT, GET_FREQ, MODEL, name_file, file_summ);
      }

      if(iter==(NMUT-1))break;

      // Change mutation model
      float fmut[20], P_sel[20]; double sum=0;
      Compute_P_sel(P_sel, P_MF_ia, P_mut_a, len_amm);
      for(a=0; a<20; a++){fmut[a]=num_aa[a]/P_sel[a]; sum+=fmut[a];}
      sum/=len_amm; for(a=0; a<20; a++)fmut[a]/=sum;
      Get_mut_par(mut_par, P_mut_a, P_cod, Q_cod, GET_FREQ,
		  fmut, len_amm, num_dna, target.name, 1);
      kCpG=mut_par[4]; tt_ratio=mut_par[5]; TWONUC=mut_par[6];
      Output_name(name_file, dir_out, target.name, TEMP, sU1, sC1, MEANFIELD,
		  MODEL, LAMBDA, OPT_LAMBDA, NEUTRAL, N_pop, mut_par);
      //sprintf(name_file, "%s_iter%d", name_file, iter+2);
    }
    fclose(file_summ);

    Print_profiles(P_MF_ia,P_mut_a,DG_MF,LAMBDA,aa_seq,len_amm,name_file,"_MF");
    Print_exchange(P_MF_ia, P_mut_a, P_cod, Q_cod, tt_ratio, TWONUC,
		   aa_seq, res, C_nat, target.sec_str, len_amm,
		   name_file, "_MF", FORMAT, EXCHANGE, MATRIX, PRINT_E);
    Print_profiles(P_WT_ia, P_mut_a, DG_WT, Lambda, aa_seq, len_amm,
		   name_file, "_WT");
    Print_exchange(P_WT_ia, P_mut_a, P_cod, Q_cod, tt_ratio, TWONUC,
		   aa_seq, res, C_nat, target.sec_str, len_amm,
		   name_file, "_WT", FORMAT, EXCHANGE, MATRIX, PRINT_E);
  }

  /******************** End Mean-field  *********************/
  if(IT_MAX==0)return(0);

  /***********************************************************

                      Simulations of evolution
 
  *************************************************************/
  // S samples are simulated each for IT_MAX iterations

  int Samples=20;

  // Accumulate samples
  double fit1_all=0, fit2_all=0, DG1_all=0, DG2_all=0,
    Enat1_all=0, Enat2_all=0, dNdS1_all=0, dNdS2_all=0,
    accept1_all=0, accept2_all=0, seqentr1_all=0, seqentr2_all=0;
  long it_all=0;
  struct load mut_load_all, trans_load_all;
  Initialize_load(&mut_load_all);
  Initialize_load(&trans_load_all);

  // Global variables valid for all independent samples
  // Number of iterations
  int it_print=10; if(it_print > IT_MAX*0.1)it_print= IT_MAX*0.1;
  int it_trans=it_print;
  int NPRINT= 1; if(it_print)NPRINT+=(IT_MAX/it_print);
  int MUTMAX= N_pop*10;  // Exhaustive search after MUTMAX unfixed mutations  

  // Amino acid distributions and entropy
  // Reference sequence
  short *aa_seq0=malloc(len_amm*sizeof(short));
  for(i=0; i<len_amm; i++)aa_seq0[i]=aa_seq[i];
  double **aa_distr=malloc(len_amm*sizeof(double *));
  double **aa_distr0=malloc(len_amm*sizeof(double *));
  double **aa_distr_all=malloc(len_amm*sizeof(double *));
  for(i=0; i<len_amm; i++){
     aa_distr[i]=malloc(20*sizeof(double));
     for(j=0; j<20; j++)aa_distr[i][j]=0;
     aa_distr0[i]=malloc(20*sizeof(double));
     for(j=0; j<20; j++)aa_distr0[i][j]=0;
     aa_distr_all[i]=malloc(20*sizeof(double));
     for(j=0; j<20; j++)aa_distr_all[i][j]=0;
  }
  float seq_entr_mut=Sequence_entropy_mut(mut_par, codon, coded_aa);

  // Nucleotide distributions at each codon position
  long **nuc_evo=malloc(len_dna*sizeof(long *));
  for(i=0; i<len_dna; i++){
    nuc_evo[i]=malloc(4*sizeof(long));
    for(j=0; j<4; j++)nuc_evo[i][j]=0;
  }

  // Averages
  float *fit_evo=malloc(NPRINT*sizeof(float));
  float *DG_evo=malloc(NPRINT*sizeof(float));
  float *Enat_evo=malloc(NPRINT*sizeof(float));

  // Thermodynamic computations
  struct REM E_mut;
  Initialize_E_REM(&E_mut, E_wt.L, E_wt.REM,
		   E_wt.T, E_wt.S_C, E_wt.S_U, FILE_STR, 0);

  // Fitness
  if(MEANFIELD)Divide_by_Pmut(P_MF_ia, P_mut_a, len_amm, 20);
  int EXHAUSTIVE=0;
  if((MEANFIELD==0)&&(N_pop > len_amm*10))EXHAUSTIVE=1;
  int INI_EXH=EXHAUSTIVE;
  float Inverse_pop=1./N_pop;

  /***********************
         MUTANT
  ************************/
  double fitness_mut;
  char dna_new;
  int nuc_mut, nuc_new, res_mut=0, aa_new;

  /***********************
          OUTPUT
  ************************/

  char fit_def[100];
  if(NEUTRAL==0){strcpy(fit_def, "F=1./(1+exp(DeltaG))");}
  else{sprintf(fit_def,"F=Theta(%.2f*DeltaG_PDB-DeltaG)",DG_THR_COEFF);}
  // FILE *file_dna=open_file(name_file, EXT_DNA, aa_seq, 2, fit_def);
  FILE *file_ave =open_file(name_file, EXT_AVE, aa_seq, 0, fit_def);
  FILE *file_stab=open_file(name_file, EXT_OUT, aa_seq, 1, fit_def);
  char name_div[200]; FILE *file_div;
  if(PRINT_TN){
    sprintf(name_div, "%s_TN_div.dat", target.name);
    file_div=fopen(name_div, "a");
    fprintf(file_div, "# Protein %s L=%d\n", target.name, len_amm);
    fprintf(file_div, "#Tajima-Nei_divergence num_subst\n");
    printf("Writing %s\n", name_div);
  }

  /****************************************************
                 S samples of evolution
  *****************************************************/

  for(int is=0; is<Samples; is++){

    // Start the trajectory

    // Number of substitutions
    int aa_subst=0, synonymous;
    int it1=0;
    long it_subst=0;
    // Mutations per one aa substitution
    long tot_mut=0, naa_mut=0, nsyn_mut=0, nsyn_subst=0;

    // Count substitutions
    long num_syn_subst=0, num_aa_subst=0;
    long num_syn_mut=0, num_aa_mut=0;

    // Average stability and fitness
    long it_sum=0;
    double E_nat_ave=0, f_ave=0, DG_ave=0;
    double E_nat_dev=0, f_dev=0, DG_dev=0;

    // Sequence entropy
    double seq_entr_sum=0, seq_entr_dev=0, entr_dev=0;
    float seq_entr;
    int n_seq_entr=0;

    for(i=0; i<len_amm; i++)aa_seq0[i]=aa_seq[i];
    for(i=0; i<len_amm; i++){
      for(j=0; j<20; j++)aa_distr[i][j]=0;
    }

    // Loads
    struct load mut_load, trans_load;
    Initialize_load(&mut_load);
    Initialize_load(&trans_load);

    // Wild-type stability
    E_wt.DeltaG=
      Compute_DG_overT_contfreq(&E_wt, aa_seq, C_nat, i_sec);    
    // Wild-type fitness
    if(NEUTRAL==0){
      fitness_wt=1./(1+exp(E_wt.DeltaG));
    }else{
      fitness_wt=1;
      DG_thr=E_wt.DeltaG*DG_THR_COEFF;
    }
    fprintf(file_stab, "#\n#\n# Evolutionary trajectory %d\n#\n#\n",is+1);
    fprintf(file_stab, "%3d  %c%c  %.3f  %.3f %7.4f",
	    0, AMIN_CODE[aa_seq[0]], AMIN_CODE[aa_seq[0]],
	    E_wt.E_nat, E_wt.DeltaG, fitness_wt);

    // Record fitness, energy and DeltaG
    int nsam=0;
    Record_evo(fit_evo, fitness_wt, DG_evo, E_wt.DeltaG,
	       Enat_evo, E_wt.E_nat, &nsam);
    
    // Mutation rates
    Ini_count(dna_seq, len_dna, count);
    Compute_rates(rate, mut_par, tt_ratio);

    // Simulations
    while(1){
    
      if(EXHAUSTIVE){
	// Extract next substitution (requires full exploration of sequence
	// neighbors, only convenient if N is large).
	Substitute_exhaustive(&res_mut, &aa_new, &nuc_mut, &nuc_new,
			      &mut_load, &trans_load, it_subst-it_trans,
			      &fitness_mut, &E_mut, &E_wt, 
			      &naa_mut, &nsyn_mut, &nsyn_subst,
			      NEUTRAL, DG_thr, N_pop, C_nat, i_sec,
			      aa_seq, dna_seq, nuc_seq, len_dna,
			      codon, coded_aa, tt_ratio, mut_par);
	tot_mut=naa_mut+nsyn_subst+nsyn_mut+1;
	dna_new=Nuc_code(nuc_new);
	if(it_subst > it_trans){
	  num_aa_mut += naa_mut;
	  num_syn_mut+= nsyn_mut;
	  num_syn_subst += nsyn_subst;
	}
	goto Update_aa;
      }

      tot_mut++;
      if(tot_mut>=MUTMAX){
	printf("WARNING, too many unfixed mutations > %d\n", MUTMAX);
	printf("NEUTRAL= %d MEANFIELD= %d N= %d T= %.2f\n",
	       NEUTRAL, MEANFIELD, N_pop, TEMP);
	if(MEANFIELD==0){
	  printf("Performing exhaustive search of mutants\n");
	  EXHAUSTIVE=1; INI_EXH=1;
	  aa_subst=1; synonymous=0;
	  continue;
	}else{
	  exit(8);
	}
      }

      // Individual mutation
      synonymous=Mutate_seq(dna_seq, len_dna, codon, coded_aa,
			    aa_seq, len_amm, mut_par, tt_ratio, count, rate,
			    &nuc_mut, &dna_new, &res_mut, &aa_new);

      if(synonymous > 0){
	nsyn_mut++;
	if(it_subst >= it_trans)num_syn_mut++; 
	//if((MEANFIELD==0)&&(RandomFloating() > Inverse_pop))continue;
	if((RandomFloating() > Inverse_pop))continue;
	nsyn_subst++;
	if(it_subst >= it_trans)num_syn_subst++;
	nuc_new=Code_nuc(dna_new);
	goto Update_dna;

      }else if(synonymous==0){

	if(it_subst >= it_trans)num_aa_mut++;

	// Compute stability and fitness
	Copy_E_REM(&E_mut, &E_wt);
	E_mut.DeltaG=
	  Mutate_DG_overT_contfreq(&E_mut, aa_seq, C_nat, i_sec,
				   res_mut, aa_new);
	if(NEUTRAL){
	  if(E_mut.DeltaG<DG_thr){fitness_mut=1;}else{fitness_mut=0;}
	}else{
	  fitness_mut=1./(1+exp(E_mut.DeltaG));
	}
	// Selection
	if(MEANFIELD==0){
	  aa_subst=Selection(fitness_mut, fitness_wt, N_pop);
	}else{
	  aa_subst=Selection(fitness_mut, fitness_wt, N_pop);
	  //aa_subst=Detailed_balance(P_MF_ia[res_mut],aa_new,aa_seq[res_mut]);
	}
	if(aa_subst<=0)continue;

      }else{ 
	continue;  // Stop codon
      }

    Update_aa:
      // A non-synonymous substitution has occurred

      // Print thermodynamic properties (also in transient)
      it_subst++;
      fprintf(file_stab, " %2ld %3ld\n", nsyn_subst, tot_mut);
      if(INI_EXH){
	fprintf(file_stab, "# Exhaustive search\n"); INI_EXH=0;
      }
      fprintf(file_stab, "%3d  %c%c  %.3f  %.3f %7.4f",
	      res_mut, AMIN_CODE[aa_seq[res_mut]], AMIN_CODE[aa_new],
	      E_mut.E_nat, E_mut.DeltaG, fitness_mut);
      fflush(file_stab);
      if((it_subst/it_print)*it_print == it_subst){
	Record_evo(fit_evo, fitness_wt, DG_evo, E_mut.DeltaG,
		   Enat_evo, E_mut.E_nat, &nsam);
      }
    
      // Statistics after the transient
      if(it_subst >= it_trans){
	// Average old amino acid sequence, weight=tot_mut
      
	// Update amino acid distribution
	for(i=0; i<len_amm; i++)aa_distr0[i][aa_seq[i]]+=tot_mut;
      
	// Averages
	E_nat_ave +=tot_mut*E_wt.E_nat;
	E_nat_dev+=tot_mut*E_wt.E_nat*E_wt.E_nat;
	DG_ave+=tot_mut*E_wt.DeltaG;
	DG_dev+=tot_mut*E_wt.DeltaG*E_wt.DeltaG;
	f_ave+=tot_mut*fitness_wt;
	f_dev+=tot_mut*fitness_wt*fitness_wt;
	it_sum+=tot_mut;
     
	num_aa_subst++; it1++;

	// Print averages
	if(it1==it_print){
	  it1=0;
	  seq_entr=Sequence_entropy(aa_distr0, len_amm)-seq_entr_mut;
	  seq_entr_sum+=seq_entr; seq_entr_dev+=seq_entr*seq_entr;
	  n_seq_entr++;
	  for(i=0; i<len_amm; i++){
	    double *aa0=aa_distr0[i];
	    for(j=0; j<20; j++){aa_distr[i][j]+=aa0[j]; aa0[j]=0;}
	  }
	  Print_ave(file_ave, it_sum, num_aa_subst, N_pop,
		    f_ave, f_dev, E_nat_ave, E_nat_dev, DG_ave, DG_dev,
		    num_syn_subst, num_aa_subst, num_syn_mut, num_aa_mut,
		    seq_entr, mut_load, trans_load);
	  if(PRINT_TN){
	    aa_seq[res_mut]=aa_new;
	    Print_TN_div(aa_seq, aa_seq0, len_amm, it_subst, file_div);
	  }
	    
	} // end print
      }

      // Update_AA:
      aa_seq[res_mut]=aa_new;
      Copy_E_REM(&E_wt, &E_mut);
      fitness_wt=fitness_mut;
      tot_mut=0; nsyn_subst=0;

    Update_dna:
      if(aa_subst ||(synonymous > 0)){
	count[nuc_new]++;
	count[nuc_seq[nuc_mut]]--;
	for(i=0; i<len_dna; i++)nuc_evo[i][nuc_seq[i]]+=nsyn_mut;
	nsyn_mut=0;
	nuc_seq[nuc_mut]=nuc_new;
	dna_seq[nuc_mut]=dna_new;
      }

      if(num_aa_subst == IT_MAX)break;
      
    }
    fprintf(file_stab, " %2ld %3ld\n", nsyn_subst, tot_mut);


    /*********************  End of simulation  **************************/

    // Loads
    /*Compute_load(&Tload_sum, &Tload_dev, &Mload_sum, &Mload_dev, &Nload,
      E_wt, C_nat, i_sec, aa_seq, fitness_wt, dna_seq, len_dna,
      codon, coded_aa);*/

    // Entropy
    for(i=0; i<len_amm; i++){
      double *aa0=aa_distr0[i];
      for(j=0; j<20; j++){aa_distr[i][j]+=aa0[j]; aa0[j]=0;}
    }
    seq_entr=Sequence_entropy(aa_distr, len_amm)-seq_entr_mut;
    for(i=0; i<len_amm; i++){
      double *aa1=aa_distr[i];
      for(j=0; j<20; j++){aa_distr_all[i][j]+=aa1[j]; aa1[j]=0;}
    }

    if(n_seq_entr > 1){
      entr_dev=seq_entr_dev-seq_entr_sum*seq_entr_sum/n_seq_entr;
      entr_dev=sqrt(entr_dev)/n_seq_entr;
    }
    
    it_sum+=tot_mut;
    float dN_dS, accept;
    Print_final(name_file, it_sum, TEMP, sU1, sC1, sC0, MEANFIELD, LAMBDA,
		N_pop, mut_par, tt_ratio,
		fit_evo, DG_evo, Enat_evo, nsam, it_print,
		&f_ave, f_dev, &E_nat_ave, E_nat_dev, &DG_ave, DG_dev,
		seq_entr, entr_dev,
		&mut_load, &trans_load,
		num_syn_subst, num_aa_subst, num_syn_mut, num_aa_mut,
		&dN_dS, &accept, nuc_evo, len_dna, is);

    // Average samples
    it_all+=it_sum;
    fit1_all+=f_ave; fit2_all+=f_ave*f_ave;
    DG1_all+=DG_ave; DG2_all+=DG_ave*DG_ave;
    Enat1_all+=E_nat_ave; Enat2_all+=E_nat_ave*E_nat_ave;
    Sum_loads(&trans_load_all, &trans_load);
    Sum_loads(&mut_load_all, &mut_load);
    dNdS1_all+=dN_dS; dNdS2_all+=dN_dS*dN_dS;
    accept1_all+=accept; accept2_all+=accept*accept;
    seqentr1_all+=seq_entr; seqentr2_all+=seq_entr*seq_entr;

  }

  /***************************************************
                     End of samples
  ***************************************************/

  FILE *file_out=open_file(name_file,"_samples.dat",aa_seq,0, fit_def);
  fprintf(file_out, "# %d samples\n", Samples);
  fprintf(file_out, "# what mean s.e.\n");
  double x, dx; int S=Samples;
  Get_mean(&x, &dx, DG1_all, DG2_all, S, S); DG1_all=x;
  fprintf(file_out, "DeltaG/T\t%.5g\t%.2g\n", x, dx);
  fprintf(file_out, "DeltaG\t%.5g\t%.2g\n", x*TEMP, dx*TEMP);
  Get_mean(&x, &dx, Enat1_all, Enat2_all, S, S);
  fprintf(file_out, "E_nat\t%.5g\t%.2g\n", x, dx);
  Get_mean(&x, &dx, fit1_all, fit2_all, S, S);
  fprintf(file_out, "Fitness\t%.5g\t%.2g\n", x, dx);

  Get_mean(&x, &dx,mut_load_all.df_ave, mut_load_all.df_dev, S, S);
  fprintf(file_out, "Mut.load.fitness\t%.5g\t%.2g\n", x, dx);
  Get_mean(&x, &dx, mut_load_all.dG_ave, mut_load_all.dG_dev, S, S);
  fprintf(file_out, "Mut.load.DeltaG\t%.5g\t%.2g\n", x, dx);

  Get_mean(&x, &dx, trans_load_all.df_ave, trans_load_all.df_dev, S, S);
  fprintf(file_out, "Trans.load.fitness\t%.5g\t%.2g\n", x, dx);
  Get_mean(&x, &dx, trans_load_all.dG_ave, trans_load_all.dG_dev, S, S);
  fprintf(file_out, "Trans.load.DeltaG\t%.5g\t%.2g\n", x, dx);

  Get_mean(&x, &dx, dNdS1_all, dNdS2_all, Samples, Samples);
  fprintf(file_out, "dN/dS\t%.5g\t%.2g\n", x, dx);
  Get_mean(&x, &dx, accept1_all, accept2_all, Samples, Samples);
  fprintf(file_out, "acceptance\t%.5g\t%.2g\n", x, dx);
  Get_mean(&x, &dx, seqentr1_all, seqentr2_all, Samples, Samples);
  fprintf(file_out, "Seq.entropy\t%.5g\t%.2g\n", x, dx);
  fclose(file_out);
  Print_profile_evo(name_file, aa_distr_all, aa_seq0, len_amm, DG1_all, it_all);

  if(PRINT_TN)printf("Writing %s\n", name_div);

  // Free memory
  free(aa_seq);
  free(aa_seq0);
  for(i=0; i<len_amm; i++){
    free(aa_distr[i]);
    free(aa_distr0[i]);
  }
  free(aa_distr);
  free(aa_distr0);
  free(dna_seq);
  free(nuc_seq);
  free(P_MF_ia);
  free(nuc_evo);
  Empty_E_REM(&E_wt);
  Empty_E_REM(&E_mut);
  return(0);
}


void Print_ave(FILE *file_ave,
	       long it_sum, long n_subst,
	       int N_pop,
	       double f_sum, double f_dev,
	       double E_sum, double E_dev,
	       double DG_sum, double DG_dev,
	       long num_syn_subst, long num_aa_subst,
	       long num_syn_mut, long num_aa_mut,
	       float seq_entr,
	       struct load mut_load, struct load trans_load)
{
//   float t_indip=n_subst/20;
   float t_indep=1.0;

  if(ini_print==0){
    fprintf(file_ave, "# Fitness (sd)     Enat  (sd)    DG  (sd)");
    fprintf(file_ave, " seq.entropy(mut)-seq_entropy(sel)     ");
    //fprintf(file_ave, " Trans.load   (sd)  Mut.load  (sd)");
    fprintf(file_ave, "  nonsyn/syn accept  N_subst\n");
    ini_print=1;
  }

  double x, dx;
  Get_mean(&x, &dx, f_sum, f_dev, it_sum, t_indep);
  fprintf(file_ave, " %.4g\t%.2g\t", x, dx);
  Get_mean(&x, &dx, E_sum, E_dev, it_sum, t_indep);
  fprintf(file_ave, " %.4g\t%.2g\t", x, dx);
  Get_mean(&x, &dx, DG_sum, DG_dev, it_sum, t_indep);
  fprintf(file_ave, " %.4g\t%.2g\t", x, dx);
  fprintf(file_ave, " %.4f\t", -seq_entr);

  fprintf(file_ave, " %.3f\t",
	  (float)(num_aa_subst*num_syn_mut)/(float)(num_syn_subst*num_aa_mut));
  fprintf(file_ave, " %.4g\t", (float)(num_syn_subst+num_aa_subst)/it_sum);
  fprintf(file_ave, " %d\n", (int)num_aa_subst);
  fflush(file_ave);
}



void Get_mean(double *x, double *dx, float sum, float dev,
	      long it_sum, float t_indep)
{
  *x=sum/it_sum;
  double x2=dev-it_sum*(*x)*(*x);
  if(dev<=0){*dx=0; return;}
  *dx=sqrt(x2/(it_sum*t_indep));
}
    
void Print_final(char *name_file, long it_sum,
		 float TEMP, float sU1, float sC1, float sC0,
		 int MEANFIELD, float Lambda, int N_pop,
		 float *mut_par, float tt_ratio,
		 //
		 float *fit_evo, float *DG_evo, 
		 float *Enat_evo, int nsam, int step,
		 double *f_ave, double f_dev, 
		 double *E_ave, double E_dev,
		 double *DG_ave, double DG_dev,
		 float seq_entr, double seq_entr_dev,
		 struct load *mut_load, struct load *trans_load,
		 long num_syn_subst, long num_aa_subst,
		 long num_syn_mut, long num_aa_mut,
		 float *dN_dS, float *accept,
		 long **nuc_evo, int len_dna, int sample)
{
  char name_out[200];
  sprintf(name_out, "%s_final.dat", name_file);
  printf("Writing %s\n", name_out);
  FILE *file_out=NULL;

  if(sample==0){
    file_out=fopen(name_out, "w");
    int L=len_dna/3;
    // Headers   
    fprintf(file_out, "# Input thermodynamic parameters:\n");
    fprintf(file_out, "L= %d residues\n", L);
    fprintf(file_out, "TEMP=\t%.2f\n", TEMP);
    fprintf(file_out, "Entropy unfolded=\t%.1f (%.2f *L)\n", sU1*L, sU1);
    fprintf(file_out, "Entropy compact=\t%.1f (%.2f + %.2f *L)\n",
	    sC0+sC1*L, sC0, sC1);
    if(SEC_STR)
      fprintf(file_out, "#Local interactions used, factor= %.3f\n", SEC_STR);
    fprintf(file_out, "# Input population parameters:\n");
    if(MEANFIELD){
      fprintf(file_out, "Meanfield model, Lambda=\t%.2f\n", LAMBDA);
    }else{
      fprintf(file_out, "Population model, Npop=\t%d\n", N_pop);
    }
    fprintf(file_out, "# Input mutation parameters:\n");
    for(int i=0; i<4; i++)
      fprintf(file_out,"%c %.3f ",NUC_CODE[i],mut_par[i]);
    fprintf(file_out, " TT_ratio= %.2f\n", tt_ratio);
  }else{
    file_out=fopen(name_out, "a");
  }
  fprintf(file_out, "#\n#\n# sample %d ", sample+1);
  fprintf(file_out, "Attempted mutations= %ld\n", it_sum);
  fprintf(file_out, "# Results:\n");

  double x, dx, tau, r;  float t_indep=num_aa_subst/10;
  Get_mean(&x, &dx, *f_ave, f_dev, it_sum, t_indep);
  fprintf(file_out, "Fitness ave=\t%.4g %.4g\n", x, dx);
  *f_ave=x;

  x=Logistic_fit(&tau, &r, fit_evo, step, nsam);
  fprintf(file_out, "Fitness, logistic fit= %.4g tau= %.4g r=%.2g\n",x,tau,r);

  Get_mean(&x, &dx, *E_ave, E_dev, it_sum, t_indep);
  fprintf(file_out, "E_nat ave=\t%.4g %.2g\n", x, dx);
  *E_ave=x;

  x=Logistic_fit(&tau, &r, Enat_evo, step, nsam);
  fprintf(file_out, "E_nat, logistic fit= %.4g tau= %.4g r=%.2g\n",x,tau,r);

  Get_mean(&x, &dx, *DG_ave, DG_dev, it_sum, t_indep);
  fprintf(file_out, "DeltaG ave=\t%.4g %.2g\n", x, dx);
  *DG_ave=x;

  x=Logistic_fit(&tau, &r, DG_evo, step, nsam);
  fprintf(file_out, "DeltaG, logistic fit= %.4g tau= %.4g r=%.2g\n",x,tau,r);

  // Synonymous, mutation load
  *dN_dS=(float)(num_aa_subst*num_syn_mut)/(float)(num_syn_subst*num_aa_mut);
  fprintf(file_out, "dN/dS=\t%.3f\n", *dN_dS);

  *accept=(float)(num_syn_subst+num_aa_subst)/it_sum;
  fprintf(file_out, "Acceptance rate=\t%.4g\n", *accept);

  fprintf(file_out, "Seq_entropy(sel)-Seq_entropy(mut)=\t%7.4g %7.4g\n",
	  seq_entr, seq_entr_dev);

  int num=mut_load->num, indep=num/30;
  Get_mean(&x, &dx, mut_load->df_ave, mut_load->df_dev, num, indep);
  fprintf(file_out, "Mut.load.fitness=\t%.5g\t%.2g\n", x, dx);
  mut_load->df_ave=x;
  Get_mean(&x, &dx, mut_load->dG_ave, mut_load->dG_dev, num, indep);
  fprintf(file_out, "Mut.load.DeltaG=\t%.5g\t%.2g\n", x, dx);
  mut_load->dG_ave=x;

  Get_mean(&x, &dx, trans_load->df_ave, trans_load->df_dev, num, indep);
  fprintf(file_out, "Trans.load.fitness=\t%.5g\t%.2g\n", x, dx);
  trans_load->df_ave=x;
  Get_mean(&x, &dx, trans_load->dG_ave, trans_load->dG_dev, num, indep);
  fprintf(file_out, "Trans.load.DeltaG=\t%.5g\t%.2g\n", x, dx);
  trans_load->dG_ave=x;

  // Nucleotide content
  fprintf(file_out, "#\n# Final nucleotide frequencies (1, 2, 3)\n");
  int i, j, n; long nuc_count[4][3];
  for(n=0; n<4; n++){
    for(j=0; j<3; j++)nuc_count[n][j]=0;
  }
  j=0;
  for(i=0; i<len_dna; i++){
    for(n=0; n<4; n++)nuc_count[n][j]+=nuc_evo[i][n];
    j++; if(j==3)j=0;
  }
  double sum=0;
  for(n=0; n<4; n++)sum+=nuc_count[n][0];
  for(n=0; n<4; n++){
    fprintf(file_out, "%c ", NUC_CODE[n]);
    for(j=0; j<3; j++){
      fprintf(file_out," %.3f", nuc_count[n][j]/sum);
    }
    fprintf(file_out, "\n");
  }
  fclose(file_out);
}


int Read_ene_new(char *file_ene, float **interactions)
{
  FILE *file_in=fopen(file_ene,"r");
  int i, j, n=0; float ene;
  char string[200], aa1[8], aa2[8];

  if(file_in==NULL){
    printf("ERROR, energy parameter file %s not found\n", file_ene);
    exit(8);
  }
  while(fgets(string, sizeof(string), file_in)!=NULL){
    sscanf(string, "%s%s%f", aa1, aa2, &ene); n++;
    i=Code_AA(aa1[0]); j=Code_AA(aa2[0]);
    interactions[i][j]=ene; interactions[j][i]=ene;
  }
  return(n);
}


void Read_ene_par(char *file_ene, float **interactions)
{
  FILE *file_in=fopen(file_ene,"r");
  int i; char string[200];
  if(file_in==NULL){
    printf("ERROR, energy parameter file %s not found\n", file_ene);
    exit(8);
  }
  
  for(i=0; i<20; i ++){
    float *MAT=interactions[i];
    fgets(string, sizeof(string), file_in);
    sscanf(string, 
	   "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f", 
	   &MAT[0],&MAT[1],&MAT[2],&MAT[3],&MAT[4],
	   &MAT[5],&MAT[6],&MAT[7],&MAT[8],&MAT[9],&MAT[10],
	   &MAT[11],&MAT[12],&MAT[13],&MAT[14],&MAT[15],
	   &MAT[16],&MAT[17],&MAT[18],&MAT[19]);
  }
  /*fgets(string, sizeof(string), matrix); // Disulfide bonds
    sscanf(string, "%f", &interactions[210]);
    interactions[210]-=interactions[label[17][17]];*/
  (void)fclose(file_in);
  return;
}

char *Read_sequence(int *len_dna, char *inputseq)
{  
  char *sequence, string[1000], *ptr;
  FILE *file_in=fopen(inputseq, "r");
  int i=0;

  if(file_in==NULL){
    printf("WARNING, sequence file %s does not exist\n", inputseq);
    return(NULL);
  }

  printf("Reading %s\n",inputseq);
  fgets(string, sizeof(string), file_in);
  sscanf(string,"%s", seq_name); seq_name[0]=' ';
  printf("DNA sequence %s", seq_name);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    ptr=string;
    while((ptr!=NULL)&&(*ptr!='\n')){
      if((*ptr=='a')||(*ptr=='A')||(*ptr=='t')||(*ptr=='T')||
	 (*ptr=='g')||(*ptr=='G')||(*ptr=='c')||(*ptr=='C')){
	(*len_dna)++;
      }else if((*ptr!=' ')&&(*ptr!='\0')){
	printf("Wrong character %d in DNA sequence %s: %c\n",
	       *len_dna, file, *ptr);
	exit(8);
      }
      ptr++;
    }
  }
  fclose(file_in);

  if(*len_dna==0)return(NULL);

  // Reading
  sequence=(char *)malloc(*len_dna *sizeof(char));
  file_in=fopen(inputseq, "r");
  fgets(string, sizeof(string), file_in);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    ptr=string;
    while((ptr!=NULL)&&(*ptr!='\n')){
      if((*ptr!=' ')&&(*ptr!='\0')){
	*ptr=Maiuscule(*ptr); sequence[i]=*ptr;
	//int i_nuc=Code_nuc(*ptr);
	i++;
      }
      ptr++;
    }
  }
  fclose(file_in);
  
  
  return(sequence);
}

FILE *open_file(char *name_file, char *ext, short *seq, int lf, char *fit_def)
{
  FILE *file_out; char name[N_CHAR];
  int i, j; float gc=0;
  for(i=0; i<4; i++)if((NUC_CODE[i]=='G')||(NUC_CODE[i]=='C'))gc+=mut_par[i];
  sprintf(name, "%s%s", name_file, ext);
  file_out=fopen(name, "w");
  printf("Writing %s\n", name);
  fprintf(file_out, "# File %s, sequence %s,  PDB %s,  length=%d\n",
	  name_file, seq_name, target.name, len_amm);
  fprintf(file_out, "# fitness: %s\n", fit_def);
  fprintf(file_out, "# %ld iterations, random seed: %ld\n", IT_MAX, iran);
  fprintf(file_out, "# Stationary frequencies: ");
  for(i=0; i<4; i++)fprintf(file_out, "%c %.3f ", NUC_CODE[i], mut_par[i]);
  fprintf(file_out, "\n# Transition-transversion ratio= %.2f\n", tt_ratio);
  fprintf(file_out, "# GC bias: %.3f\n", gc);
  fprintf(file_out, "# T= %.3f\n", TEMP);
  fprintf(file_out, "# Config. entropy per residue, unfolded: %.3f", sU1);
  fprintf(file_out, " misfolded: %.3f\n", sC1);
  if(SEC_STR)
    fprintf(file_out, "#Local interactions used, factor= %.3f\n", SEC_STR);
  fprintf(file_out, "# N_pop= %d\n#", N_pop);

  for(i=0; i<len_amm; i++){
    if(i==(i/60)*60)fprintf(file_out,"\n# ");
    fprintf(file_out,"%c", AMIN_CODE[seq[i]]);
  }
  fprintf(file_out,"\n");
  
  if(lf==1){
    fprintf(file_out,
	    "# site aa_old-aa_new  Enat DG fitness syn_subst attempts\n");
    //D(0,n)
  }else if(lf==2){
    fprintf(file_out, "# iter [GC]_3");
    for(i=0; i<3; i++){
      fprintf(file_out, " ");
      for(j=0; j<3; j++)fprintf(file_out, "%c_%d ", NUC_CODE[i], j+1);
    }
    fprintf(file_out, "\n\n");
    /*fprintf(file_out, "  Energy alpha\n");*/
  }
  
  fflush(file_out);
  return(file_out);
}

int Print_dna(char *seq, FILE *file_out, int iter)
{
  short i, j, i_nuc, nuc_count[4][3];
  float length=(float)len_amm, gc3=0;

  for(i=0; i<4; i++)
    for(j=0; j<3; j++)
      nuc_count[i][j]=0;

  j=0;
  for(i=0; i<len_dna; i++){
    i_nuc=Code_nuc(seq[i]);
    nuc_count[i_nuc][j]++; j++;
    if(j==3)j=0;
  }

  /* GC3 */
  for(i=0; i<4; i++)
    if((NUC_CODE[i]=='G')||(NUC_CODE[i]=='C'))gc3+=nuc_count[i][2];
  fprintf(file_out,"%3d  %.3f ", iter, (gc3)/length);
  for(i=0; i<3; i++){
    fprintf(file_out," ");
    for(j=0; j<3; j++){
      fprintf(file_out," %.3f", nuc_count[i][j]/length);
    }
  }
  fprintf(file_out, "\n");
  fflush(file_out);

  return(0);
}

int Get_name(char *name, char *name_seq, int N){
  int i;
  for(i=0; i<N; i++){
    if(strncmp(name_seq+i, ".dna", 4)==0){
       name[i]='\0'; 
       return(0);
    }
    name[i]=name_seq[i];
  }
  printf("ERROR in file name %s, extension .dna not found\n", name);
  exit(8);
}

void Output_name(char *file_name, char *dir_out, char *prot_name,
		 float TEMP, float sU1, float sC1, int MEANFIELD,
		 char *MODEL, float LAMBDA, int OPT_LAMBDA, //NEW
		 int NEUTRAL, int N_pop, float *mut_par)
{
  char name[400];
  if(dir_out[0]!='\0'){sprintf(name, "%s/%s", dir_out, prot_name);}
  else{sprintf(name, "%s", prot_name);}
  if(MEANFIELD){
    //NEW
    if(OPT_LAMBDA){
      sprintf(file_name, "%s_MF_LAMBDA_OPT_%s", name, MODEL);
    }else{
      sprintf(file_name, "%s_MF_LAMBDA%.2f_%s", name, LAMBDA, MODEL);
    }
    sprintf(file_name, "%s_REM%d_T%.2f_SU%.3f", file_name, REM, TEMP, sU1);
  }else if (NEUTRAL){
    sprintf(file_name, "%s_T%.2f_SU1%.2f_SC1%.2f_NEUTRAL",
	    name, TEMP, sU1, sC1);
  }else{
    sprintf(file_name, "%s_T%.2f_SU1%.2f_SC1%.2f_N%d",
	    name, TEMP, sU1, sC1, N_pop);
  }
  sprintf(file_name, "%s_GC%.2f", file_name,
	  mut_par[Code_nuc('G')]+mut_par[Code_nuc('C')]);
  if(CpG)sprintf(file_name, "%s_CpG%.0f", file_name, mut_par[4]);
}
 
float Sequence_entropy(double **aa_distr, int L){
  int i, j; float S_sum=0, S, p; double norm=0;

  for(j=0; j<20; j++)norm+=aa_distr[0][j];
  for(i=0; i<L; i++){
    S=0;
    for(j=0; j<20; j++){
      p=aa_distr[i][j]/norm; if(p)S-=p*log(p);
    }
    S_sum+=S;
  }
  return(S_sum/L);
}

void Compute_freq_codons(float *mut_par, float *freq_aa,
			 char **codon, char *coded_aa)
{
  int i, i_aa; float w;
  for(i_aa=0; i_aa<20; i_aa++)freq_aa[i_aa]=0;
  for(i=0; i<64; i++){
    w=mut_par[Code_nuc(codon[i][0])];
    w*=mut_par[Code_nuc(codon[i][1])];
    w*=mut_par[Code_nuc(codon[i][2])];
    if(coded_aa[i]=='*')continue;
    i_aa=Code_AA(coded_aa[i]);
    if(i_aa<0)continue;
    freq_aa[i_aa]+=w;
  }
}

float Sequence_entropy_mut(float *mut_par, char **codon, char *coded_aa)
{
  float freq_aa[20]; int i; float norm=0, S=0, p;

  // Calculating amino acid distribution under mutation alone
  Compute_freq_codons(mut_par, freq_aa, codon, coded_aa);

  // Compute entropy
  for(i=0; i<20; i++)norm+=freq_aa[i];
  for(i=0; i<20; i++){
    if(freq_aa[i]){p=freq_aa[i]/norm; S-=p*log(p);}
  }
  return(S);
}

void Compute_load(double *Tload_sum, double *Tload_dev,
		  double *Mload_sum, double *Mload_dev, int *Nload,
		  struct REM *E_wt, int **C_nat, int *i_sec, short *aa_seq,
		  float fitness_wt, char *dna_seq, int len_dna, 
		  char **codon, char *coded_aa)
{
  double translation_load=0, mutation_load=0;

  // Folding thermodynamics
  float fitness, DeltaG;
  struct REM E_mut;
  Initialize_E_REM(&E_mut, E_wt->L,E_wt->REM,
		   E_wt->T,E_wt->S_C,E_wt->S_U, FILE_STR, 0);

  /*
    Load=Sum_j P(nat->Seq_j)[finess(nat)-fitness(Seq_j)]
    P(nat->Seq_j) is non zero only if sequence j is one base mutation from
    the native sequence.
    Translation load: P(nat->Seq_j)=1
    Mutation load: P(nat->Seq_j)=Mutation probability,
    i.e. f(new base) times 1 if transversion, times tt_ratio if transition    
   */

  // Mutations
  int pos=-1, res_mut=0, aa_old=-1;
  char *codon_nat=dna_seq, codon_mut[3];
  int sum_trans=0, sum_mut=0;

  for(int nuc_mut=0; nuc_mut<len_dna; nuc_mut++){
    
    pos++;
    if(pos==3){
      pos=0; res_mut++; codon_nat=dna_seq+res_mut*3;
    }
    for(int j=0; j<3; j++)codon_mut[j]=codon_nat[j];
    int nuc_wt=Code_nuc(dna_seq[nuc_mut]);

    for(int base=0; base<4; base++){

      // Mutate amino acid
      if(base==nuc_wt)continue;
      codon_mut[pos]=NUC_CODE[base];
      int aa_new=Coded_aa(codon_mut, codon, coded_aa);
      if(aa_new==aa_old)continue;         // Synonymous
      if(aa_new<0){fitness=0; goto load;} // Stop codon

      // Folding thermodynamics
      Copy_E_REM(&E_mut, E_wt);
      DeltaG=
	Mutate_DG_overT_contfreq(&E_mut, aa_seq, C_nat,
				 i_sec, res_mut, aa_new);
      fitness=1./(1+exp(DeltaG));
      
    load:
      // Translation load
      sum_trans++;
      translation_load+=fitness;

      // Mutation load
      float p=mut_par[base];
      if(base==Transition(NUC_CODE[nuc_wt]))p*=tt_ratio;
      sum_mut+=p;
      mutation_load+=fitness*p;
    }
  }
  Empty_E_REM(&E_mut);

  translation_load = fitness_wt-(translation_load)/sum_trans;
  mutation_load = fitness_wt-(mutation_load)/sum_mut;

  *Tload_sum += translation_load;
  *Tload_dev += translation_load*translation_load;
  *Mload_sum += mutation_load;
  *Mload_dev += mutation_load*mutation_load;
  (*Nload)++;
}

unsigned long randomgenerator(void){
     
     unsigned long tm;
     time_t seconds;
     
     time(&seconds);
     srand((unsigned)(seconds % 65536));
     do   /* waiting time equal 1 second */
       tm= clock();
     while (tm/CLOCKS_PER_SEC < (unsigned)(1));
     
     return((unsigned long) (rand()));

}

int Selection(float fitness, float fitness_old, int N_pop)
{
  // Moran's process:
  /* P_fix = (1-exp(-Df))/(1-exp(-N*Df)) */
  if(fitness<=0)return(0);
  //if((int)N_pop==1)return(1);
  double f_ratio= fitness_old / fitness;
  double x= (1.-f_ratio)/(1.-pow(f_ratio, N_pop));
  if(VBT){
    float rand=RandomFloating();
    printf("es= %.2f x= %.2g rand= %.2g\n", f_ratio, x, rand);
    if(rand < x)return(1);
  }else{
    if(RandomFloating() < x)return(1);
  }
  return(0);
}

int Detailed_balance(float *p, int xend, int xini)
{
  if(p[xend]>=p[xini])return(1);
  if(RandomFloating() < p[xend]/p[xini])return(1);
  return(0);
}

float *Get_counts(short *seq, int L, int Naa)
{
  float *num=malloc(Naa*sizeof(float));
  int i; for(i=0; i<Naa; i++)num[i]=0;
  for(i=0; i<L; i++){
    if((seq[i]<0)||(seq[i]>=Naa)){
      printf("ERROR, wrong symbol %d at site %d\n", seq[i], i);
      exit(8);
    }
    num[seq[i]]++;
  }
  //for(i=0; i<Naa; i++)num[i]/=L;
  return(num);
 }

void Record_evo(float *fit_evo, double fitness,
		float *DG_evo,  double DG,
		float *Enat_evo, double E_nat, int *nsam)
{
  fit_evo[*nsam]=fitness;
  DG_evo[*nsam]=DG;
  Enat_evo[*nsam]=E_nat;
  (*nsam)++;
}


void WildType_Mutants(float *Lambda_opt, double *DG_ave, float **P_WT_ia,
		      float *P_mut_a, int **C_nat, int *i_sec, short *aa_seq,
		      struct REM *E_wt, char *name_file, FILE *file_summary)
{
  printf("Wild type mutants. DG_wt=%.2f REM=%d\n",
	 E_wt->DeltaG, E_wt->REM);
  *DG_ave=0;
  int L=E_wt->L;
  float **DG_mut=Allocate_mat2_f(L, 20), DG;
  struct REM E_mut;
  Initialize_E_REM(&E_mut, L, E_wt->REM,
		   E_wt->T, E_wt->S_C, E_wt->S_U, FILE_STR, 0);

  for(int res_mut=0; res_mut<L; res_mut++){
    float *DDG=DG_mut[res_mut];
    //printf("%d  ", res_mut);
    for(int aa=0; aa<20; aa++){
      if(aa==aa_seq[res_mut]){
	DDG[aa]=0;
      }else{
	Copy_E_REM(&E_mut, E_wt);
	DG=Mutate_DG_overT_contfreq(&E_mut, aa_seq, C_nat, i_sec, res_mut, aa);
	DDG[aa]=DG-E_wt->DeltaG;
      }
    }
  }
  printf("%d mutants computed\n", L*19);

  int N_Lambda=12, k; float Lambda_step=0.1, Lambda;
  float Lambda_k[N_Lambda], lik_k[N_Lambda], lik_opt=-1000;
  float **P_opt_ia=Allocate_mat2_f(L, 20);
  printf("#Lambda log_likelihood\n");
  for(k=0; k<N_Lambda; k++){
    Lambda=k*Lambda_step;
    Compute_P_WT(P_WT_ia, Lambda, DG_mut, P_mut_a, L);
    lik_k[k]=Compute_lik(P_WT_ia, L, aa_seq);
    Lambda_k[k]=Lambda;
    printf("%.3f %.3f\n", Lambda, lik_k[k]);
    if((k==0)||(lik_k[k] > lik_opt)){
      lik_opt=lik_k[k]; *Lambda_opt=Lambda;
      Copy_P(P_opt_ia, P_WT_ia, L, 20);
    }
  }
  Lambda=Find_max(lik_k, Lambda_k, N_Lambda, 0, Lambda+0.1);
  if(Lambda==0){
    Lambda=0.05;
    printf("WARNING, optimal Lambda=0, changing to %.2f\n", Lambda);
  }
  Compute_P_WT(P_WT_ia, Lambda, DG_mut, P_mut_a, L);
  float lik=Compute_lik(P_WT_ia, L, aa_seq);
  if((lik>lik_opt)||(*Lambda_opt==0)){
    lik_opt=lik; *Lambda_opt=Lambda;
  }else{
    Copy_P(P_WT_ia, P_opt_ia, L, 20);
  }
  printf("Optimal likelihood= %.3f Lambda=%.3f\n", lik_opt, *Lambda_opt);

  //printf("\n"); for(i=0; i<20; i++)printf("%.2f ", P_WT[i]); printf("\n");
  // Compute and print properties of the distribution
  struct MF_results wt_opt;
  wt_opt.lik=lik_opt;
  wt_opt.Lambda=*Lambda_opt;
  wt_opt.Tf=E_wt->Tf;
  wt_opt.h=Hydro_ave(P_WT_ia, hydro, L);
  double KLD=0, DDG=0;
  for(int i=0; i<L; i++)KLD+=KL(P_WT_ia[i], P_mut_a, 20);
  wt_opt.KL=KLD/L;
  for(int i=0; i<L; i++){
    float *P_WT=P_WT_ia[i], *ddG=DG_mut[i];
    for(int aa=0; aa<20; aa++)DDG+=P_WT[aa]*ddG[aa];
  }
  wt_opt.DG=(DDG+E_wt->DeltaG); //L
  (*DG_ave)=wt_opt.DG;

  Print_results(wt_opt, "WT ", file_summary);

  Empty_E_REM(&E_mut);
  Empty_matrix_f(DG_mut, L);
  Empty_matrix_f(P_opt_ia, L);
}

void Compute_P_WT(float **P_WT_ia, float Lambda,
		  float **DG_mut, float *P_mut_a, int L)
{
  int i, aa;
  for(i=0; i<L; i++){
    float *P_WT=P_WT_ia[i], *DDG=DG_mut[i]; double Z=0;
    for(aa=0; aa<20; aa++){
      P_WT[aa]=P_mut_a[aa]*exp(-Lambda*(DDG[aa]));
      Z+=P_WT[aa];
    }
    for(aa=0; aa<20; aa++)P_WT[aa]/=Z;
  }
}

void Print_matrix(struct protein target){
  char name_out[200];
  sprintf(name_out, "%s_Contact_matrix.cm", target.name);
  FILE *file_out=fopen(name_out, "w"); int i;
  printf("Writing contact matrix in file %s\n", name_out);
  fprintf(file_out, "# %d %d %s\n",
	  target.length, target.n_cont, target.name);
  for(i=0; i<target.length; i++){
    short *Ci=target.contact[i];
    while(*Ci>=0){
      fprintf(file_out, "%d %d\n", i, *Ci); Ci++;
    }
  }
  fclose(file_out);
}

void Initialize_load(struct load *load){
  load->df_ave=0; load->df_dev=0;
  load->dG_ave=0; load->dG_dev=0;
  load->num=0;
}

void Sum_loads(struct load *load, struct load *load1){
  load->df_ave+=load1->df_ave;
  load->df_dev+=(load1->df_ave)*(load1->df_ave);
  load->dG_ave+=load1->dG_ave;
  load->dG_dev+=(load1->dG_ave)*(load1->dG_ave);
  load->num++;
}

FILE *Open_summary(char *name_file, struct REM E_wt,
		   float *mut_par, short *aa_seq)
{
  // Print summary of results
  char name[100]; int L=E_wt.L;
  sprintf(name, "%s_summary.dat", name_file);
  printf("Writing %s\n", name);
  FILE *file_out=fopen(name, "w");
  fprintf(file_out, "# Loc.interactions: %.3f ", SEC_STR);
  fprintf(file_out, "T= %.2f SC=%.3f SU= %.3f REM= %d T_freezing= %.2f\n",
	  E_wt.T, E_wt.S_C/L, E_wt.S_U/L, E_wt.REM, E_wt.Tf);
  name[5]='\0';
  // Wild type sequence
  double h_wt=0; for(int i=0; i<L; i++)h_wt+=hydro[aa_seq[i]]; h_wt/=L;
  fprintf(file_out, "# PDB= \"%s\" DG_wt= %.1f h= %.3f\n",
	  name, E_wt.DeltaG, h_wt);
  fprintf(file_out, "#L=%3d   nuc_freq= %.3f %.3f %.3f %.3f", L,
	  mut_par[Code_nuc('G')], mut_par[Code_nuc('C')],
	  mut_par[Code_nuc('A')], mut_par[Code_nuc('T')]);
  fprintf(file_out, "   tt= %.1f CpG= %.1f twonuc= %.1g\n",
	  mut_par[4], mut_par[5], mut_par[6]);
  fprintf(file_out, "# Model Lambda likelihood KL DG Tf h\n");
  return(file_out);
}

void Print_TN_div(short *aa_seq, short *aa_seq0, int L,
		  int num_aa_subst, FILE *file_out)
{
  float SI0=0.06;
  float SI=0;
  for(int i=0; i<L; i++)if(aa_seq[i]==aa_seq0[i])SI++;
  SI/=L; if(SI<=SI0)return;
  float TN=-log((SI-SI0)/(1-SI0));
  fprintf(file_out, "%.3f\t%d\n", TN, num_aa_subst);
}
