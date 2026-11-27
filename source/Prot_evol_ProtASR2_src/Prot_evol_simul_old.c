/* Program dna_evol.
   Needs a target structure and a starting DNA sequence.
   Attempts DNA mutations with given mutational bias.
   alpha, Z, energy are computed for the mutated sequence.
   The mutation is accepted or not according to selection criteria.
   The base composition is evaluated at every time step.
*/

#define INPUT_DIR "/home/ubastolla/RESEARCH/POP_DYN/INPUT/"
// Change when installing the program!!
#define FILE_CODE_DEF "gen_code_ATGC.in"
#define FILE_ENE_DEF  "energy.in"

#include "thread_contfreq.h"
#include "gen_code.h"
#include "allocate.h"
#include "protein3.h"
#include "read_pdb.h"
#include "random3.h"           /* Generating random numbers */
#include "mutation.h"
#include "codes.h"
#include "get_para_pop_dyn.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>

int ini_print;

/* New version 11/07/2012:
   Now energy parameters and genetic codes are given in header files
   energy_BKV.h gen_code.h
   It is possible to give as input either a PDB file (extension .pdb is
   needed) or a file with a nucleotide sequence and the name of the target
   pdb, which must be contained in the defeult input file structures.in
*/

//#define FILE_IN "Pop_evol.in"
#define N_CHAR 300          // Max. length of file names
#define EXT_AVE "_ave.dat"  // Extension for output file with averages
#define EXT_DNA "_dna.dat"  // Extension for output file with dna statistics
#define EXT_OUT "_stab.dat" // Extension for output file with folding
                            // thermodynamics and fitness

int REM3=0;   // Second or third cumulant of the energy?

// Input codes
void Read_ene_par(char *, float **interactions);
int Read_ene_new(char *, float **interactions);
char *Read_sequence(int *len_dna, char *inseq);
unsigned long randomgenerator(void);

// Output
int Output_name(char *file_name, char *dir, char *name,
		float TEMP, float S0, int N_pop, float *freq_nuc);
FILE *open_file(char *, char *, short *, int, char *fit_def);
int Print_dna(char *, FILE *, int);
void Print_ave(FILE *file_ave,
	       double it_sum, float t_indip,
	       int N_pop,
	       double f_sum, double f_dev,
	       double E_sum, double E_dev,
	       double DG_sum, double DG_dev,
	       float seq_entr,
	       double num_syn_subst, double num_aa_subst,
	       double Tload_sum, double Tload_dev,
	       double Mload_sum, double Mload_dev,
	       int Nload);
void Print_mean(FILE *file_out, float sum, float dev,
		float n_sum, float t_indip);
void Print_final(char *name_file, double it_sum,
		 float TEMP, float s0,
		 int N_pop, float *freq_nuc,
		 double f_ave, double f_dev, 
		 double E_ave, double E_dev,
		 double DG_ave, double DG_dev,
		 float seq_entr, double seq_entr_dev,
		 double Tload_sum, double Tload_dev,
		 double Mload_sum, double Mload_dev,
		 int Nload,
		 double num_syn_subst, double num_aa_subst,
		 char *dna_seq, int len_dna);

// Calculations
int Selection(float fitness, float fitness_old, int N_pop);
float Sequence_entropy(int **aa_distr, int L);
float Sequence_entropy_mut(float *freq_nuc, char **codon, char *coded_aa);
void Compute_freq_codons(float *freq_nuc, float *freq_aa,
			 char **codon, char *coded_aa);
void Compute_load(double *translation_load, double *mutation_load,
		  short *aa_seq, int len_amm, char *dna_seq, int len_dna,
		  float fitness_wt, double E_nat, double E1,
		  double E2, double E23, double E3, float Conf_entropy, 
		  char **codon, char *coded_aa);

// Input parameters
static long IT_MAX=1000;
static int N_pop=100;
static float freq_nuc[4], tt_ratio=2;
static char file[N_CHAR], seq_name[N_CHAR], dir_out[N_CHAR];
static char seq_name[N_CHAR];
static int len_amm, len_dna;
static int count[4];
static float rate[4];
unsigned long iran;


main(int argc, char **argv){

  /***********************
          INPUT
  ************************/
  // Input files
  char file_pdb[N_CHAR], file_str[N_CHAR], file_seq[N_CHAR];
  char chain[]="\0\0\0";
  char fit_def[1000], name_file[N_CHAR];
  char FILE_CODE[N_CHAR], FILE_ENE[N_CHAR], FILE_STR[N_CHAR];

  // Genetic code
  // char *codon[64], coded_aa[64], name_code[200];

  /***********************
         SEQUENCES
  ************************/
  short *aa_seq=NULL, *aa_seq0=NULL;
  char *dna_seq=NULL, nuc_new;
  int nuc_mut, res_mut=0, aa_new;

  /***********************
         Wild Type
  ************************/
  float fitness_wt;

  /***********************
         MUTANT
  ************************/
  float fitness_mut;

  /***********************
          SIMULATION
  ************************/
  int i, j, aa_subst, synonimous;
  long it_mut=0, n_change=0, syn_subst=0;
  int it_print, it1=0;
  double it_trans, it_subst=0;

  /***********************
          OUTPUT
  ************************/
  FILE *file_stab, *file_dna, *file_ave, *file_out;
  // Average stability and fitness
  double it_sum=0;
  double E_nat_ave=0, f_ave=0, DG_ave=0;
  double E_nat_dev=0, f_dev=0, DG_dev=0;

  // Count substitutions
  double num_syn_subst=0, num_aa_subst=0;
  // Sequence entropy
  double seq_entr_sum=0, seq_entr_dev=0, entr_dev;
  float seq_entr, seq_entr_mut;
  int n_seq_entr=0;
  int **aa_distr, **aa_distr0;
  // Loads
  double Tload_sum=0, Tload_dev=0, Mload_sum=0, Mload_dev=0;
  int Nload=0; 

  // Input operations
  sprintf(FILE_CODE, FILE_CODE_DEF);
  sC1=0.1; sC0=0; s0=0.0; TEMP=1; // default
  Get_para(argc, argv, file_pdb, chain, file_seq, FILE_STR,
	   &N_pop, &TEMP, &s0, &sC0, &sC1, &IT_MAX, freq_nuc,
	   &tt_ratio, &REM3, dir_out);
  sprintf(file_str, "%s%s", INPUT_DIR, FILE_STR);

  // Random numbers
  iran=randomgenerator();
  InitRandom( (RANDOMTYPE)iran);

  // Compute contact matrix from PDB file 
  len_amm=Get_pdb(&target, file_pdb, chain);
  if(len_amm > 0){
    aa_seq=target.aa_seq;
    Fill_C_nat(len_amm, target.contact);
  }else{
    // or read contact matrix from precomputed matrices
    Get_target(file_str, file_pdb, &len_amm);
    if(len_amm>0){
      printf("contact matrix %s found in %s, %d residues\n",
	     file_pdb, file_str, len_amm);
    }else{
      printf("contact matrix %s not found in %s\n",file_pdb, file_str);
      exit(8);
    }
  }

  if(file_seq[0]!='\0'){
    // Read sequence from file
    char inputseq[N_CHAR];
    sprintf(inputseq, "%s%s", INPUT_DIR, file_seq);  
    if(Check_file(inputseq)==0){
      printf("WARNING, directory %s does not contain file %s\n",
	     INPUT_DIR, file_seq);
    }else{
      dna_seq=Read_sequence(&len_dna, inputseq);
      if(dna_seq==NULL)goto extract_dna;
      if((len_amm>0)&&(aa_seq!=NULL)&&
	 (Compare_amm_dna(dna_seq,len_dna,aa_seq,len_amm,codon,coded_aa)==0))
	len_dna=0;
    }
  }

 extract_dna:
  if(aa_seq==NULL){
    aa_seq=malloc(len_amm*sizeof(short));
    Translate_new(dna_seq, aa_seq, len_amm, codon, coded_aa);
  }
  if(len_dna<=0){
    printf("Estracting DNA sequence from AA sequence\n");
    dna_seq=Extract_dna(&len_dna, len_amm, aa_seq, codon, coded_aa);
  }

  // Mutation
  Ini_count(dna_seq, len_dna, count);
  Compute_rates(rate, freq_nuc, tt_ratio);

  // Number of iterations
  // IT_MAX*=len_amm;
  it_trans=len_amm; if(it_trans > IT_MAX*0.1)it_trans= IT_MAX*0.1;
  it_print=1000; if(it_print > IT_MAX*0.1)it_print= IT_MAX*0.1;

  // Reference sequence
  aa_seq0=malloc(len_amm*sizeof(short));
  for(i=0; i<len_amm; i++)aa_seq0[i]=aa_seq[i];

  // Amino acid distribution
  aa_distr=(int **)malloc(len_amm*sizeof(int *));
  for(i=0; i<len_amm; i++){
     aa_distr[i]=(int *)malloc(20*sizeof(int));
    for(j=0; j<20; j++)aa_distr[i][j]=0;
  }
  aa_distr0=(int **)malloc(len_amm*sizeof(int *));
  for(i=0; i<len_amm; i++){
    aa_distr0[i]=(int *)malloc(20*sizeof(int));
    for(j=0; j<20; j++)aa_distr0[i][j]=0;
  }
  seq_entr_mut=Sequence_entropy_mut(freq_nuc, codon, coded_aa);

  // Initialize energy calculations. Obsolete
  /*interactions=Allocate_mat2_f(20, 20);
  sprintf(file_ene, "%s%s", INPUT_DIR, FILE_ENE);
  if(strncmp(file_ene, "energy.in", 9)==0){
    Read_ene_par(file_ene, interactions);
  }else{
    Read_ene_new(file_ene, interactions);
    printf("Reading energy parameters from file %s\n", file_ene);
    }*/

  /************ Change april 2012 *******************************/
  double E_nat_wt=0.0, E1_wt=0.0, E2_wt=0.0, E23_wt=0.0, E3_wt=0.0;
  double E_nat_mut, E1_mut, E2_mut, E23_mut, E3_mut;
  float DeltaG, DeltaG_wt, T0;
  if(REM3){
    DeltaG_wt=Compute_DG_REM3(aa_seq, len_amm, &E_nat_wt, &E1_wt,
			      &E2_wt, &E23_wt, &E3_wt, file_str);
    sprintf(name_file, "%s%s_DeltaG_3.dat", dir_out, target.name);
    T0=Print_DG_REM3(E_nat_wt, E1_wt, E2_wt, E23_wt, E3_wt, name_file);
  }else{
    DeltaG_wt=
      Compute_DG_REM2(aa_seq, len_amm, &E_nat_wt, &E1_wt, &E2_wt, file_str);
    sprintf(name_file, "%s%s_DeltaG_2.dat", dir_out, target.name);
    T0=Print_DG_REM2(E_nat_wt, E1_wt, E2_wt, name_file);
  }
  printf("Temperature with smallest |DG|: T=%.2f\n", T0);
  if(DeltaG_wt > 0){
    printf("WARNING, unstable target structure! DG=%.2f at T=%.2f\n",
	   DeltaG, TEMP);
  }
  // Fitness
  fitness_wt=1./(1+exp(DeltaG_wt)); // exit(8);
  strcpy(fit_def, "F=1./(1+exp(DeltaG))");

  /******************** End change *********************/

  // Output operations (wild type sequence needed)
  Output_name(name_file, dir_out, target.name, TEMP, s0, N_pop, freq_nuc);
  file_dna=open_file(name_file, EXT_DNA, aa_seq, 2, fit_def);
  file_ave=open_file(name_file, EXT_AVE, aa_seq, 0, fit_def);
  file_stab=open_file(name_file, EXT_OUT, aa_seq, 1, fit_def);

  // Printing
  fprintf(file_stab, "%3d  %c  %.3f  %.3f %7.4f",
	  res_mut, AMIN_CODE[aa_seq[0]], E_nat_wt, DeltaG_wt, fitness_wt);
  fflush(file_stab);
  
  // Simulating evolution
  while(1){
    
    // Mutation
    it_mut++;
    synonimous=Mutate_seq(dna_seq, len_dna, codon, coded_aa,
			  aa_seq, len_amm, freq_nuc, tt_ratio, count, rate,
			  &nuc_mut, &nuc_new, &res_mut, &aa_new);

    if(synonimous > 0){
      if(it_subst >= it_trans)num_syn_subst++;
      syn_subst++;
      goto Update_dna;

    }else if(synonimous==0){

      E_nat_mut=E_nat_wt; E1_mut=E1_wt; E2_mut=E2_wt;
      if(REM3){
	E23_mut=E23_wt; E3_mut=E3_wt;
	DeltaG=Mutate_DG_REM3(aa_seq, len_amm, res_mut, aa_new, &E_nat_mut,
			      &E1_mut, &E2_mut, &E23_mut, &E3_mut);
      }else{
	DeltaG=Mutate_DG_REM2(aa_seq, len_amm, res_mut, aa_new,
			      &E_nat_mut, &E1_mut, &E2_mut);
      }
      fitness_mut=1./(1+exp(DeltaG));
      aa_subst=Selection(fitness_mut, fitness_wt, N_pop);
      if(aa_subst<=0)continue;
    }else{ 
      continue;  // Stop codon
    }

    // A substitution has occurred

    // Print thermodynamic properties (also in transient)
    it_subst++;
    fprintf(file_stab, " %2d %3d\n", syn_subst, it_mut);
    fprintf(file_stab, "%3d  %c  %.3f  %.3f %7.4f",
	    res_mut, AMIN_CODE[aa_new], E_nat_mut, DeltaG, fitness_mut);
    //fflush(file_stab);
    
    // Statistics after the transient
    if(it_subst >= it_trans){
      // Average old amino acid sequence, weight=it_mut
      
      // Update amino acid distribution
      for(i=0; i<len_amm; i++)aa_distr0[i][aa_seq[i]]+=it_mut;
      
      // Averages
      E_nat_ave +=it_mut*E_nat_wt; E_nat_dev+=it_mut*E_nat_wt*E_nat_wt;
      DG_ave+=it_mut*DeltaG_wt; DG_dev+=it_mut*DeltaG_wt*DeltaG_wt;
      f_ave+=it_mut*fitness_wt; f_dev+=it_mut*fitness_wt*fitness_wt;
      it_sum+=it_mut;
     
      num_aa_subst++; it1++;
      // Print averages
      if(it1==it_print){
	it1=0;
	/* Loads
	double translation_load, mutation_load;
	Compute_load(&translation_load, &mutation_load,
		     aa_seq, len_amm, dna_seq, len_dna,
		     fitness_wt, E_nat_wt, E1_wt, E2_wt, E23_wt, E3_wt,
		     Conf_entropy, codon, coded_aa);
	Tload_sum+=translation_load;
	Tload_dev+=translation_load*translation_load;
	Mload_sum+=mutation_load;
	Mload_dev+=mutation_load*mutation_load;
	Nload++; */
      }
      seq_entr=Sequence_entropy(aa_distr0, len_amm)-seq_entr_mut;
      for(i=0; i<len_amm; i++){
	for(j=0; j<20; j++){
	  aa_distr[i][j]+=aa_distr0[i][j]; aa_distr0[i][j]=0;
	}
      }
      seq_entr_sum+=seq_entr; seq_entr_dev+=seq_entr*seq_entr; n_seq_entr++;
      Print_ave(file_ave, it_sum, num_aa_subst, N_pop,
		f_ave, f_dev, E_nat_ave, E_nat_dev, DG_ave, DG_dev,
		seq_entr, num_syn_subst, num_aa_subst,
		Tload_sum, Tload_dev, Mload_sum, Mload_dev, Nload);
    }
    
  Update_AA:
    aa_seq[res_mut]=aa_new;
    DeltaG_wt=DeltaG;
    fitness_wt=fitness_mut;
    E_nat_wt=E_nat_mut;
    E1_wt=E1_mut; E2_wt=E2_mut;
    if(REM3){E23_wt=E23_mut; E3_wt=E3_mut;}
    it_mut=0; syn_subst=0;

  Update_dna:
    if(aa_subst ||(synonimous > 0)){
      count[Code_nuc(dna_seq[nuc_mut])]--;
      count[Code_nuc(nuc_new)]++;
      dna_seq[nuc_mut]=nuc_new;
    }

    if(num_aa_subst == IT_MAX)break;
      
  }
  fprintf(file_stab, " %2d %d\n", syn_subst, it_mut);



  /*********************** Final results    *****************************/
  // Loads
  {
    double translation_load, mutation_load;
    Compute_load(&translation_load, &mutation_load,
		 aa_seq, len_amm, dna_seq, len_dna,
		 fitness_wt, E_nat_wt, E1_wt, E2_wt, E23_wt, E3_wt,
		 Conf_entropy, codon, coded_aa);
    Tload_sum+=translation_load;
    Tload_dev+=translation_load*translation_load;
    Mload_sum+=mutation_load;
    Mload_dev+=mutation_load*mutation_load;
    Nload++;
  }
  // Entropy
  for(i=0; i<len_amm; i++){
    for(j=0; j<20; j++){
      aa_distr[i][j]+=aa_distr0[i][j]; aa_distr0[i][j]=0;
    }
  }
  entr_dev=seq_entr_dev-seq_entr_sum*seq_entr_sum/n_seq_entr;
  entr_dev=sqrt(entr_dev)/n_seq_entr;
  seq_entr=Sequence_entropy(aa_distr, len_amm)-seq_entr_mut;

  it_sum+=it_mut;
  Print_final(name_file, it_sum, TEMP, s0, N_pop, freq_nuc,
	      f_ave, f_dev, E_nat_ave, E_nat_dev, DG_ave, DG_dev,
	      seq_entr, entr_dev,
	      Tload_sum, Tload_dev, Mload_sum, Mload_dev, Nload,
	      num_syn_subst, num_aa_subst, dna_seq, len_dna);
  //Print_dna(dna_seq, file_dna, iter);
  
  free(aa_seq);
  free(aa_seq0);
  for(i=0; i<len_amm; i++){
     free(aa_distr[i]);
     free(aa_distr0[i]);
  }
  free(aa_distr);
  free(aa_distr0);
  free(dna_seq);

  return(0);
}


void Print_ave(FILE *file_ave,
	       double it_sum, float n_subst,
	       int N_pop,
	       double f_sum, double f_dev,
	       double E_sum, double E_dev,
	       double DG_sum, double DG_dev,
	       float seq_entr,
	       double num_syn_subst, double num_aa_subst,
	       double Tload_sum, double Tload_dev,
	       double Mload_sum, double Mload_dev,
	       int Nload)
{
//   float t_indip=n_subst/20;
   float t_indip=1.0;

  if(ini_print==0){
    fprintf(file_ave, "# N  fitness (sd)     E    (sd)    gap  (sd)");
    fprintf(file_ave, " seq.entropy     Trans.load   (sd)  Mut.load  (sd)");
    fprintf(file_ave, "  syn/nonsyn reject  samples\n");
    ini_print=1;
  }
  fprintf(file_ave, "%3d ", N_pop);

  Print_mean(file_ave, f_sum, f_dev, it_sum, t_indip);
  Print_mean(file_ave, E_sum, E_dev, it_sum, t_indip);
  Print_mean(file_ave, DG_sum, DG_dev, it_sum, t_indip);
  fprintf(file_ave, " %.4f ", seq_entr);
  Print_mean(file_ave, Tload_sum, Tload_dev, Nload, Nload);
  Print_mean(file_ave, Mload_sum, Mload_dev, Nload, Nload);
  fprintf(file_ave, " %.3f  %.4f ", num_syn_subst/num_aa_subst,
	  1.-(num_syn_subst+num_aa_subst)/it_sum);
  fprintf(file_ave, " %.0f\n", n_subst);
  fflush(file_ave);
}

void Print_mean(FILE *file_out, float sum, float dev,
		float n_sum, float t_indip)
{
//   fprintf(file_out, " %.4f %.4f ", sum/n_sum,
// 	  sqrt((dev-sum*sum/n_sum)/(n_sum*t_indip)));
  fprintf(file_out, " %.8f %.8f ", sum/n_sum,
	  sqrt((dev-sum*sum/n_sum)/(n_sum*t_indip)));
}
    
void Print_final(char *name_file, double it_sum,
		 float TEMP, float s0,
		 int N_pop, float *freq_nuc,
		 double f_ave, double f_dev, 
		 double E_ave, double E_dev,
		 double DG_ave, double DG_dev,
		 float seq_entr, double seq_entr_dev,
		 double Tload_sum, double Tload_dev,
		 double Mload_sum, double Mload_dev,
		 int Nload,
		 double num_syn_subst, double num_aa_subst,
		 char *dna_seq, int len_dna)
{
  FILE *file_out; char name_out[200];
  short i, j, i_nuc, nuc_count[4][3], length;
//   float t_indep=num_aa_subst/10;
  float t_indep=1.0;

  sprintf(name_out, "%s_final.dat", name_file);
  printf("Writing %s\n", name_out);
  file_out=fopen(name_out, "w");

  // Headers
  fprintf(file_out, "# TEMP S0 Npop GC tt_ratio  ");
  fprintf(file_out, " fitness (sd) E_nat (sd) DG (sd) seq_entr (sd)");
  fprintf(file_out, " Trans.load   (sd)  Mut.load  (sd)");
  fprintf(file_out, " syn/nosyn  reject samples");
  for(i=0; i<3; i++){
    fprintf(file_out, "  ");
    for(j=0; j<3; j++)fprintf(file_out, "%c%d ", NUC_CODE[i], j+1);
  }
  fprintf(file_out, "\n#\n");

  // Parameters
  fprintf(file_out, "%.2f %.3f %3d ", TEMP, s0, N_pop);
  fprintf(file_out, " %.2f", freq_nuc[Code_nuc('G')]+freq_nuc[Code_nuc('C')]);
  fprintf(file_out, " %.1f    ", tt_ratio);

  // Results
  Print_mean(file_out, f_ave, f_dev, it_sum, t_indep);
  Print_mean(file_out, E_ave, E_dev, it_sum, t_indep);
  Print_mean(file_out, DG_ave, DG_dev, it_sum, t_indep);
  fprintf(file_out, " %7.4g %7.4g ", seq_entr, seq_entr_dev);
  Print_mean(file_out, Tload_sum, Tload_dev, Nload, Nload);
  Print_mean(file_out, Mload_sum, Mload_dev, Nload, Nload);

  // Synonimous, mutation load
  fprintf(file_out, " %.3f %.4f ", num_syn_subst/num_aa_subst,
	  1.-(num_syn_subst+num_aa_subst)/it_sum);
  fprintf(file_out, " %.0f  ", it_sum);

  // base content in DNA
  for(i=0; i<4; i++)
    for(j=0; j<3; j++)
      nuc_count[i][j]=0;
  j=0; length=len_dna/3;
  for(i=0; i<len_dna; i++){
    i_nuc=Code_nuc(dna_seq[i]);
    nuc_count[i_nuc][j]++; j++;
    if(j==3)j=0;
  }
  for(i=0; i<3; i++){
    fprintf(file_out," ");
    for(j=0; j<3; j++){
      fprintf(file_out," %.3f", nuc_count[i][j]/(float)length);
    }
  }
  fprintf(file_out, "\n");
  fclose(file_out);
}


int Read_ene_new(char *file_ene, float **interactions)
{
  FILE *file_in=fopen(file_ene,"r");
  int i, j, n=0; double norm=0; float ene;
  char string[200], aa1[8], aa2[8];

  if(file_in==NULL){
    printf("ERROR, energy parameter file %s not found%s\n", file_ene);
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
  int i, j; double norm=0;
  char string[200];
  if(file_in==NULL){
    printf("ERROR, energy parameter file %s not found%s\n", file_ene);
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
  char string[1000], *ptr;
  FILE *file_in=fopen(inputseq, "r");
  char *sequence, name[N_CHAR];
  int i=0, j=0, amm, i_nuc;

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
	i_nuc=Code_nuc(*ptr); i++;
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
  sprintf(name, "%s%s\0", name_file, ext);
  file_out=fopen(name, "w");
  printf("Writing %s\n", name);
  fprintf(file_out, "# File %s, sequence %s,  PDB %s,  length=%d\n",
	  name_file, seq_name, target.name, len_amm);
  fprintf(file_out, "# fitness: %s\n", fit_def);
  fprintf(file_out, "# %ld iterations, random seed: %ld\n", IT_MAX, iran);
  fprintf(file_out, "# Stationary frequencies: ");
  for(i=0; i<4; i++)fprintf(file_out, "%c %.3f ", NUC_CODE[i], freq_nuc[i]);
  fprintf(file_out, "\n# Transition-transversion ratio= %.2f\n", tt_ratio);
  for(i=0; i<4; i++)if((NUC_CODE[i]=='G')||(NUC_CODE[i]=='C'))gc+=freq_nuc[i];
  fprintf(file_out, "# GC bias: %.3f\n", gc);
  fprintf(file_out, "# T= %.3f\n", TEMP);
  fprintf(file_out, "# Configurational entropy per residue: %.3f\n", s0);
  fprintf(file_out, "# N_pop= %d\n#", N_pop);

  for(i=0; i<len_amm; i++){
    if(i==(i/60)*60)fprintf(file_out,"\n# ");
    fprintf(file_out,"%c", AMIN_CODE[seq[i]]);
  }
  fprintf(file_out,"\n");
  
  if(lf==1){
    fprintf(file_out,
	    "# site aa  Enat DG fitness syn_subst attempt\n"); //D(0,n)
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
  printf("ERROR in file name %s (file), dot not found\n");
  exit(8);
}

int Output_name(char *file_name, char *dir_out, char *prot_name,
		float TEMP, float S0, int N_pop, float *freq_nuc)
{
  char name[400];
  if(dir_out[0]!='\0'){sprintf(name, "%s/%s", dir_out, prot_name);}
  else{sprintf(name, "%s", prot_name);}
  sprintf(file_name, "%s_T%.2f_S0%.2f_N%d_GC%.2f\0",
	  name, TEMP, s0, N_pop,
	  freq_nuc[Code_nuc('G')]+freq_nuc[Code_nuc('C')]);
}
 
float Sequence_entropy(int **aa_distr, int L){
  int i, j; float S_sum=0, S, p; float norm=0;

  for(j=0; j<20; j++)norm+=aa_distr[0][j];
  for(i=0; i<L; i++){
    S=0;
    for(j=0; j<20; j++){
      if(aa_distr[i][j]){
	p=aa_distr[i][j]/norm; S-=p*log(p);
      }
    }
    S_sum+=S;
  }
  return(S_sum/L);
}

void Compute_freq_codons(float *freq_nuc, float *freq_aa,
			 char **codon, char *coded_aa)
{
  int i, i_aa; float w;
  for(i_aa=0; i_aa<20; i_aa++)freq_aa[i_aa]=0;
  for(i=0; i<64; i++){
    w=freq_nuc[Code_nuc(codon[i][0])];
    w*=freq_nuc[Code_nuc(codon[i][1])];
    w*=freq_nuc[Code_nuc(codon[i][2])];
    i_aa=Code_AA(coded_aa[i]);
    if(i_aa<0)continue;
    freq_aa[i_aa]+=w;
  }
}

float Sequence_entropy_mut(float *freq_nuc, char **codon, char *coded_aa)
{
  float freq_aa[20]; int i; float norm=0, S=0, p;

  // Calculating amino acid distribution under mutation alone
  Compute_freq_codons(freq_nuc, freq_aa, codon, coded_aa);

  // Compute entropy
  for(i=0; i<20; i++)norm+=freq_aa[i];
  for(i=0; i<20; i++){
    if(freq_aa[i]){p=freq_aa[i]/norm; S-=p*log(p);}
  }
  return(S);
}

void Compute_load(double *translation_load, double *mutation_load,
		  short *aa_seq, int len_amm, char *dna_seq, int len_dna,
		  float fitness_wt,  double E_nat_wt, double E1_wt,
		  double E2_wt, double E23_wt, double E3_wt,
		  float Conf_entropy, char **codon, char *coded_aa)
{
  // Mutations
  int nuc_mut=0, nuc_wt, base, pos=-1, j;
  int res_mut=0, aa_new, aa_old;
  char *codon_nat=dna_seq, codon_mut[3];
  // Folding thermodynamics
  float fitness, DeltaG;
  double sum_mut=0, sum_trans=0, p;
  double H_wt=E_nat_wt-Conf_entropy;
  double E_nat, E1, E2, E23=E23, E3;

  /*
    Load=Sum_j P(nat->Seq_j)[finess(nat)-fitness(Seq_j)]
    P(nat->Seq_j) is non zero only if sequence j is one base mutation from
    the native sequence.
    Translation load: P(nat->Seq_j)=1
    Mutation load: P(nat->Seq_j)=Mutation probability,
    i.e. f(new base) times 1 if transversion, times tt_ratio if transition    
   */

  *translation_load=0; *mutation_load=0;

  for(nuc_mut=0; nuc_mut<len_dna; nuc_mut++){
    
    pos++;
    if(pos==3){
      pos=0; res_mut++; codon_nat=dna_seq+res_mut*3;
    }
    for(j=0; j<3; j++)codon_mut[j]=codon_nat[j];
    nuc_wt=Code_nuc(dna_seq[nuc_mut]);

    for(base=0; base<4; base++){

      // Mutate amino acid
      if(base==nuc_wt)continue;
      codon_mut[pos]=NUC_CODE[base];
      aa_new=Coded_aa(codon_mut, codon, coded_aa);
      if(aa_new==aa_old)continue;         // Synonymous
      if(aa_new<0){fitness=0; goto load;} // Stop codon

      // Folding thermodynamics
      E_nat=H_wt; E1=E1_wt; E2=E2_wt;
      if(REM3){
	E23=E23_wt; E3=E3_wt;
	DeltaG=Mutate_DG_REM3(aa_seq, len_amm, res_mut, aa_new,
				  &E_nat, &E1, &E2, &E23, &E3);
      }else{
	DeltaG=
	  Mutate_DG_REM2(aa_seq, len_amm, res_mut, aa_new, &E_nat, &E1, &E2);
      }
      fitness=1./(1+exp(DeltaG));
      
    load:
      // Translation load
      sum_trans++;
      *translation_load+=fitness;

      // Mutation load
      p=freq_nuc[base];
      if(base==Transition(NUC_CODE[nuc_wt]))p*=tt_ratio;
      sum_mut+=p;
      *mutation_load+=fitness*p;
    }
  }

  *translation_load = fitness_wt-(*translation_load)/sum_trans;
  *mutation_load = fitness_wt-(*mutation_load)/sum_mut;

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
  double f_ratio, x;
  if((int)N_pop==1)return(1);
  if(fitness<=0)return(0);
  f_ratio= fitness_old / fitness;
  x= (1.-f_ratio)/(1.-pow(f_ratio, N_pop));
  if(RandomFloating() < x)return(1); return(0);
}
