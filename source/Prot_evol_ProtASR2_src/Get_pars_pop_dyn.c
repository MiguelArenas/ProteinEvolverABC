#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N_CHAR 400
#include "codes.h"
#include "get_para_pop_dyn.h"

static void help(char *);
static int
Read_parameters(char *FILE_IN,
		char *file_pdb, char *chain,
		char *file_pdb_list, char *pdb_dir,
		char *file_ali, char *file_dna,	
		char *tnm, char *tnm_mut_para,
		char **file_str_mut, int *N_str_mut, int *IWT,
		int *MEANFIELD,  int *MF_COMP,
		int *OPT_REG, int *SCORE_CV, float *REG,
		char *file_str, float *TEMP, float *sU1,
		float *sC0, float *sC1, int *REM, float *SEC_STR,
		int *REMUT, int *GET_FREQ, float *freq_nuc,
		float *trans_ratio, float *kCpG, float *TWONUC,
		int *PRINT_E, char *EXCHANGE, char *MATRIX,
		int *FORMAT, int *ALL_MUTS,
		int *IT_MAX, int *Samples, int *NEUTRAL, int *N_pop, 
		int *OPT_LAMBDA, float *LAMBDA, float *DG_OPT, char *MODEL,
		char *dir_out);
static int Find_string(char *flag, char **argv, int n_arg,int m,char *string);


int Get_para(int argc, char **argv,
	     char *file_pdb, char *chain, 
	     char *file_pdb_list, char *pdb_dir,
	     char *file_ali, char *file_dna,	
	     char *tnm, char *tnm_mut_para, 
	     char **file_str_mut, int *N_str_mut, int *IWT, 
	     int *MEANFIELD,  int *MF_COMP,
	     int *OPT_REG, int *SCORE_CV, float *REG,
	     char *file_str, float *TEMP, float *sU1,
	     float *sC0, float *sC1, int *REM, float *SEC_STR,
	     int *REMUT, int *GET_FREQ, float *freq_nuc,
	     float *trans_ratio, float *kCpG, float *TWONUC,
	     int *PRINT_E, char *EXCHANGE, char *MATRIX,
	     int *FORMAT, int *ALL_MUTS,
	     int *IT_MAX, int *Samples, int *NEUTRAL, int *N_pop, 
	     int *OPT_LAMBDA, float *LAMBDA, float *DG_OPT, char *MODEL,
	     char *dir_out)
{
  char FILE_IN[N_CHAR]="", string[N_CHAR];
  int r, i, j, i_nuc; float sum=0, p;

  if(Find_string("-h", argv, argc, 2, string))help(argv[0]);
  if(argc < 2){
    printf("ERROR, either input file or options must be specified\n");
    help(argv[0]);
  }
  file_pdb[0]='\0';

  for(j=1; j<argc; j++){
    if(strncmp(argv[j], "-file", 5)==0){
      j++; strcpy(FILE_IN, argv[j]); break;
    }
  }
  if(FILE_IN[0]!='\0'){
    printf("Reading parameters for %s from %s\n", argv[0], FILE_IN);
    r=Read_parameters(FILE_IN, file_pdb, chain, file_pdb_list, pdb_dir, 
		      file_ali, file_dna,
		      tnm, tnm_mut_para, file_str_mut, N_str_mut, IWT,
		      MEANFIELD, MF_COMP, OPT_REG, SCORE_CV, REG, 
		      file_str, TEMP, sU1, sC0, sC1, REM, SEC_STR,
		      REMUT, GET_FREQ, freq_nuc, trans_ratio, kCpG, TWONUC,
		      PRINT_E, EXCHANGE, MATRIX, FORMAT, ALL_MUTS,
		      IT_MAX, Samples, NEUTRAL, N_pop, 
		      OPT_LAMBDA, LAMBDA, DG_OPT, MODEL, dir_out);
    if(r==0){
      printf("ERROR, input file %s does not exist\n", FILE_IN);
      help(argv[0]);
    }
  }


  /*********************************************/
  char w[80];
    printf("! Changing parameters with command line:\n");
  for(j=1; j<argc; j++){
    if(strncmp(argv[j], "-file", 5)==0){
      j++; strcpy(FILE_IN, argv[j]);
    }else if(strncmp(argv[j], "-pdblist", 9)==0){
      j++; sscanf(argv[j], "%s", file_pdb_list);
    }else if(strncmp(argv[j], "-pdbdir", 8)==0){
      j++; sscanf(argv[j], "%s", pdb_dir);
    }else if(strncmp(argv[j], "-pdb", 4)==0){
      j++; sscanf(argv[j], "%s", file_pdb);
    }else if(strncmp(argv[j], "-tnm", 4)==0){
      j++; sscanf(argv[j], "%s", tnm);
    }else if(strncmp(argv[j], "-mut_para", 9)==0){
      j++; sscanf(argv[j], "%s", tnm_mut_para);
    }else if(strncmp(argv[j], "-chain", 3)==0){
      j++; sscanf(argv[j], "%s", chain);
    }else if(strncmp(argv[j], "-ali", 4)==0){
      j++; sscanf(argv[j], "%s", file_ali);
    }else if(strncmp(argv[j], "-dna", 4)==0){
      strcpy(file_dna, argv[j+1]); 
      printf("Sequence file: %s\n", file_dna);

    }else if(strncmp(argv[j], "-temp", 5)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p>0){
	*TEMP=p; printf("Temperature: %3f\n", p); j++;
      }

    }else if(strncmp(argv[j], "-pop", 4)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p >=1){
	*N_pop=p; printf("Population size: %d\n", *N_pop); j++;
      }
    }else if(strncmp(argv[j], "-tt", 3)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p>1){
	*trans_ratio=p; j++;
	printf("Transition-transversion ratio: %.2f\n", *trans_ratio);
      }
    }else if(strncmp(argv[j], "-meanfield", 10)==0){
      *MEANFIELD=1;
      printf("Meanfield computation of substitution rates\n");
    }else if(strncmp(argv[j], "-model", 6)==0){
      j++; sscanf(argv[j], "%s", w);
      if((strcmp(w, "NAT")==0)||
	 (strcmp(w, "ALL")==0)||
	 (strcmp(w, "DG")==0)){
	strcpy(MODEL, w);
      }else{
	printf("WARNING, %s is not an allowed model. ", w);
	printf("Allowed models are NAT ALL DG\n");
	printf("Using default: %s\n", MODEL);
      }
    }else if(strncmp(argv[j], "-lambda", 7)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p>0){
	*LAMBDA=p; printf("Lambda parameter for meanfield: %3f\n", p); j++;
	*OPT_LAMBDA=0;
      }
    }else if(strncmp(argv[j], "-sU1", 3)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p>0){
	*sU1=p; printf("Unfolded entropy per residue: %3f\n", *sU1); j++;
      }
    }else if(strncmp(argv[j], "-it", 3)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p>10){
	*IT_MAX=p; printf("Substitutions per residue: %d\n", *IT_MAX); j++;
      }
    }else if(strncmp(argv[j], "-freq", 5)==0){
      j++;
      if(strcmp(w, "nuc")==0){*GET_FREQ=1;}
      else if(strcmp(w, "aa")==0){*GET_FREQ=2;}
      else if(strcmp(w, "input")==0){*GET_FREQ=0;}
      else{
	printf("WARNING, %s is not an allowed fit for frequencies. ", w);
	printf("Allowed fits are aa, nuc, input\n");
	printf("Using default GET_FREQ= %d\n", *GET_FREQ);
      }
    }else if(strncmp(argv[j], "-f", 2)==0){
      sscanf(argv[j]+2, "%s", w); i_nuc=Code_nuc(w[0]);
      if((i_nuc<0)||(i_nuc>3)){
	printf("WARNING, %c is not an allowed nucleotideºn", w[0]);
	continue;
      }
      sscanf(argv[j+1], "%f", &p);
      if((p>0)&&(p<1)){
	freq_nuc[i_nuc]=p; j++;
	printf("freq_nuc(%c): %f\n", NUC_CODE[i_nuc], freq_nuc[i_nuc]);
      }
    }else if(strncmp(argv[j], "-gc", 3)==0){
      sscanf(argv[j+1], "%f", &p);
      if((p>0)&&(p<1)){
	freq_nuc[Code_nuc('G')]=0.5*p;
	freq_nuc[Code_nuc('C')]=0.5*p;
	freq_nuc[Code_nuc('T')]=0.5*(1-p);
	freq_nuc[Code_nuc('A')]=0.5*(1-p);
	for(i=0; i<4; i++)
	  printf("freq_nuc(%c): %3f\n", NUC_CODE[i], freq_nuc[i]);
	j++;
      }
      /*}else if(strncmp(argv[j], "-gencode", 8)==0){
      strcpy(file_code, argv[j+1]); j++;
      printf("Genetic code: %s\n", file_code);*/
    }else if(strncmp(argv[j], "-d", 2)==0){
      strcpy(dir_out,argv[j+1]); j++;
      printf("Directory in is: %s\n",dir_out);
    }else{
      printf("WARNING, unrecognized option %s\n", argv[j]);
    }
  }

  printf("! End changing parameters\n\n");

  /*********************************************/

  // Checks
  if(file_pdb[0]=='\0' && file_pdb_list[0]=='\0'){
    printf("ERROR, PDB file not specified\n");
    help(argv[0]);
  }else if(file_pdb[0] && file_pdb_list[0]){
    printf("WARNING, PDB files input both as list and as single file "
	   "I will ignore the latter\n");
  }

  if((*N_pop < 1)&&(*NEUTRAL==0)&&(*MEANFIELD==0)){
    printf("ERROR, %d is too small population size\n", *N_pop);
    printf("Check the input file %s\n", FILE_IN);
  }

  for(i=0; i<4; i++)sum+=freq_nuc[i];
  if(sum!=1){
    printf("Warning, nucleotide frequencies are not normalized.\n");
    printf("Normalizing\n");
  }
  for(i=0; i<4; i++)freq_nuc[i]/=sum;
  return(0);
}

int Read_parameters(char *FILE_IN,
		    char *file_pdb, char *chain, 
		    char *file_pdb_list, char *pdb_dir,
		    char *file_ali, char *file_dna,	
		    char *tnm, char *tnm_mut_para,
		    char **file_str_mut, int *N_str_mut, int *IWT,
		    int *MEANFIELD,  int *MF_COMP,
		    int *OPT_REG, int *SCORE_CV, float *REG,
		    char *file_str, float *TEMP, float *sU1,
		    float *sC0, float *sC1, int *REM, float *SEC_STR,
		    int *REMUT, int *GET_FREQ, float *freq_nuc,
		    float *trans_ratio, float *kCpG, float *TWONUC,
		    int *PRINT_E, char *EXCHANGE, char *MATRIX,
		    int *FORMAT, int *ALL_MUTS,
		    int *IT_MAX, int *Samples, int *NEUTRAL, int *N_pop, 
		    int *OPT_LAMBDA, float *LAMBDA, float *DG_OPT, char *MODEL,
		    char *dir_out)
{
  FILE *file_in=fopen(FILE_IN, "r");
  char string[1000], nuc[3];
  int i_nuc, idum; float p, sum=0;

  if(file_in==NULL){
    printf("WARNING, input file %s does not exist\n", FILE_IN);
    return(0);
  }
  printf("Reading parameters in %s\n", FILE_IN);

  char READ[100];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    // Protein files
    if(strncmp(string, "PDB=", 4)==0){
      sscanf(string+4,"%s", file_pdb);
    }else if(strncmp(string, "CHAIN=", 6)==0){
      sscanf(string+6,"%s", chain);
    }else if(strncmp(string, "PDB_LIST=", 9)==0){
      sscanf(string+10,"%s", file_pdb_list);
    }else if(strncmp(string, "ALI=", 4)==0){
      sscanf(string+4,"%s", READ);
      if(strncmp(READ, "NULL", 4)==0)continue;
      strcpy(file_ali, READ);
    }else if(strncmp(string, "DNA=", 4)==0){
      sscanf(string+4,"%s", READ);
      if(strncmp(READ, "NULL", 4)==0)continue;
      strcpy(file_dna, READ);
    }else if(strncmp(string, "TNM=", 4)==0){
      sscanf(string+5,"%s", tnm);
    }else if(strncmp(string, "STR_MUT", 7)==0){
      sscanf(string+8,"%s", READ);
      if(strncmp(READ, "NULL", 4)==0)continue;
      if(*N_str_mut<4){
	strcpy(file_str_mut[*N_str_mut], READ); (*N_str_mut)++;
      }else{
	printf("WARNING, 2 files with structural effects of mutations ");
	printf("have been already read, file STR_MUT=%s will be ignored\n",
	       READ);
      }
    }else if(strncmp(string, "IWT", 3)==0){
      sscanf(string+4,"%d", &idum);
      if(idum<=0 || idum>20){
	printf("WARNING, IWT= %d is not allowed, using default %d\n",
	       idum, *IWT);
      }else{*IWT=idum;}
      // Computations
    }else if(strncmp(string, "ALL_MUTS", 8)==0){
      sscanf(string+9,"%d", ALL_MUTS);
      // Substitution matrices
    }else if(strncmp(string, "MEANFIELD=", 10)==0){
      sscanf(string+10,"%d", MEANFIELD);
    }else if(strncmp(string, "MF_COMP", 7)==0){
      sscanf(string+8,"%d", MF_COMP);
    }else if(strncmp(string, "OPT_REG", 7)==0){
      sscanf(string+8,"%d", OPT_REG);
    }else if(strncmp(string, "SCORE_CV", 8)==0){
      sscanf(string+9,"%d", SCORE_CV);
    }else if(strncmp(string, "REG", 3)==0){
      sscanf(string+4,"%f", &p);
      if((p<0)){
	printf("WARNING, regularization %.4f is not an allowed value\n",p);
	printf("Keeping default %.4f\n", *REG);
      }else{
	*REG=p;
      }
      // Exchangeability matrices
    }else if(strncmp(string, "PRINT_E", 7)==0){
      sscanf(string+8,"%d", PRINT_E);
    }else if(strncmp(string, "FORMAT=", 7)==0){
      char format[80];
      sscanf(string+7,"%s", format);
      if(strncmp(format, "PAML", 4)==0){
	*FORMAT=1;
      }else{
	printf("WARNING, non valid option FORMAT= %s\n", format);
      }
    }else if(strncmp(string, "EXCHANGE=", 9)==0){
      char format[80];
      sscanf(string+9,"%s", format);
      if(strncmp(format, "MUT", 3)==0){
	*EXCHANGE='M';
      }else if(strncmp(format, "FLUX", 4)==0){
	*EXCHANGE='F';
      }else if(strncmp(format, "RATE", 4)==0){
	*EXCHANGE='Q';
      }else if(strncmp(format, "EXCH", 4)==0){
	*EXCHANGE='E';
      }else{
	printf("WARNING, non valid option EXCHANGE= %s\n", format);
	printf("Valid options: FLUX (default), RATE, EXCH, MUT\n");
      }
    }else if(strncmp(string, "MATRIX=", 6)==0){
      char format[80];
      sscanf(string+7,"%s", format);
      if((strncmp(format, "LG", 2)!=0)&&
	 (strncmp(format, "WAG", 3)!=0)&&
	 (strncmp(format, "JTT", 3)!=0)&&
	 (strncmp(format, "OPT", 3)!=0)){
	printf("WARNING, matrix %s not defined\n", format);
	printf("Valid options are LG WAG JTT OPT\n");
	printf("Keeping default %s\n", MATRIX);
      }else{
	strcpy(MATRIX, format);
      }
      printf("MATRIX= %s\n", MATRIX);
      // Old stuff
    }else if(strncmp(string, "OPT_LAMBDA=", 11)==0){
      sscanf(string+11,"%d", OPT_LAMBDA);
    }else if(strncmp(string, "MODEL=", 6)==0){
      char w[80]; sscanf(string+6, "%s", w);
      if((strcmp(w, "NAT")!=0)&&
	 (strcmp(w, "ALL")!=0)&&
	 (strcmp(w, "DG")!=0)){
	printf("WARNING, %s is not an allowed model. ", w);
	printf("Allowed models are NAT ALL DG\n");
	printf("Using default: %s\n", MODEL);
      }else{
	strcpy(MODEL, w);
      }
    }else if(strncmp(string, "DG_OPT=", 7)==0){
      sscanf(string+7,"%f", DG_OPT);
    }else if(strncmp(string, "LAMBDA=", 7)==0){
      sscanf(string+7,"%f", LAMBDA);
      //printf("Lambda parameter for meanfield: %1f\n", *LAMBDA);
      // Amino acid frequencies
    }else if(strncmp(string, "REMUT", 5)==0){
      sscanf(string+6,"%d", REMUT);
    }else if(strncmp(string, "GET_FREQ=", 9)==0){
      sscanf(string+9,"%d", GET_FREQ);
    }else if(strncmp(string, "FREQ", 4)==0){
      sscanf(string+4,"%s%f", nuc, &p);
      i_nuc=Code_nuc(nuc[0]); freq_nuc[i_nuc]=p; sum+=p;
    }else if(strncmp(string, "TT_RATIO=", 9)==0){
      sscanf(string+9,"%f", trans_ratio);
    }else if(strncmp(string, "kCpG=", 4)==0){
      sscanf(string+5,"%f", kCpG); 
    }else if(strncmp(string, "TWONUCMUT=", 10)==0){
      sscanf(string+10,"%f", TWONUC);
      if((*TWONUC < 0)||(*TWONUC >1)){
	printf("ERROR, the parameter TWONUCMUT must be set ");
	printf("between zero and one, but it is %.3f\n", *TWONUC);
	exit(8);
      }

      // Thermodynamic parameters:
    }else if(strncmp(string, "TEMP=", 5)==0){
      sscanf(string+5,"%f", TEMP);
    }else if(strncmp(string, "SU1=", 3)==0){
      sscanf(string+4,"%f", sU1); 
    }else if(strncmp(string, "SC1=", 4)==0){
      sscanf(string+4,"%f", sC1); 
    }else if(strncmp(string, "SC0=", 4)==0){
      sscanf(string+4,"%f", sC0); 
    }else if(strncmp(string, "REM=", 4)==0){
      int dum;
      sscanf(string+4,"%d", &dum);
      if((dum<0)||(dum>3)){
	printf("WARNING, the variable REM can only be set to 0, 1, 2 or 3\n");
      }else{
	*REM=dum;
      }
    }else if(strncmp(string, "A_LOC", 5)==0){
      sscanf(string+6,"%f", SEC_STR);
    }else if(strncmp(string, "FILE_STR=", 9)==0){
      sscanf(string+9,"%s", file_str);
      // Simulations:
    }else if(strncmp(string, "TMAX=", 5)==0){
      sscanf(string+5,"%d", IT_MAX);
    }else if(strncmp(string, "Samples=", 8)==0){
      sscanf(string+8,"%d", Samples);
    }else if(strncmp(string, "NEUTRAL=", 8)==0){
      sscanf(string+8,"%d", NEUTRAL);
    }else if(strncmp(string, "NPOP=", 5)==0){
      sscanf(string+5,"%d", N_pop);
    }else{
      printf("WARNING, uninterpreted line:\n%s", string); 
    }
  }
  fclose(file_in);
  return(1);
}

int Find_string(char *flag, char **argv, int n_arg, int m, char *string){
  int i;
  for(i=0; i<n_arg; i++){
    if(strncmp(argv[i], flag, m)==0){
      i++; if(i<n_arg)strcpy(string, argv[i]);
      return(1);
    }
  }
  return(0);
}

void help(char *prog){

  printf("\n\nFORMAT of input file:\n");
  printf("#================================================================\n");
  printf("# A) Input files\n");
  printf("PDB=/data/ortizg/databases/pdb/1tre.pdb  # file_pdb\n");
  printf("CHAIN= A (for specific chain) ALL (for all chains)\n");
  printf("FILE_STR=/home/ubastolla/RESEARCH/PROT_EVOL/INPUT/structures.in\n");
  printf("# List of alternative contact matrices for misfold computations\n");
  printf("SEQ=	tpis.dna	# nucleotifde sequence (optional)\n");
  printf("ALI=	family.aln	# FASTA file with MSA (optional)\n");
  printf("STR_MUT=1treA_mut_DE.dat # Structural effects of mutations\n");
  printf("#                           computed with the program TNM\n");
  printf("#================================================================\n");
  printf("## B) Substitution models\n");
  printf("MEANFIELD= 1	     # Generate substitution models?\n");
  printf("# Subst. models are generated based on folding stability, struct.\n");
  printf("# conservation and combination of both. Parameters are amino acid\n");
  printf("# frequencies and selection parameter Lambda, which is determined\n");
  printf("# minimizing the KL divergence between model and regularized\n");
  printf("# distribution from PDB seq and input MSA.\n");
  printf("OPT_REG=1    ! Automatically determine the regularization param.\n");
  printf("REG=0.35     ! Regularization param. if OPT_REG=0\n");
  printf("SCORE_CV=1   ! Optimize REG with Cv (1) or |KL_mod-KL_reg| (0)\n");
  printf("MF_COMP=1    ! Perform (1) or omit (0) mean-field computations of\n");
  printf("# stability constrained model (slow), otherwise only wild-type\n");
  printf("# computation is performed (faster and often better performing).\n");
  printf("#================================================================\n");
  printf("# C)  Amino acid frequencies\n");
  printf("REMUT=1             # Determine a.a. freq twice, the first time\n");
  printf("# time fitting observed frequencies with a.a. frequencies alone,\n");
  printf("# the second time fitting observed frequencies with full model\n");
  printf("# that includes selection.\n");
  printf("GET_FREQ=3	      # Allowed: 0,1,2,3\n");
  printf("# 0= Use input mutation parameters\n");
  printf("# 1= Fit mutation parameters from prot sequences\n");
  printf("# 2= Combine fitted mutation model and a.a. frequencies\n");
  printf("# 3= Get background distribution from amino acid frequencies\n");
  printf("# Parameters of the mutation model if GET_FREQ=0:\n");
  printf("FREQ A 0.25	       # f(A)\n");
  printf("FREQ T 0.25	       # f(T)\n");
  printf("FREQ C 0.25	       # f(C)\n");
  printf("FREQ G 0.25	       # f(G)\n");
  printf("kCpG=5	       # Enhanced rate at CpG dinucleotides\n");
  printf("TT_RATIO=1	       # transition-transversion ratio (not CpG)\n");
  printf("TWONUCMUT=0.1	       # Ratio between 1-nuc and 2-nuc mutations\n");
  printf("#===============================================================\n");
  printf("# D) Exchangeability matrix\n");
  printf("EXCHANGE=FLUX	       # Allowed: FLUX (default), RATE, EXCH, MUT\n");
  printf("MATRIX=JTT	       # Empirical exchange matrix (JTT, WAG)\n");
  printf("#===============================================================\n");
  printf("# E) Output\n");
  printf("PRINT_E=0            # Print exchangeability matrix at all sites?\n");
  printf("FORMAT=PAML	       # Use PAML format for exchangeability matrix\n");
  printf("ALL_MUTS=1           # Predict the effect of all nucl. mut?\n");
  printf("#================================================================\n");
  printf("# F) Thermodynamic model\n");
  printf("TEMP=	0.5	       # Temperature\n");
  printf("SU1=	0.065	       # configurational entropy per res (unfold)\n");
  printf("SC1=  0.065	       # configurational entropy per res (misfold)\n");
  printf("SC0=  0.0	       # configurational entropy offset (misfold)\n");
  printf("REM=   2	       # Use 0,1,2 moments of misfolding energy\n");
  printf("A_LOC= 1             # Use secondary structure propensities?\n");
  printf("#================================================================\n");
  printf("# G) Simulations of evolution\n");
  printf("TMAX=   000		# ITMAX: # of substitutions\n");
  printf("Samples= 1		# Independent trajectories simulated\n");
  printf("NEUTRAL= 0		# ");
  printf("#1:Neutral fitness landscape 0: Fitness=1/(1+exp(DG/T))\n");
  printf("NPOP=	10		");
  printf("# effective population size (if MEANFIELD=0, NEUTRAL=0)\n");
  printf("\n\n");

  printf("USAGE %s options:\n", prog);
  printf("  -pdblist <file.pdb_list>   # List of pdb files\n");
  printf("  -pdbdir  <path of pdb files>\n");
  printf("  -pdb <file.pdb>   # Unique pdb file, mandatory if not list\n");
  printf("  -chain <chain>    # Chain identifiers (ex. AB)\n");
  printf("  -ali <file>       # file with aligned proteins (opt)\n");
  printf("  -tnm <tnm program for computing str. effect of muts. (opt)>\n");
  printf("  -mut_para <parameters for computing str. effect of muts. (opt)>\n");
  printf("  -dna  <file_dna>  # file with DNA sequece (optional)\n");
  printf("  -file <configuration file> (ex. Prot_evol.in), optional\n");
  printf("Evolution simulations:\n");
  printf("  -it  <IT_MAX>     # Number of iterations per individual and a.a\n");
  printf("  -pop  <N_pop>     # population size\n");
  printf("Mean-field computations:\n");
  printf("  -meanfield        # Meanfield computation of subst. rates\n");
  printf("  -lambda <lambda>  # Parameter for meanfield\n");
  printf("                    # If lambda not given, use optimal lambda\n");
  printf("Mutational model:\n");
  printf("  -freq <FREQ>      # Criterion to fit frequencies. Allowed:\n");
  printf("                    # nuc (fit nucleotides from prot sequence)\n");
  printf("                    # aa (fit amino acids from prot sequence)\n");
  printf("                    # input (get nucleotides from input)\n");
  printf("  -gc <GC_bias>     # frequency of nucleotides G+C\n");
  printf("  -fG <freq>        # Input frequency for each nucleotide\n");
  printf("  -tt <trans_ratio> # transition-transversion ratio (>1)\n");
  printf("Thermodynamic model:\n");
  printf("  -temp <TEMP>      # temperature\n");
  printf("  -sU1  <sU1>       # Unfolded entropy per residue\n");


  //printf("-gencode <file>     # file with genetic code\n");
  //printf("-f<n> <freq_nuc>    # frequency of nucleotide n, n=A,T,G,C\n");


  exit(8);
}
