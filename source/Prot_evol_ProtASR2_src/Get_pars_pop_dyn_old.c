#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N_CHAR 400
#include "codes.h"
#include "get_para_pop_dyn.h"

static void help();
static int
Read_parameters(char *FILE_IN, char *file_pdb, char *chain,
		char *file_seq, char *file_str, 
		int *N_pop, float *TEMP, float *s0, float *sC0,
		float *sC1, int *REM3, long *IT_MAX,
		float *freq_nuc, float *trans_ratio);
static int Find_string(char *flag, char **argv, int n_arg,int m,char *string);


int Get_para(int argc, char **argv, char *file_pdb,
	     char chain[], char *file_seq, char *file_str,
	     int *N_pop, float *TEMP, float *s0,
	     float *sC0, float *sC1, long *IT_MAX,
	     float *freq_nuc, float *trans_ratio, int *REM3,
	     char *dir_out)
{
  char FILE_IN[N_CHAR], string[N_CHAR];
  int i, j, i_nuc; float sum=0, p;

  if(Find_string("-h", argv, argc, 2, string))help();
  if(argc < 2){
    printf("ERROR, input file name must be specified\n");
    help();
  }

  strcpy(FILE_IN, argv[1]);
  Read_parameters(FILE_IN, file_pdb, chain, file_seq, file_str,
		  N_pop, TEMP, s0, sC0, sC1, REM3, IT_MAX,
		  freq_nuc, trans_ratio);

  /*********************************************/

  if(argc <= 2)goto check;
  printf("! Changing parameters with command line:\n");

  for(j=2; j<argc; j++){
    if(strncmp(argv[j], "-pop", 4)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p >=1){
	*N_pop=p; printf("Population size: %d\n", *N_pop); j++;
      }
    }else if(strncmp(argv[j], "-seq", 4)==0){
      strcpy(file_seq, argv[j+1]); 
      printf("Sequence file: %s\n", file_seq);
    }else if(strncmp(argv[j], "-tt", 3)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p>1){
	*trans_ratio=p; j++;
	printf("Transition-transversion ratio: %.2f\n", *trans_ratio);
      }
    }else if(strncmp(argv[j], "-temp", 5)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p>0){
	*TEMP=p; printf("Temperature: %3f\n", p); j++;
      }
    }else if(strncmp(argv[j], "-s0", 3)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p>0){
	*s0=p; printf("Entropy per residue: %3f\n", *s0); j++;
      }
    }else if(strncmp(argv[j], "-it", 3)==0){
      sscanf(argv[j+1], "%f", &p);
      if(p>10){
	*IT_MAX=p; printf("Substitutions per residue: %ld\n", *IT_MAX); j++;
      }
    }else if(strncmp(argv[j], "-f", 2)==0){
      sscanf(argv[j]+2, "%d", &i_nuc);
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
 check:
  if(*N_pop < 1){
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

int Read_parameters(char *FILE_IN, char *file_pdb, char *chain,
		    char *file_seq, char *file_str, 
		    int *N_pop, float *TEMP, float *s0, float *sC0,
		    float *sC1, int *REM3, long *IT_MAX,
		    float *freq_nuc, float *trans_ratio)
{
  FILE *file_in=fopen(FILE_IN, "r");
  char string[1000], nuc[3], dumm[80];
  int i_nuc; float p, sum=0;

  if(file_in==NULL){
    printf("ERROR, input file %s does not exist\n", FILE_IN);
    help();
  }
  printf("Reading parameters in %s\n", FILE_IN);

  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "PDB=", 4)==0){
      sscanf(string+4,"%s", file_pdb);
    }else if(strncmp(string, "CHAIN=", 6)==0){
      sscanf(string+6,"%s", dumm); *chain=dumm[0];
    }else if(strncmp(string, "SEQ=", 4)==0){
      sscanf(string+4,"%s", file_seq);
    }else if(strncmp(string, "NPOP=", 5)==0){
      sscanf(string+5,"%d", N_pop);
    }else if(strncmp(string, "TEMP=", 5)==0){
      sscanf(string+5,"%f", TEMP);
    }else if(strncmp(string, "S0=", 3)==0){
      sscanf(string+3,"%f", s0); 
    }else if(strncmp(string, "SC1=", 4)==0){
      sscanf(string+4,"%f", sC1); 
    }else if(strncmp(string, "SC0=", 4)==0){
      sscanf(string+4,"%f", sC0); 
    }else if(strncmp(string, "REM3=", 5)==0){
      sscanf(string+5,"%d", REM3);
    }else if(strncmp(string, "FREQ", 4)==0){
      sscanf(string+4,"%s%f", nuc, &p);
      i_nuc=Code_nuc(nuc[0]); freq_nuc[i_nuc]=p; sum+=p;
    }else if(strncmp(string, "TT_RATIO=", 9)==0){
      sscanf(string+9,"%f", trans_ratio);
    }else if(strncmp(string, "TMAX=", 5)==0){
      sscanf(string+5,"%ld", IT_MAX);
    }else if(strncmp(string, "FILE_STR=", 9)==0){
      sscanf(string+9,"%s", file_str);
    }else{
      printf("WARNING, uninterpreted line:\n%s", string); 
    }
  }
  fclose(file_in);
  return(0);
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

void help(){
  printf("USAGE Pop_evol <file name> ");
  printf("(ex. Pop_evol.in) with default parameters\n");
  printf("ARGUMENTS:\n"); 
  printf("-temp <TEMP>        # temperature\n");
  printf("-pop  <N_pop>       # population size\n");
  printf("-seq  <file_seq>    # file with DNA sequece\n");
  printf("-tt   <trans_ratio> # transition-transversion ratio (>1)\n");
  printf("-s0   <s0>          # Configurational entropy per residue\n");
  printf("-it   <IT_MAX>      # Number of iterations per individual and a.a\n");
  printf("-gc   <GC_bias>     # frequency of nucleotides G+C\n");
  printf("-gencode <file>     # file with genetic code\n");
  printf("-d <directory>      # directory name for the output files\n");
  //printf("-f<n> <freq_nuc>    # frequency of nucleotide n, n=A,T,G,C\n");

  printf("\n\nFORMAT of input file:\n");
  printf("PDB= <pdb file or pdb code in contact matrix file>\n");
  printf("FILE_STR= <file with alternative contact matrices>\n");
  printf("CHAIN= <native protein chain identifier> (default: first chain)\n");
  printf("SEQ= <file with DNA sequence of the PDB protein>\n");
  printf("NPOP= <Effective population size>\n");
  printf("TEMP= <Environmental temperature of evolution>\n");
  printf("S0=   <Unfolded conformational entropy per residue>\n");
  printf("SC1=  <Misfolded conformational entropy per residue>\n");
  printf("SC0=  <Misfolded conformational entropy, offset>\n");
  printf("REM3= 0: Only second moment of misfolded energy, 1: third moment\n");
  printf("FREQ A 0.25 <Nucleotide frequencies>\n");
  printf("FREQ A 0.25 <Nucleotide frequencies>\n");
  printf("FREQ A 0.25 <Nucleotide frequencies>\n");
  printf("FREQ A 0.25 <Nucleotide frequencies>\n");
  printf("TT_RATIO= <Transition-transversion ratio>\n");
  printf("TMAX= <Simulation length (number of mutations)>\n");
  exit(8);
}
