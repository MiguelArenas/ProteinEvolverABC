#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "codes.h"
#include "input.h"
//#define AMIN_CODE "AEQDNLGKSVRTPIMFYCWHX"
//#define NUC_CODE  "UACG"
//#define NUC_CODE "ATGC"
#define CODE_NAME "Codes.c"

/* Numeric codes */
int Code_AA(char res){
  short i; char r;
  //i=(int)res; if(i>96){r=(char)(i-32);}else{r=res;}
  if(res > 'a'){r=(char)((int)res-32);}else{r=res;}
  for(i=0; i<20; i++)if(r==AMIN_CODE[i])return(i);
  if(r=='*')return(-1); // STOP codon
  if((res=='-')||(res=='.'))return(20);
  printf("WARNING, wrong aa type %c\n", res);
  return(0);
}

int Code_nuc(char nuc){
  int i;
  //char n; i=(int)nuc; if(i>96){n=(char)(i-32);}else{n=nuc;}
  //if(n=='T')n='U';
  for(i=0; i<4; i++)if(nuc==NUC_CODE[i])return(i);
  printf("Error, nuc %c\n", nuc);
  exit(8);
}


/* One letter codes */
char Nuc_code(int nuc){
  if((nuc<0)||(nuc>3)){
    printf("Error, wrong nucleotide code %d\n", nuc); exit(8);
  }
  return(NUC_CODE[nuc]);
}

char Amin_code(int amm){
  if((amm<0)||(amm>20)){
    printf("Error, wrong amino acid code %d\n", amm); exit(8);
  }
  return(AMIN_CODE[amm]);
}

int Transition(char nuc){
  if(nuc=='A'){
    return(Code_nuc('G'));
  }else if((nuc=='T')||(nuc=='U')){
    return(Code_nuc('C'));
  }else if(nuc=='G'){
    return(Code_nuc('A'));
  }else if(nuc=='C'){
    return(Code_nuc('T'));
  }
  printf("Error in Transition(), wrong nucleotide %c\n", nuc);
  exit(8);
}

void Read_nuc_freq(char *file_nuc_freq, float *freq_nuc,
	      float *trans_ratio, float *mut_rate)
{
  FILE *file_in=
    Open_file_r(file_nuc_freq, CODE_NAME, "nucleotide frequencies");
  char string[200], nuc[3]; int i; float f, sum=0;

  printf("Reading stationary nucleotide frequencies in %s\n",file_nuc_freq);
  for(i=0; i<4; i++){
    fgets(string, sizeof(string), file_in);
    sscanf(string, "%s %f", nuc, &f);
    freq_nuc[Code_nuc(nuc[0])]=f; sum+=f;
    //printf("%c %.3f\n", nuc[0], f);
  }
  fgets(string, sizeof(string), file_in);
  sscanf(string, "%f", trans_ratio);
  //printf("Transition-transversion ratio = %.1f\n", *trans_ratio);
  fgets(string, sizeof(string), file_in);
  sscanf(string, "%f", mut_rate);
  //printf("Mutation rate = %.2g\n", *mut_rate);
  fclose(file_in);

  if(sum!=1){
    //printf("Warning, frequencies not normalized. Normalizing:\n");
    for(i=0; i<4; i++){
      freq_nuc[i]/=sum; printf("%c %.3f\n", Nuc_code(i), freq_nuc[i]);
    }
  }
}

void Read_mut_mat(char *file_mut_mat, float **nuc_mut_mat)
{
  FILE *file_in=Open_file_r(file_mut_mat, CODE_NAME, "mutation matrix");
  int line, i_nuc[4], i, j; float f[4], sum=0, rate;
  char string[200], nuc[4][3];

  printf("Reading mutation matrix in %s\n", file_mut_mat);
  fgets(string, sizeof(string), file_in);
  sscanf(string, "%s%s%s%s", nuc[0], nuc[1], nuc[2], nuc[3]);
  for(i=0; i<4; i++){
    i_nuc[i]=Code_nuc(nuc[i][0]); printf("%s ", nuc[i]);
  }
  printf("\n");
  for(line=0; line<4; line++){
    fgets(string, sizeof(string), file_in);
    sscanf(string, "%s %f%f%f%f", nuc[0], &f[0], &f[1], &f[2], &f[3]);
    j=Code_nuc(nuc[0][0]);
    if(j!=i_nuc[line]){
      printf("Error, different order in lines and columns\n"); exit(8);
    }
    printf("%s ", nuc[0]);
    for(i=0; i<4; i++){
      nuc_mut_mat[j][i_nuc[i]]=f[i]; printf("%.3f", f[i]);
    }
    printf("\n");
  }
  fgets(string, sizeof(string), file_in);
  sscanf(string, "%f", &rate);
  fclose(file_in);

  for(i=0; i<4; i++){
    sum=0;
    for(j=0; j<4; j++)if(j!=i)sum+=nuc_mut_mat[i][j];
    if(nuc_mut_mat[i][i]!=0)
      printf("Warning, wrong matrix format, M[%d][%d] expected zero\n",
	     i, i);
    for(j=0; j<4; j++)if(i!=j)nuc_mut_mat[i][j]*=rate;
    nuc_mut_mat[i][i]=1.-rate*sum;
    if(nuc_mut_mat[i][i]<0){
      printf("Error in %s, rate too large\n", CODE_NAME); exit(8);
    }
  }
}    
