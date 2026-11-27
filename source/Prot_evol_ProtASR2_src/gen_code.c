#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "codes.h"
#include "input.h"
#include "random3.h"
#include "gen_code.h"
#include "externals.h"
#define CODE_NAME "gen_code.c"

// Codons
char coded_aa[64]="FFLLLLLLIIIMVVVVSSSSPPPPTTTTAAAAYYHHQQNNKKDDEECCWRRRRSSRRGGGG***";
char *codon[64]={"TTT","TTC","TTA","TTG","CTT","CTC","CTA","CTG","ATT","ATC","ATA","ATG","GTT","GTC","GTA","GTG","TCT","TCC","TCA","TCG","CCT","CCC","CCA","CCG","ACT","ACC","ACA","ACG","GCT","GCC","GCA","GCG","TAT","TAC","CAT","CAC","CAA","CAG","AAT","AAC","AAA","AAG","GAT","GAC","GAA","GAG","TGT","TGC","TGG","CGT","CGC","CGA","CGG","AGT","AGC","AGA","AGG","GGT","GGC","GGA","GGG","TAA","TAG","TGA"};



// Internal
void Extract_triplet(char *dna_seq, short *aa_seq, int i,
		     char *coded_aa, char **codon);


/* Genetic code */
void Read_code(char *FILE_CODE, char *coded_aa, char **codon){
  FILE *file_in;
  int i=0, j; char string[100];

  file_in=Open_file_r(FILE_CODE, CODE_NAME, "gen.code");
  printf("Reading %s\n", FILE_CODE);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    for(j=0; j<3; j++)codon[i][j]=string[j];
    coded_aa[i]=string[4]; i++;
  }
  printf("%d codons read\n", i);
  fclose(file_in);
}

char *Extract_dna(int *len_dna, int len_amm, short *aa_seq,
		  char **codon, char *coded_aa)
{
  char *dna_seq=malloc(3*len_amm*sizeof(char)); int i;
  (*len_dna)=3*len_amm;

  for(i=0; i<len_amm; i++){
    Extract_triplet(dna_seq, aa_seq, i, coded_aa, codon);
  }
  if(Compare_amm_dna(dna_seq,*len_dna,aa_seq,len_amm,codon,coded_aa)==0)
    exit(8);
  return(dna_seq);
}

void Extract_triplet(char *dna_seq, short *aa_seq, int i,
		     char *coded_aa, char **codon)
{
  char amm=AMIN_CODE[aa_seq[i]];
  int jcod[64], n=0, j;
  for(j=0; j<64; j++)if(coded_aa[j]==amm){jcod[n]=j; n++;}
  int iran=(RandomFloating()*n); j=jcod[iran];
  int ki=3*i;
  dna_seq[ki]=codon[j][0];
  dna_seq[ki+1]=codon[j][1];
  dna_seq[ki+2]=codon[j][2];
}

int Translate_aa(char **aa_trans, char *dna_seq,
		 int len_dna, char **codon, char *coded_aa)
{
  // Test if aa_seq and translation of dna_seq coincide
  int len_amm=len_dna/3, L;
  *aa_trans= malloc(len_amm*sizeof(char));
  L=Translate_char(*aa_trans, len_amm, dna_seq, len_dna, codon, coded_aa);
  if(L)len_amm=L;
  return(len_amm);
}

int Translate_char(char *aa_trans, int len_amm, char *dna_seq, int len_dna,
		    char **codon, char *coded_aa)
{
  for(int i=0; i<len_amm; i++){
    int j=3*i;
    int aa=Coded_aa(&(dna_seq[j]), codon, coded_aa);
    if(aa<0){
      char ter[4];
      ter[0]=dna_seq[j]; ter[1]=dna_seq[j+1]; ter[2]=dna_seq[j+2];ter[3]='\0';
      printf("WARNING Termination codon %s found at res %d instead of %d\n",
	     ter, i, len_amm);
      return(i);
    }
    aa_trans[i]=AMIN_CODE[aa];
  }
  return(0);
}

int Compare_amm(char *aa1, int L1, char *aa2, int L2)
{
  // Seq. 1 must be included in seq. 2
  int TOL=3, nmut=0, ini=-1, end=-1;
  for(int i=0; i<L1; i++){
    if(i>=L2)break;
    if(aa1[i]!=aa2[i]){
      //printf("WARNING, different amino acid sequences, i=%d %c %c\n",
      //   i, aa1[i],aa2[i]);
      if(ini==-1)ini=i;
      end=i; nmut++;
    }
  }
  if(nmut){
    printf("%d different residues ", nmut);
    printf("first mismatch: i=%d %c %c ",ini, aa1[ini],aa2[ini]);
    printf("last mismatch: i=%d %c %c\n",end, aa1[end],aa2[end]);
  }
  if(nmut<TOL)return(1);
  return(0);
}


int  Compare_amm_dna(char *dna_seq, int len_dna, short *aa_seq, int len_amm,
		     char **codon, char *coded_aa)
{
  // Test if aa_seq and translation of dna_seq coincide
  short *aa_test=malloc(len_amm*sizeof(short));
  int i, nmut=0;
  Translate_new(dna_seq, aa_test, len_amm, codon, coded_aa);
  for(i=0; i<len_amm; i++){
    if(aa_test[i]!=aa_seq[i]){
      printf("WARNING, different amino acid sequences, i=%d %c %c",
	     i, AMIN_CODE[aa_seq[i]],AMIN_CODE[aa_test[i]] );
      printf("  Extracting triplet\n");
      Extract_triplet(dna_seq, aa_seq, i, coded_aa, codon);
      nmut++;
    }
  }
  free(aa_test);
  return(1);
}

int Translate_new(char *dna_seq, short *aa_seq, int length,
		  char **codon, char *coded_aa)
{
  short i, j=0; char ter[4];
  for(i=0; i<length; i++){
    j=3*i;
    aa_seq[i]=Coded_aa(&(dna_seq[j]), codon, coded_aa);
    //aa_seq[i]=Triplet_code(&(dna_seq[j]));
    if(i==(i/60)*60)printf("\n# ");
    printf("%c", AMIN_CODE[aa_seq[i]]);
    if(aa_seq[i]<0){
      ter[0]=dna_seq[j]; ter[1]=dna_seq[j+1]; ter[2]=dna_seq[j+2];ter[3]='\0';
      printf("\nError! Termination codon %s found at %d\n", ter, i);
      exit(8);
    }
  }
  printf("\n# ");
  return(0);
}

int Coded_aa(char *i_codon, char **codon, char *coded_aa){
  int i;
  for(i=0; i<64; i++){
    char *cod=codon[i];
    if((i_codon[0]==cod[0])&&
       (i_codon[1]==cod[1])&&
       (i_codon[2]==cod[2])){
      if(coded_aa[i]=='*')return(-1);
      return(Code_AA(coded_aa[i]));
    }
  }
  printf("ERROR, codon %c%c%c does not exist\n",
	 i_codon[0], i_codon[1], i_codon[2]);
  exit(8);
}

int Codon_num(char *cod, char **codon){
  int i_codon;
  for(i_codon=0; i_codon<64; i_codon++)
    if(strncmp(cod, codon[i_codon], 3)==0)return(i_codon);
  printf("Warning, codon %c%c%c not found\n", cod[0], cod[1], cod[2]);
  return(-1);
}

int Code_codon(char *i_codon, char **codon){
  int i;
  for(i=0; i<64; i++){
    char *cod=codon[i];
    if((i_codon[0]==cod[0])&&
       (i_codon[1]==cod[1])&&
       (i_codon[2]==cod[2]))return(i);
  }
  printf("ERROR, codon %c%c%c does not exist\n",
	 i_codon[0], i_codon[1], i_codon[2]);
  exit(8);
}

float Weight_codon(char *c, float *f){
  return(f[Code_nuc(c[0])]*f[Code_nuc(c[1])]*f[Code_nuc(c[2])]);
}

float Weight_codon_CpG(char *c, float *f){
  // The mutation rate at CG dinucleotides is enhanced by f[5],
  // possibly giving raise to TG or CA dinucleotides

  float w=f[Code_nuc(c[0])]*f[Code_nuc(c[1])]*f[Code_nuc(c[2])];
  if((strncmp(c, "CG",2)==0)||(strncmp(c+1, "CG",2)==0)){
    //return(w*(1-2*f[5]));
    return(w*f[7]);
  }else if((strncmp(c,  "CA",2)==0)||(strncmp(c,  "TG",2)==0)){
    //return(w+f[2]*f[3]*f[5]*f[Code_nuc(c[2])]);
    return(w+f[8]*f[Code_nuc(c[2])]);
  }else if((strncmp(c+1,"CA",2)==0)||(strncmp(c+1,"TG",2)==0)){
    //return(w+f[2]*f[3]*f[5]*f[Code_nuc(c[0])]);
    return(w+f[8]*f[Code_nuc(c[0])]);
  }
  return(w);
}


int Sum_aa(float *num_aa, char *file_ali){
  if(file_ali==NULL){
    printf("WARNING, there is no name of FASTA file with protein sequences\n");
    return(0);
  }
  FILE *file_in=fopen(file_ali, "r");
  if(file_in==NULL){
    printf("WARNING, FASTA file with protein sequences %s does not exist\n",
	   file_ali); return(0);
  }
  int L_MAX=50000; // Maximum number of amino acids in a sequence
  char string[L_MAX], *s; int n=0, l=0, a;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){ n++; continue;}
    s=string;
    while(*s!='\n'){
      a=Code_AA(*s); s++;
      if((a>=0)&&(a<20)){num_aa[a]++; l++;}
    }
  }
  fclose(file_in);
  printf("%d sequences with on average %.0f residues each read in file %s\n",
	 n,l/(float)n,file_ali);
  for(a=0; a<20; a++)printf("%c %.3f  ", Amin_code(a), num_aa[a]);
  printf("\n");
  return(1);
}

void Normalize_distr(float *P, int N){
  double sum=0; int i;
  for(i=0; i<N; i++)sum+=P[i];
  for(i=0; i<N; i++)P[i]/=sum;
}
