#include "alignments.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "NeedlemanWunsch.h"

static int Check_MSA(int *L, char *file_ali);

char **Read_MSA(int *N, int *L, char ***name_seq, int **selected,
		char *file_ali, float thr){

  if(file_ali==NULL)return(NULL);

  *N=Check_MSA(L, file_ali);
  char **MSA=malloc(*N*sizeof(char *));
  *name_seq= malloc(*N*sizeof(char *));
  int NCHAR=1000, n;
  for(n=0; n<*N; n++){
    MSA[n]=malloc(*L*sizeof(char));
    (*name_seq)[n]=malloc(NCHAR*sizeof(char));
  }

  FILE *file_in=fopen(file_ali, "r");
  char string[10000], *seq=NULL; n=-1;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){
      n++; seq=MSA[n];
      sscanf(string+1, "%s", (*name_seq)[n]);
    }else{
      char *s=string;
      while((*s!='\n')&&(*s!='\0')){*seq=*s; seq++; s++;}
    }
  }
  fclose(file_in);

  *selected=malloc(*N*sizeof(int)); int m=0;
  for(n=0; n<*N; n++){
    int l=0, k; for(k=0; k<*L; k++)if(MSA[n][k]!='-')l++;
    if(l > (*L)*thr){(*selected)[n]=1; m++;}
    else{(*selected)[n]=0;}
  }

  printf("%d sequences with %d letters read in file %s, %d selected\n",
	 *N,*L,file_ali, m);
  return(MSA);
}

int Check_MSA(int *L, char *file_ali)
{
  FILE *file_in=fopen(file_ali, "r");
  if(file_in==NULL){
    printf("ERROR, MSA file %s does not exist\n", file_ali); exit(8);
  }

  printf("Reading MSA file %s\n", file_ali);
  char string[10000]; int n=-1, l=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){
      if(n==0){
	*L=l;
      }else if(l!=*L){
	printf("ERROR, different sequence lengths: 0 %d %d %d\n", *L, n, l);
	exit(8);
      }
      n++; l=0;
    }else{
      char *s=string;
      while((*s!='\n')&&(*s!='\0')){l++; s++;}
    }
  }
  if(n==0){
    *L=l;
  }else if(l!=*L){
    printf("ERROR, different sequence lengths: %d %d\n", *L, l);
    exit(8);
  }
  n++;
  fclose(file_in);
  if(n==0){
    printf("ERROR, no sequence found in MSA file %s\n", file_ali);
    exit(8);
  }
  return(n);
}

int Find_seq(int *ali_seq, int *seq_L, float *seq_id,
	     char *seq_PDB, int L_pdb, char **MSA, int n_seq, int L_ali)
{
  int iseq, i, j; 
  for(i=0; i<L_pdb; i++)ali_seq[i]=-1;
  for(i=0; i<n_seq; i++){
    int l=0; for(j=0; j<L_ali; j++)if(MSA[i][j]!='-')l++; 
    seq_L[i]=l;
  }

  int SHIMAX=10, LMIN=0.75*L_pdb, GMAX=20; char *msa;
  int shift, type, idmax=0, seqopt=0, shiftopt=0, topt=0;
  for(iseq=0; iseq<n_seq; iseq++){
    msa=MSA[iseq];
    type=0; // Only gaps in sequence
    for(shift=0; shift <= SHIMAX; shift ++){
      int id=0, gap=0, j=0;
      for(i=shift; i<L_pdb; i++){
	if(j>=L_ali)break;
	while(msa[j]!=seq_PDB[i]){
	  j++; gap++; if((j>=L_ali)||(gap>GMAX))goto check_0;
	}
	id++; j++;
      }
    check_0:
      if(id >LMIN){
	j=0;
	for(i=shift; i<L_pdb; i++){
	  if(j>=L_ali)break;
	  while(msa[j]!=seq_PDB[i]){j++; if(j>=L_ali)goto end;}
	  ali_seq[i]=j; j++;
	}
	goto end;
      }
      if(id > idmax){idmax=id; seqopt=iseq; shiftopt=shift; topt=0;}
    }
    type=1; // Only gaps in structure
    for(shift=0; shift <= SHIMAX ; shift ++){
      int id=0, gap=0, i=0;
      for(j=shift; j<L_ali; j++){
	if(i>=L_pdb)break;
	while(msa[j]!=seq_PDB[i]){
	  i++; gap++; if((i>=L_pdb)||(gap>GMAX))goto check_1;
	}
	id++; i++;
      }
    check_1:
      if(id > LMIN){
	i=0;
	for(j=shift; j<L_ali; j++){
	  if(j>=L_ali)break;
	  while(msa[j]!=seq_PDB[i]){i++; if(i>=L_pdb)goto end;}
	  ali_seq[i]=j; i++;
	}
	goto end;
      }
      if(id > idmax){idmax=id; seqopt=iseq; shiftopt=shift; topt=type;}
    }
  }
  // Not found: Make complete alignment
  {
    type=-1;
    int VBS=0;  // Verbose:
    int IDE=1;  // Use identity to score alignment
    int GAP=7;  // Gap opening penalty
    char *ali1=malloc((L_pdb+L_ali)*sizeof(char));
    char *ali2=malloc((L_pdb+L_ali)*sizeof(char));
    char *seq=malloc(L_ali*sizeof(char));
    int *label=malloc(L_ali*sizeof(int));
    for(iseq=0; iseq<n_seq; iseq++){
      int kseq=iseq, nali;
      if(iseq==0){kseq=seqopt;}else if(iseq==seqopt){kseq=0;}  
      j=0; msa=MSA[kseq];
      for(i=0; i<L_ali; i++){
	if(msa[i]!='-'){seq[j]=msa[i]; label[j]=i; j++;}
      }
      int al=alignNW(seq_PDB, L_pdb, seq, seq_L[kseq],
		     VBS, IDE, GAP, ali1, ali2, &nali);
      if(al==0){
	printf("WARNING, alignment failed\n"); continue;
      }
      int i1=0, i2=0, id=0, ali=0;
      for(i=0; i<nali; i++){
	if((ali1[i]!='-')&&(ali2[i]!='-')){
	  ali++; if(ali1[i]==ali2[i])id++;
	  ali_seq[i1]=label[i2];
	}
	if(ali1[i]!='-')i1++;
	if(ali2[i]!='-')i2++;
      }
      free(ali1); free(ali2); free(seq); free(label);
      if(id > LMIN)goto end;
      if(id > idmax){idmax=id; seqopt=kseq; topt=-1;}
    }
  }

  printf("WARNING, query sequence not found in MSA\n");
  printf("query sequence:\n");
  for(i=0; i<L_pdb; i++)printf("%c", seq_PDB[i]);
  printf("\n");
  printf("Max. aligned length: %d L_min= %d seq %d shift %d type %d\n",
	 idmax, LMIN, seqopt, shiftopt, topt); 
  for(i=0; i<L_ali; i++)printf("%c",MSA[seqopt][i]);
  printf("\n");
  return(-1);
  // found

 end:
  for(i=0; i<n_seq; i++){
    int id=0, l=L_pdb;
    for(j=0; j<L_pdb; j++){
      if(ali_seq[j]<0){l--; continue;}
      if(MSA[i][ali_seq[j]]==seq_PDB[j]){id++;}
      else if(MSA[i][ali_seq[j]]=='-'){l--;}
    }
    seq_id[i]=(float)id/(float)l;
    l=0; for(j=0; j<L_ali; j++)if(MSA[i][j]!='-')l++; 
    seq_L[i]=l;
  }

  return(iseq);
}
