#include "coord.h"
#include "protein3.h"
#include "read_str_mut.h"
#include "codes.h"
#include "allocate.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static int Get_rank(int *code_AA, char *string);
static int Get_rank2(int *code_AA, char *string);
void f_sort(float d[], int n, int *i_rank);

int Read_name_str(char *STR_MUT_TYPE, char *file_str_mut)
{
  int all=0;
  if(file_str_mut[0]=='\0')return(0);
  char *s=file_str_mut; // is file name mut_DE or mut_RMSD?
  while(*s!='\0'){
    if(strncmp(s, "mut", 3)==0){
      int j=7; 
      if(strncmp(s+4, "DE", 2)==0){strcpy(STR_MUT_TYPE, "DE"); j=7;}
      else if(strncmp(s+4, "RMSD", 4)==0){strcpy(STR_MUT_TYPE, "RMSD"); j=9;}
      else{strcpy(STR_MUT_TYPE, "UNK");
	printf("WARNING, unknown file name %s\n", file_str_mut);
      }
      if(strncmp(s+j, "all", 3)==0){
	all=1; strcat(STR_MUT_TYPE, "2");
      }
      break;
    }
    s++;
  }
  printf("Type of mutational effects: %s allpair= %d\n", STR_MUT_TYPE, all);
  return(all);
}


int Read_str_mut(float **Str_mut, char *file_str_mut,
		 struct protein target, int *res_index, int IWT)
{
  if(file_str_mut[0]=='\0')return(0);
  FILE *file_in=fopen(file_str_mut, "r");
  if(file_in==NULL){
    printf("WARNING, file %s with structural mutations does not exist\n",
	   file_str_mut); return(0);
  }
  char string[9000]; 
  int wrong_thr=10, k=0, wrong=0, L=target.L_PDB;
  short *aseq=target.aa_seq;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(k>=L){
      printf("WARNING, line %s in file %s > number of residues %d\n",
	     string, file_str_mut, L); wrong++; k++; continue;
    }
    int i=res_index[k];
    if(i<0 || i>= target.length){
      printf("WARNING, sequence index of res %d = %d not in [0,%d]\n",
	     k, i, target.length-1); wrong++;
    }else if(string[0]!=AMIN_CODE[aseq[i]]){
      printf("WARNING, amino acid %c%d in file %s different from %c ",
	     string[0], k+1, file_str_mut, AMIN_CODE[aseq[i]]);
      printf("in target sequence\n"); wrong++;
    }
    k++;
  }
  fclose(file_in);

  if(k!=L){
    printf("WARNING, file %s does contains %d lines but %d residues in prot\n",
	   file_str_mut, k, L); wrong+=abs(k-L);
  }
  if(wrong){
    printf("WARNING, there were %d errors in file %s\n", wrong, file_str_mut);
    if(wrong>wrong_thr){printf("str mut not read\n"); return(0);}
  }

  printf("Reading predicted structural effects of mutations in file %s\n",
	 file_str_mut);
  file_in=fopen(file_str_mut, "r");
  int len_amm=target.length;
  int num[len_amm]; for(k=0; k<len_amm; k++)num[k]=0;
  k=0; char dumm[10]; int code_AA[20], ncont, a;
  for(a=0; a<20; a++)code_AA[a]=a;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'){
      if(strncmp(string, "#Mut", 4)==0){
	if(Get_rank(code_AA, string)<0)
	  for(a=0; a<20; a++)code_AA[a]=a;
      }
      continue;
    }
    if(k>=L)goto next;
    int i=res_index[k];
    if(i<0 || i>= len_amm)goto next;
    float *m=Str_mut[i]; num[i]++;
    float s[20];
    sscanf(string,"%s%d %f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",
	   dumm, &ncont, s, s+1, s+2, s+3, s+4, s+5, s+6, s+7, s+8, s+9,
	   s+10, s+11, s+12, s+13, s+14, s+15, s+16, s+17, s+18, s+19);
    for(a=0; a<20; a++)m[code_AA[a]]+=s[a];
  next:
    k++;
  }
  fclose(file_in);

  int Na1=20-1, i_rank[20]; // Set deformation of wild-type
  for(int i=0; i<len_amm; i++){
    float *Str_i=Str_mut[i], Y; int wt=aseq[i];
    if(num[i]>1)for(a=0; a<20; a++)Str_mut[i][a]/=num[i];
    f_sort(Str_i, 20, i_rank);
    Y=0; for(a=0; a<IWT; a++)Y+=Str_i[i_rank[Na1-a]]; Str_i[wt]=Y/IWT;
  }

  return(1);
}

int Read_str_mut_all(float ***Str_mut, char *file_str_mut,
		     struct protein target, int *res_index)
{
  if(file_str_mut[0]=='\0')return(0);
  FILE *file_in=fopen(file_str_mut, "r");
  if(file_in==NULL){
    printf("WARNING, file %s with structural mutations does not exist\n",
	   file_str_mut); return(0);
  }
  char string[9000];
  int wrong_thr=10, k=0, wrong=0, L=target.L_PDB;
  short *aseq=target.aa_seq;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(k>=L){
      printf("WARNING, line %s in file %s > number of residues %d\n",
	     string, file_str_mut, L); wrong++; k++; continue;
    }
    int i=res_index[k];
    if(i<0 || i>= target.length){
      printf("WARNING, sequence index of res %d = %d not in [0,%d]\n",
	     k, i, target.length-1); wrong++;
    }else if(string[0]!=AMIN_CODE[aseq[i]]){
      printf("WARNING, amino acid %c%d in file %s different from %c ",
	     string[0], k+1, file_str_mut, AMIN_CODE[aseq[i]]);
      printf("in target sequence\n"); wrong++;
    }
    k++;
  }
  fclose(file_in);

  if(k!=L){
    printf("WARNING, file %s does contains %d lines but %d residues in prot\n",
	   file_str_mut, k, L); wrong+=abs(k-L);
  }
  if(wrong){
    printf("WARNING, there were %d errors in file %s\n", wrong, file_str_mut);
    if(wrong>wrong_thr){printf("str mut not read\n"); return(0);}
  }

  printf("Reading predicted structural effects of mutations in file %s\n",
	 file_str_mut);
  file_in=fopen(file_str_mut, "r");
  int len_amm=target.length;
  int num[len_amm]; for(k=0; k<len_amm; k++)num[k]=0;
  int code_AA[20], a, b; for(a=0; a<20; a++)code_AA[a]=a;
  k=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'){
      if(strncmp(string, "#pos", 4)==0){
	if(Get_rank2(code_AA, string)<0)
	  for(a=0; a<20; a++)code_AA[a]=a;
      }
      continue;
    }
    if(k>=L)goto next;
    int i=res_index[k];
    if(i<0 || i>= len_amm)goto next;
    float d[20][20];
    char *s=string; int a=0, b=1;
    //printf("%s",string);
    while(*s!='\t' && *s!=' ' && *s!='\n'){s++;} // word 
    while(*s=='\t' || *s==' '){s++;} // space
    while(*s!='\n'){
      sscanf(s,"%f", &(d[a][b]));
      //printf("%d-%d %.3g ",a,b, d[a][b]);
      b++; if(b==20){a++; b=a+1; if(a==19)break;} 
      while(*s!='\t' && *s!=' ' && *s!='\n'){s++;} // word 
      while(*s=='\t' || *s==' '){s++;} // space
    }
    if(i==0)printf("last read: Df[%d][%d]\n", a-1, b-1);
    float **Dfi=Str_mut[i]; num[i]++;
    for(a=0; a<20; a++){
      int aa=code_AA[a];
      float *Dfia=Dfi[aa];
      for(b=a+1; b<20; b++){
	int bb=code_AA[b];
	Dfia[bb]=d[a][b];
	Dfi[bb][aa]=-d[a][b];
      }
    }
  next:
    k++;
  }
  fclose(file_in);
  for(int i=0; i<len_amm; i++){
    if(num[i]>0){
      float **Dfi=Str_mut[i];
      for(a=0; a<20; a++){
	Dfi[a][a]=0;
	for(b=a+1; b<20; b++){
	  Dfi[a][b]/=num[i];
	  Dfi[b][a]=-Dfi[a][b];
	}
      }
    }
  }
  return(1);
}

int Mean_DDG(float **exp_ia, float ***Df, float **P_ia, int L)
{
  int i, a, b;
  for(i=0; i<L; i++){
    float  *ei=exp_ia[i], *P=P_ia[i], **Dfi=Df[i];
    for(b=0; b<20; b++){
      double sum=0;
      for(a=0; a<20; a++){sum+=P[a]*Dfi[a][b];}
      ei[b]=sum;
    }
  }
  return(0);
}

void WT_DDG(float **exp_ia, float ***Df, short *aa_seq, int L)
{
  for(int i=0; i<L; i++){
    int a=aa_seq[i];
    float *ei=exp_ia[i], *Dfi=Df[i][a], min=1000;
    for(int b=0; b<20; b++){
      if(a!=b){ei[b]=Dfi[b]; if(ei[b]<min)min=ei[b];}
    }
    ei[a]=min;
  }
}

int Get_rank(int *code_AA, char *string)
{
  char *AA[20]; int a;
  for(a=0; a<20; a++){
    AA[a]=malloc(10*sizeof(char)); strcpy(AA[a], "");
  }
  char dumm[10], dumm2[10]; 
  sscanf(string,"%s%s %s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
	 dumm, dumm2,
	 AA[0], AA[1], AA[2], AA[3], AA[4], AA[5], AA[6], AA[7], AA[8], AA[9],
	 AA[10], AA[11], AA[12], AA[13], AA[14],
	 AA[15], AA[16], AA[17], AA[18], AA[19]);
  printf("Amino acid order: ");
  for(a=0; a<20; a++)printf("%s ",AA[a]);
  printf("\n");
  for(a=0; a<20; a++){
    int k; for(k=0; k<20; k++)if(AA[a][0]==AMIN_CODE[k])break;
    if(k==20){
      printf("WARNING, amino acid %d = %s not found in string %s\n",
	     a, AA[a], string); return(-1);
    }
    code_AA[a]=k;
  }
  return(0);
}

int Get_rank2(int *code_AA, char *string)
{
  char *AA[20]; int a; 
  for(a=0; a<20; a++){
    AA[a]=malloc(10*sizeof(char)); strcpy(AA[a], "");
  }
  char dumm[10]; 
  sscanf(string,"%s %s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
	 dumm,
	 AA[1], AA[2], AA[3], AA[4], AA[5], AA[6], AA[7], AA[8], AA[9],
	 AA[10], AA[11], AA[12], AA[13], AA[14],
	 AA[15], AA[16], AA[17], AA[18], AA[19]);
  // AE	AQ AD AN AL AG AK AS AV	AR AT AP AI AM AF AY AC AW AH
  AA[0][1]=AA[1][0]; AA[0][2]='\0';
  //for(a=1; a<20; a++){AA[a][0]=AA[a][1]; AA[a][1]='\0';}
  printf("Amino acid order: ");
  for(a=0; a<20; a++)printf("%c ",AA[a][1]);
  printf("\n");
  for(a=0; a<20; a++){
    int k; for(k=0; k<20; k++)if(AA[a][1]==AMIN_CODE[k])break;
    if(k==20){
      printf("WARNING, amino acid %d = %c not found in string %s\n",
	     a, AA[a][1], AMIN_CODE); return(-1);
    }
    code_AA[a]=k;
  }
  return(0);
}

void f_sort(float d[], int n, int *i_rank)
{
  // Sort the vector d from large to small.
  // Returns i_rank[i]= index of object at rank i
  int k,j,i, jmax, not_ranked[n];
  for (i=0; i<n; i++)not_ranked[i]=1;
  for (k=0; k<n; k++) {
    for(i=0; i<n; i++)if(not_ranked[i])break;
    float d_max=d[jmax=i];
    for (j=i+1;j<n; j++){
      if(not_ranked[j]&&(d[j] >= d_max))d_max=d[jmax=j];
    }
    i_rank[k]=jmax; not_ranked[jmax]=0;
  }
}

