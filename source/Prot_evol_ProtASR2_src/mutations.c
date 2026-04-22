#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "coord.h"
#include "mut_del.h"
#include "codes.h"

static int Assign_mutation(struct mutation *mut, char *string);
static int Assign_deletion(struct deletion *del, char *string);
static int Count_words(char *string);
int Find_site(char amm, int rmut, char chain, struct residue *res, int N_res);

int Read_mut(struct mutation **mutations, char *FILE_MUT,
	     struct deletion **deletions, int *Ndel)
{
  // Mutation can be input as file or as a single string, for instance S100A
  FILE *file_in=fopen(FILE_MUT, "r");
  if(file_in == NULL){
    printf("WARNING, mutation file %s does not exist\n", FILE_MUT);
    return(0);
  }

  int nd=0, nm=0;
  char string[200];
  printf("Reading mutation file %s\n", FILE_MUT);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "DEL", 3)==0){nd++;}
    else{nm++;}
  }
  fclose(file_in);

  if(nm)*mutations=malloc(nm*sizeof(struct mutation));
  if(nd)*deletions=malloc(nd*sizeof(struct mutation));
  file_in=fopen(FILE_MUT, "r"); nd=0; nm=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "DEL", 3)==0){
      if(Assign_deletion((*deletions)+nd, string))nd++;
    }else{
      if(Assign_mutation((*mutations)+nm, string))nm++;
    }
  }
  fclose(file_in);
  *Ndel=nd;
  return(nm);
}

int Assign_mutation(struct mutation *mut, char *string)
{
  int m=Count_words(string);
  if(m==0){
    printf("WARNING, mutation string %s does not contain mutations\n",
	   string); return(0);
  }
  printf("Mutant %s%d mutations\n", string, m);

  mut->amm1_mut=malloc(m*sizeof(char));
  mut->amm2_mut=malloc(m*sizeof(char));
  mut->chain=malloc(m*sizeof(char));
  mut->imut=malloc(m*sizeof(int));
  mut->rmut=malloc(m*sizeof(int));
  mut->nmut=m; m=0;
  mut->interface=0;

  char *s=string, *s0;
  // Format: A22G
  while(1){
    while((*s==' ')||(*s=='\t'))s++;
    if((*s=='\n')||(*s=='\r'))break;
    // Assign amm1
    mut->amm1_mut[m]=*s; s++; s0=s;
    // Assign amm2 and delete it
    while((*s!=' ')&&(*s!='\0')&&(*s!='\n')&&(*s!='\r'))s++;
    s--; mut->amm2_mut[m]=*s; *s=' ';
    // Assign res
    sscanf(s0, "%d", &(mut->rmut[m]));
    // Assign chain
    mut->chain[m]='*';
    // Assign name
    if(m)strcat(mut->name, "_");
    char nam[20];
    sprintf(nam, "%c%d%c", mut->amm1_mut[m], mut->rmut[m], mut->amm2_mut[m]);
    strcat(mut->name, nam);
    m++;
    // Next word
    while((*s!=' ')&&(*s!='\n')&&(*s!='\t'))s++;
  }
  printf("Mutant %s %d mutations m=%d\n", mut->name, mut->nmut, m);
  return(1);
}

int Assign_deletion(struct deletion *del, char *string)
{
  int m=Count_words(string)-1;
  if(m==0){
    printf("WARNING, deletion string %s does not contain deletions\n",
	   string); return(0);
  }

  del->aa1=malloc(m*sizeof(char));
  del->aa2=malloc(m*sizeof(char));
  del->chain=malloc(m*sizeof(char));
  //del->i1=malloc(m*sizeof(int));
  //del->i2=malloc(m*sizeof(int));
  del->res1=malloc(m*sizeof(int));
  del->res2=malloc(m*sizeof(int));
  del->interface=0;
  int Ldel=0;

  char *s=string, ss[10]; m=0;
  // Skip the record DEL
  while((*s!=' ')&&(*s!='\n')&&(*s!='\t'))s++;
  while((*s==' ')||(*s=='\t'))s++;
  while(1){
    del->aa1[m]=*s; s++;
    while((*s!='-')&&(*s!='\n')&&(*s!='\r'))s++;
    if(*s=='-')*s=' ';
    sscanf(ss, "%d", &(del->res1[m])); s++;
    del->aa2[m]=*s; s++;
    while((*s!=' ')&&(*s!='\n')&&(*s!='\r'))s++;
    sscanf(ss, "%d", &(del->res2[m]));
    Ldel+=(del->res2[m]-del->res1[m]);
    del->chain[m]='*';
    while((*s==' ')||(*s=='\t'))s++;
    /*if((*s!='\n')&&((*(s+1)==' ')||(*(s+1)=='\t')||(*(s+1)=='\n'))){
      del->chain[m]=*s; s++; while((*s==' ')||(*s=='\t'))s++;
      }*/
    if(m)strcat(del->name, "_");
    char nam[20];
    sprintf(nam, "%c%d-%c%d",
	    del->aa1[m], del->res1[m], del->aa2[m], del->res2[m]);
    strcat(del->name, nam);
    m++; if((*s=='\n')||(*s=='\r'))break;
  }
  del->ndel=m;
  del->Ldel=Ldel;
  return(1);
}

int Count_words(char *string){
  int m=0;
  char *s=string;
  printf("Counting words ");
  while(1){
    while((*s==' ')||(*s=='\t'))s++;
    if((*s=='\n')||(*s=='\0')||(*s=='\r'))break;
    printf("%d: %c ", m+1, *s);
    while((*s!=' ')&&(*s!='\n')&&(*s!='\t'))s++;
    m++;
  }
  printf("\n");
  return(m);
}

int Construct_mut(short *aa_seq, struct mutation *mut,
		  struct residue *res, int L)
{	
  if(mut->nmut==0)return(-1);
  printf("Mutation %s\n", mut->name);
  int k, imut;
  for(k=0; k<mut->nmut; k++){
    imut=Find_site(mut->amm1_mut[k],mut->rmut[k],mut->chain[k], res, L);
    if(imut<0)return(-1);
    aa_seq[imut]=Code_AA(mut->amm2_mut[k]);
  }
  return(0);
}

int Construct_del(short *aa_seq, struct deletion *del,
		  struct residue *res, int L)
{	
  if(del->ndel==0)return(-1);
  int k, i1, i2, j;
  for(k=0; k<del->ndel; k++){
    i1=Find_site(del->aa1[k], del->res1[k], del->chain[k], res, L);
    i2=Find_site(del->aa2[k], del->res2[k], del->chain[k], res, L);
    if((i1<0)||(i2<0))return(-1);
    for(j=i1; j<=i2; j++)aa_seq[j]=-1; //20;
  }
  return(0);
}

int Find_site(char amm, int rmut, char chain, struct residue *res, int N_res)
{
  int k, l, resnum;

  for(k=0; k<N_res; k++){
    sscanf(res[k].pdbres, "%d", &resnum);
    if((resnum==rmut)&&((chain=='*')||(chain==res[k].chain))){
      if(res[k].amm==amm){
	return(k);
      }else{
	printf("WARNING, mutation from %c%d but A[%d]=%c, order in PDB: %d\n",
	       amm, rmut, rmut, res[k].amm, k);
	printf("Sequence:\n");
	for(l=0; l<=k; l++)printf("%c", res[l].amm);
	printf("  %d-%d\n", resnum, rmut);
	for(l=k+1; l<N_res; l++)printf("%c", res[l].amm);
	printf("  %s-%s\n",res[rmut+1].pdbres, res[N_res-1].pdbres);
	//return(-1);
	return(k);
      }
    }
  }
  printf("ERROR, mutated residue %d chain %c is not present\n",
	 rmut, chain);
  return(-1);
}


