#include "codes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_ali.h"

static int Align_nogap(int *ali_col, char *aa_seq, int L_PDB,
		       char *ali, int L_ali, int nalimin);

int Read_ali(char ***msa_seq, int *L_ali, char ***name_seq, char *file_ali)
{
  if(file_ali==NULL){
    printf("WARNING, there is no name of FASTA file with alignment\n");
    return(0);
  }
  FILE *file_in=fopen(file_ali, "r");
  if(file_in==NULL){
    printf("WARNING, FASTA file with alignment %s does not exist\n",
	   file_ali); return(0);
  }

  int n=0, l=0; char string[10000], *s;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){n++; continue;}
    else if(n==1){
      s=string; while((*s!='\n')&&(*s!='\r')&&(*s!='\0')){l++; s++;}
    }
  }
  fclose(file_in);
  printf("%d sequences of length %d read in alignment %s\n",n,l,file_ali);

  // Allocate
  int N_seq=n; char *ali_n; *L_ali=l;
  *msa_seq=malloc(N_seq*sizeof(char **));
  *name_seq=malloc(N_seq*sizeof(char *));
  for(n=0; n<N_seq; n++){
    (*msa_seq)[n]=malloc(*L_ali*sizeof(char));
    (*name_seq)[n]=malloc(40*sizeof(char));
  }
  
  file_in=fopen(file_ali, "r"); n=-1;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){
      if((n>=0)&&(l!=(*L_ali))){
	printf("ERROR, variable number of columns in alignment\n");
	printf("Sequence %d has %d symbols instead of %d\n",n,l,*L_ali);
	printf("Just read: %s\n", string);
	for(n=0; n<N_seq; n++)free((*msa_seq)[n]);
	free(msa_seq); *L_ali=0;
	return(0);
      }
      n++; ali_n=(*msa_seq)[n]; l=0;
      sscanf(string+1, "%s", (*name_seq)[n]);
      continue;
    }
    s=string; 
    while((*s!='\n')&&(*s!='\r')&&(*s!='\0')){
      if(l>=*L_ali){
	printf("ERROR reading alignment, sequence %d has %d > %d symbols\n",
	       n, l, *L_ali); printf("Just read: %s\n", string);
	for(n=0; n<N_seq; n++){free((*msa_seq)[n]);}
	free(msa_seq); *L_ali=0;
	return(0);
      }
      *ali_n=*s; s++; l++; ali_n++;
    }
  }
  fclose(file_in);
  printf("Read %d sequences with L_ali %d\n", N_seq,*L_ali);

  // Remove columns with all gaps
  int allgaps=0, i, inew=0;
  for(i=0; i<(*L_ali); i++){
    int res=0;
    for(n=0; n<N_seq; n++){if((*msa_seq)[n][i]!='-'){res=1; break;}}
    if(res==0){allgaps++; continue;}
    if(inew!=i){
      for(n=0; n<N_seq; n++){(*msa_seq)[n][inew]=(*msa_seq)[n][i];}
    }
    inew++;
  }
  if(allgaps){for(n=0; n<N_seq; n++){(*msa_seq)[n][inew]='\0';}}
  (*L_ali)-=allgaps;
  if(inew!=(*L_ali)){
    printf("ERROR removing gaps from MSA, allgaps=%d but Lali=%d != %d\n",
	   allgaps, inew, *L_ali); exit(8);
  }
  if(allgaps){
    // Print MSA without gaps
    printf("Printing MSA without columns with all gaps\n");
    FILE *file_out=fopen(file_ali, "w");
    for(n=0; n<N_seq; n++){
      fprintf(file_out, ">%s\n", (*name_seq)[n]);
      ali_n=(*msa_seq)[n];
      for(i=0; i<(*L_ali); i++){
	fprintf(file_out, "%c",*ali_n); ali_n++;
      }
      fprintf(file_out, "\n");
    }
    fclose(file_out);
  }
  return(N_seq);		  
}

int Align_ali_PDB(short **ali_seq, float *seq_id, float *Seq_id_ave, 
		  char **msa_seq, char **name_seq, int N_seq, int L_ali,
		  short *i_seq, int L_seq_PDB, int L_str_PDB,
		  char *pdbname, int *ali_PDB)
{
  int npdb=-1, n, i;
  // Align to PDB sequence
  int ali_col[L_seq_PDB], gapmax=10;
  char aa_seq[L_seq_PDB]; int nalimax=0, nalimin=L_str_PDB-gapmax;
  for(i=0; i<L_seq_PDB; i++)aa_seq[i]=AMIN_CODE[i_seq[i]];
  for(n=0; n<N_seq; n++){
    int nali=
      Align_nogap(ali_col, aa_seq, L_seq_PDB, msa_seq[n], L_ali, nalimin);
    if(nali > nalimax)nalimax=nali;
    if(nali >= nalimin){break;}
  }
  if(n==N_seq){
    printf("WARNING, the PDB sequence could not be aligned to any seq\n");
    printf("Maximum alignment length: %d\n", nalimax);
    return(-1);
  }
  npdb=n;
  for(i=0; i<L_ali; i++)ali_PDB[i]=-1;
  printf("PDB sequence aligned to sequence %d of %d nali= %d\n",
	 npdb,N_seq,nalimax);
  printf("Seq PDB: ");
  for(i=0; i<L_seq_PDB; i++)printf("%c",aa_seq[i]);
  printf("\n");
  printf("Seq MSA: ");
  for(i=0; i<L_seq_PDB; i++){
    int c=ali_col[i];
    if(c<0){printf("-");}
    else{
      ali_PDB[c]=i;
      printf("%c", msa_seq[npdb][c]);
    }
  }
  printf("\n");
 
  for(i=0; i<L_seq_PDB; i++){
    int c=ali_col[i];
    for(n=0; n<N_seq; n++){
      char *ali=msa_seq[n];
      if((c>=0)&&(ali[c]!='-')&&(ali[c]!='.')){
	int a=Code_AA(ali[c]);
	if((a<0)||(a>19)){printf("ERROR wrong AA %c\n", ali[c]); a=-1;}
	ali_seq[n][i]=a;
	//printf("%c", ali[c]);
      }else{
	ali_seq[n][i]=-1;
      }
    }
  }


  for(int j=0; j<L_seq_PDB; j++){
    if(ali_seq[0][j]>=0){printf("%c",AMIN_CODE[ali_seq[0][j]]);}
    else{printf("-");}
  }
  printf("\n");
  
  // Mean sequence identity with PDB
  *Seq_id_ave=0;
  for(n=0; n<N_seq; n++){
    if(n==npdb){seq_id[n]=1; continue;}
    int s=0;
    for(i=0; i<L_ali; i++){
      if(msa_seq[n][i]==msa_seq[npdb][i] && msa_seq[n][i]!='-')s++;
    } 
    seq_id[n]=s/(float)L_seq_PDB;
    *Seq_id_ave+=s;
  }
  *Seq_id_ave/=(L_seq_PDB*(N_seq-1));
  return(npdb);
}

int Align_nogap(int *ali_col, char *aa_seq, int L_PDB,
		char *ali, int L_ali, int nalimin)
{
  char *a1, *a2; int L1, i1, L2, i2;
  int DBG=0;
  int nali=0,  nali2=0, gap=0, start=0;

  a1=ali; L1=L_ali; gap=0;
  a2=aa_seq; L2=L_PDB; i2=0; nali2=0;
  for(i1=0; i1<L_PDB; i1++)ali_col[i1]=-1;
  for(i1=0; i1<L1; i1++){
    if(a1[i1]!=a2[i2]){
      if(a1[i1]=='-')continue;
      if(nali2==1){ali_col[i2-1]=-1; nali2=0;}
      if(a2[i2]==a1[i1+1] && gap==0){gap=1; continue;}
      if(DBG && nali>100)printf("%c %c %d %d\n",a2[i2],a1[i1],i2,nali2);
      while((a2[i2]!=a1[i1])&&(i2<L2))i2++;
      if(i2==L2)break;
    }
    if(a1[i1]==a2[i2]){
      ali_col[i2]=i1; gap=0;
      int j1=i1-1;
      if(i2 && ali_col[i2-1]<j1 && a2[i2-1]==a1[j1]){
	ali_col[i2-1]=j1;
	for(int j2=i2-2; j2>=0; j2--)
	  if(ali_col[j2]==j1){ali_col[j2]=-1; break;}
      }
      i2++; nali2++;
    }else if(gap){
      break;
    }
  }
  if(DBG)printf("nali2= %d\n", nali2);
  if(nali2>=nalimin)return(nali2);

  a1=aa_seq; L1=L_PDB;
  a2=ali; L2=L_ali; i2=0;
  if(DBG){
    for(i1=0; i1<L_ali; i1++)printf("%c", ali[i1]);
    printf("\n");
  }
  for(i1=0; i1<L_PDB; i1++)ali_col[i1]=-1;
  for(i1=0; i1<L1; i1++){
    while(a2[i2]=='-' && i2<L2)i2++;
    if(i2==L2)break;
    if(a1[i1]!=a2[i2]){
      if(start==0)continue; // Initial part of PDB sequence, it may be 
                            // disordered and some programs do not read it
      if(nali==1){ali_col[i1-1]=-1; nali=0;}
      if(a2[i2]==a1[i1+1] && gap==0){gap=1; continue;}
      if(DBG && nali>100)printf("%c %c %d %d\n", a1[i1],a2[i2],i1,nali);
      while((a2[i2]!=a1[i1])&&(i2<L2))i2++;
      if(i2==L2)break;
    }
    if(start==0)start=1;
    if(a1[i1]==a2[i2]){
      ali_col[i1]=i2; i2++; nali++; if(gap)gap=0;
    }else if(gap){
      break;
    }
  }
  if(DBG)printf("nali= %d\n", nali);
  if(nali>=nalimin)return(nali);
  
  if(nali>nali2){return(nali);}
  else{return(nali2);}
}

