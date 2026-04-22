#include "allocate.h"
#include "protein3.h"
#include "coord.h"
#include "read_pdb.h"
#include "contact_matrix.h"
#include "Sec_str_all.h"
#include "codes.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

static void GetPdbId(char *pdb_file_in, char *pdbid);

// Contact matrices
float cont_thr_a=8, cont_thr_a2;
float cont_thr_b=8, cont_thr_b2;
float cont_thr_c=4.5, cont_thr_c2;
int init_map=0;
static int Contact_alpha(struct residue res_i, struct residue res_j);
static int Contact_beta (struct residue res_i, struct residue res_j);
static int Contact(struct residue res_i, struct residue res_j);
static atom *Find_atom(struct residue seq, char code[2]);
static short **Store_map(int N_res, struct residue *seq, short *num_cont,
			 int N_cont,struct contact *Contact_list, int l_cont);
void Get_seqres(char *string, char *seqres, int *n_seqres);
int Cluster_chains(char **seq, int *nchain_out, int *res_index,
		   int **ini_seq, int **n_seq,
		   struct residue *res, char *file,
		   char *chain_to_read, int nchain, int nres);

// Routines:
int Get_pdb(struct protein *prot,
	    struct residue **res, struct res_short **res_short,
	    char *file, char *chain, int **res_index)
{
  int natoms, natoms2, ANISOU, nmr, i;

  prot->length=0;
  int nres=Read_coord(file, &nmr, NULL, NULL, chain, &ANISOU, &natoms);
 
  if(nres<=0){
    printf("WARNING, pdb file %s, no residue found in chain %s\n",
	   file, chain); return(0);
  }

  *res=malloc(nres*sizeof(struct residue));
  atom *atom_read=malloc(natoms*sizeof(atom));
  int nres2=Read_coord(file, &nmr, *res, atom_read, chain, &ANISOU, &natoms2);

  if((nres2!=nres)||(natoms2!=natoms)){
    printf("ERROR, inconsistent number of residues/atoms in %s\n", file);
    printf("Before: %d %d Now: %d %d\n", nres, natoms, nres2, natoms2);
    exit(8);
  }

  // File name
  char pdbid[200];
  GetPdbId(file, pdbid); if(chain[0]==' ')chain[0]='_';
  printf("PDB %s chain %s nres=%d\n",pdbid,chain,nres);
  sprintf(prot->name, "%s%s", pdbid, chain);

  // Count chains
  for(i=0; i<sizeof(chain); i++){if(chain[i]=='\0')break;}
  int nchain=i;
  if(nchain==0){
    printf("WARNING, 0 chains found in %s setting to 1\n", chain);
    nchain=1;
  }
  printf("prot %s has %d chain conformations %s\n",
	 prot->name, nchain, chain);

  // Cluster chains according to seqres. Output: residue indexes
  // seq contains the nres1 SEQRES redisues
  // res_index contains the correspondence btw nres PDB residues and SEQRES
  char *seq=NULL; int *ini_seq=NULL, *n_seq=NULL;
  *res_index=malloc(nres*sizeof(int));
  int nc, nres1=Cluster_chains(&seq, &nc,*res_index, &ini_seq, &n_seq,
			       *res, file, chain, nchain, nres);
  printf("prot %s has %d over %d chains different by sequence ",
	 prot->name, nc, nchain);
  printf("and %d over %d different residues\n", nres1, nres);

  // Set chains
  prot->nchain=nc;
  prot->ini_chain=malloc(nc*sizeof(int));
  prot->len_chain=malloc(nc*sizeof(int));
  for(int ic=0; ic<nc; ic++){
    prot->ini_chain[ic]=ini_seq[ic];
    prot->len_chain[ic]=n_seq[ic];
  }

  // Representative sequence (one per cluster of identical chains)
  prot->L_PDB=nres;
  prot->length=nres1;
  short *aa_seq=malloc(nres1*sizeof(short));
  (*res_short)=malloc(nres1*sizeof(res_short));
  prot->aa_seq=aa_seq;
  for(i=0; i<nres1; i++){
    aa_seq[i]=Code_AA(seq[i]);
    (*res_short)[i].seq=seq[i];
    (*res_short)[i].aa=aa_seq[i];
    (*res_short)[i].ncont=0;
    sprintf((*res_short)[i].pdbres, "  -1");
  }

  // Residue index
  for(i=0; i<nres; i++){
    int k=(*res_index)[i];
    if(k<0 || k>=nres1){
      printf("ERROR, index %d of res %d not allowed, max= %d\n",
	     k, i, nres1-1); exit(8);
    }
    sprintf((*res_short)[k].pdbres, (*res)[i].pdbres); 
  }

  // Contact matrix and number of contacts
  int n_cont;
  short **contact=Contact_matrix(*res, nres, &n_cont, 'c', IJ_MIN);
  int **Cont_map=Allocate_mat2_i(nres1, nres1);
  prot->n_cont=0;
  for(i=0; i<nres; i++){
    int k=(*res_index)[i];
    short *j=contact[i];
    while(*j>=0){
      int kj=(*res_index)[*j]; j++;
      if(kj>=0){
	(*res_short)[k].ncont++; (*res_short)[kj].ncont++; 
	Cont_map[k][kj]=1; Cont_map[kj][k]=1; prot->n_cont++;
      }
    }
  }
  prot->contact=malloc(nres1*sizeof(short *));
  for(i=0; i<nres1; i++){
    prot->contact[i]=malloc(nres1*sizeof(short));
    short *c=prot->contact[i];
    for(int j=i+1; j<nres1; j++)if(Cont_map[i][j]){*c=j; c++;}
    *c=-1;
  }
  for(i=0; i<nres1;i++)free(Cont_map[i]); free(Cont_map);
  for(i=0; i<nres; i++)free(contact[i]);  free(contact);
  prot->cont_list=Contact2Contlist(prot->contact,nres1,prot->n_cont);

  // Secondary structure
  prot->sec_str=malloc(nres1*sizeof(char));
  prot->i_sec=malloc(nres1*sizeof(int));
  for(i=0; i<nres1; i++){
    (*res_short)[i].sec_str='D';
    prot->sec_str[i]='D'; prot->i_sec[i]=0;
  }
  Read_secondary_structure(*res, file, chain, nres);
  struct residue *r=*res;
  for(int i=0; i<nres; i++){
    int k=(*res_index)[i];
    if(r->sec_str != ' '){
      prot->sec_str[k]=r->sec_str;
      prot->i_sec[k]=r->i_sec;
    }else if(prot->sec_str[k]!='H' && prot->sec_str[k]!='E'){
      prot->sec_str[k]='C';
      prot->i_sec[k]=r->i_sec;
    }
    (*res_short)[k].sec_str=prot->sec_str[k];
    r++;
  }

  if(1){
    printf("Sequence:            ");
    for(i=0; i<nres1; i++)printf("%c", AMIN_CODE[prot->aa_seq[i]]);
    printf("\n");
    printf("Secondary structure: ");
    for(i=0; i<nres1; i++)printf("%c", prot->sec_str[i]);
    printf("\n");
    printf("16 states:           ");
    for(i=0; i<nres1; i++)printf("%c", SEC_EL[prot->i_sec[i]]);
    printf("\n");
  }
  return(nres1);
}


void GetPdbId(char *pdb_file_in, char *pdbid){
     /* This subroutine pretends to get the
        PDB id from a pdb file name, and ressembles
	quite a lot my "old-and-dirty-Perl" days */
	
  int start=0, end=0, i,j;
  for(i=strlen(pdb_file_in)-1;i>=0;i--){
    char *ptr=pdb_file_in+i;
    if (*ptr=='.'){
      end=i-1;
    }else if((*ptr=='/')||(*ptr=='\\')){
      start=i+1; break;
    }
  }
  j=0;
  for (i=start;i<=end;i++){
    pdbid[j]=pdb_file_in[i];
    j++;
  }
  pdbid[j]='\0';
}
/***************************************************************************/
/*                           Contact matrices                              */
/***************************************************************************/

short **Contact_matrix(struct residue *seq, int N_res, int *N_cont,
		       char l_cont, int ij_min)
{
  short i_res, j_res, contact=0, num_cont[L_MAX], **Cont_map;
  struct contact Contact_list[L_MAX*40], *cont;
  struct residue *seq_i;

  if(init_map==0){
    cont_thr_a2=cont_thr_a*cont_thr_a;
    cont_thr_b2=cont_thr_b*cont_thr_b;
    cont_thr_c2=cont_thr_c*cont_thr_c;
    init_map++;
  }

  (*N_cont)=0;
  cont=Contact_list;
  for(i_res=0; i_res<N_res; i_res++){
    num_cont[i_res]=0; seq_i=seq+i_res;
    for(j_res=i_res+ij_min; j_res< N_res; j_res++){
      if(l_cont=='c'){
	contact=Contact(*seq_i, seq[j_res]);
      }else if(l_cont=='b'){
	contact=Contact_beta(*seq_i, seq[j_res]);
      }else if(l_cont=='a'){
	contact=Contact_alpha(*seq_i, seq[j_res]);
      }
      if(contact){
	cont->res1=i_res; cont->res2=j_res; cont++;
	(*N_cont)++; num_cont[i_res]++; 
      }
    }
  }

  Cont_map=Store_map(N_res, seq, num_cont, *N_cont, Contact_list, 1);

  float cont_thr=cont_thr_c;
  if(l_cont=='a'){cont_thr=cont_thr_a;}
  else if(l_cont=='b'){cont_thr=cont_thr_b;}
  printf("Contact type %c, thr= %.2f, %d native contacts\n",
	 l_cont, cont_thr, *N_cont);

  return(Cont_map);
}

int Contact(struct residue res_i, struct residue res_j){
  atom *atom1=res_i.atom_ptr, *atom2;
  float dx, dy, dz, *r1, *r2; int i, j;

  //printf("%c %c\n", res_i.amm, res_j.amm);
  for(i=0; i<res_i.n_atom; i++){
    r1=atom1->r;
    atom2=res_j.atom_ptr;
    for(j=0; j<res_j.n_atom; j++){
      r2=atom2->r;
      dx=(*(r1)  -*(r2));   if(fabs(dx)>cont_thr_c) goto new;
      dy=(*(r1+1)-*(r2+1)); if(fabs(dy)>cont_thr_c) goto new;
      dz=(*(r1+2)-*(r2+2)); if(fabs(dz)>cont_thr_c) goto new;
      if((dx*dx+dy*dy+dz*dz)<=cont_thr_c2) return(1);
    new:
      atom2++;
    }
    atom1++;
  }
  return(0);
}


int Contact_beta(struct residue res_i, struct residue res_j){
  atom *atom1=Find_atom(res_i, "CB"), *atom2=Find_atom(res_j, "CB");
  float dx, dy, dz;

  if((atom1==NULL)||(atom2==NULL))return(0);
  float *r1=atom1->r, *r2=atom2->r;
   dx=(*(r1)  -*(r2));   if(fabs(dx)>cont_thr_b) return(0);
   dy=(*(r1+1)-*(r2+1)); if(fabs(dy)>cont_thr_b) return(0);
   dz=(*(r1+2)-*(r2+2)); if(fabs(dz)>cont_thr_b) return(0);
  if((dx*dx+dy*dy+dz*dz)<=cont_thr_b2)return(1);
  return(0);
}


int Contact_alpha(struct residue res_i, struct residue res_j)
{
  atom *atom1=Find_atom(res_i, "CA"), *atom2=Find_atom(res_j, "CA");
  float dx, dy, dz;
  float *r1=atom1->r, *r2=atom2->r;
   dx=(*(r1)  -*(r2));   if(fabs(dx)>cont_thr_a) return(0);
   dy=(*(r1+1)-*(r2+1)); if(fabs(dy)>cont_thr_a) return(0);
   dz=(*(r1+2)-*(r2+2)); if(fabs(dz)>cont_thr_a) return(0);
  if((atom1==NULL)||(atom2==NULL))return(0);

  if((dx*dx+dy*dy+dz*dz)<=cont_thr_a2)return(1);
  return(0);
}

short **Store_map(int N_res, struct residue *seq, short *num_cont,
		  int N_cont,struct contact *Contact_list, int l_cont)
{
  short **Cont_map=malloc(N_res*sizeof(short *));
  int i_res, j_res, i_cont;
  struct contact *cont;

  for(i_res=0; i_res< N_res; i_res++){
    Cont_map[i_res]=malloc((num_cont[i_res]+1)*sizeof(short));
    Cont_map[i_res][num_cont[i_res]]=-1;
    num_cont[i_res]=0;
    if(l_cont)seq[i_res].n_cont=0;
  }

  cont=Contact_list;
  for(i_cont=0; i_cont<N_cont; i_cont++){
    i_res=cont->res1; j_res=cont->res2; cont++;
    Cont_map[i_res][num_cont[i_res]]=j_res;
    num_cont[i_res]++;
    if(l_cont){
      seq[i_res].n_cont++; seq[j_res].n_cont++;
    }
  }
  return(Cont_map);
}

atom *Find_atom(struct residue seq, char code[2])
{
  atom *atom1=seq.atom_ptr; int i;
  for(i=0; i<seq.n_atom; i++){
    if(strncmp(atom1->name, code, 2)==0)return(atom1);
    atom1++;
  }
  return(NULL);
}

struct contact *Contact2Contlist(short **contact, int nres, int ncont)
{
  struct contact *cont_list=malloc(ncont*sizeof(struct contact));
  struct contact *cont=cont_list; int i;
  for(i=0; i<nres; i++){
    short *Ci=contact[i];
    while(*Ci >= 0){
      cont->res1=i; cont->res2=*Ci; cont++; Ci++;
    }
  }
  return(cont_list);
}

int Cluster_chains(char **seq, int *nchain_out, int *res_index,
		   int **ini_seq, int **n_seq,
		   struct residue *res, char *file,
		   char *chain_to_read, int nchain, int nres)
{
  // Read SEQRES
  FILE *file_in=fopen(file, "r");
  if(file_in==NULL){
    printf("\nERROR reading SEQRES, file %s not found\n", file);
    exit(8);
  }
  char string[1000]; int i, nmax=nres*2;
  int n_seqres[nchain]; char *seqres[nchain], *seq3[nchain];
  for(i=0; i<nchain; i++){
    n_seqres[i]=0;
    seq3[i]=malloc(3*nmax*sizeof(char));
    seqres[i]=malloc(nmax*sizeof(char));
  }
  int n_exo=0, NEXO=400; char *res_exo[NEXO], *res_std[NEXO];

  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string,"SEQRES", 6)==0){
      char chain=string[11]; int ichain=-1;
      for(i=0; i<nchain; i++)if(chain==chain_to_read[i]){ichain=i; break;}
      if(ichain<0)continue;
      Get_seqres(string, seq3[ichain], n_seqres+ichain);
    }else if(strncmp(string,"MODRES", 6)==0){
      Read_modres(res_exo, res_std, string, &n_exo);
    }else if(strncmp(string,"ATOM", 4)==0){
      break;
    }
  }
  fclose(file_in);
  int error=0;
  for(i=0; i<nchain; i++){
    if(n_seqres[i]==0 || n_seqres[i]>=nmax){
      printf("ERROR reading SEQRES of chain %c, nres=%d max.all.= %d\n",
	     chain_to_read[i], n_seqres[i], nmax); error=1;
    }
    char *s1=seqres[i], *s3=seq3[i];
    for(int k=0; k<n_seqres[i]; k++){
      *s1=Code_3_1(s3);
      if(*s1=='X')*s1=Het_res(s3, res_exo, res_std, n_exo);
      s1++; s3+=3;
    }
  }
  if(error)exit(8);

  // Cluster chains according to SEQRES
  int cluster[nchain]; for(i=0; i<nchain; i++)cluster[i]=i;
  int nc=nchain;
  for(i=0; i<nchain; i++){
    if(cluster[i]<i)continue; // chain already clustered
    char ci=chain_to_read[i];
    for(int j=i+1; j<nchain; j++){
      if(cluster[j]<j)continue; // chain already clustered
      if(n_seqres[i]!=n_seqres[j])continue;
      char cj=chain_to_read[j]; int d=0;
      for(int k=0; k<n_seqres[i]; k++)if(seqres[i][k]!=seqres[j][k])d++;
      printf("Chains %c and %c have %d differences\n", ci, cj, d);
      if(d==0){
	printf("Joining %c and %c\n", ci, cj); nc--;
	if(cluster[i]<cluster[j]){cluster[j]=cluster[i];}
	else{cluster[i]=cluster[j];}
      }
    }
  }
  printf("%d different chains found:\n", nc);
  int irep[nc]; for(i=0; i<nc; i++)irep[i]=-1;
  int cl_new[nchain]; for(i=0; i<nchain; i++)cl_new[i]=-1;
  int ini=0;
  *ini_seq=malloc(nc*sizeof(int));
  *n_seq=malloc(nc*sizeof(int));
  for(int cl=0; cl<nc; cl++){
    int cmin=nchain;
    for(i=0; i<nchain; i++){
      if(cl_new[i]<0 && cluster[i]<cmin){cmin=cluster[i]; irep[cl]=i;}
    }
    if(cmin==nchain){
      printf("ERROR, representative of cluster %d / %d not found\n",cl,nc);
      exit(8);
    }
    for(i=0; i<nchain; i++)if(cluster[i]==cmin)cl_new[i]=cl;
    (*ini_seq)[cl]=ini; (*n_seq)[cl]=n_seqres[irep[cl]];
    ini+=(*n_seq)[cl];
    printf("seq.cl. %d rep: %d ini: %d len: %d\n",
	   cl, irep[cl], (*ini_seq)[cl], (*n_seq)[cl]);
  }
  error=0;
  for(i=0; i<nchain; i++){
    if(cl_new[i]<0){
      printf("ERROR, chain %d has no assigned cluster out of %d\n",i,nc);
      error++;
    }
  }
  if(error)exit(8);
  *nchain_out=nc;
  int nres1=ini;

  *seq=malloc(nres1*sizeof(char)); char *s=*seq;
  for(i=0; i<nc; i++){
    char *sr=seqres[irep[i]]; int n=n_seqres[irep[i]];
    for(int j=0; j<n; j++){*s=*sr; s++; sr++;}
  }

  // Assign SEQRES index to residues in PDB.
  printf("Assigning index to %d residues in %d chains\n", nres, nchain);
  for(i=0; i<nres; i++)res_index[i]=-1;
  int j, ic, n, k; ini=0; //int ini_chain[nchain];
  char old_chain='*', *aa1;
  struct residue *r=res;
  for(i=0; i<nres; i++){
    if(r->chain!=old_chain){
      old_chain=r->chain; ic=-1;
      for(k=0; k<nchain; k++)if(r->chain==chain_to_read[k]){ic=k; break;}
      if(ic<0){printf("ERROR chain %c not found\n",r->chain); exit(8);}
      //ini_chain[ic]=i;
      int iclus=cl_new[ic]; ic=irep[iclus];
      j=0; aa1=seqres[ic]; n=n_seqres[ic]; ini=(*ini_seq)[iclus];
    }
    while(*aa1!=r->amm && j<n){aa1++; j++;} //&& *aa1!='X' 
    if(*aa1==r->amm)res_index[i]=ini+j;
    r++;
  }
  error=0;
  for(i=0; i<nres; i++){
    if(res_index[i]<0){
      printf("ERROR, residue %c%d chain %c not assigned\n",
	     res[i].amm, i, res[i].chain); error++;
    }
  }
  if(error){
    printf("%d not assigned residues, exiting\n", error); exit(8);
  }
  return(nres1);
}

void Get_seqres(char *string, char *seq3, int *n_seqres)
{
  char *ptr=string+19; 

  /* If amino acid chain, aa=1 */
  char *seq=seq3+3*(*n_seqres);
  for(int k=0; k<13; k++){
    for(int j=0; j<3; j++){*seq+=*ptr; seq++; ptr++;}
    (*n_seqres)++; ptr++; if(*ptr==' ')return;
  }

}

