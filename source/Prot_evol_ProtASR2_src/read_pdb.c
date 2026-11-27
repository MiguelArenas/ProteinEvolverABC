#include "allocate.h"
#include "protein3.h"
#include "coord.h"
#include "read_pdb.h"
#include "contact_matrix.h"
#include "Sec_str_all.h"
#include "codes.h"
#include "externals.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

int N_disulf, *Disulf_res1, *Disulf_res2;
float *Disulf_d;

static void GetPdbId(char *pdb_file_in, char *pdbid);

// Contact matrices
float cont_thr_a=8, cont_thr_a2;
float cont_thr_b=8, cont_thr_b2;
float cont_thr_c=4.5, cont_thr_c2;
float disulf_thr=3.0, disulf_thr_2;
int init_map=0;
static int Contact_alpha(struct residue res_i, struct residue res_j);
static int Contact_beta (struct residue res_i, struct residue res_j);
static int Contact(struct residue res_i, struct residue res_j);
static float Disulfide(struct residue *res_i, struct residue *res_j);
static atom *Find_atom(struct residue seq, char code[2]);
static int Get_chain(char chain, char *chains, int nchain);
short **Contact_list(struct residue *seq, int N_res, int *N_cont,
		     char l_cont, int ij_min);
short **Store_map(int N_res, struct residue *seq, short *num_cont,
		  int N_cont,struct contact *Contact_list, int l_cont);
void Get_seqres(char *string, char *seqres, int *n_seqres);
int Cluster_chains(char **seq, int *nchain_out,
		   int *clust, int *num_ch, int *res_index,
		   int **ini_seq, int **n_seq,
		   struct residue *res, char *file,
		   char *chain_to_read, int nchain, int nres);

// Routines:
int Get_pdb(struct protein *prot,
	    struct residue **res, struct res_short **res_short,
	    char *file, char *chain, int **res_index, int nmod)
{
  int natoms, natoms2, ANISOU, nmr, i;

  prot[0].length=0;
  int nres=Read_coord(file, &nmr, NULL, NULL, chain, &ANISOU, &natoms, 0);
  // The last argument of kmod is the model to read; if 0 it reads the
  // first structure even if the record MODEL is not present

  if(nres<=0){
    printf("WARNING, pdb file %s, no residue found in chain %s\n",
	   file, chain); return(0);
  }
  printf("Chains read: %s\n", chain);

  *res=malloc(1.1*nres*sizeof(struct residue));
  atom *atom_read=malloc(1.1*natoms*sizeof(atom));
  int nres2=
    Read_coord(file, &nmr, *res, atom_read, chain, &ANISOU, &natoms2, 0);

  if((nres2!=nres)||(natoms2!=natoms)){
    printf("ERROR, inconsistent number of residues/atoms in %s %s\n",
	   file, chain);
    printf("Before: %d %d Now: %d %d\n", nres, natoms, nres2, natoms2);
    exit(8);
  }

  // File name
  char pdbid[20];
  GetPdbId(file, pdbid); if(chain[0]==' ')chain[0]='_';
  printf("PDB %s chain %s nres=%d\n",pdbid,chain,nres);
  sprintf(prot[0].name, "%s%s", pdbid, chain);

  // Count chains
  i=0; while(chain[i]!='\0')i++;
  int nchain=i;
  if(nchain==0){
    printf("WARNING: 0 chains found in %s, setting to 1\n", chain);
    nchain=1;
  }
  printf("prot %s is a complex of %d chains %s\n",
	 prot[0].name, nchain, chain);

  /* Cluster chains according to seqres. Output:
     nc clusters with total of nres_diff different a.a.
     cluster i contains num_ch chains
     chain i belongs to cluster clust[i]
     seq contains the nres1 SEQRES redisues
     res_index contains the correspondence btw nres structured
     PDB residues and nres_diff SEQRES residues
     ini_aa_cl[i]: index of first aa in cluster i (1-nc)
     n_aa_cl[i]:  number of aa in cluster i (1-nc)
     ini_aa_ch[i]: index of first aa in chain i (1-nchain)
     n_aa_ch[i]:  number of aa in chain i (1-nchain)
  */
  char *seq=NULL; int *ini_aa_cl=NULL, *n_aa_cl=NULL;
  int clust[nchain], num_ch[nchain];
  *res_index=malloc(nres*sizeof(int));
  int nc, nres_diff=
    Cluster_chains(&seq, &nc, clust, num_ch, *res_index, &ini_aa_cl, &n_aa_cl,
		   *res, file, chain, nchain, nres);
  int nres_tot=n_aa_cl[0]*num_ch[0]; // Number of AA in all SEQRES records
  int CLUSTER=1, ic;
  int num_chain=num_ch[0];
  for(ic=1; ic<nc; ic++){
    //if(num_ch[ic]!=num_ch[0])CLUSTER=0;
    nres_tot += n_aa_cl[ic]*num_ch[ic];
    if(num_ch[i]>num_chain)num_chain=num_ch[ic];
  }
  printf("prot %s has %d over %d chains different by sequence with "
	 "%d different residues, and %d over %d ordered residues\n",
	 prot[0].name, nc, nchain, nres_diff, nres, nres_tot);

  /*int n_aa_ch[nchain], ini_aa_ch[nchain], ini=0; 
  for(i=0; i<nchain; i++){
    n_aa_ch[i]=n_aa_cl[clust[i]];
    ini_aa_ch[i]=ini; ini+=n_aa_ch[i];
    //printf("%d c=%d l=%d\n", i, clust[i], n_aa_ch[i]);
    }*/

  int nres1=0, nch1=0;
  if(CLUSTER==0){
    nres1=nres_tot; nch1=nchain; num_chain=1;
    //nres1=nres_diff; nch1=nc; num_chain=1;
    printf("Not all clusters have the same number of chains, do not cluster\n");
    //for(i=0; i<nchain; i++)clust[i]=i;
  }else{
    nres1=nres_diff; nch1=nc; num_chain=num_ch[0]; 
    printf("All %d clusters have %d chains, cluster them\n", nc, num_chain);
  }

  // Order chains
  prot[0].nchain=nch1;
  prot[0].num_chain=num_chain;
  prot[0].chain_all=nchain;
  prot[0].ini_chain=malloc(nch1*sizeof(int));
  prot[0].len_chain=malloc(nch1*sizeof(int));
  // For each residue in the structure, seq_index represents the seq position
  { 
    int k=0, l;
    for(ic=0; ic<nc; ic++){
      l=n_aa_cl[ic];
      prot[0].ini_chain[ic]=k;
      prot[0].len_chain[ic]=l;
      k+=l;
    }
    if(k != nres_diff){
      printf("ERROR clustering %d chains, %d != %d residues\n",nch1,k,nres1);
      for(int i=0; i<nchain; i++)
	printf("Chain %d cl= %d L= %d\n", i, clust[i], n_aa_cl[clust[i]]);
      exit(8);
    }
  }

  // Align all seqs to cluster representatives
  int all_to_clust[nres_tot]; int k=0;
  for(ic=0; ic<nchain; ic++){
    int c=clust[ic], l=n_aa_cl[c], kj=ini_aa_cl[c];
    for(int j=0; j<l; j++){
      all_to_clust[k]=kj; k++; kj++;
    }
    if(k> nres_tot || kj>nres_diff){
      printf("ERROR aligning chain %d, %d > %d total residues or "
	     "%d > %d different residues\n", ic,k,nres_tot,kj,nres_diff);
      exit(8);
    }
  }
  /*printf("all_to_clust= ");
  for(i=0; i<nres_tot; i++){printf(" %d", all_to_clust[i]);}
  printf("\n");*/
  
  // Redirect res_index to cluster representative
  if(CLUSTER){
    for(i=0; i<nres; i++){
      int k=(*res_index)[i];
      (*res_index)[i]=all_to_clust[k];
      if((*res_index)[i] >= nres1){
	printf("ERROR, residue %d has index %d >= %d k= %d\n",
	       i, (*res_index)[i], nres1, k); exit(8);
      }
    }
  }

  prot[0].L_PDB=nres;   // Structured residues
  prot[0].length=nres1; // Residues in record SEQRES of repr. seq.
  // WARNING: nres1 can be larger than nres (if disordered residues)
  // or smaller (if several identical chains)
  // Representative sequence (one per cluster of identical chains)
  short *aa_seq=malloc(nres1*sizeof(short));
  (*res_short)=malloc(nres1*sizeof(struct res_short));
  prot[0].aa_seq=aa_seq;
  for(i=0; i<nres1; i++){
    if(CLUSTER){k=i;}
    else{k=all_to_clust[i];}
    (*res_short)[i].seq=seq[k];
    aa_seq[i]=Code_AA(seq[k]);
    (*res_short)[i].aa=aa_seq[i];
    (*res_short)[i].str_index=-1;
    (*res_short)[i].sec_str='D';
    (*res_short)[i].n_cont=0;
    sprintf((*res_short)[i].pdbres, "  -1");
  }

  // Residue index
  for(i=0; i<nres; i++){
    int k=(*res_index)[i];
    if(k<0 || k>=nres1){
      printf("ERROR, index %d of res %d not allowed, max= %d\n",
	     k, i, nres1-1); exit(8);
    }
    sprintf((*res_short)[k].pdbres, "%s", (*res)[i].pdbres);
    (*res_short)[k].str_index=i;
  }
  Read_secondary_structure(*res, file, chain, nres);

  // Read other models
  for(int jmod=0; jmod<nmod; jmod++){

    struct protein *model=prot+jmod;
    if(jmod){
      Read_coord(file, &nmr, *res, atom_read, chain, &ANISOU, &natoms2, jmod);
      model->aa_seq=NULL;
    }
    model->contact=NULL;
    model->cont_list=NULL;
    model->sec_str=NULL;
    model->i_sec=NULL;
    model->length=nres1;
    model->L_PDB=nres;

    // Contact matrix and number of contacts
    int n_cont[nres1], CMAX=30*num_chain;
    
    model->contact=malloc(nres1*sizeof(short *));
    for(i=0; i<nres1; i++){
      model->contact[i]=malloc(CMAX*sizeof(short));
      for(int j=0; j<CMAX; j++)model->contact[i][j]=-1;
      n_cont[i]=0;
    }
    short **contact=
      Contact_list(*res, nres, &model->n_cont, 'c', IJ_MIN);

    // Inter-chain contacts
    struct inter_chain *inter_chain=NULL;
    int N_inter_chain_max=0;
    model->N_inter_chain=0;
    if(nchain > 1){
      N_inter_chain_max=0.25*model->n_cont;
      inter_chain=malloc(N_inter_chain_max*sizeof(struct inter_chain));
      model->inter_chain=inter_chain;
    }else{
      model->inter_chain=NULL;
    }

    for(i=0; i<nres; i++){
      int k=(*res_index)[i];
      if(k<0 || k>=nres1)continue;
      short *j=contact[i];
      while(*j>=0 && *j<nres){
	int kj=(*res_index)[*j];
	if(kj>=0 && kj<nres1){
	  model->contact[k][n_cont[k]]=kj; n_cont[k]++;
	  (*res_short)[k].n_cont++; (*res_short)[kj].n_cont++;
	  // Inter-chain contact
	  if((*res)[i].chain != (*res)[*j].chain){
	    if(model->N_inter_chain >= N_inter_chain_max){
	      printf("ERROR Interface %d%c-%d%c, > %d inter_chain contacts "
		     "for %d chains\n", i, (*res)[i].chain,
		     *j, (*res)[*j].chain, N_inter_chain_max, nchain);
	      exit(8);
	    }
	    struct inter_chain *store=inter_chain+model->N_inter_chain;
	    (model->N_inter_chain)++;
	    store->res1=k;
	    store->chain1=(*res)[i].chain;
	    store->ichain1=Get_chain(store->chain1, chain, nchain);
	    store->res2=kj;
	    store->chain2=(*res)[*j].chain;
	    store->ichain2=Get_chain(store->chain2, chain, nchain);
	  }
	  j++;
	}
      }
    }
    if(nchain > 1){
      printf("%d inter chain contacts\n", model->N_inter_chain);
    }

    for(i=0; i<nres1; i++)model->contact[i][n_cont[i]]=-1;
    model->cont_list=Contact2Contlist(model->contact,nres1,model->n_cont);
    if(num_chain!=1){
      for(i=0; i<nres1; i++)(*res_short)[i].n_cont/=num_chain;
    }


    // Disulfide bonds
    float Dis_thr=3.5;
    int N_dis_max=nres/10+1;
    N_disulf=0;
    if(Disulf_res1)free(Disulf_res1);
    Disulf_res1=malloc(N_dis_max*sizeof(int));
    if(Disulf_res2)free(Disulf_res2);
    Disulf_res2=malloc(N_dis_max*sizeof(int));
    if(Disulf_d)free(Disulf_d);
    Disulf_d=malloc(N_dis_max*sizeof(float));
    for(i=0; i<nres; i++){
      int k=(*res_index)[i];
      if(k<0 || k>=nres1 || seq[k]!='C')continue;
      short *j=contact[i];
      while(*j>=0 && *j<nres){
	int kj=(*res_index)[*j]; j++;
	if(kj>=0 && kj<nres1 && seq[kj]=='C'){
	  float d=Disulfide(*res+k, *res+kj);
	  if(d>0 && d<Dis_thr){
	    if(N_disulf==N_dis_max){
	      printf("ERROR, too many disulfides: %d for %d long protein\n",
		     N_disulf, nres); exit(8);
	    }
	    Disulf_res1[N_disulf]=k;
	    Disulf_res2[N_disulf]=kj;
	    Disulf_d[N_disulf]=d;
	    N_disulf++;
	  }
	}
      }
    }
    printf("%d disulfide bonds found\n", N_disulf);

    // Empty contact
    for(i=0; i<nres; i++)free(contact[i]);
    free(contact);

    // Secondary structure
    model->sec_str=malloc(nres1*sizeof(char));
    model->i_sec=malloc(nres1*sizeof(int));
    for(i=0; i<nres1; i++){
      model->sec_str[i]='D'; model->i_sec[i]=0;
    }
    struct residue *r=*res;
    for(int i=0; i<nres; i++){
      int k=(*res_index)[i];
      if(k<0 || k>=nres1){r++; continue;}
      if(r->sec_str != ' '){
	model->sec_str[k]=r->sec_str;
	model->i_sec[k]=r->i_sec;
      }else{// if(prot->sec_str[k]!='H' && prot->sec_str[k]!='E'){
	model->sec_str[k]='C';
	model->i_sec[k]=0; //r->i_sec;
      }
      (*res_short)[k].sec_str=model->sec_str[k];
      r++;
    }
    
  } // end models


  if(1){
    printf("Sequence:            ");
    for(i=0; i<nres1; i++)printf("%c", AMIN_CODE[prot[0].aa_seq[i]]);
    printf("\n");
    printf("Secondary structure: ");
    for(i=0; i<nres1; i++)printf("%c", prot[0].sec_str[i]);
    printf("\n");
    printf("16 states:           ");
    for(i=0; i<nres1; i++)printf("%c", SEC_EL[prot[0].i_sec[i]]);
    printf("\n");
    //printf("i_sec:               ");
    //for(i=0; i<nres1; i++)printf("%d", prot[0].i_sec[i]);
    //printf("\n");
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

short **Contact_list(struct residue *seq, int N_res, int *N_cont,
		     char l_cont, int ij_min)
{
  int CMAX=50;
  short **Cont_map=malloc(N_res*sizeof(short *));

  if(init_map==0){
    cont_thr_a2=cont_thr_a*cont_thr_a;
    cont_thr_b2=cont_thr_b*cont_thr_b;
    cont_thr_c2=cont_thr_c*cont_thr_c;
    disulf_thr_2=disulf_thr*disulf_thr;
    init_map++;
  }

  (*N_cont)=0;
  int num_cont[N_res];
  for(int i_res=0; i_res<N_res; i_res++)seq[i_res].n_cont=0;
  for(int i_res=0; i_res<N_res; i_res++){
    Cont_map[i_res]=malloc(CMAX*sizeof(short));
    if(Cont_map[i_res]==NULL){
      printf("Out of memory in Cont_map %d\n", i_res); exit(8);
    }
    num_cont[i_res]=0;
    struct residue *seq_i=seq+i_res; int contact=0;
    for(int j_res=i_res+ij_min; j_res< N_res; j_res++){
      if(l_cont=='a'){
	contact=Contact_alpha(*seq_i, seq[j_res]);
      }else if(l_cont=='b'){
	contact=Contact_beta(*seq_i, seq[j_res]);
      }else{
	contact=Contact(*seq_i, seq[j_res]);
      }
      if(contact){
	Cont_map[i_res][num_cont[i_res]]=j_res;
	num_cont[i_res]++; (*N_cont)++;
	seq[i_res].n_cont++; seq[j_res].n_cont++;
      }
    }
    if(num_cont[i_res]>=CMAX){
      printf("ERROR, %d contacts for residue %i, max is %d\n",
	     num_cont[i_res], i_res, CMAX); exit(8);
    }
  }
  for(int i_res=0; i_res<N_res; i_res++){
    Cont_map[i_res][num_cont[i_res]]=-1;
    /*printf("%3d:  ", i_res);
    for(int j=0; j<num_cont[i_res]; j++)printf("%d ",Cont_map[i_res][j]);
    printf("\n");*/
  }

  cont_thr=cont_thr_c;
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

float Disulfide(struct residue *res_i, struct residue *res_j)
{
  // Only cysteine residues
  if(res_i->amm != 'C' || res_j->amm != 'C')return(-1);
  atom *atom1=res_i->atom_ptr, *atom2;
  float dx, dy, dz, *r1, *r2; int i, j;

  //printf("%c %c\n", res_i->amm, res_j->amm);
  for(i=0; i<res_i->n_atom; i++){
    if(strncmp(atom1->name, "SG", 2)) goto new_1;
    r1=atom1->r;
    atom2=res_j->atom_ptr;
    for(j=0; j<res_j->n_atom; j++){
      if(strncmp(atom2->name, "SG", 2)) goto new_2;
      r2=atom2->r;
      dx=(*(r1)  -*(r2));   if(fabs(dx)>disulf_thr) goto new_2;
      dy=(*(r1+1)-*(r2+1)); if(fabs(dy)>disulf_thr) goto new_2;
      dz=(*(r1+2)-*(r2+2)); if(fabs(dz)>disulf_thr) goto new_2;
      float d2=dx*dx+dy*dy+dz*dz;
      if(d2<=disulf_thr_2)return(sqrt(d2));
    new_2:
      atom2++;
    }
  new_1:
    atom1++;
  }
  return(-1);
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
    if(Cont_map[i_res]==NULL){
      printf("Out of memory in Store_map %d\n", i_res); exit(8);
    }
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
  if(cont_list==NULL){
    printf("Out of memory in Contact2Contlist %d %d\n", nres, ncont);
    exit(8);
  }
  for(i=0; i<nres; i++){
    short *Ci=contact[i];
    while(*Ci >= 0){
      cont->res1=i; cont->res2=*Ci; cont++; Ci++;
    }
  }
  return(cont_list);
}

int Cluster_chains(char **seq, int *nchain_out,
		   int *clust_new, int *num_ch, int *res_index,
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
  char string[1000]; int i, nmax=nres*3;
  int n_seqres[nchain]; char *seqres[nchain], *seq3[nchain];
  for(i=0; i<nchain; i++){
    n_seqres[i]=0;
    seq3[i]=malloc(3*nmax*sizeof(char));
    seqres[i]=malloc(nmax*sizeof(char));
    if(seq3[i]==NULL || seqres[i]==NULL){
      printf("Out of memory in Cluster_chains seq chain %d %d\n", i, nmax);
      exit(8);
    }
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
      if(*s1=='X'){
	printf("(residue %d of %d chain %d)\n", k, n_seqres[i], i);
	*s1=Het_res(s3, res_exo, res_std, n_exo);
      }
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
      printf("Chains %c and %c have %d differences, ", ci, cj, d);
      if(d==0){
	printf("joining them\n"); nc--;
	if(cluster[i]<cluster[j]){cluster[j]=cluster[i];}
	else{cluster[i]=cluster[j];}
      }else{
	printf("not joining them\n");
      }
    }
  }
  printf("%d different chains found:\n", nc);
  int irep[nc]; for(i=0; i<nc; i++)irep[i]=-1;
  for(i=0; i<nchain; i++){clust_new[i]=-1; num_ch[i]=0;}
  int ini=0;
  *ini_seq=malloc(nchain*sizeof(int));
  *n_seq=malloc(nchain*sizeof(int));
  for(int cl=0; cl<nc; cl++){
    int cmin=nchain;
    for(i=0; i<nchain; i++){
      if(clust_new[i]<0 && cluster[i]<cmin){cmin=cluster[i]; irep[cl]=i;}
    }
    if(cmin==nchain){
      printf("ERROR, representative of cluster %d / %d not found\n",cl,nc);
      exit(8);
    }
    for(i=0; i<nchain; i++){
      if(cluster[i]==cmin){clust_new[i]=cl; num_ch[cl]++;}
    }
    (*ini_seq)[cl]=ini; (*n_seq)[cl]=n_seqres[irep[cl]];
    ini+=(*n_seq)[cl];
    printf("seq.cl. %d %d chains ini: %d len: %d rep: %d\n",
	   cl, num_ch[cl], (*ini_seq)[cl], (*n_seq)[cl], irep[cl]);
  }
  error=0;
  for(i=0; i<nchain; i++){
    if(clust_new[i]<0){
      printf("ERROR, chain %d has no assigned cluster out of %d\n",i,nc);
      error++;
    }
  }
  if(error)exit(8);
  *nchain_out=nc;
  int nres1=ini;

  *seq=malloc(nres1*sizeof(char)); char *s=*seq;
  if(*seq ==NULL){
    printf("Out of memory in Cluster_chains seq %d\n", nres1);
    exit(8);
  }
  for(i=0; i<nc; i++){
    char *sr=seqres[irep[i]]; int n=n_seqres[irep[i]];
    for(int j=0; j<n; j++){*s=*sr; s++; sr++;}
  }

  // Assign SEQRES index to residues in PDB.
  printf("Assigning index to %d structured residues in %d chains\n",
	 nres, nchain);
  for(i=0; i<nres; i++)res_index[i]=-1;
  int j=0, ic, n=n_seqres[0], k; ini=0; //int ini_chain[nchain];
  char old_chain='*', *sc=seqres[0];
  struct residue *r=res;
  for(i=0; i<nres; i++){
    if(r->chain!=old_chain){
      old_chain=r->chain; ic=-1;
      for(k=0; k<nchain; k++)if(r->chain==chain_to_read[k]){ic=k; break;}
      if(ic<0){printf("ERROR chain %c not found\n",r->chain); exit(8);}
      if(0){ // res_index points to chain cluster
	int iclus=clust_new[ic]; ic=irep[iclus]; ini=(*ini_seq)[iclus];
      }else{ // res_index points to all chains
	ini=0; for(int c=0; c<ic; c++)ini+=n_seqres[c];
      }
      j=0; sc=seqres[ic]; n=n_seqres[ic];
    }
    // int ires; sscanf(r->pdbres, "%d", &ires);
    while(sc[j] !=r->amm && j<n){j++;} //&& *aa1!='X' 
    if(sc[j]==r->amm){res_index[i]=ini+j; j++;}
    r++;
  }
  printf("Largest index: %d\n", res_index[nres-1]);
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
  for(i=0; i<nchain; i++){
    free(seqres[i]); free(seq3[i]);
  }

  return(nres1);
}

void Get_seqres(char *string, char *seq3, int *n_seqres)
{
  char *ptr=string+19; 

  /* If amino acid chain, aa=1 */
  char *seq=seq3+3*(*n_seqres);
  for(int k=0; k<13; k++){
    for(int j=0; j<3; j++){*seq=*ptr; seq++; ptr++;}
    (*n_seqres)++; ptr++; if(*ptr==' ')return;
  }

}

 int Get_chain(char chain, char *chains, int nchain){
   for(int i=0; i<nchain; i++)if(chains[i]==chain)return(i);
   printf("ERROR, chain %c not found in %s\n", chain, chains);
   exit(8); return(-1);
 }
