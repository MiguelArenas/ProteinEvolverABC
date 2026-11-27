#include "coord.h"
#include "protein3.h"
#include "contact_matrix.h"
#include "externals.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define DUMM 1000000

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

short **Contact_matrix(struct residue *seq, int N_res, int *N_cont,
		       char l_cont, int IJ_MIN)
{
  short i_res, j_res, contact=0, num_cont[L_MAX], **Cont_map;
  struct contact Contact_list[L_MAX*40], *cont;

  if(init_map==0){
    cont_thr_a2=cont_thr_a*cont_thr_a;
    cont_thr_b2=cont_thr_b*cont_thr_b;
    cont_thr_c2=cont_thr_c*cont_thr_c;
    init_map++;
  }

  (*N_cont)=0;
  cont=Contact_list;

  for(i_res=0; i_res<N_res; i_res++){
    num_cont[i_res]=0;
    for(j_res=i_res+IJ_MIN; j_res< N_res; j_res++){
      if(l_cont=='c'){
	contact=Contact(seq[i_res], seq[j_res]);
      }else if(l_cont=='b'){
	contact=Contact_beta(seq[i_res], seq[j_res]);
      }else if(l_cont=='a'){
	contact=Contact_alpha(seq[i_res], seq[j_res]);
      }
      if(contact){
	cont->res1=i_res; cont->res2=j_res; cont++;
	(*N_cont)++; num_cont[i_res]++; 
      }
    }
  }

  Cont_map=Store_map(N_res, seq, num_cont, *N_cont, Contact_list, 1);

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

  for(i=0; i<res_i.n_atom; i++){
    r1=atom1->r;
    atom2=res_j.atom_ptr;
    for(j=0; j<res_j.n_atom; j++){
      r2=atom2->r;
      dx=((*r1)  -(*r2));   if(fabs(dx)>cont_thr_c) goto new;
      dy=((*r1+1)-(*r2+1)); if(fabs(dy)>cont_thr_c) goto new;
      dz=((*r1+2)-(*r2+2)); if(fabs(dz)>cont_thr_c) goto new;
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
   dx=((*r1)-(*r2));     if(fabs(dx)>cont_thr_b) return(0);
   dy=((*r1+1)-(*r2+1)); if(fabs(dy)>cont_thr_b) return(0);
   dz=((*r1+2)-(*r2+2)); if(fabs(dz)>cont_thr_b) return(0);
  if((dx*dx+dy*dy+dz*dz)<=cont_thr_b2)return(1);
  return(0);
}


int Contact_alpha(struct residue res_i, struct residue res_j)
{
  atom *atom1=Find_atom(res_i, "CA"), *atom2=Find_atom(res_j, "CA");
  float dx, dy, dz;
  float *r1=atom1->r, *r2=atom2->r;
   dx=((*r1)-(*r2));     if(fabs(dx)>cont_thr_a) return(0);
   dy=((*r1+1)-(*r2+1)); if(fabs(dy)>cont_thr_a) return(0);
   dz=((*r1+2)-(*r2+2)); if(fabs(dz)>cont_thr_a) return(0);
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
