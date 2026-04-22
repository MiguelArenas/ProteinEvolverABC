#include "coord.h"
#include "Sec_str_all.h"
#include "Sec_str_16_NR.h"
#include "Sec_str_comp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"
#include "codes.h"
#include "allocate.h"
#include <math.h>

int print_loc=1;
double *E_loc_misf=NULL;
double *Set_E_loc_misf();

/* 
int N_SS; // Number of sec.str. types 
char SEC_EL[N_SS];
double **Propensity; // Energy parameters
*/

/*
  Read secondary structure
*/
struct secondary{
  char type, chain;
  char ini_res[6];
  char end_res[6];
};

char SS_label[16]="CTS01234HcXABEYZ";
float f_SS[16]={0.1014,0.1005,0.0796,0.0248,0.0396,0.0370,0.0338,0.0285,0.2036,0.0378,0.0081,0.0306,0.0438,0.1549,0.0489,0.0271};
double E_loc_0[20];
float **E_loc_over_T=NULL;

int Read_sec_str(struct secondary *sec_str, int *N_sec_str, char *string,
		 char chain, char *chain_to_read, char type);
void Assign_sec_str_16(struct residue *res, int nres);
int SS_code(char sec_str);
int Find_sec_str(int *i1, int *i2, struct residue *res, int nres,
		 int ires, struct secondary *sec);
int Find_res(struct residue *res, int nres, int kres, char *pdbres, char chain);

void Combine_propensity(int i1, int i2);
#define SEC_STR_MAX 100000

int Read_secondary_structure(struct residue *res, char *file,
			     char *chain_to_read, int nres)
{
  FILE *file_in=fopen(file, "r");
  if(file_in==NULL){
    printf("ERROR, file %s not found\n", file); exit(8);
  }
  char string[200], chain;
  int N_sec_str=0;
  struct secondary sec_str[SEC_STR_MAX];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string,"HELIX ", 6)==0){
      chain=string[19];
      Read_sec_str(sec_str, &N_sec_str,  string, chain, chain_to_read, 'H');
    }else if(strncmp(string,"TURN ", 5)==0){
      chain=string[19];
      Read_sec_str(sec_str, &N_sec_str,  string, chain, chain_to_read, 'T');
    }else if(strncmp(string,"SHEET ", 6)==0){
      chain=string[21];
      Read_sec_str(sec_str, &N_sec_str,  string, chain, chain_to_read, 'E');
    }else if(strncmp(string,"ATOM", 4)==0){
      break;
    }
  }
  fclose(file_in);
  if(N_sec_str > SEC_STR_MAX){
    printf("ERROR, too many secondary structure elements %d max: %d\n",
	   N_sec_str, SEC_STR_MAX); exit(8);
  }

  int i, i1, i2, ires=0, isec;
  for(i=0; i<nres; i++)res[i].sec_str=' ';
  struct secondary *sec=sec_str;
  for(isec=0; isec<N_sec_str; isec++){
    int found=Find_sec_str(&i1, &i2, res, nres, ires, sec);
    if(found){
      for(i=i1; i<=i2; i++)res[i].sec_str=sec->type;
      ires=i2;
    }else{
      printf("WARNING Sec. str. %c%d (%s - %s %c) not found\n",
	     sec->type, isec, sec->ini_res, sec->end_res, sec->chain);
    }
    sec++;
  }
  Assign_sec_str_16(res, nres);
  return(1);
}

void Assign_sec_str_16(struct residue *res, int nres)
{
  /* Helix = 0(1234H...c)X
     Strand= A (BE.. Y)Z
     Coil= C T S
  */

  struct residue *r=res; char ss_old='C';
  for(int ires=0; ires<nres; ires++){
    /*if((r->i_aa<0)||(r->i_aa==20)){
      r->c_sec='C'; r++; continue;
      }*/
    if((r->sec_str=='G')||(r->sec_str=='H')||(r->sec_str=='I')){
      r->c_sec='H';
      if(ss_old!='H'){
	if(ss_old=='0'){r->c_sec='1';}
	else if(ss_old=='1'){r->c_sec='2';}
	else if(ss_old=='2'){r->c_sec='3';}
	else if(ss_old=='3'){r->c_sec='4';}
	//else if(ss_old=='3'){r->c_sec='H';}
	else if(ss_old=='4'){r->c_sec='H';}
	else{r->c_sec='0';}
	//else{r->c_sec='1'; if(ss_old=='C')(r-1)->c_sec='0';}
      }
    }else if((r->sec_str=='E')||(r->sec_str=='B')){
      r->c_sec='E';
      if((ss_old!='E')&&(ss_old!='B')){
	r->c_sec='B';
	if(ss_old=='C')(r-1)->c_sec='A';
      }
    }else if((r->sec_str=='T')||(r->sec_str=='S')){
      r->c_sec=r->sec_str;
    }else{
      r->c_sec='C';
      if(ss_old=='H'){
	/*(r-1)->c_sec='c';
	  r->c_sec='X';*/
	(r-1)->c_sec='X';
      }else if(ss_old=='E'){
	(r-1)->c_sec='Y';
	//r->c_sec='Z';
      }
    }
    ss_old=r->c_sec;
    r++;
  }
  r=res;
  for(int ires=0; ires<nres; ires++){
    r->i_sec=SS_code(r->c_sec); r++;
  }
}

int Find_sec_str(int *i1, int *i2, struct residue *res, int nres,
		 int ires, struct secondary *sec)
{
  int jres, ini, end;
  for(jres=ires; jres<nres; jres++)if(res[jres].chain==sec->chain)goto fres;
  for(jres=0; jres<ires; jres++)if(res[jres].chain==sec->chain)goto fres;
  printf("Sec. str. chain %c not found\n", sec->chain);
  return(0); // not found
 fres:
  ini=Find_res(res, nres, jres, sec->ini_res, sec->chain);
  if(ini<0){
    printf("Initial res. %s %c not found\n", sec->ini_res, sec->chain);
    return(0); // not found
  }
  end=Find_res(res, nres, ini, sec->end_res, sec->chain);
  if(end<0){
    printf("Final res. %s %c not found\n", sec->end_res, sec->chain);
    return(0); // not found
  }
  if(ini<end){*i1=ini; *i2=end;}else{*i1=end; *i2=ini;}
  return(1);
}

int Find_res(struct residue *res, int nres, int kres, char *pdbres, char chain)
{
  int jres=kres;
  while(res[jres].chain==chain && jres<nres){
    if(strcmp(res[jres].pdbres, pdbres)==0)return(jres);
    jres++;
  }
  jres=kres;
  while(res[jres].chain==chain && jres>=0){
    if(strcmp(res[jres].pdbres, pdbres)==0)return(jres);
    jres--;
  }
  return(-1);
}


int Read_sec_str(struct secondary *sec_str, int *N_sec_str, char *string,
		 char chain, char *chain_to_read, char type)
{
  if(*chain_to_read!='*'){
    for(int i=0; i<sizeof(chain_to_read); i++)
      if(chain==chain_to_read[i])goto read;
    return(0);
  }

 read:
  sec_str[*N_sec_str].type=type;
  char pdbres[6]; pdbres[5]='\0';
  /* initial residue */
  if(type=='H'){
    pdbres[0]=string[21];  pdbres[1]=string[22];  pdbres[2]=string[23]; 
    pdbres[3]=string[24];  pdbres[4]=string[25];
  }else if(type=='E'){
    pdbres[0]=string[22];  pdbres[1]=string[23];  pdbres[2]=string[24]; 
    pdbres[3]=string[25];  pdbres[4]=string[26];
  }else if(type=='T'){
    pdbres[0]=string[20];  pdbres[1]=string[21];  pdbres[2]=string[22]; 
    pdbres[3]=string[23];  pdbres[4]=string[24];
  }
  strcpy(sec_str[*N_sec_str].ini_res, pdbres);
  if(type=='T'){
    strcpy(sec_str[*N_sec_str].end_res, pdbres);
    sec_str[*N_sec_str].chain=chain;
    (*N_sec_str)++;
  }
  /* final residue */
  if((type=='H')||(type=='E')){
    pdbres[0]=string[33];  pdbres[1]=string[34];  pdbres[2]=string[35]; 
    pdbres[3]=string[36];  pdbres[4]=string[37];
  }else if(type=='T'){
    pdbres[0]=string[31];  pdbres[1]=string[32];  pdbres[2]=string[33]; 
    pdbres[3]=string[34];  pdbres[4]=string[35];
    strcpy(sec_str[*N_sec_str].ini_res, pdbres);
  }    
  strcpy(sec_str[*N_sec_str].end_res, pdbres);
  sec_str[*N_sec_str].chain=chain;
  (*N_sec_str)++;
  return(1);
}

/*
  Sec_str_code
 */
int SS_code(char sec_str) {
  int i;
  for(i=0; i<N_SS; i++)
    if(sec_str==SEC_EL[i])return(i);

  // Equivalences

  // Alpha helices
  if(sec_str=='4'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='H')return(i);
  }else if(sec_str=='c'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='H')return(i);
  }else if(sec_str=='X'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='H')return(i);
  }else if(sec_str=='1'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='H')return(i);
  }else if(sec_str=='2'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='H')return(i);
  }else if(sec_str=='3'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='H')return(i);

    // No distinction along strands
  }else if(sec_str=='A'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='C')return(i);
  }else if(sec_str=='B'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='E')return(i);
  }else if(sec_str=='Y'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='E')return(i);
  }else if(sec_str=='Z'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='C')return(i);

  // Parallel-antiparallel beta
  }else if(sec_str=='a'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='A')return(i);
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='C')return(i);
  }else if(sec_str=='b'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='B')return(i);
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='e')return(i);
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='E')return(i);
  }else if(sec_str=='e'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='E')return(i);
  }else if(sec_str=='y'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='Y')return(i);
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='e')return(i);
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='E')return(i);
  }else if(sec_str=='z'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='Z')return(i);
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='C')return(i);

    // Other types
  }else if(sec_str=='T'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='C')return(i);
  }else if(sec_str=='S'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='C')return(i);
  }else if(sec_str=='W'){
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='y')return(i);
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='Y')return(i);
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='e')return(i);
    for(i=0; i<N_SS; i++)if(SEC_EL[i]=='E')return(i);
  }

  printf("WARNING, sec.str=%c not defined", sec_str); //exit(8);
  printf(" allowed: %s. ", SEC_EL);
  return(0);
}

/*
  Compute energy
*/
void Initialize_E_loc(float T, float A_LOC){

  Combine_propensity(0, 1);  // T->C
  Combine_propensity(0, 2);  // S->C
  //Combine_propensity(0, 15); // Z->C
  //Combine_propensity(8, 7);  // 4->H

  int ia, is;

  if(E_loc_over_T)Empty_matrix_f(E_loc_over_T, N_SS);
  E_loc_over_T=malloc(N_SS*sizeof(float *));
  for(ia=0; ia<20; ia++)E_loc_0[ia]=0;
  float c=A_LOC/T;
  for(is=0; is<N_SS; is++){
    float *E=malloc(20*sizeof(float));
    E_loc_over_T[is]=E;
    for(ia=0; ia<20; ia++){
      E[ia]=c*E_loc[is][ia];
      E_loc_0[ia]+=f_SS[is]*exp(-E[ia]);
    }
  }
  for(ia=0; ia<20; ia++)E_loc_0[ia]=-log(E_loc_0[ia]);

  for(is=0; is<N_SS; is++){
    for(ia=0; ia<20; ia++){
      E_loc_over_T[is][ia]-=E_loc_0[ia];
    }
  }

  char nameout[80]="E_loc_0.txt";
  FILE *file_out=fopen(nameout, "r");

  if(file_out==NULL){
    file_out=fopen(nameout, "w");
    printf("Writing %s\n", nameout);
    fprintf(file_out, "# A_LOC= %.3f T= %.2f\n", A_LOC, T);
    for(ia=0; ia<20; ia++)fprintf(file_out, "   %c   ", Amin_code(ia));
    fprintf(file_out, "\n");
    for(ia=0; ia<20; ia++)fprintf(file_out, "%.2g ", E_loc_0[ia]);
    fprintf(file_out, "\n");
    fclose(file_out);
  }

  if(0){
    printf("Writing local interactions\n");
    for(is=0; is<N_SS; is++){
      printf("%c: ", SEC_EL[is]);
      for(ia=0; ia<20; ia++)printf(" %.4f", E_loc_over_T[is][ia]);
      printf("\n");
    }
  }

}

double Compute_E_loc(int *i_sec, short *aa_seq, int L){
  double E=0;
  for(int i=0; i<L; i++){
    if((aa_seq[i]>=0)&&(aa_seq[i]<20)){
      E+=E_loc_over_T[i_sec[i]][aa_seq[i]];
    }
  }
  if(1){
    int pos[N_SS], np=0; float ee[N_SS], ns[N_SS];
    char pos_res[N_SS][L];
    for(int i=0; i<N_SS; i++){pos[i]=0; ee[i]=0; ns[i]=0;}
    for(int i=0; i<L; i++){
      if((aa_seq[i]>=0)&&(aa_seq[i]<20)){
	int is=i_sec[i];
	float e=E_loc_over_T[is][aa_seq[i]];
	ee[is]+=e; ns[is]++;
	if(e<0){continue;}
	pos_res[is][pos[is]]=Amin_code(aa_seq[i]);
	pos[is]++; 
	/*printf(" %.2f (%c %c %d)\n",
	  e, SEC_EL[is], Amin_code(aa_seq[i]), i);*/
      }
    }
    if(print_loc){
      print_loc=0;
      FILE *file_out=fopen("Local_interactions.dat", "w");
      fprintf(file_out, "type pos E\n");
      for(int i=0; i<N_SS; i++){
	if(ns[i]){
	  np+=pos[i];
	fprintf(file_out, "%c %.3f %.3f\t",
		SEC_EL[i], pos[i]/ns[i], ee[i]/ns[i]);
	for(int j=0; j<pos[i]; j++)
	  fprintf(file_out, " %c", pos_res[i][j]);
	fprintf(file_out, "\n");
      }
      }
      fprintf(file_out, "E_loc= %.1f %.3f positive values\n",
	      E, np/(float)L);
      fclose(file_out); //exit(8);
    }
  }
  return(E);
}

double Mutate_E_loc_misf(int aa_old, int aa_new){
  if(E_loc_misf==NULL)E_loc_misf=Set_E_loc_misf();
  return(E_loc_misf[aa_new]-E_loc_misf[aa_old]);
}

double *Set_E_loc_misf(){
  E_loc_misf=malloc(21*sizeof(double));
  /*// Boltzmann average of local interactions
  for(int a=0; a<20; a++){
    double E_misf=0, Z=0;
    for(int j=0; j<N_SS; j++){
      float E=E_loc_over_T[j][a], w=exp(-E);
      E_misf+=w*E; Z+=w;
    }
    E_loc_misf[a]=E_misf/Z;
    }*/
  for(int a=0; a<20; a++){
    E_loc_misf[a]=E_loc_over_T[0][a]; // iss=0 => C structure
  }
  E_loc_misf[20]=0;
  return(E_loc_misf);
}

double Compute_E_loc_misf(short *aa_seq, int L){


  if(E_loc_misf==NULL)E_loc_misf=Set_E_loc_misf();

  double E=0;
  for(int i=0; i<L; i++){
    if((aa_seq[i]>=0)&&(aa_seq[i]<20)){
	E+=E_loc_misf[aa_seq[i]];
    }
  }
  return(E);
}

double Delta_E_loc(int i_sec, short aa_old, short aa_new){
  return(E_loc_over_T[i_sec][aa_new]-E_loc_over_T[i_sec][aa_old]);
}

void Combine_propensity(int i1, int i2){
  float f1=f_SS[i1], f2=f_SS[i2], f12=f1+f2;
  double *E1=E_loc[i1], *E2=E_loc[i2];
  for(int a=0; a<20; a++){
    double y=(f1*exp(-E1[a])+f2*exp(-E2[a]))/f12;
    E1[a]=-log(y);
  }
  f_SS[i1]=f12; f_SS[i2]=0;
}

/*
  Read_propensities
*/

float **Read_propensity(char *name_in)
{
  /*
    NEW format of the file:
    List of AA ### AA_1 ... AA_20
    Data: N_SS lines for each SS; 20 columns, L(SS,AA)
  */

  FILE *file_in=fopen(name_in, "r");
  if(file_in==NULL){
    printf("ERROR, cannot find file %s with propensities\n", name_in);
    exit(8);
  }

  char string[1000];
  N_SS=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]!='#')N_SS++;
  }
  fclose(file_in);
  printf("Reading %d sec str types in %s\n", N_SS, name_in);

  float **Propensity=Allocate_mat2_f(N_SS, 20);
  char *s;  float loc;
  int i_aa, i, n_ss=0, l_aa[20];
  file_in=fopen(name_in, "r");
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'){
      if((string[1]!='#')||(string[2]!='#'))continue;
      s=string+4;
      for(i=0; i<20; i++){
	while((*s)==' '){s++;} l_aa[i]=Code_AA(*s); s++;
      }
    }else{
      s=string+3; while(*s!='(')s++; SEC_EL[n_ss]=*(s+1); s+=3;
      for(i=0; i<20; i++){
	i_aa=l_aa[i];
	loc=Read_column(&s, string, i_aa);
	Propensity[n_ss][i_aa]=loc; //exp(-loc);
      }
      n_ss++;
    }
    if(n_ss==N_SS)break;
  }
  fclose(file_in);
  for(i=0; i<N_SS; i++){printf(" %c", SEC_EL[i]);} printf("\n\n");

  /*
    printf("### ");
    for(i_aa=0; i_aa<20; i_aa++)printf("        %c", Amin_code(i_aa));
    printf("\n");
    for(i=0; i<N_SS; i++){
    printf("Loc(%c) ", SEC_EL[i]);
    for(i_aa=0; i_aa<20; i_aa++)printf(" %7.4f", -log(Propensity[i][i_aa]));
    printf("\n");
  }
  */
  return(Propensity);
}
