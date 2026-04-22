#include "coord.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "protein3.h"
#include "read_pdb.h"
#include "externals.h"

#define AMIN_CODE "AEQDNLGKSVRTPIMFYCWHX"
#define ASTPATH    "/data/ortizg/databases/astral_40/"
#define PDBPATH    "/data/ortizg/databases/pdb/"
#define PDBEXT     ".pdb" /*.ent.Z*/
#define PDBCAT     "zcat"  /* zcat */
#define PDBTMP     "pdb.dat"

static int n_atom;
static char *res_exo[400], *res_std[400];

static int Next_residue(int *N_res, int *start, struct residue *seq,
			atom *first_atom, short *i_atom,
			char *res_type_old, int *res_num_old, char *icode_old,
			char *chain_old, int *hetatm_old, int n_exo,
			char *res_type, int res_num, char icode, char chain,
			int hetatm);
static short Write_residue(char *res_type_old, atom *first_atom, short i_atom,
			   struct residue *pt, int nres, int n_exo, int res_num,
			   char icode, char chain, int hetatm);

static int Code_AA(char res);
int Get_compression(char *pdb_name);

//int N_ATOM_MAX=5000000; 

int Read_coord(char *pdb_name, int *nmr, struct residue *seq, atom *atoms,
	       char *chain_to_read, int *ANISOU, int *natoms, int kmod)
{
  int N_res=0, n_exo=0, Compression=0, start=0;
  FILE *file_in; int read_all=0;
  char string[200], command[300];

  atom *atom_ptr, atom_read[100], *first_atom=NULL;
  if(atoms){atom_ptr=atoms;}else{atom_ptr=atom_read;}
  int read_model=0; if(kmod==0)read_model=1;

  short i_atom=0, alternative=0;
  int hetatm=0, hetatm_old=0;
  int i, res_num, res_num_old=10000;
  char altloc, altloc_sel=' ';
  char chain_old='#', old_chain='#', chain='Z';
  char res_type[5], res_type_old[5], icode=0, icode_old;
  float x, y, z;
  char file_name[500];
  int nchain=1, ic, ichain=0, ini_chain=0;
  n_atom=0;

  /* Open file */
  //printf("Reading %s\n", pdb_name);
  Compression=Get_compression(pdb_name); 
  if(Compression){
    sprintf(command, "%s %s > %s\n", PDBCAT, pdb_name, PDBTMP);
    system(command); strcpy(file_name, PDBTMP);
  }else{
    sprintf(file_name, "%s", pdb_name);
  }
  file_in=fopen(file_name, "r");
  if(file_in==NULL){
    printf("WARNING, file %s not found\n", file_name); return(0);
  }
  //printf("Reading %s\n", file_name);

  // Count chains to read
  printf("Reading %s chains %s\n", pdb_name, chain_to_read);
  if(strncmp(chain_to_read,"ALL",3)==0){
    read_all=1; ichain=-1;
    printf("Reading all chains\n");
  }else{
    read_all=0; i=0;
    while((chain_to_read[i]!='\0')&&(chain_to_read[i]!=' ')){
      printf("%c", chain_to_read[i]); i++;
    }
    nchain=i;
    chain_to_read[i]='\0';
    if(nchain==0)nchain=1;
    printf(" %d chains to read: ", nchain);
    for(i=0; i<nchain; i++)printf("%c", chain_to_read[i]);
    printf("\n");
  }

  *nmr=0; int imod=-1;
  strcpy(res_type_old,"xxx");
  while(fgets(string, sizeof(string), file_in)!=NULL){

    if(strncmp(string,"EXPDTA", 6)==0){
      if(strncmp(string+10, "NMR", 3)==0)*nmr=1;
      continue;
      /* NMR structure */
    }else if(strncmp(string,"MODRES", 6)==0){
      Read_modres(res_exo, res_std, string, &n_exo);
      continue;
    }else if(strncmp(string,"MODEL", 5)==0){
      imod++; if(imod==kmod)read_model=1;
      continue;
    }else if(strncmp(string,"ENDMDL", 6)==0 && read_model){
      break;                     
    }else if(strncmp(string,"ATOM", 4)==0){
      /* Standard residue or DNA basis */
      if(read_model==0)continue;
      hetatm=0; if(ini_chain==0)ini_chain=1;
    }else if(strncmp(string,"HETATM", 6)==0){
      /* Cofactor or exotic residue */
      if(read_model==0)continue;
      if(ini_chain==0)continue;
      hetatm=1;
    }else if(strncmp(string,"ANISOU", 6)==0){
      // Anisotropic structure factor
      if(read_model==0)continue;
      if(atoms==NULL)continue;
      if(*ANISOU==0)*ANISOU=1;
      sscanf(string+28, "%f %f %f %f %f %f",
	     &(atom_ptr->anisou[0][0]), &(atom_ptr->anisou[1][1]),
	     &(atom_ptr->anisou[2][2]), &(atom_ptr->anisou[0][1]),
	     &(atom_ptr->anisou[0][2]), &(atom_ptr->anisou[1][2]));
      for(int i=0; i<3; i++){
	for(int j=i; j<3; j++){
	  atom_ptr->anisou[i][j]*=0.0001;
	  if(j!=i)atom_ptr->anisou[j][i]=atom_ptr->anisou[i][j];
	}
      }
      // Check
      //for(i=0; i<3; i++)B+=aniso[i][i]; B*=26.319; // 8pi^2/2
      //printf("B= %.2f %.2f %d\n", atom_ptr->B_factor, B, n_atom);
      continue;
    }else if((strncmp(string,"TER",3)==0)&&(N_res>0)){
      Next_residue(&N_res, &start, seq, first_atom, &i_atom,
		   res_type_old, &res_num_old, &icode_old,
		   &chain_old, &hetatm_old, n_exo,
		   res_type, res_num, icode, chain, hetatm);
      ini_chain=0;
      continue;
      /* end model */

      /*}else if(strncmp(string,"HELIX ", 6)==0){
	Read_sec_str(string, chain, 'H'); continue;
	}else if(strncmp(string,"SHEET ", 6)==0){
	Read_sec_str(string, chain, 'E'); continue;
	}else if(strncmp(string,"TURN ", 5)==0){
	Read_sec_str(string, chain, 'T'); continue;*/

    }else{
      continue;
    }

    chain=string[21];
    if(read_all){
       if(chain!=old_chain){
	 old_chain=chain; ichain++; chain_to_read[ichain]=chain;
      }
    }else{
     if((*chain_to_read==' ')||(*chain_to_read=='\0'))*chain_to_read=chain;
      for(ic=0; ic<nchain; ic++)
	if(chain==chain_to_read[ic]){ichain=ic; goto read;}
      continue;
    }

  read:
    /* Read atom name */
    if((string[13]=='H')||(string[12]=='H') ||
       (string[13]=='D') ||(string[12]=='D'))continue;

    /* Read residue; check if water molecule */
    res_type[0]=string[17]; res_type[1]=string[18]; res_type[2]=string[19];
    res_type[3]='\0';
    if((hetatm==1)&&((strncmp(res_type,"HOH",3)==0)||
		     (strncmp(res_type,"DOD",3)==0)))continue;
    
    icode=string[26]; string[26]=' ';


    /* Read coordinates */
    sscanf(string+22,"%d %f %f %f", &res_num, &x, &y, &z);
    
    /* Check if alternative conformation */
    if((string[72]=='A')&&(string[73]=='L')&&(string[74]=='T')&&
       (string[75]!='1')&&(string[75]!=' '))continue;

    altloc=string[16];
    if(altloc!=' '){
      if(altloc_sel==' ')altloc_sel=altloc;
      if(altloc!=altloc_sel)continue;
    }

    if((icode!=icode_old)&&(res_num==res_num_old)){
      if(alternative==1){
	continue;
      }else{
	atom *atom_old=first_atom;
	float dx, dy, dz;
	dx=x-atom_old->r[0]; dy=y-atom_old->r[1]; dz=z-atom_old->r[2];
	if((dx*dx+dy*dy+dz*dz)<.5){
	  alternative=1; continue;  
	}
      }         
    }

    /* New residue */
    if((res_num!=res_num_old)||(icode!=icode_old)||
       (strncmp(res_type, res_type_old, 3)!=0)){
      Next_residue(&N_res, &start, seq, first_atom, &i_atom,
		   res_type_old, &res_num_old, &icode_old,
		   &chain_old, &hetatm_old, n_exo,
		   res_type, res_num, icode, chain, hetatm);
    }

    if(atoms){
      atom_ptr=atoms+n_atom;
    }else{
      atom_ptr=atom_read+i_atom;
    }
    if(i_atom==0)first_atom=atom_ptr;

    atom_ptr->r[0] = x;
    atom_ptr->r[1] = y;
    atom_ptr->r[2] = z;
    if(string[12]!=' '){ // Hydrogen atoms
      for(ic=0; ic<4; ic++)atom_ptr->name[ic]=string[12+ic];
    }else{
      for(ic=0; ic<3; ic++)atom_ptr->name[ic]=string[13+ic];
      atom_ptr->name[3]=' ';
    }
    sscanf(string+56, "%f", &atom_ptr->occupancy);
    sscanf(string+60, "%f", &atom_ptr->B_factor);
    atom_ptr->chain=ichain;
    i_atom++; n_atom++;
  }
  Next_residue(&N_res, &start, seq, first_atom, &i_atom,
	       res_type_old, &res_num_old, &icode_old,
	       &chain_old, &hetatm_old, n_exo,
	       res_type, res_num, icode, chain, hetatm);
  if(N_res >= L_MAX){
    printf("\n ERROR, more than %d residues found\n", L_MAX); exit(8);
  }
  fclose(file_in);
  if(Verbose)printf("%3d residues\n", N_res);
  if(Compression){
    sprintf(command, "rm -f %s\n", PDBTMP); system(command);
  }
  *natoms=n_atom;
  if(read_all){nchain=ichain+1;}
  chain_to_read[nchain]='\0';
  //printf("%d chains: %s\n",nchain, chain_to_read);
  printf("%d residues\n", N_res);
  return(N_res);
}

char Het_res(char *res_type_old, char **res_exo, char **res_std, int n_exo)
{
  for(int i=n_exo-1; i>=0; i--){
      if(strncmp(res_type_old,res_exo[i],3)==0){
	return(Code_3_1(res_std[i]));
      }
  }
  return('X');
}

static short Write_residue(char *res_type_old, atom *first_atom, short i_atom,
			   struct residue *seq, int nres, int n_exo,
			   int res_num, char icode, char chain, int hetatm)
{
  short het=1, exo;
  char amm;

  /* Check amino acid type */
  amm=Code_3_1(res_type_old);
  if(amm != 'X'){
    // Standard residue
    if(hetatm){het=1; goto discard;} // Standard res. and HETATM => cofactor
    het=0; exo=0;
  }else{
    // Non-standard residue
    printf("res= %d (%d) het=%d %d atoms chain= %c\n",
	   res_num, nres, hetatm, i_atom, chain);
    exo=1;
    amm=Het_res(res_type_old, res_exo, res_std, n_exo);
    if(amm!='X')het=0;
  }

  /* Check backbone */
  if(het && (i_atom >=3)){
    // If backbone atoms exist: Modified residue
    int i_N=0, i_CA=0, i_C=0, i;
    atom *atom_ptr=first_atom;
    for(i=0; i<i_atom; i++){
      if(strncmp(atom_ptr->name, "N ", 2)==0){
	i_N=1;
      }else if(strncmp(atom_ptr->name, "CA", 2)==0){
	i_CA=1;
      }else if(strncmp(atom_ptr->name, "C ", 2)==0){
	i_C=1;
      }
      if(i_N && i_CA && i_C){het=0; break;}
      atom_ptr++;
    }
  }

  if(het==0){
    if(seq){
      struct residue *ptr_tmp=seq+nres;
      ptr_tmp->atom_ptr=first_atom;
      ptr_tmp->n_atom=i_atom;
      ptr_tmp->amm=amm; ptr_tmp->exo=exo;
      ptr_tmp->i_aa=Code_AA(amm);
      if(ptr_tmp->i_aa<0){
	printf("Unknown residue %s %d%c\n",
	       res_type_old, res_num, icode);
	ptr_tmp->i_aa=0;
      }
      ptr_tmp->chain=chain;
      sprintf(ptr_tmp->pdbres, "%4d%c", res_num, icode);
    }
    return(het);
  }

  //printf("%s %d%c  %d %d\n", res_type_old, res_num, icode, hetatm, het);

 discard:

  if(het)
    printf("Group %s %d %c %c ch. %c (%d atoms) not a residue hetatm=%d\n",
	   res_type_old, res_num, icode, amm, chain, i_atom, hetatm);

  return(het);
}

int Next_residue(int *N_res, int *start, struct residue *seq,
		 atom *first_atom, short *i_atom,
		 char *res_type_old, int *res_num_old, char *icode_old,
		 char *chain_old, int *hetatm_old, int n_exo,
		 char *res_type, int res_num, char icode, char chain,
		 int hetatm)
{
  if((*start)==0){
    *start=1;
  }else if(*i_atom){
    int het=Write_residue(res_type_old, first_atom, *i_atom, seq, *N_res, n_exo,
			  *res_num_old, *icode_old, *chain_old, *hetatm_old);
    if(het==0){(*N_res)++;}else{n_atom-=(*i_atom);} (*i_atom)=0;
  }
  strcpy(res_type_old,res_type); *icode_old=icode;
  *res_num_old=res_num; *chain_old=chain; *hetatm_old=hetatm;
  
  return(0);
}

char Code_3_1(char *res){

  if(strncmp(res,"ALA",3)==0){return('A');
  }else if(strncmp(res,"GLU",3)==0){return('E');
  }else if(strncmp(res,"GLN",3)==0){return('Q');
  }else if(strncmp(res,"ASP",3)==0){return('D');
  }else if(strncmp(res,"ASN",3)==0){return('N');
  }else if(strncmp(res,"LEU",3)==0){return('L');
  }else if(strncmp(res,"GLY",3)==0){return('G');
  }else if(strncmp(res,"LYS",3)==0){return('K');
  }else if(strncmp(res,"SER",3)==0){return('S');
  }else if(strncmp(res,"VAL",3)==0){return('V');
  }else if(strncmp(res,"ARG",3)==0){return('R');
  }else if(strncmp(res,"THR",3)==0){return('T');
  }else if(strncmp(res,"PRO",3)==0){return('P');
  }else if(strncmp(res,"ILE",3)==0){return('I');
  }else if(strncmp(res,"MET",3)==0){return('M');
  }else if(strncmp(res,"PHE",3)==0){return('F');
  }else if(strncmp(res,"TYR",3)==0){return('Y');
  }else if(strncmp(res,"CYS",3)==0){return('C');
  }else if(strncmp(res,"TRP",3)==0){return('W');
  }else if(strncmp(res,"HIS",3)==0){return('H');
  }else if(strncmp(res,"HIE",3)==0){return('H');
  }else if(strncmp(res,"HID",3)==0){return('H');
  }else if(strncmp(res,"HIP",3)==0){return('H');
  }else if(strncmp(res,"ASX",3)==0){return('N');
  }else if(strncmp(res,"GLX",3)==0){return('Q');
  }else{
    printf("WARNING, a.a. %c%c%c not known\n", *res,*(res+1),*(res+2));
    return('X');
  }
}

int Code_AA(char res){
  short i; char r;
  i=(int)res; if(i>96){r=(char)(i-32);}else{r=res;}
  for(i=0; i<20; i++)if(r==AMIN_CODE[i])return(i);
  if(res=='X')return(0);
  if((res!='-')&&(res!='.')&&(res!='*'))
    printf("Warning, wrong aa type %c\n", res);
  return(-1);
}

int Get_compression(char *pdb_name){
  char *tmp=pdb_name;
  while(*tmp!='\0'){
    if((*tmp=='.')&&(*(tmp+1)=='g')&&(*(tmp+2)=='z'))return(1);
    tmp++;
  }
  return(0);
}

int Count_models_PDB(char *pdb_name){

  int Compression=Get_compression(pdb_name);
  char string[200], command[200], file_name[500];
  if(Compression){
    sprintf(command, "%s %s > %s\n", PDBCAT, pdb_name, PDBTMP);
    system(command); strcpy(file_name, PDBTMP);
  }else{
    sprintf(file_name, "%s", pdb_name);
  }
  FILE *file_in=fopen(file_name, "r");
  if(file_in==NULL){
    printf("\nWARNING, file %s not found\n", file_name); return(-1);
  }
  int N_model=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string, "MODEL", 5)==0)N_model++;
  }
  fclose(file_in);
  if(N_model==0)N_model=1;
  return(N_model);
}

int Read_modres(char **res_exo, char **res_std, char *string, int *n_exo)
{
  res_exo[*n_exo]=malloc(3*sizeof(char));
  res_std[*n_exo]=malloc(3*sizeof(char));
  res_exo[*n_exo][0]=string[12]; res_std[*n_exo][0]=string[24];
  res_exo[*n_exo][1]=string[13]; res_std[*n_exo][1]=string[25];
  res_exo[*n_exo][2]=string[14]; res_std[*n_exo][2]=string[26];
  for(int j=0; j<*n_exo; j++){
    if(strncmp(res_exo[j],res_exo[*n_exo],3)==0)return(0);
  }
  (*n_exo)++; return(0);
} 
