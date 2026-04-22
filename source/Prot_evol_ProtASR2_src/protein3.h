#define IJ_MIN 3  // Minimum value of |i-j| for contacts 
#define N_PROT_MAX 10000
#define N_PARAM 211
#define N_STR 100

/******************************

 *          Functions         *
 
 ******************************/

extern int Threading(short *, float *, float *, float, float *, int);
extern int Ini_thread(char *, int);

/*****************************

 *         Variables         *
 
 *****************************/

struct contact{
  short res1;
  short res2;
};

struct inter_chain{
  int res1;
  int res2;
  int ichain1;
  int ichain2;
  char chain1;
  char chain2;
};

struct protein{
  char name[120];
  // based on conformation
  int L_PDB;  // number of residues with PDB coord
  int nchain; // stored chains
  int num_chain; // each stored chain consists of nchain_num chains
  int chain_all; // nchain *num_chain
  int *ini_chain;
  int *len_chain;
  short *aa_seq;
  char *sec_str;
  int *i_sec;
  // based on seqres
  int length; // length of unique sequence
  char *seqres;
  int nc_seqres;
  int *ini_seqres;
  int *len_seqres;
  // Correspondence
  int *res_index;
  // contacts
  int n_cont;
  short **contact;
  struct contact *cont_list;
  struct inter_chain *inter_chain;
  int N_inter_chain;
};

struct state{
  short first;
  float energy;
  float overlap;
  float alpha;
  short n_cont;
  short **cont;
  struct protein *prot_ptr;
};

extern struct protein prot[N_PROT_MAX], target;
//extern float **interactions;









