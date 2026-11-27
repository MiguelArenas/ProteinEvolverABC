struct mutation{
  int nmut;
  char *amm1_mut;
  char *amm2_mut;
  int *rmut;
  int *imut;
  char *chain;
  char name[100];
  int interface;
};
struct deletion{
  int ndel;
  int Ldel;
  char *aa1; // Amino acid at first deleted position
  char *aa2;
  int *res1; // PDB numbering of residues
  int *res2;
  int *i1;
  int *i2;  // Order in PDB file
  char *chain;
  char name[100];
  int interface;
};
int Read_mut(struct mutation **mutations, char *FILE_MUT,
	     struct deletion **deletions, int *Ndel);
int Construct_mut(short *aa_seq, struct mutation *mut,
		  struct residue *res, int L);
int Construct_del(short *aa_seq, struct deletion *del,
		  struct residue *res, int L);
