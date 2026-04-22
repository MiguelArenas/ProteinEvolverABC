
#define L_MAX 40000       /* Maximal length of a chain */


typedef struct {
  //float x,y,z;
  float r[3];
  float B_factor;
  char name[5];
  int last_rot_axe; // last axis counted from protein start
  int res;
  char aa;
  int i_next, i_num;
  float anisou[3][3];
  float occupancy;
  int chain;
} atom;

struct axe{
  float v[3];                    // axe versor
  float offset[3];               // orign of the axis
  double global_shift[3];        // glob_shift[k]=<r-offset[k]>_i>k X axe[k]
  double global_rot[3];          // 
  double local_rot[3];
  double local_shift[3];         // glob_shift-
  float shift[3];
  char type;
  int res;
  int atom1, atom2;              // atoms that define the axe
  int first_rot_atom;            // first and last atoms moved by the axe
  int last_rot_atom;
  int first_rot_kin;             // first reference atom moved by the axe
  int last_rot_kin;
  int chain;
};

struct interaction{
  int i1, i2; // atom2 > atom1
  float r2, rmin, sec_der;
  int type;
};

struct residue{
  atom *atom_ptr;
  int n_atom;
  char chain;
  char pdbres[6];
  char amm;
  short i_aa;
  short exo;
  char sec_str;
  char c_sec;
  int i_sec;
  short n_cont;
  float PE;
  float asa;
  float B_factor;
  float B_NM;
};

struct res_short{
  char seq;
  int str_index;
  int aa;
  char pdbres[6];
  int n_cont;
  char sec_str;
};

struct chain{
  int ini_atom;
  int natoms;
  int ini_axe;
  int mainaxes;
  int ini_side;
  int nsides;
  int ini_res;
  int nres;
  int nref;
  int ini_cart;
  int ncart;
  char label;
  int *alignres;
};

#define FILE_MSG   "message.out"


//extern int Verbose;
//extern char cont_type;
//extern float cont_thr, cont_thr2;


short **Compute_map(int N_res, struct residue *seq, int *N_cont);
void Average_B(struct residue *seq, int N);
int Read_coord(char *pdb_name, int *nmr, struct residue *seq, atom *atoms,
	       char *chain_to_read, int *ANISOU, int *natoms, int kmod);
