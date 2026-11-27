
struct REM{
  float DeltaG;
  float G_misf;
  float G_unf;
  float Tf;
  double E_loc;
  double E_loc_misf;
  double E_nat;
  double *c1U1;
  double E1;
  double E2;
  double E2cont1, E2cont2;
  double E2site1, E2site2;
  double E3;
  double E3cont1, E3cont2, E3cont3;
  double E3site1, E3site2;
  double E22, E23;
  //double E321, E332, E342;
  //double e342, E343_2, E343_3;
  //double E353, E363;
  //double *c332U2;
  int L;     // Total length, including unstructured residues
  int L_str; // Only structured residues
  int REM;
  float T; // Temperature
  float S_C; // Entropy of compact conformations
  float S_U; // Entropy of unfolded conformations
};

void Compute_E_nat(struct REM *E,
		   short *aa_seq, float **C_nat,
		   int *i_sec, char *c_sec);
void Mutate_Enat(struct REM *E, short *aa_seq, float **C_nat,
		 int *is, int res_mut, short aa_new);
void Initialize_E_REM(struct REM *E, int L, int REM,
		      float TEMP, float S_C, float S_U,
		      char *file_str);
void Empty_E_REM(struct REM *E);
void Copy_E_REM(struct REM *E2, struct REM *E1);
float Mutate_DG_overT_contfreq(struct REM *E, short *aa_seq,
			      float **C_nat, int *i_sec, char *c_sec,
			      int res_mut, short aa_new);
float Compute_hydro(short *aa_seq, int L);

float Compute_DG_overT_contfreq(struct REM *E,
				short *aa_seq, float **C_nat,
				int *i_sec, char *c_sec,
				int Verbose);
float Compute_DG_overT_threading(struct REM *E,
				 short *aa_seq, float **C_nat,
				 int *i_sec, char *c_sec);
float Print_DG_contfreq(struct REM *E, char *name);
void Test_contfreq(struct REM *E, short *aa_seq, float **C_nat,
		   int *i_sec, char *c_sec, char *nameout);
float **Fill_C_nat(int length, short **contact, int num_chain);
float **Get_target(char *file_str, char *name_tar, int *len_tar);
float Compute_Tfreeze(struct REM *E);

void Set_to_zero(struct REM *E);
float G_misfold(struct REM *E);
float DeltaG(float E_nat, float G_misf, float S_U);

//extern int DTAIL_MAX;
//extern int LEN_MAX;
//extern int LEN;
//extern float *Cont_L;
//extern float *CNc1_L;
//extern float *CNc2_L;
//extern float *Cnc1_L;
//extern float **nc_nc_L;

//extern char AA_code[21];
//extern float Econt[21][21];
//extern float hydro[21];
//extern char FILE_STR[200];
//extern float **Econt_T, T_ratio;
//extern float SEC_STR;
//extern float **E_loc_over_T;

//float conf_entropy;
//extern double *Cont_freq, NC1, NC2, NC3;
//extern double *C221_ij, *C232_i, C242;
//extern double *C321_ij, *C342_ij, *C343_3_i, *C353_i, C363;
//extern double **C332_ij, **C343_22_ij;
//extern double *M221_ij, *M232_i, M242;
//extern double *M321_ij, *M3X2_ij, *M332_i, M3X3;
//extern double *nc1i, *nc2i, *CNc;
//extern double nc2_sum, nc12_sum;

