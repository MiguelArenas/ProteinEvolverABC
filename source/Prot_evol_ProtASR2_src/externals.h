// Mutation parameters
//extern char coded_aa[64], *codon[64]; 
#define MUTPAR 10
extern float mut_par[MUTPAR];
extern int NPAR;
// contact definition
//extern char cont_type;
extern float cont_thr, cont_thr2;
// contact statistics
extern int DTAIL_MAX;
extern int LEN_MAX;
extern int LEN;
extern float *Cont_L;
extern float *Cont_1L;
extern float *Cont_2L;
extern float *Cont_3L;
extern float *C2_L;
extern float *C3_L;
extern float *Cnc1_L;
extern float **nc_nc_L;
extern float **nc_nc_Nc_L;
extern float NC1, NC2;
extern float *nc_L, *nc2_L;
//extern double nc2_sum, nc12_sum;
//extern double *nc1i=NULL, *nc2i=NULL;
extern int tail(int i, int L);

extern char FILE_STR[200];
// Contact energy function
extern char AA_code[21];
//extern float Econt[21][21];
extern float **Econt_T, T_ratio;
extern float **E_loc_over_T;
extern float hydro[21];
extern float SEC_STR;
//extern char SEC_EL[16];
extern float TEMP;
extern float sC1, sC0, SC, s0; //, Conf_entropy, K2Thr
// general
extern int Verbose;
// AA distr
extern float Entr_ave, Entr_reg;
extern int *i_sec;

// Disulfide bonds
extern int N_disulf, *Disulf_res1, *Disulf_res2;
extern float *Disulf_d;

extern int N_disulf_seq;
extern int Len_disulf_seq;
extern float Log_len_disulf_seq;

// Secondary structure
extern int N_helix, N_strand;
extern int N_helix_seq, N_strand_seq;
