//void Compute_exchange_mut(float **exchange, float **rate,
//			  float *P_cod, float **Q_cod);
void Get_mut_par(float *mut_par, float *P_mut, float *P_cod, float **Q_cod,
		 int GET_FREQ, float *num_aa, int L_aa, float *num_dna, 
		 char *name_file, int iter);
void Normalize_distr(float *P, int N);
int Code_codon(char *i_codon, char **codon);
float Weight_codon(char *c, float *f);
float Weight_codon_CpG(char *c, float *f);
int Sum_aa(float *num_aa, char *file_ali);
int Translate_aa(char **aa_trans, char *dna_seq,
		 int len_dna, char **codon, char *coded_aa);
int Translate_char(char *aa_trans, int len_amm, char *dna_seq, int len_dna,
		   char **codon, char *coded_aa);
int Compare_amm(char *aa_trans, int L, char *aa_seq, int len_amm);

extern char coded_aa[64], *codon[64]; 
//#define MUTPAR 10
//extern float mut_par[MUTPAR];
//extern int NPAR;
