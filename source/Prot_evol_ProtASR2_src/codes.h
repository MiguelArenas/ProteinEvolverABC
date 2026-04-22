#define NUC_CODE "ATGC"
#define AMIN_CODE "AEQDNLGKSVRTPIMFYCWH-"

void Read_code(char *FILE_CODE, char *coded_aa, char **codon);
int Code_AA(char res), Code_nuc(char nuc);
int Codon_num(char *cod, char **codon);
char Amin_code(int), Nuc_code(int);
void Read_nuc_freq(char *file_nuc_freq, float *freq_nuc,
		   float *trans_ratio, float *mut_rate);
void Read_mut_mat(char *file_mut_mat, float **nuc_mut_mat);
int Transition(char nuc);
int Coded_aa(char *i_codon, char **codon, char *coded_aa);
char *Extract_dna(int *len_dna, int len_amm, short *aa_seq,
		  char **codon, char *coded_aa);
int  Compare_amm_dna(char *dna_seq, int len_dna, short *aa_seq, int len_amm,
		     char **codon, char *coded_aa);
int Translate_new(char *dna_seq, short *aa_seq, int length,
		  char **codon, char *coded_aa);
