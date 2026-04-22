#define NUC_CHAR "ATGC"
#define AMIN_CODE "AEQDNLGKSVRTPIMFYCWH"

int Compute_rates(float *rate, float *freq_nuc, float trans_ratio);
int Get_name(char *name, char *name_seq, int N);
int Mutate_seq(char *dna_seq, char **codon, char *coded_aa,
	       short *aa_seq, short *aa_seq_new,
	       int *nuc_mut, char *nuc_new, int *res_mut, int *aa_new);
