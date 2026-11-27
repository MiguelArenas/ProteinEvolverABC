void Ini_count(char *dna_seq, int len_dna, int *count);
int Mutate_seq(char *dna_seq, int len_dna, char **codon, char *coded_aa,
	       short *aa_seq, int len_nat,
	       float *freq_nuc, float tt_ratio, int *count, float *rate, 
	       int *nuc_mut, char *nuc_new, int *res_mut, int *aa_new);
int Compute_rates(float *rate, float *freq_nuc, float trans_ratio);
