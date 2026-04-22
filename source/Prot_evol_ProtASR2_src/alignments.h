char **Read_MSA(int *N, int *L, char ***name_seq, int **selected,
		char *file_ali, float thr);
int Find_seq(int *ali_seq, int *seq_L, float *seq_id,
	     char *seq_PDB, int len, char **MSA, int n_seq, int L_ali);
