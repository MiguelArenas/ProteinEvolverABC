int Read_name_str(char *STR_MUT_TYPE, char *file_str_mut);
int Read_str_mut(float **Str_mut, char *file_str_mut,
		 struct protein target, int *res_index, int IWT);
int Read_str_mut_all(float ***Df, char *file_str_mut,
		     struct protein target, int *res_index);
int Mean_DDG(float **Str_mut, float ***Df, float **P_ia, int len_amm);
void WT_DDG(float **exp_ia, float ***DD, short *aa_seq, int len_amm);
