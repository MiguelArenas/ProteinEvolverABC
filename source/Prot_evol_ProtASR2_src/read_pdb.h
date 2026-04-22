int Get_pdb(struct protein *prot,
	    struct residue **res, struct res_short **res_short,
	    char *file, char *chain, int **res_index, int nmod);
int Count_models_PDB(char *file_pdb);
int Read_modres(char **res_exo, char **res_std, char *string, int *n_exo);
char Het_res(char *res_type_old, char **res_exo, char **res_std, int n_exo);
char Code_3_1(char *res);
