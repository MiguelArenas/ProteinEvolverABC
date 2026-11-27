struct load{
  float df;
  float dG;
  double df_ave;
  double df_dev;
  double dG_ave;
  double dG_dev;
  int num;
};
int Substitute_exhaustive(int *Res_mut, int *AA_new,
			  int *Nuc_mut, int *Nuc_new,
			  struct load *mut_load,
			  struct load *trans_load,
			  int it_load, double *fitness_wt,
			  struct REM *E_mut, struct REM *E_wt,
			  long *naa_mut, long *nsyn_mut, long *syn_subst,
			  int NEUTRAL, float DG_thr, int N_pop,
			  float **C_nat, int *i_sec, char *c_sec,
			  short *aa_seq,
			  char *dna_seq, short *nuc_seq, int len_dna, 
			  char **codon, char *coded_aa,
			  float tt_ratio, float *mut_par);
