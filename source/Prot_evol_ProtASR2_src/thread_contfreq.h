int Fill_C_nat(int length, short **contact);
float Compute_DG_REM3(short *aa_seq, int L, double *E_nat,
		      double *E1_misf, double *E2_misf,
		      double *E23_misf, double *E3_misf,
		      char *file_str);
float Mutate_DG_REM3(short *aa_seq, int L, int res_mut, short aa_new,
			 double *E_nat, double *E1, double *E2,
			 double *E23, double *E3);
float Compute_DG_REM2(short *aa_seq, int L, double *E_nat,
		      double *E1_misf, double *E2_misf, char *file_str);
float Mutate_DG_REM2(short *aa_seq, int L, int res_mut, short aa_new,
		     double *E_nat, double *E1, double *E2);
int Get_target(char *file_str, char *name_tar, int *len_tar);
float Print_DG_REM3(double E_nat, double E1, double E2,
		   double E23, double E3, char *name);
float Print_DG_REM2(double E_nat, double E1, double E2, char *name);

extern float TEMP;
extern float sC1, sC0, SC, s0, Conf_entropy, K2Thr, TEMP2;
extern int LEN;

