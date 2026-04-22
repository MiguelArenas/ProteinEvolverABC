struct MF_results{
  float Lambda[2];  // Selection parameters
  float score;   // score
  float KL_mut;  // KL from mutational distribution
  float KL_reg;  // KL divergence of model from regularized frequencies
  float KL_mod;  // KL divergence of reg.freq.from model
  float lik_reg; // Likelihood of reg frequencies wrt model
  float lik_MSA; // Likelihood of MSA wrt model
  float DG;      // Average DeltaG over sequence distribution 
  float Tf;      // Average Tf (?)
  float h;       // Protein average of site-specific average hydro
  float entropy; // Protein average of site-specific average entropy
  float rate;    // Protein average of site-specific average rate
  int conv;
};

//extern float Entr_ave, Entr_reg;

int Fixed_Lambda(struct MF_results *MF_res, float **P_MF_ia, float *P_mut_a,
		 float **C_nat, int *i_sec, char *c_sec,
		 float **f_reg_ia, float **f_msa_ia,
		 float *wi, struct REM E_wt, float Lambda, char *name_file,
		 int L, int Naa);
int Optimize_Lambda(struct MF_results *MF_res,
		    float **P_MF_ia, float *P_mut_a,
		    float **C_nat, int *i_sec, char *c_sec,
		    float **f_reg_ia, float **f_msa_ia, float *wi,
		    short *aa_seq, struct REM E_wt,
		    float DG_OPT, int GET_FREQ, char *MODEL,
		    char *name_file, FILE *file_summary, int L, int Naa);
int meanfield(float **P_MF_ia, struct MF_results *res,
	      float *P_mut_a, float **P_ini_ia,
	      float **C_nat, int *i_sec, char *c_sec,
	      float Lambda, struct REM E_wt,
	      int L, int Naa, int ALL);
void Divide_by_Pmut(float **P_MF_ia, float *P_mut_a, int L, int Naa);
float Compute_lik(float **P_MF_ia, float **freq_ia, int L, int Naa);
void Compute_score(struct MF_results *MF_res, float **P_MF_ia,
		   float **f_reg_ia, float **f_msa_ia, float *wi,
		   int L, int Naa);
void Copy_P(float **P1_ia, float **P2_ia, int L, int Naa);
void Compute_P_sel(float *P_sel, float **P_MF_ia, float *P_mut_a,
		   int L, int Naa);
float Find_max(float *y, float *x, int n, float MIN, float MAX);
void Test_distr(struct MF_results *MF_res, float **P_MF_ia,
		float **f_reg_ia, float **f_msa_ia, float *wi,
		float **C_nat, int *i_sec, char *c_sec,
		struct REM E_wt, int L, int Naa);
float Average_entropy(float **P_MF_ia, int L, int Naa);
float Hydro_ave(float **P_ia, float *hydro, int L);
void Print_results(struct MF_results res, char *what, FILE *file_out);
float KL(float *P, float *Q, int n); // Kullback-Leibler divergence
