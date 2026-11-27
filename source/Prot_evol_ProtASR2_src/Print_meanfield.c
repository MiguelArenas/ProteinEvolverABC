#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "coord.h"
#include "Print_meanfield.h"
#include "codes.h"
#include "output.h"
#include "gen_code.h"
#include "allocate.h"
//#include "wag.h"
//#include "jtt.h"
//#include "lg.h"
#include "REM.h"
#include "externals.h"
//

extern int PRINT_ALL_EXCH;
extern float LG_matrix[20][20], f_LG[20];
extern float JTT[20][20], fJTT[20];
extern float WAG[20][20], fWAG[20];
extern char AA_WAG[20];

int AVE_FLUX=0; // 1: Flux with average frequencies equal to average flux
                // 0: Rate with average frequencies equal to average rate

char AA_string[]="ACDEFGHIKLMNPQRSTVWY";
char AA_PAML[]="ARNDCQEGHILKMFPSTWYV";
//#define AMIN_CODE "AEQDNLGKSVRTPIMFYCWHX"

float hscale[]={0.1366,-0.0484,0.0325,-0.1233,-0.0345,0.4251,-0.0464,-0.0101,-0.0432,0.4083,0.0363,0.0589,0.0019,0.4172,0.1747,0.4076,0.3167,0.2745,0.2362,0.0549};

static void Change_AA_order(int *iaa, char *AA);
//static void Compute_Q_sel(float **Q_sel, float *P_MF);
static void Print_matrix(float **exchange, int *iaa, float norm,
			 FILE *file_mat, char *AA_string,
			 float *p, int FORMAT);
void Compute_rate_matrix(float **rate_matrix, float *P_aa,
			 float *P_cod, float **Q_cod);
void Compute_exchange_mut_sym(float **S_mut, float **rate_matrix, float *P_aa);
float **Empirical_exchangeability(float *f_emp, char *MATRIX);
float **Rescale_exchangeability(float **exch_emp, float *f_emp,
				char EXCHANGE, float **Pair_MF, float *F_MF);
float **Exch_flux_abs(float **exch_emp, float *f_emp,
		      float **Pair_MF, float *F_MF);
void Normalize_exchange(float **exch, float *f, int n);
// Output: site-specific rate and exchangeability matrix (may be NULL)
// Either compute exch_HB_flux: pointer to it, exchange_glob, f_glob
// Or exch_HB_flux=NULL, exchange_glob=exch_HB_flux, f_glob=NULL
void Exch_Halpern(float *rate_HB, float **exchange_HB_i,
		  float **exchange_HB_flux, // output
		  float **P_MF_ia, float *P_mut, int L, // Input
		  float **exchange_glob, float *f_glob);
float **Pair_matrix(float *F_MF, float **P_MF_ia, int L, float *P_mut_a,
		    int *iaa, char *name_file, char *AA_string, int FORMAT);


/*static void Compute_exchange_mut_old(float **S_mut,
				     float *mut_par,
				     float tt_ratio,float TWONUC); */
/* static void Update_F(float *F_mut, float *S, char *cod2, float *w,
		     int a, int c, int j, int n, int nt, float tt_ratio); */

double Mean_hydro(float *P, float *hydro, int n);
void Sum_matrix(float *ncont, int **Cnat, int L);
float Corr_coeff(float *xx, float *yy, int n);
float Compute_rate_abs(float *p, float **exchange, int n);
float Compute_rate_Q(float *p, float **rate, int n);
float Compute_rate_sel(float *p, float *p_mut, float **rate, int n);
float Corr_vM(float *p, float **exch, int n);

/******************* Codes ****************************/
int Print_profiles(float **P_MF_ia, char *tag, double DG_ave, float Lambda,
		   float *P_mut_a, short *aa_seq, int L, char *MATRIX,
		   char *nameout, int PRINT_GLOB, float *wi, FILE *out)
{
  int iaa[20], i, a;
  Change_AA_order(iaa, AA_PAML); //, AA_string

  char name_file[100]; sprintf(name_file, "%s_%s", nameout, tag);
  char name[200]; sprintf(name, "%s_AA_profiles.txt", name_file);
  FILE *file_out=Output_file(name_file, "AA_profiles", "txt");
  fprintf(out, "Printing site-specific amino acid frequencies in %s\n", name); 
 
  fprintf(file_out, "# Predicted site specific amino acid frequencies ");
  fprintf(file_out, "from model %s, Lambda= %.3f\n", tag, Lambda);
  fprintf(file_out, "# ave(DeltaG/T)= %.2f\n", DG_ave);
  fprintf(file_out, "# Flux model: %s\n", MATRIX);
  // corresponding to Npop= %.3g, Lambda*exp(-DG_ave)

  // Print amino acid names
  fprintf(file_out, "#pos AA_PDB");
  fprintf(file_out, "    %c", AA_PAML[0]); //AA_string[0]
  for(a=1; a<20; a++)fprintf(file_out, "     %c", AA_PAML[a]); //AA_string[a]
  fprintf(file_out, " Entropy\n");
  // Print mutational distribution
  fprintf(file_out, "#P_MUT P=");
  double S=0;
  for(a=0; a<20; a++){
    int ia=iaa[a]; float p=P_mut_a[ia]; if(p)S+=p*log(p);
    fprintf(file_out, " %.3f", p);
  }
  fprintf(file_out, "   %.2f\n", -S);

  // Print mean-field distribution
  int ii=1;
  for(i=0; i<L; i++){
    if(wi[i]==0)continue;
    fprintf(file_out, "%3d %c ", ii, Amin_code(aa_seq[i])); ii++;
    S=0;
    for(a=0; a<20; a++){
      int ia=iaa[a]; float p=P_MF_ia[i][ia]; if(p)S+=p*log(p);
      fprintf(file_out, " %.3f", p);
    }
    fprintf(file_out, "   %.2f\n", -S);
  }
  fclose(file_out);

  // Print global profile
  if(PRINT_GLOB){
    sprintf(name, "%s_AA_profile_global.txt", name_file);
    fprintf(out, "Printing global amino acid frequencies in %s\n", name);
    file_out=Output_file(name_file, "AA_profile_global", "txt");
    fprintf(file_out, "#AA P_obs P_MF P_mut diff_MF diff_mut\n");
    double Psum[20], Pobs[20]; int ll=L;
    for(a=0; a<20; a++){Psum[a]=0; Pobs[a]=0;}
    for(i=0; i<L; i++){
      if(aa_seq[i]>=0 && aa_seq[i]<20){Pobs[aa_seq[i]]++;}
      else{ll--;}
      for(a=0; a<20; a++)Psum[a]+=P_MF_ia[i][a];
    }
    for(a=0; a<20; a++){Psum[a]/=L; Pobs[a]/=ll;}
    for(a=0; a<20; a++){
      int ia=iaa[a];
      fprintf(file_out, "%c %.2f %.2f %.2f  %5.2f %5.2f\n", 
	      AA_string[a], Pobs[ia]*100, Psum[ia]*100, P_mut_a[ia]*100,
	      (Pobs[ia]-Psum[ia])*100, (Pobs[ia]-P_mut_a[ia])*100);
    }
  }

  return(0);
}

int Print_profile_evo(char *name, double **P_ia, short *aa_seq, int L,
		      double DG_ave, long it_sum)
{
  FILE *file_out=Output_file(name, "AA_profiles_evo", "txt");
  int *iaa=malloc(20*sizeof(int)), i, a;
  Change_AA_order(iaa, AA_string);

  fprintf(file_out, "# Evolutionary distribution, %ld substitutions\n",
	  it_sum);
  fprintf(file_out, "# ave(DeltaG/T)= %.2f\n", DG_ave/it_sum);
  fprintf(file_out, "#pos AA");
  for(a=0; a<20; a++)fprintf(file_out, "    %c", AA_string[a]);
  fprintf(file_out, " Entropy\n");

  // Print distribution
  for(i=0; i<L; i++){
    fprintf(file_out, "%3d %c ", i+1, Amin_code(aa_seq[i]));
    double S=0, norm=0;
    for(a=0; a<20; a++)norm+=P_ia[i][a];
    for(a=0; a<20; a++){
      int ia=iaa[a]; float p=P_ia[i][ia]/norm; if(p)S+=p*log(p);
      fprintf(file_out, " %.3f", p);
    }
    fprintf(file_out, "   %.2f\n", -S);
  }
  fclose(file_out);
  free(iaa);
  return(0);
}

/* int Print_subst_rate(float **P_MF_ia, float *P_mut_a,
		     float *P_cod, float **Q_cod,
		     float tt_ratio, short *aa_seq,
		     int **Cnat, char *sec_str, int L,
		     char *name_file, float TWONUC)
{
  int iaa[20], i;
  Change_AA_order(iaa, AA_string);
  float **exchange=Allocate_mat2_f(20, 20);
  float ncont[L];
  Sum_matrix(ncont, Cnat, L);
  float entropy[L], Rate[L], hydro_i[L], hydro_ave[L];

  FILE *file_out=Output_file(name_file, "rate_profile", "dat");
  fprintf(file_out, "#SITE AA sec.str. ncont rate entropy hydro ave_hydro\n");
  float P_mut[20], **rate_mut=Allocate_mat2_f(20, 20);
  Compute_rate_matrix(rate_mut, P_mut, P_cod, Q_cod);
  Compute_exchange_mut_sym(exchange, rate_mut, P_mut);
  for(i=0; i<L; i++){
    float *p=P_MF_ia[i];
    Rate[i]=Compute_rate_sel(p, P_mut, exchange, 20);
    entropy[i]=Entropy(P_MF_ia[i], 20);
    hydro_i[i]=hscale[aa_seq[i]];
    hydro_ave[i]=Mean_hydro(P_MF_ia[i], hscale, 20);
    fprintf(file_out, "%3d %c %c %2.0f   %5.3f  %5.2f %5.2f %5.2f\n",
	    i+1, Amin_code(aa_seq[i]), sec_str[i],  ncont[i], Rate[i],
	    entropy[i], hydro_i[i], hydro_ave[i]);
  }
  float r=Corr_coeff(ncont, hydro_i, L);
  fprintf(file_out, "# Corr(ncont, hydro)=     %6.3f %d sites\n", r, L);
  r=Corr_coeff(ncont, hydro_ave, L);
  fprintf(file_out, "# Corr(ncont, ave_hydro)= %6.3f %d sites\n", r, L);
  r=Corr_coeff(ncont, entropy, L);
  fprintf(file_out, "# Corr(ncont, entropy)=   %6.3f %d sites\n", r, L);
  r=Corr_coeff(ncont, Rate, L);
  fprintf(file_out, "# Corr(ncont, rate)=      %6.3f %d sites\n", r, L);
  r=Corr_coeff(entropy, Rate, L);
  fprintf(file_out, "# Corr(rate, entropy)=    %6.3f %d sites\n", r, L);

  fclose(file_out);
  Empty_matrix_f(exchange, 20);
  Empty_matrix_f(rate_mut, 20);
  return(0);
}
*/

int Print_exchange(float **P_MF_ia, char *TAG, struct res_short *res, int L,
		   float *P_mut_a, float *P_cod, float **Q_cod,
		   float tt_ratio, float TWONUC,
		   char *nameout, int FORMAT, char EXCHANGE, char *MATRIX,
		   int PRINT_EXCH, int PRINT_GLOB, int PRINT_MUT,
		   float *wi, FILE *out)
{
  // Amino acid order
  int iaa[20], i, a, b;
  if(FORMAT==0){Change_AA_order(iaa, AA_string);}
  else{Change_AA_order(iaa, AA_PAML);}
  char name_file[300]; sprintf(name_file, "%s_%s", nameout, TAG);

  // Type of exchangeability model
  char EXCH=EXCHANGE, exc_type[200], namexc[20], type[100];
 //float **exchange_global;
  if(EXCH=='M'){
    strcpy(namexc, "MUT");
    sprintf(type, "exchangeability_HB_MUT");
    sprintf(exc_type,
	    "# Exchangeability: genetic code, tt_ratio=%.2f TWONUC=%.3f\n",
	    tt_ratio, TWONUC); 
     //exchange_global=exchange_mut;
  }else if(EXCH=='E'){
    strcpy(namexc, "EXCHANGE");
    sprintf(type, "exchangeability_HB_%s_EMP", MATRIX);
    sprintf(exc_type, "# Exchangeability: %s\n", MATRIX);
    //exchange_global=exchange_emp;
  }else{
    if(EXCH!='F'){
      printf("WARNING, exchangeability model %c not defined\n", EXCH);
      EXCH='F';printf("Using default %c\n", EXCH);
    }
    sprintf(exc_type,"# Exchangeability such that Mean FLUX as %s model\n",
	    MATRIX);
    //exchange_global=exchange_HB_flux;
  }
  fprintf(out, "%s", exc_type);
  
  // Empirical model
  float f_emp[20], **exchange_emp=Empirical_exchangeability(f_emp, MATRIX);
  Normalize_exchange(exchange_emp, P_mut_a, 20);
  float Rate_HB_emp[L];

  // Global flux model
  // Selection matrix sum_i P_i[a]*P_i[b] for computing flux
  float F_MF[20], **Pair_MF=Pair_matrix(F_MF, P_MF_ia, L, P_mut_a, iaa,
					name_file, AA_string, FORMAT);
  float **exchange_flux=Exch_flux_abs(exchange_emp, f_emp, Pair_MF, F_MF);
  Normalize_exchange(exchange_flux, P_mut_a, 20);

  // Mutation model
  float R_mut[20], Rate_HB_mut[L]; for(i=0; i<L; i++)Rate_HB_mut[i]=0;
  float norm_rate;

  FILE *file_mat; char name[500];

  // FILE <>_rate_mut.dat
  // Print exchangeability matrices of mutation model
  if(PRINT_MUT){
    // Compute rates of mutation model
    sprintf(name, "%s_rate_mut.dat", nameout);
    fprintf(out, "Printing Mutational exchangeability ");
    fprintf(out, "matrices in %s\n", name);
    file_mat=fopen(name, "w"); printf("Writing %s\n", name);
    fprintf(file_mat,"#ab E^mut(a,b) ba  E^mut(b,a) E_flux(%s)\n", MATRIX);

    float **rate_mut=Allocate_mat2_f(20,20), P_mut[20];
    Compute_rate_matrix(rate_mut, P_mut, P_cod, Q_cod);

    float **exchange_mut=Allocate_mat2_f(20, 20);
    Compute_exchange_mut_sym(exchange_mut, rate_mut, P_mut);
    Normalize_exchange(exchange_mut, P_mut, 20);

    Exch_Halpern(Rate_HB_mut, NULL, NULL,
		 P_MF_ia, P_mut, L, exchange_mut, P_mut);

    // Halpern-Bruno rate obtained with the _mut exchangeability matrix
    for(a=0; a<20; a++){
      R_mut[a]=-rate_mut[a][a];
      for(b=a+1; b<20; b++){
	fprintf(file_mat, "%c%c %.2f %c%c %.2f %.2f\n", 
		Amin_code(a), Amin_code(b), rate_mut[a][b]/P_mut[b],
		Amin_code(b), Amin_code(a), rate_mut[b][a]/P_mut[a],
		exchange_flux[b][a]);
      }
    }
    fclose(file_mat);

    Empty_matrix_f(rate_mut, 20);
    Empty_matrix_f(exchange_mut, 20);
  }

  //// Flux model
  float Rate_HB_flux[L], **exchange_HB_flux=NULL;
  // Print exchangeability matrix in <>_exchangeability_sites_LG_FLUX.txt
  // Print rates in rate_profile.dat
  if(PRINT_EXCH){

    exchange_HB_flux=Allocate_mat2_f(20,20);
    Exch_Halpern(Rate_HB_emp, NULL, exchange_HB_flux,
		 P_MF_ia, P_mut_a, L, exchange_emp, f_emp);
    Normalize_exchange(exchange_HB_flux, P_mut_a, 20);
    norm_rate=Compute_rate_abs(F_MF, exchange_HB_flux, 20);

    // Print global exchangebility matrix
    if(1){
      sprintf(type, "exchangeability_global_%s_%s", MATRIX, "FLUX");
      file_mat=Output_file(name_file, type, "txt");
      sprintf(name, "%s_%s.txt", name_file, type);
      fprintf(out, "Printing global exchangeability matrix of type ");
      fprintf(out, "%s %c in %s\n", MATRIX, EXCH, name);
      fprintf(file_mat, exc_type);
      if(FORMAT==0){
	fprintf(file_mat, "AA ");
	for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_string[a]);
      }else{
	for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
      }
      fprintf(file_mat, "\n");
      Print_matrix(exchange_HB_flux, iaa, norm_rate, file_mat,
		   AA_string, P_mut_a, FORMAT);
      if(FORMAT){
	fprintf(file_mat, "\n// end of data. The rest are notes\n");
	for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
	fprintf(file_mat, "\n");
      }
      fclose(file_mat);
    }

    float **exchange_HB_i=Allocate_mat2_f(20,20);

    // Print site-specific exchange matrices
    if(PRINT_ALL_EXCH){
      sprintf(type, "exchangeability_sites_%s_%s", MATRIX, "FLUX");
      file_mat=Output_file(name_file, type, "txt");
      sprintf(name, "%s_%s.txt", name_file, type);
      fprintf(out, "Printing site-specific exchangeability matrices of type ");
      fprintf(out, "%s %c in %s\n", MATRIX, EXCH, name);
      fprintf(file_mat, exc_type);
      if(FORMAT==0){
	fprintf(file_mat, "AA ");
	for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_string[a]);
      }else{
	for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
      }
      fprintf(file_mat, "\n");
    }

    int ii=1;
    for(i=0; i<L; i++){
      if(wi[i]==0){continue;}
      Exch_Halpern(Rate_HB_flux+i, exchange_HB_i, NULL,
		   P_MF_ia+i, P_mut_a, 1, exchange_HB_flux, NULL);
      if(PRINT_ALL_EXCH){
	fprintf(file_mat, "SITE %3d %c\n", ii, res[i].seq); ii++;
	Print_matrix(exchange_HB_i, iaa, norm_rate, file_mat,
		     AA_string, P_MF_ia[i], FORMAT);
      }
    }

    if(PRINT_ALL_EXCH){
      if(FORMAT){
	fprintf(file_mat, "\n// end of data. The rest are notes\n");
	for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
	fprintf(file_mat, "\n");
      }
      fclose(file_mat);
    }
    Empty_matrix_f(exchange_HB_i, 20);

    sprintf(name, "%s_rate_profile.dat", name_file);
    fprintf(out, "Printing rates of all positions in %s\n", name);
    FILE *file_out=Output_file(name_file, "rate_profile", "dat");
    fprintf(file_out, "# Protein %s L= %d\n", nameout, L);
    fprintf(file_out, "# Local interactions coefficient=%.3f\n", SEC_STR);
    fprintf(file_out, "# Frequencies rate entropy corr(f_a,sum_b E_ab)\n");
    float p[20]; for(a=0; a<20; a++)p[a]=0.05;
    fprintf(file_out, "# Equal       %.3f %.3f %6.3f\n",
	    Compute_rate_abs(p, exchange_flux, 20),
	    log(20.),
	    Corr_vM(p, exchange_flux, 20));
    fprintf(file_out, "# mean-field  %.3f %.3f %6.3f\n",
	    Compute_rate_abs(F_MF, exchange_flux, 20),
	    Entropy(F_MF, 20),
	    Corr_vM(F_MF, exchange_flux, 20));
    fprintf(file_out, "# %s_%s     %.3f %.3f %6.3f\n", MATRIX, namexc,
	    Compute_rate_abs(f_emp, exchange_flux, 20),
	    Entropy(f_emp, 20),
	    Corr_vM(f_emp, exchange_flux, 20));
    fprintf(file_out, "#SITE AA sec.str. ncont  entropy ave_hydro");
    fprintf(file_out, " rate_HB_mut rate_HB_emp rate_HB_flux");
    //fprintf(file_out, " rate_abs_mut rate_abs_emp rate_abs_flux");
    fprintf(file_out, "\n");
    
    
    struct res_short *rs=res;
    float entropy[L], hydro_i[L], hydro_ave[L];
    double hydro_tot=0;
    for(i=0; i<L; i++){
      if(wi[i]==0){continue;}
      float *P_ia=P_MF_ia[i];
      entropy[i]=Entropy(P_ia, 20);
      hydro_ave[i]=Mean_hydro(P_ia, hscale, 20);
      hydro_i[i]=hscale[res[i].aa];
      hydro_tot+=hydro_ave[i];
      fprintf(file_out, "%s %c %c %d %.3f %5.2f",
	      rs->pdbres, rs->seq, rs->sec_str, rs->n_cont,
	      entropy[i], hydro_ave[i]);
      rs++;
      fprintf(file_out, " %.4f %.4f %.4f\n",
	      Rate_HB_mut[i], Rate_HB_emp[i], Rate_HB_flux[i]);
      //Rate_abs_mut[i]=Compute_rate_Q(P_ia, rate_mut, 20);
      //Rate_abs_mut[i]=Compute_rate_sel(P_ia, P_mut, exchange_mut, 20);
      /*Rate_abs_mut[i]=Compute_rate_abs(P_ia, exchange_mut, 20);
	Rate_abs_emp[i]=Compute_rate_abs(P_ia, exchange_emp, 20);
	Rate_abs_flux[i]=Compute_rate_abs(P_ia, exchange_flux, 20);
	fprintf(file_out, "  %.4f %.4f %.4f",
	Rate_abs_mut[i], Rate_abs_emp[i], Rate_abs_flux[i]);*/
    }
    fprintf(file_out, "# %d sites\n", L);
    fprintf(file_out, "# Mean averaged hydrophobicity= %.4f\n", hydro_tot/L);
    float r, ncont[L]; for(i=0; i<L; i++)ncont[i]=res[i].n_cont;
    r=Corr_coeff(ncont, hydro_i, L);
    fprintf(file_out, "# Corr(ncont, hydro)=     %6.3f\n", r);
    r=Corr_coeff(ncont, hydro_ave, L);
    fprintf(file_out, "# Corr(ncont, ave_hydro)= %6.3f\n", r);
    r=Corr_coeff(ncont, entropy, L);
    fprintf(file_out, "# Corr(ncont, entropy)=   %6.3f\n", r);
    r=Corr_coeff(ncont, Rate_HB_flux, L);
    fprintf(file_out, "# Corr(ncont, rate_HB_flux)=  %6.3f\n", r);
    r=Corr_coeff(ncont, Rate_HB_emp, L);
    fprintf(file_out, "# Corr(ncont, rate_HB_emp)=  %6.3f\n", r);
    r=Corr_coeff(ncont, Rate_HB_mut, L);
    fprintf(file_out, "# Corr(ncont, rate_HB_mut)=  %6.3f\n", r);
    r=Corr_coeff(entropy, Rate_HB_flux, L);
    fprintf(file_out, "# Corr(entropy, rate_HB_flux)=%6.3f\n", r);
    r=Corr_coeff(ncont, Rate_HB_emp, L);
    fprintf(file_out, "# Corr(entropy, rate_HB_emp)=%6.3f\n", r);
    r=Corr_coeff(entropy, Rate_HB_mut, L);
    fprintf(file_out, "# Corr(entropy, rate_HB_mut)=%6.3f\n", r);
    fclose(file_out);
  }

  // FILE exchangeability_glob_%s_FLUX
  // Global exchangeability matrix with average frequencies F_MF
  if(PRINT_GLOB){
    sprintf(type, "exchangeability_glob_%s_%s", MATRIX, "FLUX");
    sprintf(name, "%s_%s.txt", name_file, type);
    fprintf(out,
	    "Printing global exchangeability matrix of type %s FLUX in %s\n",
	    MATRIX, name);
    file_mat=Output_file(name_file, type, "txt");
    fprintf(file_mat,"# %s: %s matrix\n", "FLUX", MATRIX);


    float **exchange_ave=Allocate_mat2_f(20,20);
    for(a=0; a<20; a++){
      for(b=0; b<a; b++){
	exchange_ave[a][b]= exchange_emp[a][b]*f_emp[a]*f_emp[b]
	  /(F_MF[a]*F_MF[b]);
	exchange_ave[b][a]=exchange_ave[a][b];
      }
    }
    Print_matrix(exchange_ave,iaa,norm_rate,file_mat,AA_string,F_MF,FORMAT);
    if(exchange_ave)Empty_matrix_f(exchange_ave, 20);
    if(FORMAT){
      fprintf(file_mat,"// end of data. The rest are notes\n");
      for(a=0; a<20; a++)fprintf(file_mat, "%c\t", AA_PAML[a]);
      fprintf(file_mat, "\n");
      for(a=0; a<20; a++)fprintf(file_mat, "%d\t", iaa[a]);
      fprintf(file_mat, "\nL= %d\n", L);
      fprintf(file_mat, "P_mut:\n");
      for(a=0; a<20; a++)fprintf(file_mat, "%.4f\t", P_mut_a[iaa[a]]);
      fprintf(file_mat, "\n");
    }
    fclose(file_mat);
  }

  // FILE <>_AArates.dat
  // Print rates per amino acid in matrix <>_AArates.dat
  if(PRINT_GLOB){
    sprintf(name, "%s_AArates.dat", name_file);
    file_mat=fopen(name, "w"); printf("Writing %s\n", name);
    fprintf(out, "Printing rates for all amino acids in %s\n", name);
    fprintf(file_mat, "#(s)=rescaled to reduce selection\n");
    fprintf(file_mat,"#AA Q^mut(a,a) sum_b(E^emp(s)_ab/19) ");
    fprintf(file_mat,"sum_b(E^emp(s)_ab*f_b) sum_b(E^emp_ab*f_b) hydro\n");
    float R_emp_s[20], Rf_emp_s[20], Rf_emp[20];
    for(a=0; a<20; a++){
      double R=0; for(b=0; b<20; b++)if(b!=a)R+=exchange_flux[a][b];
      R_emp_s[a]=R/19;
      R=0; for(b=0; b<20; b++)if(b!=a)R+=exchange_flux[a][b]*f_emp[b];
      Rf_emp_s[a]=R;
      R=0; for(b=0; b<20; b++)if(b!=a)R+=exchange_emp[a][b]*f_emp[b];
      Rf_emp[a]=R;
      fprintf(file_mat, "%c %.3f %.3f %.3f  %.3f %.2f\n",
	      Amin_code(a), R_mut[a], R_emp_s[a], Rf_emp_s[a],
	      Rf_emp[a], hscale[a]);
    }
    fprintf(file_mat, "# Corr(R_mut,R_emp(s))= %.3f\n",
	    Corr_coeff(R_mut, R_emp_s, 20));
    fprintf(file_mat, "# Corr(R_mut,Rf_emp(s))= %.3f\n",
	    Corr_coeff(R_mut, Rf_emp_s, 20));
    fprintf(file_mat, "# Corr(R_mut,Rf_emp)= %.3f\n",
	    Corr_coeff(R_mut, Rf_emp, 20));
    fprintf(file_mat, "# Corr(R_mut, h)=    %.3f\n",
	    Corr_coeff(R_mut, hscale, 20));
    fprintf(file_mat, "# Corr(R_emp(s), h)=    %.3f\n",
	    Corr_coeff(R_emp_s, hscale, 20));
    fprintf(file_mat, "# Corr(Rf_emp(s), h)=   %.3f\n",
	    Corr_coeff(Rf_emp_s, hscale, 20));
    fprintf(file_mat, "# Corr(Rf_emp, h)=   %.3f\n",
	    Corr_coeff(Rf_emp, hscale, 20));
    fclose(file_mat);
  }
  
  if(exchange_emp)Empty_matrix_f(exchange_emp, 20);
  if(exchange_flux)Empty_matrix_f(exchange_flux, 20);
  if(exchange_HB_flux)Empty_matrix_f(exchange_HB_flux, 20);
  if(Pair_MF)Empty_matrix_f(Pair_MF, 20);
  return(0);
}

float Compute_rate_abs(float *p, float **exchange, int n){
  // F= sum_{a!=b} P_a P_b E_ab
  double F=0; int a, b;
  for(a=1; a<n; a++){
    double Fa=0; for(b=0; b<a; b++)Fa += exchange[a][b]*p[b];
    F+=p[a]*Fa;
  }
  return(2*F);
}

float Compute_rate_Q(float *p, float **rate, int n){
  // R= sum_{a!=b} P_a Q_ab = -sum_a P_a Q_aa
  double R=0;
  for(int a=0; a<n; a++)R+=p[a]*rate[a][a];
  return(-R);
}

float Compute_rate_sel(float *p, float *p_mut, float **exchange, int n)
{
  double R=0; int a, b; float P_mix[n];
  for(a=0; a<n; a++)P_mix[a]=sqrt(p[a]*p_mut[a]);
  for(a=0; a<n; a++){
    double Ra=0;
    for(b=0; b<n; b++){
      if(b!=a)Ra+=exchange[a][b]*P_mix[b];
    }
    R+=P_mix[a]*Ra;
  }
  return(R);
}


float Corr_vM(float *p, float **exch, int n){
  float *e=malloc(n*sizeof(float)); int a, b;
  for(a=0; a<n; a++){
    e[a]=0; for(b=0; b<n; b++)if(b!=a)e[a]+=exch[a][b];
  }
  return(Corr_coeff(p, e, n));
}

void Change_AA_order(int *iaa, char *AA){
  int a; for(a=0; a<20; a++)iaa[a]=Code_AA(AA[a]);
}



/* void Compute_exchange_mut_old(float **S_mut,
			      float *mut_par, float tt_ratio, float TWONUC)
{
  // S(a,b)=sum_{c\in C(a) c'\in C(b)} w_c*Q(c,c')/(W_a*W_b)
   int a, b, c;
  for(a=0; a<20; a++)for(b=0; b<20; b++)S_mut[a][b]=0;
  double *norm=malloc(20*sizeof(double));
  for(a=0; a<20; a++)norm[a]=0;
  float w[64];
  for(c=0; c<64; c++)w[c]=Weight_codon(codon[c], mut_par);
  float r=1; if(TWONUC>0)r=TWONUC; // rate

  for(c=0; c<64; c++){
    if(coded_aa[c]=='*')continue; // Stop codon
    a=Code_AA(coded_aa[c]); norm[a]+=w[c];
    //printf("%2d a=%c b= ", c, Amin_code(a));
    float *Fa=S_mut[a];
    int j1, n1; char cod2[3], *cod=codon[c];
    for(j1=0; j1<3; j1++)cod2[j1]=cod[j1];
    for(j1=0; j1<3; j1++){
      char nj=cod[j1]; int nt=Transition(nj);
      for(n1=0; n1<4; n1++){
	if(Nuc_code(n1)==nj)continue;
	// cod2[j]=n
	float S=r;
	Update_F(Fa,&S,cod2,w,a,c,j1,n1,nt,tt_ratio,codon,coded_aa);
	if(TWONUC<=0)continue;
	int j2, n2, j3, n3;
	for(j2=j1+1; j2<3; j2++){
	  int nt2=Transition(cod[j2]);
	  for(n2=0; n2<4; n2++){
	    if(Nuc_code(n2)==cod[j2])continue;
	    float S2=S*r;
	    Update_F(Fa,&S2,cod2,w,a,c,j2,n2,nt2,tt_ratio,codon,coded_aa);
	    for(j3=j2+1; j3<3; j3++){
	      int nt3=Transition(cod[j3]);
	      for(n3=0; n3<4; n3++){
		if(Nuc_code(n3)==cod[j3])continue;
		float S3=S2*r;
		Update_F(Fa,&S3,cod2,w,a,c,
			 j3,n3,nt3,tt_ratio,codon,coded_aa);
	      } // End n3
	      cod2[j3]=cod[j3];
	    } // End j3
	  } // End n2
	  cod2[j2]=cod[j2];
	} // End j2
      } // End n1
      cod2[j1]=cod[j1];
    } // End j
  } // End c loop
  
  // Impose symmetry
  double sum=0;
  for(a=0; a<20; a++)sum+=norm[a];
  for(a=0; a<20; a++)norm[a]/=sum;
  for(a=0; a<20; a++){
    for(b=0; b<a; b++){
      S_mut[a][b]=(S_mut[a][b]+S_mut[b][a])/(2*norm[a]*norm[b]);
      S_mut[b][a]=S_mut[a][b];
    }
  }

  // Impose normalization Q_aa=-sum_b Q_ab
  //~ for(a=0; a<20; a++){
  //~   double sum=0;
  //~   for(b=0; b<20; b++){
  //~     if(a!=b)sum+=S_mut[a][b]*norm[b];
  //~   }
  //~   S_mut[a][a]=-sum/norm[a];
  //~   }
}*/

/* void Update_F(float *F_mut, float *S, char *cod2, float *w,
	      int a, int c, int j, int n, int nt, float tt_ratio)
{
  cod2[j]=Nuc_code(n);
  int c2=Code_codon(cod2, codon);
  if(coded_aa[c2]=='*')return; // Stop codon
  int b=Code_AA(coded_aa[c2]); if(b==a)return;
  if(n==nt)(*S)*=tt_ratio;
  F_mut[b]+=w[c]*w[c2]*(*S);
  //printf("%c ", Amin_code(b));
} */

/*  
void Compute_Q_sel(float **Q_sel, float *P_MF)
{
  int a, b;
  for(a=0; a<20; a++){
    double sum=0;
    for(b=0; b<20; b++){
      if(b==a)continue;
      if(P_MF[b]>=P_MF[a]){Q_sel[a][b]=1;}
      else{Q_sel[a][b]=P_MF[b]/P_MF[a];}
      sum+=Q_sel[a][b];
    }
    Q_sel[a][a]=-sum;
  }
}
*/

void Sum_matrix(float *ncont, int **Cnat, int L){
  int i, j; float sum;
  for(i=0; i<L; i++){
    sum=0; for(j=0; j<L; j++)sum+=Cnat[i][j];
    ncont[i]=sum;
  }
}

float Average_entropy(float **P_MF_ia, int L){
  double entropy=0;
  for(int i=0; i<L; i++)entropy+=Entropy(P_MF_ia[i], 20);
  return(entropy/=L);
}

float Entropy(float *PP, int n){
  double S=0, norm=0; float *p=PP; int i; 
  for(i=0; i<n; i++){
    if(*p){S+=(*p)*log(*p); norm+=(*p);} p++;
  }
  if(norm==0)return(-1); //Columns without observations
  S = S/norm -log(norm);
  return(-S);
}

float Corr_coeff(float *xx, float *yy, int n){
  double x1=0, x2=0, y1=0, y2=0, xy=0;
  int i; float *x=xx, *y=yy;
  for(i=0; i<n; i++){
    x1 += *x; x2+= (*x)*(*x);
    y1 += *y; y2+= (*y)*(*y);
    xy += (*x)*(*y); x++; y++;
  }
  float r=(n*xy-y1*x1)/sqrt((n*x2-x1*x1)*(n*y2-y1*y1));
  return(r);
}

double Mean_hydro(float *P, float *hydro, int n){
  double ave=0; int a;
  for(a=0; a<n; a++)ave+=hydro[a]*P[a];
  return(ave);
}

void Print_matrix(float **exchange, int *iaa, float norm,
		  FILE *file_mat, char *AA_string, float *p, int FORMAT)
{
  int a,b;
  if(FORMAT==0){
    for(a=0; a<20; a++){
      fprintf(file_mat, "%c", AA_string[a]);
      float *e=exchange[iaa[a]];
      for(b=0; b<20; b++)fprintf(file_mat, "\t%.4f", e[iaa[b]]/norm);
      fprintf(file_mat, "\n");
    }
  }else{
    for(a=1; a<20; a++){
      float *e=exchange[iaa[a]];
      for(b=0; b<a; b++)fprintf(file_mat, "%.4f\t", e[iaa[b]]/norm);
      fprintf(file_mat, "\n");
    }
    fprintf(file_mat, "\n");
    for(a=0; a<20; a++)fprintf(file_mat, "%.4f\t", p[iaa[a]]);
    fprintf(file_mat, "\n");
  }
}

void Compute_exchange_mut_sym(float **S_mut, float **rate_matrix, float *P_aa)
{
  // S(a,b)=sum_{c\in C(a) c'\in C(b)} w_c*Q(c,c')/(W_a*W_b)
  // From flux to exchangeability; symmetrize
  for(int a=0; a<20; a++){
    for(int b=0; b<a; b++){
      S_mut[a][b]=(P_aa[a]*rate_matrix[a][b]+P_aa[b]*rate_matrix[b][a])
	/(2*P_aa[a]*P_aa[b]);
      S_mut[b][a]=S_mut[a][b];
    }
    S_mut[a][a]=rate_matrix[a][a]/P_aa[a];
  }
}
  
void Compute_rate_matrix(float **rate_matrix, float *P_aa,
			 float *P_cod, float **Q_cod)
{
  int a, b, c, d;

  /*float P_cod[64], **Q_cod=Allocate_mat2_f(64, 64); int CpG=1; 
  Compute_P_mut(P_aa, P_cod, Q_cod, mut_par, codon, coded_aa, NULL, NULL);*/
  
  for(a=0; a<20; a++){
    P_aa[a]=0;
    for(b=0; b<20; b++)rate_matrix[a][b]=0;
  }
  for(c=0; c<64; c++){
    if(coded_aa[c]=='*')continue;
    a=Code_AA(coded_aa[c]);
    float p=P_cod[c], *Q=Q_cod[c];
    P_aa[a]+=p;
    for(d=0; d<64; d++){
      if(coded_aa[d]=='*')continue;
      b=Code_AA(coded_aa[d]);
      rate_matrix[a][b]+=p*Q[d];
    }
  } 
  double sum=0; // sum_a!=b P(a)Q(a,b) 
  for(a=0; a<20; a++){
    for(b=0; b<20; b++){
      if(a==b)continue;
      sum+=rate_matrix[a][b];
      rate_matrix[a][b]/=P_aa[a];
    }
  }
  for(a=0; a<20; a++){
    double R=0;
    for(b=0; b<20; b++){
      if(a==b)continue;
      rate_matrix[a][b]/=sum;
      R+=rate_matrix[a][b];
    }
    rate_matrix[a][a]=-R;
  }
}


float **Empirical_exchangeability(float *f_emp, char *MATRIX)
{
  float **exch_emp=Allocate_mat2_f(20, 20);
  
  int a, b;
  // Choose empirical exchangeability matrix
  int iwag[20]; float *f, *exch[20];
  if(strncmp(MATRIX, "WAG", 3)==0){
    f=fWAG; for(a=0; a<20; a++)exch[a]=WAG[a];
  }else if(strncmp(MATRIX, "JTT", 3)==0){
    strcpy(MATRIX, "JTT");
    f=fJTT; for(a=0; a<20; a++)exch[a]=JTT[a];
  }else{
    strcpy(MATRIX, "LG");
    f=f_LG; for(a=0; a<20; a++)exch[a]=LG_matrix[a];
  }
  Change_AA_order(iwag, AA_WAG);
  for(a=0; a<20; a++){
    int ia=iwag[a]; f_emp[ia]=f[a];
    for(b=0; b<a; b++){
      int ib=iwag[b];
      exch_emp[ia][ib]=exch[a][b];
      exch_emp[ib][ia]=exch[a][b];
    }
  }
  return(exch_emp);
}

float Flux_exch_MSA(float **f_msa, int L, float sum_msa,
		    float exch_emp[20][20], float fe[20])
{
  // Symmetrize
  int a, b, i;
  float exch[20][20];
  for(a=0; a<20; a++){
    for(b=0; b<a; b++){
      exch[a][b]=exch_emp[a][b];
      exch[b][a]=exch_emp[a][b];
    }
  }

  // Normalize exchangebility matrix
  float rate=0;
  for(a=0; a<20; a++){
    float Q_aa=0;
    for(b=0; b<20; b++){if(b!=a)Q_aa+=exch[a][b]*fe[b];}
    exch[a][a]=-Q_aa/fe[a];
    rate+=fe[a]*Q_aa;
  }
  for(a=0; a<20; a++){
    for(b=0; b<20; b++){exch[a][b]/=rate;}
  }

  // Choose empirical exchangeability matrix
  double Lik=0;
  int iwag[20];
  Change_AA_order(iwag, AA_WAG);
  for(i=0; i<L; i++){
    float *f=f_msa[i];
    for(a=0; a<20; a++){
      int ia=iwag[a]; if(f[ia]==0)continue;
      for(b=0; b<a; b++){
	int ib=iwag[b]; if(f[ib]==0)continue;
	Lik+=f[ia]*f[ib]*exch[a][b];
      }
    }
  }
  return(Lik/(sum_msa*sum_msa));
}

void Normalize_exchange(float **exch, float *f, int n)
{
  double sum=0; int a, b;
  for(a=1; a<n; a++){
    double Ra=0; for(b=0; b<a; b++)Ra+=exch[a][b]*f[b];
    sum+=f[a]*Ra;
    exch[a][a]=-Ra/f[a];
  }
  sum*=2;
  for(a=1; a<20; a++){
    for(b=0; b<a; b++){
      exch[a][b]/=sum; exch[b][a]=exch[a][b];
    }
    exch[a][a]/=sum;
  }
}

float **Rescale_exchangeability(float **exch_emp, float *f_emp,
				char EXCHANGE, float **Pair_MF, float *F_MF)
{
  // Rescale exchangeability
  float **exchange_emp=Allocate_mat2_f(20,20); int a, b;
  if(EXCHANGE=='E'){
    // Same exchangeability as in empirical model
    for(a=0; a<20; a++){
      for(b=0; b<a; b++){
	exchange_emp[a][b]=exch_emp[a][b];
      }
    }
  }else if(EXCHANGE=='Q'){
    // Impose that average rate (Q=Ef) is the same as in empirical model
    for(a=0; a<20; a++){
      for(b=0; b<a; b++){
	exchange_emp[a][b]=exch_emp[a][b]*
	  (f_emp[a]+f_emp[b])/(F_MF[a]+F_MF[b]);
      }
    }
  }else{ //if(EXCHANGE=='F'){
    // Impose that average flux (F) is as in empirical model
    for(a=0; a<20; a++){
      for(b=0; b<a; b++){
	exchange_emp[a][b]= exch_emp[a][b]*f_emp[a]*f_emp[b]/Pair_MF[a][b];
      }
    }
  }
  // Symmetrize
  for(a=0; a<20; a++)
    for(b=a+1; b<20; b++)exchange_emp[a][b]=exchange_emp[b][a];
  return(exchange_emp);
}

float **Exch_flux_abs(float **exch_emp, float *f_emp,
		      float **Pair_MF, float *F_MF)
{
  // Rescale exchangeability
  float **exchange_flux=Allocate_mat2_f(20,20); int a, b;
  // Impose that average flux (F) is as in empirical model
  for(a=0; a<20; a++){
    for(b=0; b<a; b++){
	exchange_flux[a][b]= exch_emp[a][b]*f_emp[a]*f_emp[b]/Pair_MF[a][b];
	exchange_flux[b][a]=exchange_flux[a][b];
    }
  }
  return(exchange_flux);
}

float **Pair_matrix(float *F_MF, float **P_MF_ia, int L, float *P_mut_a,
		    int *iaa, char *name_file, char *AA_string, int FORMAT)
{
  float **Pair_MF=Allocate_mat2_f(20,20);
  double sum=0; int a, b, i;
  for(a=0; a<20; a++){
    double P=0; for(i=0; i<L; i++)P+=P_MF_ia[i][a];
    F_MF[a]=P; sum+=P;
    for(b=0; b<a; b++){
      P=0; for(i=0; i<L; i++)P+=P_MF_ia[i][a]*P_MF_ia[i][b];
      P/=L; Pair_MF[a][b]=P; Pair_MF[b][a]=P;
    }
  }
  for(a=0; a<20; a++)F_MF[a]/=sum;
  if(0){
    // Print selection matrix
    FILE *file_out=Output_file(name_file, "selection", "txt");
    fprintf(file_out, "L=%d\n", L);
    Print_matrix(Pair_MF, iaa, 1, file_out, AA_string, P_mut_a, FORMAT);
    for(a=0; a<20; a++)fprintf(file_out, "%c\t", Amin_code(iaa[a]));
    fprintf(file_out, "\n"); fclose(file_out);
  }
  return(Pair_MF);
}


void Exch_Halpern(float *rate_HB, float **exchange_HB,
		  float **exchange_HB_flux, // output
		  float **P_MF_ia, float *P_mut, int L,  // Input
		  float **exchange_emp, float *f_emp)
{
  float P_MIN=0.001, log_MIN=log(P_MIN);
  int i,a,b;
  if(exchange_HB_flux){
    for(a=0; a<20; a++)for(b=0; b<20; b++)exchange_HB_flux[a][b]=0;
  }
  for(i=0; i<L; i++){
    if(rate_HB)rate_HB[i]=0;
    float *Pi=P_MF_ia[i];
    float P_sel[20], log_P_sel[20];
    for(a=0; a<20; a++){
      P_sel[a]=Pi[a]/P_mut[a];
      if(P_sel[a]>P_MIN){log_P_sel[a]=log(P_sel[a]);}
      else{log_P_sel[a]=log_MIN;}
    }
    for(a=1; a<20; a++){
      double R_a=0, Fix;
      // Fixation probability from a to b
      for(b=0; b<a; b++){
	if(P_sel[a]!=P_sel[b]){
	  Fix=(log_P_sel[a]-log_P_sel[b])/(P_sel[a]-P_sel[b]);
	}else{
	  Fix=1;
	}
	double e_fix=exchange_emp[a][b]*Fix;
	R_a+=Pi[b]*e_fix;
	if(exchange_HB){
	  //float **exchange_HB_i=exchange_HB[i];
	  exchange_HB[a][b]=e_fix;
	  exchange_HB[b][a]=e_fix;
	}
	if(exchange_HB_flux)
	  exchange_HB_flux[a][b]+=Pi[a]*Pi[b]*Fix;
      }
      if(rate_HB){rate_HB[i]+=Pi[a]*R_a;}
    }
    if(rate_HB)rate_HB[i]*=2;
  }
  if(exchange_HB_flux){
    for(a=0; a<20; a++){
      for(b=0; b<a; b++){
	exchange_HB_flux[a][b]=
	  L*f_emp[a]*f_emp[b]*exchange_emp[a][b]/exchange_HB_flux[a][b];
      }
    }
  }
  return;
}
