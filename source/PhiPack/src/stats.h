#ifndef STATS
#define STATS

void identity_permutation(int *permutation, int n);


void sample_permutation(int *permutation, int n);


/* Print out histogram of incompatibilities */
void inc_stats(int max_size,int num_chars,cbool verbose,inc_type **inc_matrix,int **counts, FILE *logfile);

double PHI(inc_type **inc_matrix,int *num_states,int *permutation, int num_chars, int k);


double NSS(inc_type **inc_matrix, int *permutation, int num_chars);


void get_f_and_g(inc_type **inc_matrix,int *num_states,int num_chars, double *f_values,double *g_values);

void write_diag(const char *file_name,int k,double small_variance, double mean, inc_type** inc_matrix,int *num_states, int num_chars);

void analytic_mean_variance(inc_type** inc_matrix, int* num_states, int num_chars,int k_val ,double *mean, double *varnce) ;
#endif
