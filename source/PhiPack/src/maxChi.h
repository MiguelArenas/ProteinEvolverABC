#include "global.h"
#ifndef MAXCHI
#define MAXCHI

void maxchi(double windowScale,cbool verbose,FILE *logfile,fileType inFile,align_type **alignment, site *site_desc, int *site_states, char **taxa_names,int num_taxa, int num_sites, int num_trials, int *emp_MAXCHI);

void all_pairs_maxchi( align_type **alignment, int num_taxa, int num_sites, int windowsize, int *permutation, site *site_desc,double *maxvalue,char **taxa_names,cbool report_breakpoint,FILE *logfile);

#endif
