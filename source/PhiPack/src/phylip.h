#ifndef PHYLIP
#define PHYLIP

/* char skip_all_space(FILE *in_file); */


/* char skip_newlines(FILE *in_file); */

/* char skip_non_newline_space(FILE *in_file); */



/* void read_strict_name(FILE *in_file,char *name, int capacity); */




/* void read_relaxed_name(FILE *in_file,char *name, int capacity); */

/* int read_seq_bit(FILE *in_file, int base_limit, align_type* states,int index, int *num_bases); */



void read_phylip(const char* file_name, cbool strict, char ***taxa_names, align_type ***alignment, int *num_taxa, int *num_sites);

void write_phylip(FILE* cur_stream, cbool strict, char **taxa_names,align_type **alignment, int num_taxa,int num_sites);

#endif
