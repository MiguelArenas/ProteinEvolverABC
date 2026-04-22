#ifndef FASTA
#define FASTA

void write_fasta(FILE* cur_stream,  char **taxa_names,align_type **alignment, int num_taxa,int num_sites);

void read_fasta(const char* file_name,  char ***taxa_names, align_type ***alignment, int *num_taxa, int *num_sites);

#endif
