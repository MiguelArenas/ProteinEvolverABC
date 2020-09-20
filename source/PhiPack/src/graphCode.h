#ifndef GRAPH
#define GRAPH

void create_incompat_pic(int num_states,int max_inc,inc_type **inc_matrix, cbool embellish, int num_chars);

void create_profile_pic(int num_sites,int k,int num_inf, double *break_values,site *site_desc);

#endif
