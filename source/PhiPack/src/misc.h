#include<stdlib.h>
#include<stdio.h>
#ifndef MISC
#define MISC

char *strsub(char *s1, const char *s2, int start);

char skip_all_space(FILE *in_file);

char skip_non_newline(FILE *in_file);


char skip_newlines(FILE *in_file);

char skip_non_newline_space(FILE *in_file);

#endif
