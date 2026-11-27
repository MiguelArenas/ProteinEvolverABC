#include <stdio.h>
char Maiuscule(char);
char Minuscule(char);
int Check_file(char *name);
FILE *Open_file_r(char *name, char *code_name, char *message);
float Read_column(char **s, char *string, int i);
void Read_hydro_list(char *file_hydro, char *hydro_name, float *hydro);
