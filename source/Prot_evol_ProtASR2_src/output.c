#include <stdio.h>
#include <string.h>
#include "output.h"

FILE *Output_file(char *name_file, char *what, char *ext)
{
  char nameout[200];
  sprintf(nameout, "%s_%s.%s", name_file, what, ext);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);
  return(file_out);
}
