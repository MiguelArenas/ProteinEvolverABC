#include "misc.h"
#include <ctype.h>


void error(const char * message) 
{
  fprintf(stderr,"\nERROR: %s\n",message);
  fprintf(stderr,"Exiting...\n");
  
  exit(1);
  
}


char *strsub(char *s1, const char *s2, int start)
{
  char *p;
  p=s1;
  s2=s2+start;
  while( *s2 != '\0')
    {
      *p=*s2;
      p++;
      s2++;
    }
  *p='\0';
  return p;
  
}

char skip_all_space(FILE *in_file)
{
  char ch;
  
  ch=fgetc(in_file);
  while((ch != EOF) && isspace(ch))
    ch=fgetc(in_file);
  
  if(ch!=EOF)
    ungetc(ch,in_file);
  
  return ch;
}


char skip_newlines(FILE *in_file)
{
  char ch;
  
  ch=fgetc(in_file);
  while((ch != EOF) && (ch == '\n'))
    ch=fgetc(in_file);
  
  if(ch!=EOF)
    ungetc(ch,in_file);
  
  return ch;

}

char skip_non_newline(FILE *in_file)
{
  char ch;
  
  ch=fgetc(in_file);
  while((ch != EOF) && (ch != '\n'))
    ch=fgetc(in_file);
  
  if(ch!=EOF)
    ungetc(ch,in_file);
  
  return ch;
}

char skip_non_newline_space(FILE *in_file)
{
  char ch;
  
  ch=fgetc(in_file);
  while((ch != EOF) && (ch != '\n') && isspace(ch))
    ch=fgetc(in_file);
  
  if(ch!=EOF)
    ungetc(ch,in_file);
  
  return ch;
}
