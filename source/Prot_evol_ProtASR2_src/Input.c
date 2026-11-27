/* Input */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "input.h"
#include "codes.h"
#define CODE_NAME "Input.c"


char Minuscule(char s){
  int i=(int)s;
  if((i>=65)&&(i<=90)){return((char)(i+32));}else{return(s);}
}

char Maiuscule(char s){
  int i=(int)s;
  if(i>96){return((char)(i-32));}else{return(s);}
}

FILE *Open_file_r(char *name, char *code_name, char *message){
  FILE *file_in=fopen(name, "r");
  if(file_in==NULL){
    printf("ERROR in %s, input file %s not found (%s)\n",
	   code_name, name, message); exit(8);
  }
  return(file_in);
}

int Check_file(char *name){
  FILE *file_in=fopen(name, "r");
  if(file_in!=NULL){
    fclose(file_in); return(1);
  }
  return(0);
}

/* Reading hydro scales */
void Read_hydro_list(char *file_hydro, char *hydro_name, float *hydro)
{
  FILE *file_in=Open_file_r(file_hydro, CODE_NAME, "hydro scales ");
  char string[200], name[20], amm; int i, ia; float h;
  printf("Reading hydrophobicity scale in %s\n", file_hydro);

  while(fgets(string, sizeof(string), file_in)!=NULL){
    sscanf(string, "%s", name);
    if(strncmp(name, hydro_name, 3)==0){
      printf("Scale: %s\n", hydro_name);
      for(i=0; i<20; i++){
	fgets(string, sizeof(string), file_in);
	sscanf(string, "%c%f", &amm, &h);
	ia=Code_AA(amm); hydro[ia]=h;
	if(ia<0){
	  printf("Error in %s, wrong aa code %c. Goodbye\n",
		 file_hydro, amm); exit(8);
	}
      }
      fclose(file_in);
      return;
    }else{
      for(i=0; i<21; i++)fgets(string, sizeof(string), file_in);
    }
  }

  /* Normalizing the scale */
  {
    float ave=0, dev=0, h;
    for(i=0; i<20; i++){h=hydro[i]; ave+=h; dev+=h*h;}
    ave/=20; dev=sqrt(dev/20.-ave*ave);
    for(i=0; i<20; i++)hydro[i]/=dev;
  }

  printf("Error, scale %s not found in %s\n", hydro_name, file_hydro);
  exit(8);
}

float Read_column(char **s, char *string, int i){
  /* Reads the first column and updates the char pointer *s */
  float data; int read=0;
  sscanf(*s, "%f", &data);
  while(**s!='\n'){
    if(**s!=' '){
      if(read==0)read=1;
    }else{
      if(read)return(data);
    }
    (*s)++;
  }
  if(read)return(data);
  printf("WARNING, column %d not read in line:\n%s", i, string);
  return(0);
}
