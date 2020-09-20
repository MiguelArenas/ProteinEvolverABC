
#ifndef GLOBAL
#define GLOBAL

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>

#define PHYLIP_SIZE 10
#define MAX_SIZE 1000

/* For reading in data... */
#define CHAR_START  33
#define CHAR_END  126
//#define MAX_STATES  94

#define MAX_STATE 100


/* Two ways to treating missing data - either as no data or as an extra state */
#define GLOBAL_GAP  '-'

#define GLOBAL_MISSING  '?'
#define MISSING_STATE 95

#define IGNORE_STATE 101
extern const char DNA_ALPHA[];
extern const int DNA_ALPHA_SIZE;

extern const char DNA_MISSING[];
extern const int DNA_MISSING_SIZE;

extern const char DNA_AMBIG[];
extern const int DNA_AMBIG_SIZE;

extern const char AA_ALPHA[] ;
extern const int AA_ALPHA_SIZE ;
extern const char AA_MISSING[];
extern const int  AA_MISSING_SIZE;
extern const char AA_AMBIG[];
extern const int AA_AMBIG_SIZE;



//#define MISSING_STATE 101


typedef enum {DNA,AA,OTHER} alignmentClass;

typedef enum {FALSE, TRUE} cbool;
typedef enum {constant, polymorphic} polystatus;
typedef enum {uninformative, informative} infstatus;
typedef enum {nongapped,gapped} gapstatus;

typedef enum {strict,relaxed,fasta} fileType;

typedef struct {
  polystatus poly;
  infstatus inf;
  gapstatus gap;
  int       orig_index;
  int       num_states;
  int       num_missing;
} site;


typedef int inc_type;
typedef char align_type;


cbool memberOf(const char *set, const int num, char ch);

cbool validState(alignmentClass alignKind, char ch);

cbool missing_ambig_State(alignmentClass alignKind, char ch);

cbool valid_gap(alignmentClass alignKind, char ch);
#endif
