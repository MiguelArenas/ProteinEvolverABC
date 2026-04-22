#ifndef PROFILE_OPTIONS
#define PROFILE_OPTIONS

typedef struct{

    
  alignmentClass alignKind;
  fileType inFile;

  char seqFile[MAX_SIZE+1];

  cbool doPerm;
  int ntrials;
  cbool otherStats;
  
  double winSize;
  int k;

  cbool verbose;
  
  int breakWindow;
  int stepSize;
 

} options;


#endif
