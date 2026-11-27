#ifndef OPTIONS
#define OPTIONS

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
  
  cbool printMatrix;
  int graphType;
  cbool embellish;
  
  cbool printBreakpoints;
  int breakWindow;
  
 

} options;


#endif
