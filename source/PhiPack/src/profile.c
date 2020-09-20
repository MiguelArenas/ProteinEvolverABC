/*  
   Copyright (c)2005, Trevor Bruen
   All rights reserved.                          

   Any feedback is very welcome.
   email: trevor@mcb.mcgill.ca
*/

#include <time.h>

#include "global.h"
#include "profile_options.h"
#include "misc.h"
#include  "pairScore.h"


#include "normal.h"
#include "stats.h"
#include "mem.h"
#include "phylip.h"
#include "fasta.h"
#include "seqManip.h"


void print_usage()
{
  
  fprintf(stderr, "Input file not specified!!! (use -f|-r|-s)\n");
  
  fprintf(stderr, "Usage: Phi   [-f|-s|-r] Filename  [-t] AlignmentType\n"); 
  fprintf(stderr, "             [-w #] [-k #] [- v] [-n #] [-m #]\n\n");
  
  
  fprintf(stderr, "Options:\n"); 
  fprintf(stderr, "  -f: Filename = FASTA format\n");
  fprintf(stderr, "  -s: Filename = Strict phylip file\n");
  fprintf(stderr, "  -r: Filename = Relaxed phylip file\n");
  fprintf(stderr, "  -t: AlignmentType = D|A|O where D=DNA\n"); 
  fprintf(stderr, "                      A=AA and O=OTHER [default DNA]\n");
  
  fprintf(stderr, "  -w: # = Change default window size for PHI [default w = 100]\n");
  
 
  
  fprintf(stderr, "  -v: Verbose [default = FALSE]\n");
  
  fprintf(stderr, "  -n: # = Change scanning size for regions to test [default n = 1000]\n");
  fprintf(stderr, "  -m: # = Step size for profile method [default m = 25]\n");
  
  
  
  

}


void get_params(int argc, char**argv, options *opt) 
{ 
  char *cur,ch,nextch;
  char temp[MAX_SIZE+1];
  int i;
  cbool inFileFound=FALSE;
  
  opt->doPerm=FALSE; 
  opt->winSize=100;
  opt->k=0;
  opt->ntrials=1000;
  opt->alignKind=DNA;
  opt->breakWindow=1000;
  opt->stepSize=25;
 
  opt->verbose=FALSE;
  
  if(argc < 3)
    {
      print_usage();
      exit(1);
    }
  
  
  for(i=1;i<argc;i++)
    {
      cur=argv[i];
      ch=cur[0];
      
      //cout<<"Have "<<cur<<" and "<<ch<<endl;
      if(ch == '-')
	{
	  ch=toupper(cur[1]);
	  switch(ch)
	    {
	  
	      /* Type of sequence */
	    case 'T':
	      i++;
	      cur=argv[i];
	      nextch=toupper(cur[0]);
	      if(nextch == 'D')
		opt->alignKind=DNA;
	      else if(nextch == 'A')
		opt->alignKind=AA;
	      else if(nextch == 'O')
		opt->alignKind=OTHER;
	      else
		opt->alignKind=OTHER;
	      //Error
	      break;
	      
	      /* Fasta file */ 
	    case 'F':
	      if(!inFileFound)
		{
		  inFileFound=TRUE;
		  opt->inFile=fasta;
		  i++;
		  strsub(opt->seqFile,argv[i],0);
		  break;

		}
	      else
		{
		  //Error
		}
	      break;
	      
	      /* Strict phylip */
	    case 'S':
	      if(!inFileFound)
		{
		  inFileFound=TRUE;
		  opt->inFile=strict;
		  i++;
		  strsub(opt->seqFile,argv[i],0);
		  break;

		}
	      else
		{
		  //Error
		}
	      break;
	      
	      /* Relaxed phylip */
	    case 'R':
	       if(!inFileFound)
		{
		  inFileFound=TRUE;
		  opt->inFile=relaxed;
		  i++;
		  strsub(opt->seqFile,argv[i],0);
		  break;

		}
	      else
		{
		  //Error
		}
	       break;
	    
	    case 'V':
	      opt->verbose=TRUE;
	      break;
	      

	      
	      /* Or window size */
	    case 'W':
	      i++;
	      strsub(temp,argv[i],0);
	      opt->winSize=atof(temp);
	      break;
	   
	   
	      /* Window size for profile method */
	    case 'N':
	      i++;
	      strsub(temp,argv[i],0);
	      opt->breakWindow=atoi(temp);
	      break;
	      
	      /* Step size for profile method*/
	    case 'M':
	      i++;
	      strsub(temp,argv[i],0);
	      opt->stepSize=atoi(temp);
	      break;

	    }
	}
      else
	{
	  //Error
	}
    }
  // Error if no file
  if(!(inFileFound))
    {
      print_usage();
      exit(1);
    }
  
} 
void print_vals(FILE *logfile,cbool print_val_a,cbool print_val_b,double val_a, double val_b)
{
  if(print_val_a)
     {
       fprintf(stdout,"%4.2le          ",val_a);
       fprintf(logfile,"%4.2le          ",val_a);
     }
  else
    {
      fprintf(stdout,"     --          ");
      fprintf(logfile,"     --          ");
    }
  if(print_val_b)
    {
      fprintf(stdout,"%4.2le\n",val_b);
      fprintf(logfile,"%4.2le\n",val_b);
    }
  else
    {
      fprintf(stdout,"--\n");
         fprintf(logfile,"--\n");
     }
 
}


int main(int argc, char* argv[])
{
  align_type **alignment;
  char **taxa_names;
  int num_taxa,num_sites;
  int num_inf;
  int pair_inc;
  double divg,val;

  /*Parameters */
  options opt;
  
  /* Number of states at each character and description of each site */
  int  *site_states;
  site *site_desc;

  /* Counts of incompatibilities */
  int *counts,max_state,max_inc;


  /* Alignment for informative states */
  align_type **inf_alignment;
  int *inf_states;
  site *inf_site_desc;

  /* Alignment by characters */
  align_type **char_alignment;
  
  /* Matrix of incompatibilities */
  int i,j;
  inc_type **inc_matrix;
  
  /* Log file */
  FILE *logfile;

  /* Stuff for profiling... */
  int big_counter,k;
    double *values;
  int num_tests;
  FILE *vals_file;
  align_type **new_alignment;
 

  /* Values of statistics */
  cbool valid_normal_approx;
  double orig_PHI,cur_PHI;
  double sum_PHI=0,sum_sq_PHI=0,obs_mean=0,obs_varnce=0;
  double mean=0.0,variance=0.0;
  double difference,normal_p_val=0.0;
  int emp_PHI=0;
 
  int *permutation;

 
  /* Seed random number generator */
  srand(time(NULL));
     
  get_params(argc,argv,&opt);
  
  /* Open log file */

  logfile=ffopen("Profile.log","w");
 
  
  fprintf(stdout,"Reading sequence file %s\n",opt.seqFile);
  fprintf(logfile,"Reading sequence file %s\n",opt.seqFile);
  
  switch(opt.inFile)
    {
    case strict:
      read_phylip(opt.seqFile,TRUE, &taxa_names,&alignment,&num_taxa,&num_sites);
      break;
    case relaxed:
      read_phylip(opt.seqFile,FALSE, &taxa_names,&alignment,&num_taxa,&num_sites);
      break;
    case fasta:
      read_fasta(opt.seqFile, &taxa_names,&alignment,&num_taxa,&num_sites );
      break;
    }
  
  fprintf(logfile,"Allocating space for %d taxa and %d sites\n",num_taxa,num_sites);
  
  if(num_taxa<=0 || num_sites<=0)
    error("No sequences or sitesin alignment");
  
  
  /* Validate alignment type */
  validate_alignment(alignment,opt.alignKind,num_taxa,num_sites);
  fprintf(stdout,"Alignment looks like a valid ");
  fprintf(logfile,"Alignment looks like a valid ");
  switch(opt.alignKind)
    {
    case DNA:
      fprintf(stdout,"DNA alignment.\n");
      fprintf(logfile,"DNA alignment.\n");
      break;
    case AA:
      fprintf(stdout,"AA alignment.\n");
      fprintf(logfile,"AA alignment.\n");
      break;
    case OTHER:
      fprintf(stdout,"OTHER alignment.\n");
      fprintf(logfile,"OTHER alignment.\n");
      break;
    }


  /* Big loop here? */
  num_tests=(int)((num_sites-opt.breakWindow)/opt.stepSize)+1;
  values=(double *)mmalloc(num_tests * sizeof(double));
  
    
  fprintf(stdout,"\nPHI is calculated with w as %3.0lf\n\n",opt.winSize);
  fprintf(logfile,"\nPHI is calculated with w as %3.0lf\n\n",opt.winSize);
  
  
  fprintf(stdout,"Checking for recombination using a window size of %d and a step size of %d\n\n",opt.breakWindow,opt.stepSize);
  fprintf(logfile,"Checking for recombination using a window size of %d and a step size of %d\n\n",opt.breakWindow,opt.stepSize);
  
  //  fprintf(stderr,"Allocated total before is %d\n",(int)alloc_size());
  for(big_counter=0;big_counter<num_tests;big_counter++)
    {
   
      fprintf(stdout,"Checking for recombination at %d between %d and %d\n",big_counter*opt.stepSize+opt.breakWindow/2,big_counter*opt.stepSize,big_counter*opt.stepSize+opt.breakWindow);
      fprintf(logfile,"Checking for recombination at %d between %d and %d\n",big_counter*opt.stepSize+opt.breakWindow/2,big_counter*opt.stepSize,big_counter*opt.stepSize+opt.breakWindow);
      
      extract(alignment,num_taxa,num_sites,big_counter*opt.stepSize,big_counter*opt.stepSize+opt.breakWindow,&new_alignment);
      
      // fprintf(stderr,"After extract is %d\n",(int)alloc_size());

      /* Initialize k*/
      opt.k=0;
      
      get_seq_div(new_alignment, opt.alignKind, num_taxa, opt.breakWindow, &divg);
      if(opt.verbose)
	{
	  fprintf(stdout,"Estimated diversity is (pairwise deletion - ignoring missing/ambig): %4.1lf%%\n",(divg*100.0));
	  fprintf(logfile,"Estimated diversity is (pairwise deletion - ignoring missing/ambig): %4.1lf%%\n",(divg*100.0));
	}
      
      /* Get informative sites */
      find_states(new_alignment, opt.alignKind,FALSE,num_taxa, opt.breakWindow, &site_desc,  &site_states);

      //  fprintf(stderr,"After states is %d\n",(int)alloc_size());
      
      
      get_informative(new_alignment, site_desc, site_states, num_taxa, opt.breakWindow, &num_inf, &inf_alignment, &inf_states, &inf_site_desc);
      
      
      if(opt.verbose)
	{
	  fprintf(stdout,"Found %d informative sites.\n",num_inf);
	  fprintf(logfile,"Found %d informative sites.\n",num_inf);
  	}
      
      //   fprintf(stderr,"After informative is %d\n",(int)alloc_size());
      
      
   /* Allocate possibly sparsely - big matrix sequential stride */
      inc_matrix=(inc_type **)mmalloc(num_inf * sizeof(inc_type *) );
      for(i=0;i<num_inf;i++)
	inc_matrix[i]=(inc_type *)mmalloc( num_inf * sizeof(inc_type));
      
      
      // fprintf(stderr,"After matrix is %d\n",(int)alloc_size());
   
      /* Reorder by character */
      reorder_chars(inf_alignment, opt.alignKind,FALSE,num_inf, num_taxa, &char_alignment);
      
      //  fprintf(stderr,"After character is %d\n",(int)alloc_size());
      switch(opt.alignKind)
	{
	case DNA:
	  max_state=4;
	  break;
	case AA:
	  max_state=20;
	  break;
	case OTHER:
	  max_state=MAX_STATE;
	  break;
	default:
	  max_state=0;
	}
      max_inc=max_state*max_state-2*max_state+1;
      
      //  fprintf(stdout,"Calculating all pairwise incompatibilities...\n");
      /* Now get incompatibilities between all pairs... */
      for(i=0;i<num_inf;i++)
	{
	  /* Fix / remove this... */
	  if(i % 100 == 0)
	    {
	      //if(i != 0)
		//	fprintf(stdout,"\b\b\b\b\b\b");
	      // else
	      //	fprintf(stdout,"Done: ");
	      
	      val=((double)i)*((double)num_inf)-((double)i)*((double)(i+1))/2.0;
	      val=val/(((double)num_inf)*((double)num_inf-1)/2);
	      val=val*100.0;
	      //    fprintf(stdout,"%5.1f%%",val);
	      
	      //  fflush(stdout);
	    }
	  
	  for(j=i;j<num_inf;j++)
	    {
	      if(i == j)
		inc_matrix[i][j]=(inc_type)0;
	      else
		{
		
		  pair_inc=pair_score(char_alignment, inf_states,i, j, num_inf, num_taxa);
		
		  inc_matrix[i][j]=(inc_type)pair_inc;
		  inc_matrix[j][i]=(inc_type)pair_inc;
		}
	    }
	}
      
      //  fprintf(stdout,"\b\b\b\b\b\b");
      //  fprintf(stdout,"%5.1f%%",100.0);
      //  fflush(stdout);
      //  fprintf(stdout,"\n\n");
      
      inc_stats(max_inc,num_inf,opt.verbose,inc_matrix,&counts,logfile);
      // fprintf(stderr,"After inc_stats is %d\n",(int)alloc_size());
      if(opt.k == 0)
	{
	  val=((double)num_inf)/((double)opt.breakWindow);
	  val=val*opt.winSize;
	  
	  opt.k=(int)(val+0.5);
	  
	  if(opt.k == 0)
	    opt.k++;
	}
      else
	{
	  val=((double)opt.breakWindow)/((double)num_inf);
	  val=val*opt.k;
	  opt.winSize=val;
	}
      
      if(opt.verbose)
	{
	  fprintf(stdout, "Using a window size of %3.0lf with k as %d\n",opt.winSize,opt.k);
	  fprintf(logfile, "Using a window size of %3.0lf with k as %d\n",opt.winSize,opt.k);
	}
      
      permutation=(int *)mmalloc(num_inf * sizeof(int));
      // fprintf(stderr,"After permutation is %d\n",(int)alloc_size());

      identity_permutation(permutation,num_inf);
      orig_PHI=PHI(inc_matrix,inf_states,permutation,num_inf,opt.k);
      
      if(num_inf <= 2*opt.k)
	{
	  valid_normal_approx=FALSE;
	  fprintf(stdout, "Too few informative sites to use normal approximation.\nTry doing a permutation test or increasing alignment length\nCan also try decreasing windowsize.\n\n");
	  fprintf(logfile, "Too few informative sites to use normal approximation.\nTry doing a permutation test or increasing alignment length\nCan also try decreasing windowsize.\n\n");
	}
      else
	{
	  valid_normal_approx=TRUE;
	  //	  fprintf(stdout,"\nCalculating analytical mean and variance\n");
	  //	  fprintf(logfile,"\nCalculating analytical mean and variance\n");
	  
	  /* Fix argument with 1 */
	  analytic_mean_variance(inc_matrix,inf_states,num_inf,opt.k,&mean,&variance);
	  
	  difference=mean-orig_PHI;
	  normal_p_val=normal_01_cdf((-difference)/sqrt(variance));
	  
	  values[big_counter]=normal_01_cdf((-difference)/sqrt(variance));

	  if(opt.verbose)
	    {
	      fprintf(stdout,"\nZ-value for PHI is: %le\n",((-difference)/sqrt(variance)));
	      fprintf(logfile,"\nZ-value for PHI is: %le\n",((-difference)/sqrt(variance)));
	    }
	  
	}
      
      if(opt.doPerm)
	{
	  if(opt.k >= num_inf)
	    {
	      error("Too few informative sites to test significance.  Try decreasing windowsize or increasing alignment length\n\n");
	    }
	  //	  fprintf(stdout,"\nDoing permutation test for PHI\n");
	  //	  fprintf(logfile,"\nDoing permutation test for PHI\n");
	  for(i=0;i<opt.ntrials;i++)
	    {
	      sample_permutation(permutation,num_inf);
	      cur_PHI=PHI(inc_matrix,inf_states,permutation,num_inf,opt.k);
	      if(cur_PHI <= orig_PHI)
		emp_PHI++;
	      sum_PHI=sum_PHI+cur_PHI;
	      sum_sq_PHI=sum_sq_PHI+cur_PHI*cur_PHI;
	      
	    }
	  obs_mean=sum_PHI/((double)opt.ntrials);
	  obs_varnce=(sum_sq_PHI-opt.ntrials*obs_mean*obs_mean)/(opt.ntrials-1);
	  
	}
      
      
      if(opt.verbose)
	{
	  
	  fprintf(stdout,"\n                      PHI Values\n");
	  fprintf(stdout,"                      ----------\n");
	  
	  fprintf(logfile,"\n                      PHI Values\n");
	  fprintf(logfile,"                      ----------\n");
	  
	  if(!opt.doPerm)
	    i=0;
	  else
	    i=opt.ntrials;
	  
	  fprintf(stdout,"              Analytical    (%d) Permutations\n\n",i);
	  fprintf(logfile,"              Analytical    (%d) Permutations\n\n",i);
	  
	  fprintf(stdout,"Mean:          ");
	  fprintf(logfile,"Mean:          ");
	  
	  print_vals(logfile,valid_normal_approx,opt.doPerm,mean,obs_mean);
	  
	  fprintf(stdout,"Variance:      ");
	  fprintf(logfile,"Variance:      ");
	  
	  print_vals(logfile,valid_normal_approx,opt.doPerm,variance,obs_varnce);
       
	  fprintf(stdout,"Observed:      %4.2le          %4.2le   \n\n",orig_PHI,orig_PHI);
	  fprintf(logfile,"Observed:      %4.2le          %4.2le   \n\n",orig_PHI,orig_PHI);
	
      
      
	  fprintf(stdout,"\n     **p-Value(s)**     \n");
	  fprintf(stdout,"       ----------\n\n");
      
      
	  fprintf(logfile,"\n       p-Value(s)\n");
	  fprintf(logfile,"       ----------\n\n");
	  
	}
   
      if(valid_normal_approx)
	{
	  fprintf(stdout,"PHI (Normal):        %4.2le\n",normal_p_val);
	  fprintf(logfile,"PHI (Normal):        %4.2le\n",normal_p_val);
	}
      else
	{
	  fprintf(stdout,"PHI (Normal):        --\n");
	  fprintf(logfile,"PHI (Normal):        --\n");
	}
    
      //  fprintf(stderr,"Total allocated is %d\n",(int)alloc_size());

       mfree(permutation,sizeof(int)*num_inf);
       //   fprintf(stderr,"Freed permutation %d\n",(int)alloc_size());
       
       mfree(counts,sizeof(int)*max_inc);
       //  fprintf(stderr,"Freed inc_stats %d\n",(int)alloc_size());

       for(k=0;k<num_inf;k++)
	 {
	   mfree(char_alignment[k],sizeof(align_type) * num_taxa);
	 }
       mfree(char_alignment,sizeof(align_type *) * num_inf);
       
       // fprintf(stderr,"Freed char alignment %d\n",(int)alloc_size());
     
       for(k=0;k<num_inf;k++)
	 {
	   mfree(inc_matrix[k],sizeof(inc_type) * num_inf);
	 }
        mfree(inc_matrix,sizeof(inc_type *) * num_inf);
       

	// fprintf(stderr,"Freed matrix alignment %d\n",(int)alloc_size());
       
       for(k=0;k<num_taxa;k++)
	 {
	   // fprintf(stderr,"Freeing taxa %d\n",k);
	   mfree(inf_alignment[k],sizeof(align_type) * num_inf);
	   
	 }
       // fprintf(stderr,"Freeing inf_alignment\n");
       mfree(inf_alignment,sizeof(align_type *) * num_taxa);
       
       // fprintf(stderr,"Freed inf_alignment\n");

       mfree(inf_site_desc,sizeof(site) * num_inf);
       //    fprintf(stderr,"Freed inf_site\n");
       mfree(inf_states,sizeof(int)*num_inf);
       //  fprintf(stderr,"Freed inf_states\n");
       // fprintf(stderr,"Freed informative alignment %d\n",(int)alloc_size());
       
       mfree(site_desc,sizeof(site) * opt.breakWindow);
       mfree(site_states,sizeof(int) * opt.breakWindow);
       
       // fprintf(stderr,"Freed states%d\n",(int)alloc_size());
       
       for(k=0;k<num_taxa;k++)
	 {
	   mfree(new_alignment[k],sizeof(align_type) * opt.breakWindow);
	 }
       mfree(new_alignment,sizeof(align_type *) * num_taxa);
       // fprintf(stderr,"Freed extract%d\n",(int)alloc_size());
       //  fprintf(stdout,"Total allocated after free is %d\n",(int)alloc_size());
      
    }

  fprintf(stdout,"Number of tests performed is %d\n\n",num_tests);
  fprintf(logfile,"Number of tests performed is %d\n\n",num_tests);
  
  vals_file=ffopen("Profile.csv","w");
  for(k=0;k<num_tests;k++)
    {
      fprintf(vals_file,"%d, %8.7le\n",k*opt.stepSize+(opt.breakWindow)/2,values[k]);
    }
  fclose(vals_file);

  fprintf(stdout,"\n");
  fprintf(logfile,"\n");
  fclose(logfile);
  return 0;
  
}
