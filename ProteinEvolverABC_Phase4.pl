###
### ProteinEvolverABC_Phase4.pl
###
### Script to compute SS from simulated data
###


use strict;
#use warnings;
use File::Copy;
use File::Basename;
use Cwd;

##############################################
##############################################
# Credits
system("clear");
print "\n****************************************************************************************************************";
print "\nProteinEvolverABC_Phase4.pl";
print "\nProteinEvolverABC Phase 4 does:";
print "\n- Compute the summary statistics for the simulated data";
print "\n****************************************************************************************************************\n\n";

my $Subs_init = "";
my $Subs_dip = "default";
my $SSCPE = 0; # If 0, SSCPE is not run. If 1, SSCPE runs // Here this also works for ProtASR2
my $PDB_init = "";
my $PDB_dip = "default";
my $Chain_init = "";
my $Chain_dip = "default";

my $SeqPDB_init = "";
my $SeqPDB_dip = "default";
my $SSSeqPDB_init = "";
my $SSSeqPDB_dip = "default";

my $SSPDB_init = "";
my $SSPDB_dip = "default";
my $SSChain_init = "";
my $SSChain_dip = "default";

my $NumberOfSS = 0;


##############################################
##############################################
# Loading input file, "Settings"
my $file = $ARGV[0];
unless (open(FROM,$file))
	{
	print STDERR "Cannot open file \"$file\"\n\n";
	exit;
	}
#print "> Input file uploaded: $file \n\n";	


##############################################
##############################################
# Directories
my $maindir = dirname($0);
#print "> ProteinEvolverABC directory detected: $maindir \n\n";
my $setsdir = dirname($file);
#print "> Settings directory detected: $setsdir \n\n";
system("chmod ugo+rwx \"$setsdir\"/*");

if (-e "$setsdir/SSsequences*")  
	{
	system ("rm \"$setsdir\"/SSsequences*"); # at this point no simulated SS should be here			
	}


##############################################
##############################################
# Save sequences?
open(FROM,$file);
my $save_simulations = -1;
my $save_simulations_init = "";
my $indels = -1;
my $indels_init = "";
my $show_info = -1;
my $show_info_init = "";

while (<FROM>) 
	{
	if ($_ =~ /^#/)
		{
		;
		}
	else
		{
        
        # Indels consideration
        if ($_ =~ /^*Indels=/)
            {
            ($indels_init, $indels) = split(/=/, $_);
            $indels =~ s/\n//g; # remove the end of line
            if ($indels == 0)
                {
                print " Indels (gaps) are ignored \n";
                }
            elsif ($indels == 1)
                {
                print " Indels (gaps) are considered as a new state \n";
                }
            else
                {
                print "\n\nERROR!. \"Indels\" must be 0 (ignoring indels) or 1 (indels considered as a new state)($$indels)\n\n";
                print "Type CTRL+C to abort the execution \n";
                my $HereError = <STDIN>;
                chop($HereError);
                exit;
                }
            }
        
            
		# Save simulations
		if ($_ =~ /^*SaveSimulations=/) 
			{
			($save_simulations_init, $save_simulations) = split(/=/, $_);
			$save_simulations =~ s/\n//g; # remove the end of line
			if ($save_simulations == 0)
				{
				print " Simulated data will not be saved \n";	
				}
			elsif ($save_simulations == 1)
				{
				print " Simulated data will be saved (requires space in the disk) \n";	
				}
			else
				{
				print "\n\nERROR!. \"SaveSimulations\" must be 0 (data is not saved) or 1 (data is saved) ($save_simulations)\n\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;	
				}
            }
            
        
        # Show running information (simulations and SSs) on the screen. 0 (No), 1 (Yes, but it slows down the running time)
        if ($_ =~ /^*ShowInformationScreen=/)
            {
            ($show_info_init, $show_info) = split(/=/, $_);
            $show_info =~ s/\n//g; # remove the end of line
            if ($show_info == 0)
                {
                print " Running information (simulations and summary statistics) is not shown on the screen \n";
                }
            elsif ($show_info == 1)
                {
                print " Running information (simulations and summary statistics) is shown on the screen (slows down the running time) \n";
                }
            else
                {
                print "\n\nERROR!. \"ShowInformationScreen\" must be 0 (Running information -simulations and summary statistics- is not shown) or 1 ( Running information (simulations and summary statistics) is shown)($show_info)\n\n";
                print "Type CTRL+C to abort the execution \n";
                my $HereError = <STDIN>;
                chop($HereError);
                exit;
                }
                
            print "\n";
            }
            
        
            ##############
                             
                             
                     # PDB file, for SS
                     if ($_ =~ /^SSpdb=/)
                         {
                         ($SSPDB_init, $SSPDB_dip) = split(/=/, $_);
                         $SSPDB_dip =~ s/\n//g; # remove the end of line
                                 
                         print "  PDB file for summary statistics: $SSPDB_dip\n";

                        if ($SSPDB_dip ne "default")
                             {
                             if (-e "$setsdir/$SSPDB_dip")
                                 {
                                 print "    $setsdir/$SSPDB_dip was detected for summary statistics... Ok!";
                                 }
                            else
                                 {
                                print "\nERROR!. There is not a file \"$setsdir/$SSPDB_dip\" in this directory \n";
                                print "Type CTRL+C to abort the execution \n";
                                my $error = <STDIN>;
                                chop($error);
                                exit;
                                 }
                             }
                         }
                             
                             
                     # Chain of the PDB file, for SS
                     if ($_ =~ /^SSchain=/)
                         {
                         ($SSChain_init, $SSChain_dip) = split(/=/, $_);
                         $SSChain_dip =~ s/\n//g; # remove the end of line
                                         
                         print "  Chain for summary statistics: $SSChain_dip\n";
                         }
                             
                    # Sequence of the PDB file, required for SSCPE, for SS
                    if ($_ =~ /^SSpdbsequence=/)
                        {
                            ($SSSeqPDB_init, $SSSeqPDB_dip) = split(/=/, $_);
                            $SSSeqPDB_dip =~ s/\n//g; # remove the end of line
                            
                            print "  File with the PDB sequence for summary statistics: $SeqPDB_dip\n";
                            ## $NumberOfArguments = $NumberOfArguments + 1; # add arguments -@
                                
                            if ($SSPDB_dip ne "default" && $SSSeqPDB_dip eq "default")
                                {
                                print "\n\nERROR!. The calculation of SS based on a PDB protein structure ($Subs_dip) requires a file with the PDB sequence\n\n";
                                print "Type CTRL+C to abort the execution \n";
                                my $HereError = <STDIN>;
                                chop($HereError);
                                exit;
                            }
                        }
 
		}
	}
close(FROM);



##############################################
##############################################
# Compute SS
print "> Computing summary statistics from the simulated data \"$setsdir/sequences*\".. \n";

if (-e "$setsdir/ProteinEvolverABC_Phase4.sh")
	{
	print "  \"$setsdir/ProteinEvolverABC_Phase4.sh\" was detected... Ok!\n";
	}
else  
	{
	print "\nERROR!. There is not a file \"$setsdir/ProteinEvolverABC_Phase4.sh\"\n";
	print "Type CTRL+C to abort the execution \n";
	my $error = <STDIN>;
	chop($error);
	exit;
	}
		
if (-e "$setsdir/ArgumentsSS.txt")  
	{
	print "  \"$setsdir/ArgumentsSS.txt\" was detected... Ok!\n";			
	}
else  
	{
	print "\nERROR!. There is not a file \"$setsdir/ArgumentsSS.txt\"\n";
	print "Type CTRL+C to abort the execution \n";
	my $error = <STDIN>;
	chop($error);
	exit;
	}

system("chmod ugo+rwx \"$setsdir\"/*");
system ("\"$setsdir/ProteinEvolverABC_Phase4.sh\"");




#if (glob("*.ss") )  {
#	system ("mv *.ss \"$setsdir\"");
#}


##############################################
##############################################
# Concatenate SS # ELENA CORRECTED TO INSERT THE NEW SS
print "\n\n> Concatenating summary statistics from each simulated data .. \n";
open (OUTPUT_FILE,'>',"$setsdir/SummaryStatistics.txt");

my @files_ss = <$setsdir/SSsequences*.ss>;

my $NumReplicate_ss = 0;
foreach my $file_ss (@files_ss) # all ss files
	{	
	#print "\n\n Reading $file_ss \n";		
		
	unless (open(FROM,$file_ss))
		{
		print STDERR "Cannot open file \"$file_ss\"\n\n";
		}
	
	while (<FROM>) # each line
		{
		print OUTPUT_FILE "$_\n";
		}
	close (FROM);
	#system ("rm $file_ss");
	}

close (OUTPUT_FILE);
system ("mv \"$setsdir\"/SummaryStatistics.txt \"$setsdir\"/SSsimulations.ss");
print "\n> Summary statistics placed in the file \"$setsdir/SSsimulations.ss\"\n\n";
system("chmod ugo+rwx \"$setsdir\"/*");


# cleaning from SSCPE 
if ($SSPDB_dip ne "default") # SS from the protein structure
    {
    #system ("rm sequences*_pdb.fas.log");  # already cleaned
    
    system ("rm *.mut_DE.dat");
    system ("rm *.mut_prof_RMSD.dat");
    system ("rm *.mut_RMSD.dat");
    system ("rm *.SSCPE.summary.dat");
    system ("rm *.tnm.summary.dat");

    system ("rm -f *_mutation.dat");
    system ("rm -f *hydrophobicity.dat");
    system ("rm -f *entropy.dat");
        
    system ("rm -rf *.Modes_*");
        
    system ("rm Prot_evol3");
    system ("rm -rf tnm");
    system ("rm structures.in");
    }
    
#system ("rm SSsequences0*.ss");
system ("find . -name 'SSsequences0*.ss' -delete");

##############################################
##############################################
# checking files
print "> Checking output files and preparing files for the ABC estimation phase..\n\n";

if ($SSPDB_dip eq "default") # SS from the protein structure are not calculated
    {
    $NumberOfSS = 16;
    }
else # SS from the protein structure are calculated
    {
    $NumberOfSS = 20;
    }


system ("perl \"$maindir/scripts/Prepare_Files_for_ABC.2.pl\" \"$file\" \"$NumberOfSS\"");


# Zip simulated data
if ($save_simulations == 1)
	{
	print "> Compress simulated data to \"$setsdir/SimulatedData.tar.gz\"\n";
	system ("tar czfv \"$setsdir/SimulatedData.tar.gz\" \"$setsdir\"/sequences*");
	system ("rm \"$setsdir\"/sequences*");
	}
	
# save other outputs
if ( -e "$setsdir/OtherOutputs")
	{
	}
else
	{
	system ("mkdir \"$setsdir/OtherOutputs\"");
	}
system ("mv \"$setsdir/ProteinEvolverABC_Phase3.sh\" \"$setsdir/OtherOutputs\"");
system ("mv \"$setsdir/ProteinEvolverABC_Phase4.sh\" \"$setsdir/OtherOutputs\"");
system ("mv \"$setsdir/ArgumentsSS.txt\" \"$setsdir/OtherOutputs\"");
if (-e "$setsdir/SimulatedData.tar.gz")  
	{
	system ("mv \"$setsdir/SimulatedData.tar.gz\" \"$setsdir/OtherOutputs\"");			
	}
system ("mv \"$setsdir/ProteinEvolverABC_arguments.txt\" \"$setsdir/OtherOutputs\"");
system ("mv \"$setsdir/ProteinEvolverABC_Pop_arguments.txt\" \"$setsdir/OtherOutputs\"");
system ("mv \"$setsdir/SSsimulations.ss\" \"$setsdir/OtherOutputs\"");





##############################################
##############################################
#### END ####

print "> Successful!\n";
exit;
	
	


