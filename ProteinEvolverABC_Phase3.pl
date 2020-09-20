###
### ProteinEvolverABC_Phase3.pl, vs 2.0
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
print "\nProteinEvolverABC_Phase3.pl";
print "\nProteinEvolverABC Phase 3 does:";
print "\n- Compute the summary statistics for the simulated data";
print "\n****************************************************************************************************************\n\n";


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
            
            
            
            
		}
	}
close(FROM);



##############################################
##############################################
# Compute SS
print "> Computing summary statistics from the simulated data \"$setsdir/sequences*\".. \n";

if (-e "$setsdir/ProteinEvolverABC_Phase3.sh")  
	{
	print "  \"$setsdir/ProteinEvolverABC_Phase3.sh\" was detected... Ok!\n";			
	}
else  
	{
	print "\nERROR!. There is not a file \"$setsdir/ProteinEvolverABC_Phase3.sh\"\n";
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
system ("\"$setsdir/ProteinEvolverABC_Phase3.sh\"");

#if (glob("*.ss") )  {
#	system ("mv *.ss \"$setsdir\"");
#}


##############################################
##############################################
# Concatenate SS
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
	system ("rm $file_ss");
	}

close (OUTPUT_FILE);
system ("mv \"$setsdir\"/SummaryStatistics.txt \"$setsdir\"/SSsimulations.ss");
print "\n> Summary statistics placed in the file \"$setsdir/SSsimulations.ss\"\n\n";
system("chmod ugo+rwx \"$setsdir\"/*");


##############################################
##############################################
# checking files
print "> Checking output files and preparing files for the ABC estimation phase..\n\n";
system ("perl \"$maindir/scripts/Prepare_Files_for_ABC.2.pl\" \"$file\"");


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
system ("mv \"$setsdir/ProteinEvolverABC_Phase2.sh\" \"$setsdir/OtherOutputs\"");
system ("mv \"$setsdir/ProteinEvolverABC_Phase3.sh\" \"$setsdir/OtherOutputs\"");
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
	
	


