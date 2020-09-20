###
### ProteinEvolverABC_Phase2.pl
###
### Script for reading the Settings file and make the R script to run the ABC analysis for ProteinEvolverProtABC1.2.0
###


use strict;
#use warnings;
use File::Basename;
use File::Copy;
use Cwd;
use Config;
use FindBin;


##############################################
##############################################
# Credits
system("clear");
print "\n****************************************************************************************************************";
print "\nProteinEvolverABC_Phase2.pl";
print "\nProteinEvolverABC Phase 2 does:";
print "\n- Read Settings file";
print "\n- Perform simulations";
print "\n****************************************************************************************************************\n\n";



# Operative system
my $OperativeSystem = "";
$OperativeSystem = $^O;
$OperativeSystem =~ s/\n//g; # remove the end of line
if ($OperativeSystem eq "darwin")
	{
	print "> OS detected: $OperativeSystem (Mac OS X):\n  $Config{osname} \n  $Config{archname}\n\n";	
	}
if ($OperativeSystem eq "linux")
	{
	print "> OS detected: $OperativeSystem:\n  $Config{osname} \n  $Config{archname}\n\n";	
	}
if ($OperativeSystem eq "MSWin32")
	{
	print "> OS detected: $OperativeSystem (Windows):\n  $Config{osname} \n  $Config{archname}\n";
	print "Warning!. The execution of the second step will require a Linux or Mac OS\n";
	print "Type CTRL+C to abort the execution \n\n";
	my $error = <STDIN>;
	chop($error);	
	}


##############################################
##############################################
# Loading input file, "Settings"
my $file = $ARGV[0];
unless (open(FROM,$file))  {
	print STDERR "Cannot open file \"$file\"\n\n";
	exit;
}
#print "> Input file uploaded: $file \n\n";	

my $indels = -1;
my $indels_init = "";
my $numberOfProcessors = 0;
my $numberOfProcessors_init = "";
my $NumberSimulations = 0;
my $number_simulations_init = "";
my $save_simulations = -1;
my $save_simulations_init = "";
my $show_info = -1;
my $show_info_init = "";

open(FROM,$file);
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
                print "  Indels (gaps) are ignored \n";
                }
            elsif ($indels == 1)
                {
                print "  Indels (gaps) are considered as a new state \n";
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
        
		# Number of available processors
		if ($_ =~ /^NumberOfProcessors=/) 
			{
			($numberOfProcessors_init, $numberOfProcessors) = split(/=/, $_);
			$numberOfProcessors =~ s/\n//g; # remove the end of line
			my $thisN = 0;
			$thisN = $numberOfProcessors;
			if ($thisN < 1)
				{
				print "\n\nERROR!. Number of processors must be at least 1 ($thisN)\n\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;
				}
			print "  Number of processors to run the simulations: $numberOfProcessors";
			
			if ($OperativeSystem ne "linux" && $thisN > 1)
				{
				#print "\n   Although the parallelization is only available on Linux OS";
				;
				}
			print "\n";
			}	
		
		# Number of simulations
		if ($_ =~ /^*NumberOfSimulations=/) 
			{
			($number_simulations_init, $NumberSimulations) = split(/=/, $_);
			$NumberSimulations =~ s/\n//g; # remove the end of line
			
			print "  Number of simulations: $NumberSimulations \n";
			}
			
		# Save simulations
		if ($_ =~ /^*SaveSimulations=/) 
			{
			($save_simulations_init, $save_simulations) = split(/=/, $_);
			$save_simulations =~ s/\n//g; # remove the end of line
			if ($save_simulations == 0)
				{
				print "  Simulated data will not be saved \n";	
				}
			elsif ($save_simulations == 1)
				{
				print "  Simulated data will be saved (requires space in the disk) \n";	
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
                print "  Running information (simulations and summary statistics) is not shown on the screen \n";
                }
            elsif ($show_info == 1)
                {
                print "  Running information (simulations and summary statistics) is shown on the screen (slows down the running time) \n";
                }
            else
                {
                print "\n\nERROR!. \"ShowInformationScreen\" must be 0 (Running information -simulations and summary statistics- is not shown) or 1 ( Running information (simulations and summary statistics) is shown ($show_info)\n\n";
                print "Type CTRL+C to abort the execution \n";
                my $HereError = <STDIN>;
                chop($HereError);
                exit;
                }
                
            print "\n";
            }
        
            
			
		}
	}

		

##############################################
##############################################
# Directories
my $maindir = dirname($0);
#print "> ProteinEvolverABC directory detected: $maindir \n\n";
my $setsdir = dirname($file);
#print "> Settings directory detected: $setsdir \n\n";


close(FROM);
system("chmod ugo+rwx \"$setsdir\"/*");


##############################################
##############################################
# Perform simulations
print "> Performing simulations.. \n";
if (-e "$setsdir/ProteinEvolverABC_Phase2.sh")  {
	print "  \"$setsdir/ProteinEvolverABC_Phase2.sh\" was detected... Ok!\n";			
}
else  {
	print "\nERROR!. There is not a file \"$setsdir/ProteinEvolverABC_Phase2.sh\"\n";
	print "Type CTRL+C to abort the execution \n";
	my $error = <STDIN>;
	chop($error);
	exit;
}	
if (-e "$setsdir/ProteinEvolverABC_arguments.txt")  {
	print "  \"$setsdir/ProteinEvolverABC_arguments.txt\" was detected... Ok!\n";			
}
else  {
	print "\nERROR!. There is not a file \"$setsdir/ProteinEvolverABC_arguments.txt\"\n";
	print "Type CTRL+C to abort the execution \n";
	my $error = <STDIN>;
	chop($error);
	exit;
}
system ("\"$setsdir/ProteinEvolverABC_Phase2.sh\"");

#if ( glob("sequences*") )  {
#	system ("mv sequences* \"$setsdir\"");
#}


# Generate ProteinEvolverABC_Phase3.sh and ArgumentsSS.txt for the SS (next step)
print "\n\n> Generating files \"ProteinEvolverABC_Phase3.sh\" and \"ArgumentsSS.txt\" for the phase 3... ";
if ($OperativeSystem eq "linux" && $numberOfProcessors > 1)
	{
	open (OUT_FILE,'>ProteinEvolverABC_Phase3.sh');
    
    if ($indels == 1)
        {
        print OUT_FILE "more \"$setsdir/ArgumentsSS.txt\" | xargs -n 3 -P $numberOfProcessors \"$maindir/scripts/Run_ComputeSSfromSimulatedData.sh\" ";
        }
    else
        {
        print OUT_FILE "more \"$setsdir/ArgumentsSS.txt\" | xargs -n 3 -P $numberOfProcessors \"$maindir/scripts/Run_ComputeSSfromSimulatedData_NoIndels.sh\" ";
        }
	}	
else
	{
	open (OUT_FILE,'>ProteinEvolverABC_Phase3.sh');
    
    if ($indels == 1)
        {
        print OUT_FILE "more \"$setsdir/ArgumentsSS.txt\" | xargs -n 3 \"$maindir/scripts/Run_ComputeSSfromSimulatedData.sh\" ";
        }
    else
        {
        print OUT_FILE "more \"$setsdir/ArgumentsSS.txt\" | xargs -n 3 \"$maindir/scripts/Run_ComputeSSfromSimulatedData_NoIndels.sh\" ";
        }
	}
close (OUT_FILE);

open (OUT_FILE_2,'>ArgumentsSS.txt');

my $numberSim = 0;
for ($numberSim = 1; $numberSim <= $NumberSimulations; $numberSim++)
	{
	if ($numberSim < 10)
		{
		print OUT_FILE_2 "$setsdir/sequences00000000$numberSim $save_simulations $show_info\n";
		}
	if ($numberSim > 9 && $numberSim < 100)
		{
		print OUT_FILE_2 "$setsdir/sequences0000000$numberSim $save_simulations $show_info\n";
		}
	if ($numberSim > 99 && $numberSim < 1000)
		{
		print OUT_FILE_2 "$setsdir/sequences000000$numberSim $save_simulations $show_info\n";
		}
	if ($numberSim > 999 && $numberSim < 10000)
		{
		print OUT_FILE_2 "$setsdir/sequences00000$numberSim $save_simulations $show_info\n";
		}
	if ($numberSim > 9999 && $numberSim < 100000)
		{
		print OUT_FILE_2 "$setsdir/sequences0000$numberSim $save_simulations $show_info\n";
		}
	if ($numberSim > 99999 && $numberSim < 1000000)
		{
		print OUT_FILE_2 "$setsdir/sequences000$numberSim $save_simulations $show_info\n";
		}
	if ($numberSim > 999999 && $numberSim < 10000000)
		{
		print OUT_FILE_2 "$setsdir/sequences00$numberSim $save_simulations $show_info\n";
		}
	if ($numberSim > 9999999 && $numberSim < 100000000)
		{
		print OUT_FILE_2 "$setsdir/sequences0$numberSim $save_simulations $show_info\n";
		}
	if ($numberSim > 99999999 && $numberSim < 1000000000)
		{
		print OUT_FILE_2 "$setsdir/sequences$numberSim $save_simulations $show_info\n";
		}
	if ($numberSim > 999999999)
		{
		print "Warning!. Too many simulations, check ProteinEvolverABC_Phase2.pl\n";
		print "Type CTRL+C to abort the execution \n\n";
		my $error = <STDIN>;
		chop($error);
		}
	}
close (OUT_FILE_2);
system ("chmod ugo+rwx ProteinEvolverABC_Phase3.sh");
system ("chmod ugo+rwx ArgumentsSS.txt");
system ("rm -r Results"); # remove Results folder made by ProteinEvolver

#system ("mv ProteinEvolverABC_Phase3.sh \"$setsdir\"");
#system ("mv ArgumentsSS.txt \"$setsdir\"");
#print "Ok! \n\n";


##############################################
##############################################
#### END ####

print "\n\n> Successful!\n\n";
exit;

