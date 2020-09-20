###
### ProteinEvolverABC_Phase1.pl, vs 2.0
###
### Script for reading the Settings file and make the R script to apply the prior distributions
###


use strict;
#use warnings;
use File::Copy;
use File::Basename;
use Cwd;
use Config;
use FindBin;

##############################################
##############################################
# Credits
system("clear");
print "\n****************************************************************************************************************";
print "\nProteinEvolverABC_Phase1.pl";
print "\n- Read Settings file";
print "\n- Compute the summary statistics for the target data";
print "\n- Analyze the evolutionary scenario";
print "\n- Compute the prior distributions (requires R with libraries: lattice, MCMCpack, ape, graphics)";
print "\n****************************************************************************************************************\n\n";


##############################################
##############################################
# Loading input file, "Settings"
my $file = $ARGV[0];
unless (open(FROM,$file))
	{
	print STDERR "Cannot open input file \"$file\"\n\n";
	exit;
	}
print "> Input file uploaded: $file \n\n";	


##############################################
##############################################
# Directories
my $maindir = dirname($0);
print "> ProteinEvolverABC directory detected: $maindir \n\n";
my $setsdir = dirname($file);
print "> Settings directory detected: $setsdir \n\n";

system("chmod ugo+rwx \"$maindir\"/*");
system("chmod ugo+rwx \"$maindir\"/source/*");
system("chmod ugo+rwx \"$maindir\"/scripts/*");
system("chmod ugo+rwx \"$setsdir\"/*");


##############################################
##############################################
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
# open output file and write header
print "> R output file created: \"$setsdir/MakePriorsProteinEvolverABC.r\" ...\n\n";
open (FILE_1,'>',"$setsdir/MakePriorsProteinEvolverABC.r");
print FILE_1 "# R script to apply priors and make executable files \"ProteinEvolverABC_Phase2.sh\" and \"ProteinEvolverABC_arguments.txt\" for ProteinEvolverProtABC1.2.0\n";
print FILE_1 "# By Miguel Arenas\n";
print FILE_1 "# 2020\n";
print FILE_1 "# \n\n";
print FILE_1 "## Distributions included ##\n";
print FILE_1 "# fix			# A value is fixed	(integer or non integer)			# v <- 2 \n";
print FILE_1 "# sample		# Sample (integer): lowest highest 						# v <- sample(0:23,1,replace=T)\n";
print FILE_1 "# unif			# Uniform (non integer): lowest highest; t truncated	# v <- runif(1,0.1,0.9) // v <- runif(1,0.1,0.9) + t\n";
print FILE_1 "# norm			# Normal: mean, sd; t truncated							# v <- rnorm(1,0.9,0.1) // v <- rnorm(1,0.9,0.1) + t\n";
print FILE_1 "# exp			# Exponential: rate; t truncated						# v <- rexp(1,0.4) // v <- rexp(1,0.4) + t\n";
print FILE_1 "# gamma			# Gamma: shape, rate (1/scale); t truncated				# v <- rgamma(1,0.3,0.4) // v <- rgamma(1,0.3,0.4) + t\n";
print FILE_1 "# beta			# Beta: shape1, shape2 (1/scale); t truncated			# v <- rbeta(1,0.3,0.4) // v <- rbeta(1,0.3,0.4) + t\n";
print FILE_1 "# dirichlet		# Dirichlet: alpha (vector)								# v <- rdirichlet(1, c(1,1,1,1))\n";
print FILE_1 "\n";
print FILE_1 "# Load libraries\n";
print FILE_1 "suppressPackageStartupMessages(library(lattice))\n";
print FILE_1 "suppressPackageStartupMessages(library(MCMCpack))\n";
print FILE_1 "suppressPackageStartupMessages(library(ape))\n";
print FILE_1 "suppressPackageStartupMessages(library(graphics))\n";
print FILE_1 "\n";
print FILE_1 "unlink (\"$setsdir/ProteinEvolverABC_arguments.txt\", recursive = FALSE)\n";
print FILE_1 "unlink (\"$setsdir/ProteinEvolverABC_Pop_arguments.txt\", recursive = FALSE)\n";
print FILE_1 "unlink (\"$setsdir/ProteinEvolverABC_Phase2.sh\", recursive = FALSE)\n";
#print FILE_1 "unlink (\"$setsdir/TrueSimulatedValues.txt\", recursive = FALSE)\n";


print FILE_1 "##### GENERAL SETTINGS #####\n";

# variables
my $TotalNumberLines = 0;
my $no = 0;
my $counterLines = 0;
my $number_simulations = 0;
my $number_simulations_init = "";
my $indels = -1;
my $indels_init = "";
my $save_simulations = -1;
my $save_simulations_init = "";
my $show_info = -1;
my $show_info_init = "";
my $hap_dip = 0;
my $hap_dip_init = "";
my $samplesize = 0;
my $popsize = 0;

my $distribution_popsize = "default";
my $distribution_popsize_init = "";
my $distribution_popsize_all = "";
my $value1popsize = 0;
my $value2popsize = 0;

my $growthrate = "default";
my $growthrateV = "default";
my $growthrateYes = 0;
my $distribution_GrowthRate = "default";
my $distribution_GrowthRate_init = "";
my $distribution_GrowthRate_all = "";
my $value1GrowthRate = 0;
my $value2GrowthRate = 0;
my $value3GrowthRate = 0;
my $value4GrowthRate = 0;
my $value5GrowthRate = 0;
my $IsThere_GrowthRate = 0;

my $demogperiods = "default";
my $demogperiods_init = "default";
my $demogperiodsDemesNcte = "default";
my $demogperiodsDemesNcte_init = "default";
my $demogperiodsDemesNvar = "default";
my $demogperiodsDemesNvar_init = "default";

my $datedtips = "default";
my $datedtips_init = "";

my $gtime = 0;
my $gtime_init ="";
my $gtime_all ="";
my $distribution_gtime = "default";
my $value1gtime = 0;
my $value2gtime = 0;


my $Migrationmodel = "default";
my $Migrationmodel_init = "";
my $MigrationmodelYes = 0;

my $migrationrate = "default";
my $migrationrate_init = "";
my $migrationrateYes = 0;

my $convdemes = "default";
my $convdemes_init = "";
my $convdemesYes = 0;

my $distribution_homogRec = "default";
my $distribution_homogRec_init = "";
my $distribution_homogRec_all = "";
my $value1Hrec = 0;
my $value2Hrec = 0;
my $value3Hrec = 0;
my $value4Hrec = 0;
my $value5Hrec = 0;

my $distribution_outgroup = "default";
my $distribution_outgroup_init = "";
my $distribution_outgroup_all = "";
my $value1Houtgroup = 0;
my $value2Houtgroup = 0;
my $value3Houtgroup = 0;
my $value4Houtgroup = 0;
my $value5Houtgroup = 0;


my $sequences = "default";
my $sequencesYes = 0;
my $optseq = 0;
my $optseqYes = 0;
my $mrcagmrca = "default";
my $trees = "default";
my $times = "default";
my $network = "default";
my $breakpoints = "default";
my $printomegas = "default";
my $noisy = "default";
my $noisyV = 0;

my $distribution_subsRate = "default";
my $distribution_subsRate_init = "default";
my $distribution_subsRate_all = "default";
my $value1subsR = 0;
my $value2subsR = 0;
my $value3subsR = 0;
my $value4subsR = 0;
my $value5subsR = 0;

my $Subs_init = "";
my $Subs_dip = "default";

my $distribution_FreqsAA_init = "default";
my $distribution_FreqsAA_all = "default";
my $distribution_FreqsAA1 = "default";
my $value1FreqsAA = 0;
my $value2FreqsAA = 0;
my $value3FreqsAA = 0;
my $value4FreqsAA = 0;
my $value5FreqsAA = 0;
my $value6FreqsAA = 0;
my $value7FreqsAA = 0;
my $value8FreqsAA = 0;
my $value9FreqsAA = 0;
my $value10FreqsAA = 0;
my $value11FreqsAA = 0;
my $value12FreqsAA = 0;
my $value13FreqsAA = 0;
my $value14FreqsAA = 0;
my $value15FreqsAA = 0;
my $value16FreqsAA = 0;
my $value17FreqsAA = 0;
my $value18FreqsAA = 0;
my $value19FreqsAA = 0;
my $value20FreqsAA = 0;

my $distribution_AArateHetSites_init = "";
my $distribution_AArateHetSites_all = "";
my $distribution_AArateHetSites = "default";
my $value1AArateHetSites = 0;
my $value2AArateHetSites = 0;
my $value3AArateHetSites = 0;
my $value4AArateHetSites = 0;
my $value5AArateHetSites = 0;
my $IsThere_AArateHetSites = 0;

my $distribution_AApinv_init = "default";
my $distribution_AApinv_all = "default";
my $distribution_AApinv = "default";
my $value1AAPinv = 0;
my $value2AAPinv = 0;
my $value3AAPinv = 0;
my $value4AAPinv = 0;
my $value5AAPinv = 0;
my $IsThere_AAPinv = 0;


my $MaxRho = 1;


my $NameOfPhylipFile = "";
my $NameOfPhylipFile_init = "";
my $nuc_number = -1;
my $aa_number = -1;


my $numberOfProcessors_init = "default";
my $numberOfProcessors = 1; # default
my $NumberOfArguments = 0;



##############################################
##############################################
# Reading Settings from input file
print "> Reading Settings from input main file ... \n";
while (<FROM>) # for each line
	{
	$TotalNumberLines++;
	}
close (FROM);

open(FROM,$file);
while (<FROM>) 
	{
	$counterLines++;

	if ($_ =~ /^#/)
		{
		;
		}
	else
		{
		
	# Detect Target alignment input file
	if ($_ =~ /^*NameOfPhylipFile/)
		{
		($NameOfPhylipFile_init, $NameOfPhylipFile) = split(/=/, $_);
		$NameOfPhylipFile =~ s/\n//g; # remove the end of line
		print "  Target alignment file: $NameOfPhylipFile\n";
			
		if (-e "$setsdir/$NameOfPhylipFile")
			{
			print "    $setsdir/$NameOfPhylipFile was detected... Ok!";
			}
		else
			{
			print "\nERROR!. There is not a file \"$setsdir/$NameOfPhylipFile\" in this directory \n";
			print "Type CTRL+C to abort the execution \n";
			my $error = <STDIN>;
			chop($error);
			exit;
			}	
		}
			
	# Number of simulations
	if ($_ =~ /^*NumberOfSimulations=/) 
		{
		($number_simulations_init, $number_simulations) = split(/=/, $_);
		$number_simulations =~ s/\n//g; # remove the end of line
		print "\n  Number of simulations: $number_simulations\n";
		}

            
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
			#print "\n   Although the parallelization is only available on Linux and Mac";
			;
			}
		print "\n";
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
            print "\n\nERROR!. \"ShowInformationScreen\" must be 0 (Running information -simulations and summary statistics- is not shown) or 1 ( Running information (simulations and summary statistics) is shown)($show_info)\n\n";
            print "Type CTRL+C to abort the execution \n";
            my $HereError = <STDIN>;
            chop($HereError);
            exit;
            }
        }
            
            
	
	### Demographic Settings
	# Haploid/diploid
	if ($_ =~ /^*Haploid\/Diploid=/) 
		{
		($hap_dip_init, $hap_dip) = split(/=/, $_);
		$hap_dip =~ s/\n//g; # remove the end of line
		if ($hap_dip < 1 || $hap_dip > 2)
			{
			print "\n\nERROR!. Haploid/Diploid must be 1 or 2 ($hap_dip)\n\n";
			print "Type CTRL+C to abort the execution \n";
			my $HereError = <STDIN>;
			chop($HereError);
			exit;
			}
		print "  Haploid/diploid: $hap_dip\n";
		
		if ($hap_dip == 1)
			{
			$MaxRho	= $MaxRho * 2;
			}
		if ($hap_dip == 2)
			{
			$MaxRho	= $MaxRho * 4;
			}	
		}
	
	# Population size
	if ($_ =~ /^*PopulationSize=/)
		{	
		my $ThisLineHere = $_;
		($distribution_popsize_init, $distribution_popsize_all) = split(/=/, $_);
		$distribution_popsize_all =~ s/\n//g; # remove the end of line	
			
		#($distribution_popsize, $value1popsize, $value2popsize) = split(/ /, $distribution_popsize_all); # this is only available now
		$distribution_popsize =~ s/\n//g;
		$value1popsize =~ s/\n//g;
		$value2popsize =~ s/\n//g;
		$distribution_popsize = "fix"; # this is only available now
		$value1popsize = $distribution_popsize_all; # this is only available now
		
		if ($distribution_popsize eq "fix")
			{
			print "  Population size according to: $distribution_popsize $value1popsize\n";
			$MaxRho = $MaxRho * $value1popsize;
			if ($value1popsize < 2)
				{
				print "\n\nERROR!. The population size is too small ($value1popsize)\n\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;
				}	
			}
		elsif ($distribution_popsize eq "uniform")
			{
			print "  Population size according to: $distribution_popsize $value1popsize $value2popsize\n";
			$MaxRho = $MaxRho * $value2popsize;
			if ($value1popsize < 2 || $value2popsize < 2 )
				{
				print "\n\nERROR!. The population size is too small ($value1popsize)\n\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;
				}
			}
		else
			{
			print "\n\nERROR!. Population size must be specified (unif # # or fix #)\n";
			print "Type CTRL+C to abort the execution \n";
			my $HereError = <STDIN>;
			chop($HereError);
			exit;
			}	
		$NumberOfArguments = $NumberOfArguments + 2; # add arguments -e1000 2	
		}
	
	# GrowthRate
	if ($_ =~ /^GrowthRate=/) 
		{
		if ($_ =~ /^\s*$/) # blank line
			{
			#print "   GrowthRate: not specified\n";
			}
		else
			{
			my $ThisLineHere = $_;
			($distribution_GrowthRate_init, $distribution_GrowthRate_all) = split(/=/, $_);
			$distribution_GrowthRate_all =~ s/\n//g; # remove the end of line		
				
			($distribution_GrowthRate, $value1GrowthRate, $value2GrowthRate, $value3GrowthRate, $value4GrowthRate, $value5GrowthRate) = split(/ /, $distribution_GrowthRate_all); # [fix, unif, norm, exp, gamma, beta]
			$distribution_GrowthRate =~ s/\n//g;
			$value1GrowthRate =~ s/\n//g;
			$value2GrowthRate =~ s/\n//g;
			$value3GrowthRate =~ s/\n//g;
			$value4GrowthRate =~ s/\n//g;
			$value5GrowthRate =~ s/\n//g;
					
			if ($distribution_GrowthRate eq "fix")
				{
				$NumberOfArguments = $NumberOfArguments + 2; # add arguments -g0 2.3e-03
				print "  Growth rate according to: $distribution_GrowthRate $value1GrowthRate\n";
				$IsThere_GrowthRate++;
				}
			elsif ($distribution_GrowthRate eq "uniform")
				{
				$NumberOfArguments = $NumberOfArguments + 2; # add arguments -g0 2.3e-03	
				print "  Growth rate according to: $distribution_GrowthRate $value1GrowthRate $value2GrowthRate\n";
				$IsThere_GrowthRate++;
				}
			elsif ($distribution_GrowthRate eq "normal" || $distribution_GrowthRate eq "gamma" || $distribution_GrowthRate eq "beta")
				{
				$NumberOfArguments = $NumberOfArguments + 2; # add arguments -g0 2.3e-03
				if ($value3GrowthRate eq "t")
					{
					print "  Growth rate according to: $distribution_GrowthRate $value1GrowthRate $value2GrowthRate truncated ($value3GrowthRate) at $value4GrowthRate and $value5GrowthRate\n";
					if ($value4GrowthRate > $value5GrowthRate)
						{
						print "\n     Error in Growth rate. Uncorrectly truncated (lowest > highest)!\n\n";
						print "Type CTRL+C to abort the execution \n";
						my $HereError = <STDIN>;
						chop($HereError);
						exit;	
						}
					}
				else
					{
					print "  Growth rate according to: $distribution_GrowthRate $value1GrowthRate $value2GrowthRate\n";
					}
				$IsThere_GrowthRate++;
				}
			elsif ($distribution_GrowthRate eq "exponential") 
				{
				$NumberOfArguments = $NumberOfArguments + 2; # add arguments -g0 2.3e-03
				if ($value2GrowthRate eq "t")
					{
					print "  Growth rate according to: $distribution_GrowthRate $value1GrowthRate truncated ($value2GrowthRate) at $value3GrowthRate and $value4GrowthRate\n";
					if ($value3GrowthRate > $value4GrowthRate)
						{
						print "\n     Error in Growth rate. Uncorrectly truncated (lowest > highest)!\n\n";
						print "Type CTRL+C to abort the execution \n";
						my $HereError = <STDIN>;
						chop($HereError);
						exit;	
						}
					}
				else
					{
					print "  Growth rate according to: $distribution_GrowthRate $value1GrowthRate\n";
					}
				$IsThere_GrowthRate++;	
				}
			else
				{
				#print "\nWarning!. Growth rate was not specified ($distribution_GrowthRate) (fix #; unif # #; norm # # t # #; exp # t # #; gamma # # t # #; beta # # t # #)\n";
				}	
			}	
						
		}

		
	# Demographic periods
	if ($_ =~ /^DemographicPeriods=/) 
		{
		if ($_ =~ /^\s*$/) # blank line
			{
			#print "   Demographic periods: not specified\n";
			}
		else
			{
			($demogperiods_init, $demogperiods) = split(/=/, $_);
			$demogperiods =~ s/\n//g; # remove the end of line								
			my $ss = $demogperiods;
			$ss =~ s/g//g;
			print "  Demographic periods: $ss\n";
			if ($value1GrowthRate > 0)
				{
				print "\n\nERROR!. Choose demographic periods or growth rate, not both!\n\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;
				}
			
			my @partitions_demogperiods = split(/\s+/, $ss);
			my $numberOfelements_demogperiods = scalar @partitions_demogperiods;			
			#print "		numberOfelements_demogperiods = $numberOfelements_demogperiods \n";
			$NumberOfArguments = $NumberOfArguments + $numberOfelements_demogperiods; # add arguments -g1 numbers
			}
		}		
		
	# Demographic periods in populations/species tree with constant N
	if ($_ =~ /^DemographicsSpeciesTree_ConstantPops=/) 
		{
		if ($_ =~ /^\s*$/) # [ or blank line
			{
			#print "   Populations/species tree with constant population sizes: not specified\n";
			}
		else
			{
			($demogperiodsDemesNcte_init, $demogperiodsDemesNcte) = split(/=/, $_);	
			$demogperiodsDemesNcte =~ s/\n//g;
			my $ss = $demogperiodsDemesNcte;
			$ss =~ s/g//g;
			print "   Demographic periods with constant population sizes: $ss\n";
			if (($value1GrowthRate > 0) || ($demogperiods ne "default"))
				{
				print "\n\nERROR!. Choose only one demographic option.\n\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;
				}
			
			my @partitions_demogperiodsNcte = split(/\s+/, $ss);
			my $numberOfelements_demogperiodsNcte = scalar @partitions_demogperiodsNcte;			
			$NumberOfArguments = $NumberOfArguments + $numberOfelements_demogperiodsNcte; # add arguments -g2 numbers
			}
		}
		
	# Demographic periods in populations/species tree with variable N
	if ($_ =~ /^DemographicsSpeciesTree_VariablePops=/) 
		{
		if ($_ =~ /^\s*$/) # [ or blank line
			{
			#print "   Populations/species tree with variable population sizes: not specified\n";
			}
		else
			{
			($demogperiodsDemesNvar_init, $demogperiodsDemesNvar) = split(/=/, $_);
			$demogperiodsDemesNvar =~ s/\n//g;
			my $ss = $demogperiodsDemesNvar;
			$ss =~ s/g//g;
			print "   Demographic periods with variable population sizes: $ss\n";
			if (($value1GrowthRate > 0) || ($demogperiods ne "default") || ($demogperiodsDemesNcte ne "default"))
				{
				print "\n\nERROR!. Choose only one demographic option.\n\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;
				}
				
			my @partitions_demogperiodsNvar = split(/\s+/, $ss);
			my $numberOfelements_demogperiodsNvar = scalar @partitions_demogperiodsNvar;			
			$NumberOfArguments = $NumberOfArguments + $numberOfelements_demogperiodsNvar; # add arguments -g3 numbers
			}
		}


	# Dated tips
	if ($_ =~ /^DatedTips=/)
		{
		if ($_ =~ /^\s*$/) # [ or blank line
			{
			#print "   Dated tips: not specified\n";
			}
		else
			{
			($datedtips_init, $datedtips) = split(/=/, $_);
			$datedtips =~ s/\n//g;
			my $ss = $datedtips;
			$ss =~ s/=//g;
			print "  Dated tips: $ss\n";
			
			my @partitions_datedtips = split(/\s+/, $ss);
			my $numberOfelements_datedtips = scalar @partitions_datedtips;			
			$NumberOfArguments = $NumberOfArguments + $numberOfelements_datedtips; # add arguments -= numbers
			}
		}
	
	# Generation time
	if ($_ =~ /^GenerationTime=/)
		{		
		my $ThisLineHere = $_;
		($gtime_init, $gtime_all) = split(/=/, $_);
		$gtime_all =~ s/\n//g; # remove the end of line	
			
		($distribution_gtime, $value1gtime, $value2gtime) = split(/ /, $gtime_all);
		$distribution_gtime =~ s/\n//g;
		$value1gtime =~ s/\n//g;
		$value2gtime =~ s/\n//g;
		if ($distribution_gtime eq "fix")
			{
			$NumberOfArguments = $NumberOfArguments + 1; # add arguments -/
			print "  Generation time according to: $distribution_gtime $value1gtime\n";
			if ($value1popsize < 1)
				{
				print "\n\nERROR!. Generation time is too small ($value1gtime)\n\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;
				}	
			}
		elsif ($distribution_gtime eq "uniform")
			{
			$NumberOfArguments = $NumberOfArguments + 1; # add arguments -/
			print "  Generation time according to: $distribution_gtime $value1gtime $value2gtime\n";
			if ($value1gtime < 1 || $value2gtime < 1 )
				{
				print "\n\nERROR!. Generation time is too small ($value1gtime)\n\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;
				}
			}
		else
			{
			print "\n\nERROR!. Generation time must be correctly specified (unif # # or fix #)\n";
			print "Type CTRL+C to abort the execution \n";
			my $HereError = <STDIN>;
			chop($HereError);
			exit;
			}	
		}	
	
	
	# Migration model
	if ($_ =~ /^MigrationModel=/) 
		{
		if ($_ =~ /^\s*$/)
			{
			#print "   Migration model: not specified\n";
			}
		else
			{
			($Migrationmodel_init, $Migrationmodel) = split(/=/, $_);			
			$Migrationmodel =~ s/\n//g;
			my $ss = $Migrationmodel;
			$ss =~ s/q//g;
			print "  Migration model: $ss\n";
			$MigrationmodelYes++;
			
			if (($demogperiodsDemesNcte ne "default") || ($demogperiodsDemesNvar ne "default"))
				{
				print "\n\nERROR!. Demographic in demes requires a migration model.\n\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;
				}
				
			my @partitions_Migrationmodel = split(/\s+/, $ss);
			my $numberOfelements_Migrationmodel = scalar @partitions_Migrationmodel;			
			$NumberOfArguments = $NumberOfArguments + $numberOfelements_Migrationmodel; # add arguments -q numbers
			}
		}

	# Migration rate
	if ($_ =~ /^MigrationRate=/) 
		{
		if ($_ =~ /^\s*$/)
			{
			#print "   Migration rate: not specified\n";
			}
		else
			{
			($migrationrate_init, $migrationrate) = split(/=/, $_);			
			$migrationrate =~ s/\n//g;
			my $ss = $migrationrate;
			$ss =~ s/t//g;
			print "  Migration rate: $ss\n";
			$migrationrateYes++;
			
			if (($demogperiodsDemesNcte ne "default") || ($demogperiodsDemesNvar ne "default"))
				{
				print "\n\nERROR!. Demographic in demes requires a migration rate.\n\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;
				}	
				
			my @partitions_migrationrate = split(/\s+/, $ss);
			my $numberOfelements_migrationrate = scalar @partitions_migrationrate;			
			$NumberOfArguments = $NumberOfArguments + $numberOfelements_migrationrate; # add arguments -t numbers		
			}
		}

	# Species/populationsTree - Convergence of demes
	if ($_ =~ /^Species\/populationsTree=/) 
		{
		if ($_ =~ /^\s*$/)
			{
			#print "   Convergence of demes: not specified\n";
			}
		else
			{
			($convdemes_init, $convdemes) = split(/=/, $_);			
			$convdemes =~ s/\n//g;
			my $ss = $convdemes;
			$ss =~ s/%//g;
			print "  Convergence of demes: $ss\n";
			$convdemesYes++;
			if (($MigrationmodelYes == 0) || ($migrationrateYes == 0))
				{
				print "\n\nERROR!. Species/populationsTree requires a migration model and migration rate.\n\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;
				}
				
			my @partitions_convdemes = split(/\s+/, $ss);
			my $numberOfelements_convdemes = scalar @partitions_convdemes;			
			$NumberOfArguments = $NumberOfArguments + $numberOfelements_convdemes; # add arguments -% numbers
			}
		}


	# Recombination
	if ($_ =~ /^*RecombinationRate=/) 
		{
		if ($_ =~ /^\s*$/) # [ or blank line
			{
			#print "   Recombination rate: not specified\n";			
			}
		else
			{
			($distribution_homogRec_init, $distribution_homogRec_all) = split(/=/, $_);			
			
			($distribution_homogRec, $value1Hrec, $value2Hrec, $value3Hrec, $value4Hrec, $value5Hrec) = split(/ /, $distribution_homogRec_all); # [fix, unif, norm, exp, gamma, beta]
			$distribution_homogRec =~ s/\n//g;
			$value1Hrec =~ s/\n//g;
			$value2Hrec =~ s/\n//g;
			$value3Hrec =~ s/\n//g;
			$value4Hrec =~ s/\n//g;
			$value5Hrec =~ s/\n//g;
				
			if ($distribution_homogRec eq "fix")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -r
				print "  Recombination rate according to: $distribution_homogRec $value1Hrec\n";
				$MaxRho = $MaxRho * $value1Hrec;
				}
			elsif ($distribution_homogRec eq "uniform")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -r
				print "  Recombination rate according to: $distribution_homogRec $value1Hrec $value2Hrec\n";
				$MaxRho = $MaxRho * $value2Hrec;
				}
			elsif ($distribution_homogRec eq "normal" || $distribution_homogRec eq "gamma" || $distribution_homogRec eq "beta")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -r
				if ($value3Hrec eq "t")
					{
					print "  Recombination rate according to: $distribution_homogRec $value1Hrec $value2Hrec truncated ($value3Hrec) at $value4Hrec and $value5Hrec\n";
					$MaxRho = $MaxRho * $value5Hrec;
					}
				else
					{
					print "  Recombination rate according to: $distribution_homogRec $value1Hrec $value2Hrec\n";
					$MaxRho = 0;
					}
					
				if ($value4Hrec > $value5Hrec)
					{
					print "\n     Error in recombination rate. Uncorrectly truncated (lowest > highest)!\n\n";
					print "Type CTRL+C to abort the execution \n";
					my $HereError = <STDIN>;
					chop($HereError);
					exit;	
					}
				}
			elsif ($distribution_homogRec eq "exponential") 
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -r
				if ($value2Hrec eq "t")
					{
					print "  Recombination rate according to: $distribution_homogRec $value1Hrec truncated ($value2Hrec) at $value3Hrec and $value4Hrec\n";
					$MaxRho = $MaxRho * $value4Hrec;
					}
				else
					{
					print "  Recombination rate according to: $distribution_homogRec $value1Hrec\n";
					$MaxRho = 0;
					}	
				
				if ($value3Hrec > $value4Hrec)
					{
					print "\n     Error in recombination rate. Uncorrectly truncated (lowest > highest)!\n\n";
					print "Type CTRL+C to abort the execution \n";
					my $HereError = <STDIN>;
					chop($HereError);
					exit;	
					}
					
				}
			else
				{
				#print "\n\nWarning!. Recombination rate was not specified ($distribution_homogRec) (fix #; unif # #; norm # # t # #; exp # t # #; gamma # # t # #; beta # # t # #)\n";
				}	
			}	
		}
	
	
	# Outgroup
	if ($_ =~ /^Outgroup=/) 
		{
		if ($_ =~ /^\s*$/) # [ or blank line
			{
			#print "   Outgroup: not specified\n";			
			}
		else
			{
			($distribution_outgroup_init, $distribution_outgroup_all) = split(/=/, $_);			
							
			($distribution_outgroup, $value1Houtgroup, $value2Houtgroup, $value3Houtgroup, $value4Houtgroup, $value5Houtgroup) = split(/ /, $distribution_outgroup_all); # [fix, unif, norm, exp, gamma, beta]
			$distribution_outgroup =~ s/\n//g;
			$value1Houtgroup =~ s/\n//g;
			$value2Houtgroup =~ s/\n//g;
			$value3Houtgroup =~ s/\n//g;
			$value4Houtgroup =~ s/\n//g;
			$value5Houtgroup =~ s/\n//g;
		
			if ($distribution_outgroup eq "fix")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -o
				print "  Outgroup according to: $distribution_outgroup $value1Houtgroup\n";
				}
			elsif ($distribution_outgroup eq "uniform")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -o
				print "  Outgroup according to: $distribution_outgroup $value1Houtgroup $value2Houtgroup\n";
				}
			elsif ($distribution_outgroup eq "normal" || $distribution_outgroup eq "gamma" || $distribution_outgroup eq "beta")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -o
				if ($value3Houtgroup eq "t")
					{
					print "  Outgroup according to: $distribution_outgroup $value1Houtgroup $value2Houtgroup truncated ($value3Houtgroup) at $value4Houtgroup and $value5Houtgroup\n";
					}
				else
					{
					print "  Outgroup according to: $distribution_outgroup $value1Houtgroup $value2Houtgroup\n";
					}
					
				if ($value4Houtgroup > $value5Houtgroup)
					{
					print "\n     Error in outgroup. Uncorrectly truncated (lowest > highest)!\n\n";
					print "Type CTRL+C to abort the execution \n";
					my $HereError = <STDIN>;
					chop($HereError);
					exit;	
					}
				}
			elsif ($distribution_outgroup eq "exponential") 
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -o
				if ($value2Houtgroup eq "t")
					{
					print "  Outgroup according to: $distribution_outgroup $value1Houtgroup truncated ($value2Houtgroup) at $value3Houtgroup and $value4Houtgroup\n";
					}
				else
					{
					print "  Outgroup according to: $distribution_outgroup $value1Houtgroup\n";
					}	
					
				if ($value3Houtgroup > $value4Houtgroup)
					{
					print "\n     Error in outgroup. Uncorrectly truncated (lowest > highest)!\n\n";
					print "Type CTRL+C to abort the execution \n";
					my $HereError = <STDIN>;
					chop($HereError);
					exit;	
					}
					
				}
			else
				{
				#print "\n\nWarning!. Outgroup was not specified ($distribution_outgroup) (fix #; unif # #; norm # # t # #; exp # t # #; gamma # # t # #; beta # # t # #)\n";
				}	
			}	
		}	
							
	# Print sumMRCA/GMRCA
	#$mrcagmrca = "$";
		
	# Print trees
	#$trees = "trees";

	# Print times
	#$times = "times";
	
	# Print network
	#$network = "NetworkFile";
			
	# Print breakpoints
	#$breakpoints = "breakpoints";
	
	# Noisy
	$noisyV = 1;	
	
	
	### Amino acid substitution model
	# substitution rate
	if ($_ =~ /^*SubstitutionRate=/) 
		{
		if ($_ =~ /^\s*$/) # [ or blank line
			{
			#print "   Substitution rate: not specified\n";
			}
		else
			{
			($distribution_subsRate_init, $distribution_subsRate_all) = split(/=/, $_);			
			
			
			($distribution_subsRate, $value1subsR, $value2subsR, $value3subsR, $value4subsR, $value5subsR) = split(/ /, $distribution_subsRate_all); # [fix, unif, norm, exp, gamma, beta]
			$distribution_subsRate =~ s/\n//g;
			$value1subsR =~ s/\n//g;
			$value2subsR =~ s/\n//g;
			$value3subsR =~ s/\n//g;
			$value4subsR =~ s/\n//g;
			$value5subsR =~ s/\n//g;
				
			if ($distribution_subsRate eq "fix")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -u
				print "  Substitution rate according to: $distribution_subsRate $value1subsR\n";
				}
			elsif ($distribution_subsRate eq "uniform")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -u
				print "  Substitution rate according to: $distribution_subsRate $value1subsR $value2subsR\n";
				}
			elsif ($distribution_subsRate eq "normal" || $distribution_subsRate eq "gamma" || $distribution_subsRate eq "beta")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -u
				if ($value3subsR eq "t")
					{
					print "  Substitution rate according to: $distribution_subsRate $value1subsR $value2subsR truncated ($value3subsR) at $value4subsR and $value5subsR\n";
					}
				else
					{
					print "  Substitution rate according to: $distribution_subsRate $value1subsR $value2subsR\n";
					}
					
				if ($value4subsR > $value5subsR)
					{
					print "\n     Error in substitution rate. Uncorrectly truncated (lowest > highest)!\n\n";
					print "Type CTRL+C to abort the execution \n";
					my $HereError = <STDIN>;
					chop($HereError);
					exit;	
					}
				}
			elsif ($distribution_subsRate eq "exponential") 
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -u
				if ($value2subsR eq "t")
					{
					print "  Substitution rate according to: $distribution_subsRate $value1subsR truncated ($value2subsR) at $value3subsR and $value4subsR\n";
					}
				else
					{
					print "  Substitution rate according to: $distribution_subsRate $value1subsR\n";
					}	
					
				if ($value3subsR > $value4subsR)
					{
					print "\n     Error in substitution rate. Uncorrectly truncated (lowest > highest)!\n\n";
					print "Type CTRL+C to abort the execution \n";
					my $HereError = <STDIN>;
					chop($HereError);
					exit;	
					}	
				}
			else
				{
				print "\n\nERROR!. Substitution rate was not specified ($distribution_subsRate) (fix #; unif # #; norm # # t # #; exp # t # #; gamma # # t # #; beta # # t # #)\n";
				print "Type CTRL+C to abort the execution \n";
				my $HereError = <STDIN>;
				chop($HereError);
				exit;
				}	
						
			}	
						
		}

    
            
    # Substituion model of amino acid evolution (i.e., Blosum62, CpRev, Dayhoff, DayhoffDCMUT, HIVb, HIVw, JTT, JonesDCMUT, LG, Mtart, Mtmam, Mtrev24, RtRev, VT, WAG, UserEAAM)
    if ($_ =~ /^*SubstitutionModel=/)
        {
        ($Subs_init, $Subs_dip) = split(/=/, $_);
        $Subs_dip =~ s/\n//g; # remove the end of line
        if ($Subs_dip ne "Blosum62" && $Subs_dip ne "CpRev" && $Subs_dip ne "Dayhoff" && $Subs_dip ne "DayhoffDCMUT" && $Subs_dip ne "HIVb" && $Subs_dip ne "HIVw" && $Subs_dip ne "JTT" && $Subs_dip ne "JonesDCMUT" && $Subs_dip ne "LG" && $Subs_dip ne "Mtart" && $Subs_dip ne "Mtmam" && $Subs_dip ne "Mtrev24" && $Subs_dip ne "RtRev" && $Subs_dip ne "VT" && $Subs_dip ne "WAG" && $Subs_dip ne "UserEAAM")
            {
            print "\n\nERROR!. Substitution model must be one of the following: Blosum62, CpRev, Dayhoff, DayhoffDCMUT, HIVb, HIVw, JTT, JonesDCMUT, LG, Mtart, Mtmam, Mtrev24, RtRev, VT, WAG, UserEAAM ($Subs_dip)\n\n";
            print "Type CTRL+C to abort the execution \n";
            my $HereError = <STDIN>;
            chop($HereError);
            exit;
            }
        print "  Substitution model: $Subs_dip\n";
        $NumberOfArguments = $NumberOfArguments + 1; # add arguments -@

        }
            
            

	# Amino acid frequencies. fix or dirichlet. By default equally distributed frequencies. i.e., dirichlet 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
	if ($_ =~ /^*AminoacidFrequencies=/)
		{
		($distribution_FreqsAA_init, $distribution_FreqsAA_all) = split(/=/, $_);			
	
		($distribution_FreqsAA1, $value1FreqsAA, $value2FreqsAA, $value3FreqsAA, $value4FreqsAA, $value5FreqsAA, $value6FreqsAA, $value7FreqsAA, $value8FreqsAA, $value9FreqsAA, $value10FreqsAA, $value11FreqsAA, $value12FreqsAA, $value13FreqsAA, $value14FreqsAA, $value15FreqsAA, $value16FreqsAA, $value17FreqsAA, $value18FreqsAA, $value19FreqsAA, $value20FreqsAA) = split(/ /, $distribution_FreqsAA_all);
		$distribution_FreqsAA1 =~ s/\n//g;
		$value1FreqsAA =~ s/\n//g;
		$value2FreqsAA =~ s/\n//g;
		$value3FreqsAA =~ s/\n//g;
		$value4FreqsAA =~ s/\n//g;
		$value5FreqsAA =~ s/\n//g;
		$value6FreqsAA =~ s/\n//g;
		$value7FreqsAA =~ s/\n//g;
		$value8FreqsAA =~ s/\n//g;
		$value9FreqsAA =~ s/\n//g;
		$value10FreqsAA =~ s/\n//g;
		$value11FreqsAA =~ s/\n//g;
		$value12FreqsAA =~ s/\n//g;
        $value13FreqsAA =~ s/\n//g;
        $value14FreqsAA =~ s/\n//g;
        $value15FreqsAA =~ s/\n//g;
        $value16FreqsAA =~ s/\n//g;
        $value17FreqsAA =~ s/\n//g;
        $value18FreqsAA =~ s/\n//g;
        $value19FreqsAA =~ s/\n//g;
        $value20FreqsAA =~ s/\n//g;
		
		if ($distribution_FreqsAA1 eq "fix" || $distribution_FreqsAA1 eq "dirichlet")
			{
			print "  Amino acid frequencies (1x20): $distribution_FreqsAA1 $value1FreqsAA $value2FreqsAA $value3FreqsAA $value4FreqsAA $value5FreqsAA $value6FreqsAA $value7FreqsAA $value8FreqsAA $value9FreqsAA $value10FreqsAA $value11FreqsAA $value12FreqsAA $value13FreqsAA $value14FreqsAA $value15FreqsAA $value16FreqsAA $value17FreqsAA $value18FreqsAA $value19FreqsAA $value20FreqsAA\n";
			}
		else
			{
			#print "   Warning!. Amino acid frequencies (1x20) were not specified - ($distribution_FreqsAA1), by default equal frequencies. However maybe you use CAT frequencies.\n";
			$distribution_FreqsAA1 = "fix";
			$value1FreqsAA = 0.05;
			$value2FreqsAA = 0.05;
			$value3FreqsAA = 0.05;
			$value4FreqsAA = 0.05;
			$value5FreqsAA = 0.05;
			$value6FreqsAA = 0.05;
			$value7FreqsAA = 0.05;
			$value8FreqsAA = 0.05;
			$value9FreqsAA = 0.05;
			$value10FreqsAA = 0.05;
			$value11FreqsAA = 0.05;
			$value12FreqsAA = 0.05;
            $value13FreqsAA = 0.05;
            $value14FreqsAA = 0.05;
            $value15FreqsAA = 0.05;
            $value16FreqsAA = 0.05;
            $value17FreqsAA = 0.05;
            $value18FreqsAA = 0.05;
            $value19FreqsAA = 0.05;
            $value20FreqsAA = 0.05;
			}
		
		}
				
	
	# rate of heterogeneity among sites (AA model)
	if ($_ =~ /^RateHetSites=/) 
		{
		
		if ($_ =~ /^\s*$/) # [ or blank line
			{
			#print "    Rate of heterogeneity among sites (+G): not specified\n";
			}
		else
			{
			($distribution_AArateHetSites_init, $distribution_AArateHetSites_all) = split(/=/, $_);
				
			($distribution_AArateHetSites, $value1AArateHetSites, $value2AArateHetSites, $value3AArateHetSites, $value4AArateHetSites, $value5AArateHetSites) = split(/ /, $distribution_AArateHetSites_all); # [fix, unif, norm, exp, gamma, beta]
			$distribution_AArateHetSites =~ s/\n//g;
			$value1AArateHetSites =~ s/\n//g;
			$value2AArateHetSites =~ s/\n//g;
			$value3AArateHetSites =~ s/\n//g;
			$value4AArateHetSites =~ s/\n//g;
			$value5AArateHetSites =~ s/\n//g;
		
			if ($distribution_AArateHetSites eq "fix")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -a#
				print "  Rate of heterogeneity among sites (+G) according to: $distribution_AArateHetSites $value1AArateHetSites\n";
				$IsThere_AArateHetSites++;
				}
			elsif ($distribution_AArateHetSites eq "uniform")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -a#
				print "  Rate of heterogeneity among sites (+G) according to: $distribution_AArateHetSites $value1AArateHetSites $value2AArateHetSites\n";
				$IsThere_AArateHetSites++;
				}
			elsif ($distribution_AArateHetSites eq "normal" || $distribution_AArateHetSites eq "gamma" || $distribution_AArateHetSites eq "beta")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -a#
				if ($value3AArateHetSites eq "t")
					{
					print "  Rate of heterogeneity among sites (+Gsites) according to: $distribution_AArateHetSites $value1AArateHetSites $value2AArateHetSites truncated ($value3AArateHetSites) at $value4AArateHetSites and $value5AArateHetSites\n";
					if ($value4AArateHetSites > $value5AArateHetSites)
						{
						print "\n      Error in rate of heterogeneity among sites (+G). Uncorrectly truncated (lowest > highest)!\n\n";
						print "Type CTRL+C to abort the execution \n";
						my $HereError = <STDIN>;
						chop($HereError);
						exit;	
						}
					}
				else
					{
					print "  Rate of heterogeneity among sites (+G) according to: $distribution_AArateHetSites $value1AArateHetSites $value2AArateHetSites\n";
					}
				$IsThere_AArateHetSites++;
				}
			elsif ($distribution_AArateHetSites eq "exponential") 
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -a#
				if ($value2AArateHetSites eq "t")
					{
					print "  Rate of heterogeneity among sites (+G) according to: $distribution_AArateHetSites $value1AArateHetSites truncated ($value2AArateHetSites) at $value3AArateHetSites and $value4AArateHetSites\n";
					if ($value3AArateHetSites > $value4AArateHetSites)
						{
						print "\n      Error in rate of heterogeneity among sites (+G). Uncorrectly truncated (lowest > highest)!\n\n";
						print "Type CTRL+C to abort the execution \n";
						my $HereError = <STDIN>;
						chop($HereError);
						exit;	
						}
					}
				else
					{
					print "  Rate of heterogeneity among sites (+G) according to: $distribution_AArateHetSites $value1AArateHetSites\n";
					}
				$IsThere_AArateHetSites++;	
				}
			else
				{
				#print "\n Warning!. Rate of heterogeneity among sites (+G) was not specified ($distribution_AArateHetSites) (fix #; unif # #; norm # # t # #; exp # t # #; gamma # # t # #; beta # # t # #)\n";
				}	
			}	
		}


	# proportion of invariable sites (AA model)
	if ($_ =~ /^PropInvSites=/) 
		{
		if ($_ =~ /^\s*$/) # blank line
			{
			#print "   Proportion of invariable sites (+I) (AA model): not specified\n";
			}
		else
			{
			($distribution_AApinv_init, $distribution_AApinv_all) = split(/=/, $_);	
				
			($distribution_AApinv, $value1AAPinv, $value2AAPinv, $value3AAPinv, $value4AAPinv, $value5AAPinv) = split(/ /, $distribution_AApinv_all); # [fix, unif, norm, exp, gamma, beta]
			$distribution_AApinv =~ s/\n//g;
			$value1AAPinv =~ s/\n//g;
			$value2AAPinv =~ s/\n//g;
			$value3AAPinv =~ s/\n//g;
			$value4AAPinv =~ s/\n//g;
			$value5AAPinv =~ s/\n//g;
				
			if ($distribution_AApinv eq "fix")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -i
				print "  Proportion of invariable sites (+I) according to: $distribution_AApinv $value1AAPinv\n";
				$IsThere_AAPinv++;
				}
			elsif ($distribution_AApinv eq "uniform")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -i
				print "  Proportion of invariable sites (+I) according to: $distribution_AApinv $value1AAPinv $value2AAPinv\n";
				$IsThere_AAPinv++;
				}
			elsif ($distribution_AApinv eq "normal" || $distribution_AApinv eq "gamma" || $distribution_AApinv eq "beta")
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -i
				if ($value3AAPinv eq "t")
					{
					print "  Proportion of invariable sites (+I) according to: $distribution_AApinv $value1AAPinv $value2AAPinv truncated ($value3AAPinv) at $value4AAPinv and $value5AAPinv\n";
					if ($value4AAPinv > $value5AAPinv)
						{
						print "\n     Error in proportion of invariable sites (+I). Uncorrectly truncated (lowest > highest)!\n\n";
						print "Type CTRL+C to abort the execution \n";
						my $HereError = <STDIN>;
						chop($HereError);
						exit;	
						}

					}
				else
					{
					print "  Proportion of invariable sites (+I) according to: $distribution_AApinv $value1AAPinv $value2AAPinv\n";
					}
				$IsThere_AAPinv++;
				}
			elsif ($distribution_AApinv eq "exponential") 
				{
				$NumberOfArguments = $NumberOfArguments + 1; # add arguments -i
				if ($value2AAPinv eq "t")
					{
					print "  Proportion of invariable sites (+I) according to: $distribution_AApinv $value1AAPinv truncated ($value2AAPinv) at $value3AAPinv and $value4AAPinv\n";
					if ($value3AAPinv > $value4AAPinv)
						{
						print "\n     Error in proportion of invariable sites (+I). Uncorrectly truncated (lowest > highest)!\n\n";
						print "Type CTRL+C to abort the execution \n";
						my $HereError = <STDIN>;
						chop($HereError);
						exit;
						}
					}
				else
					{
					print "  Proportion of invariable sites (+I) according to: $distribution_AApinv $value1AAPinv\n";
					}
				$IsThere_AAPinv++;	
				}
			else
				{
				#print "\nWarning!. Proportion of invariable sites (+I) was not specified ($distribution_AApinv) (fix #; unif # #; norm # # t # #; exp # t # #; gamma # # t # #; beta # # t # #)\n";
				}	
			}	
						
		}

	

		
		
	} # end #
	


	} # end of lines FROM

close (FROM);






##############################################
##############################################
# Reading Settings from data file
print "\n> Reading input protein sequences file ... \n";
print " Reading input file \"$setsdir/$NameOfPhylipFile\"... \n";
unless (open(FROM_targetAln,"$setsdir/$NameOfPhylipFile"))
	{
	print "\nERROR!. Cannot open the target alignment file: \"$setsdir/$NameOfPhylipFile\"\n\n";
	print "Type CTRL+C to abort the execution \n";
	my $error = <STDIN>;
	chop($error);
	exit;
	}
my $TotalNumberLines_target = 0;
while (<FROM_targetAln>) # for each line
	{
	$TotalNumberLines_target++;	
	if ($TotalNumberLines_target == 1)
		{
		my $num3 = -1;
		my $num4 = -1;	
		my @seq1 = split(/\s+/, $_);
		#($num3, $num4) = split(/ /, $_);	
		$num3 = @seq1[0];
		$num4 = @seq1[1];
		$num4 =~ s/\n//g; # remove the end of line
						
		$samplesize = $num3;
		$nuc_number = $num4;
		
		print " $NameOfPhylipFile with sample size $samplesize and length $nuc_number amino acids";
		
		$aa_number = $nuc_number;
						
		}	
	}
$MaxRho = $MaxRho * $nuc_number;
close (FROM_targetAln);		



##############################################
##############################################
# Compute SS from input data
print "\n\n> Computing summary statistics from input protein sequences file .. \n";
if ($indels == 1)
    {
    system ("perl \"$maindir/scripts/ComputeProtSSfromRealData.pl\" \"$file\" \"$show_info\"");
    }
else
    {
    system ("perl \"$maindir/scripts/ComputeProtSSfromRealData_NoIndels.pl\" \"$file\" \"$show_info\"");
    }



##############################################
##############################################
# Writing R script

$counterLines = 0;
print "\n> Writing R output file \"$setsdir/MakePriorsProteinEvolverABC.r\" ... ";

# Checking and writing values
if ($number_simulations < 1)
	{
	print "\n\nERROR!!!. The number of replicates must be higher than 0 ($number_simulations)\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
	}
if ($samplesize < 1)
	{
	print "\n\nERROR!!!. The sample size is too small ($samplesize)\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
	}
if ($nuc_number < 1)
	{
	print "\n\nERROR!!!. Sequence length is too small ($nuc_number)\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
	}
if ($hap_dip != 1 && $hap_dip != 2)
	{
	print "\n\nERROR!!!. Haploid/diploid value is incorrect ($hap_dip)\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
	}
if ($distribution_popsize ne "uniform" && $distribution_popsize ne "fix")
	{
	print "\n\nERROR!!!. Population size is incorrect ($distribution_popsize $value1popsize $value2popsize)\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
	} 	
if ($distribution_homogRec ne "fix" && $distribution_homogRec ne "uniform" && $distribution_homogRec ne "gamma" && $distribution_homogRec ne "beta" && $distribution_homogRec ne "normal" && $distribution_homogRec ne "exponential")	
	{
	print "\n\nERROR!!!. Recombination rate is incorrect (distribution $distribution_homogRec?)\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
	} 
if ($distribution_subsRate ne "fix" && $distribution_subsRate ne "uniform" && $distribution_subsRate ne "gamma" && $distribution_subsRate ne "beta" && $distribution_subsRate ne "normal" && $distribution_subsRate ne "exponential")	
	{
	print "\n\nERROR!!!. Substitution rate is incorrect (distribution $distribution_subsRate?)\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
	}
if ($MaxRho > 200)
	{
	print "\n\nWarning!. Rho (population recombination rate*) can be too high with these settings: $MaxRho, it should not be higher than 150 or the simulator might fail due to the simulation of an almost infinite ancestral recombination graph (ARG)\n";
	print "Rho can be computed as: Rho=4*N*r*L, where N is the population size, r is the recombination rate per site (amino acid), L is the sequence length in amino acids\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	}
#print "\n\nMax Rho (rec without \"norm, gamma, beta, exp\" not truncated): $MaxRho\n";
	

# Other arguments
$NumberOfArguments = $NumberOfArguments + 1; # add arguments -n1
$NumberOfArguments = $NumberOfArguments + 2; # add arguments -s #
$NumberOfArguments = $NumberOfArguments + 21; # add arguments -f20 #
$NumberOfArguments = $NumberOfArguments + 1; # add arguments -bsequences
$NumberOfArguments = $NumberOfArguments + 3; # add arguments -c1 1 0
$NumberOfArguments = $NumberOfArguments + 1; # add arguments -y1 or -y0
$NumberOfArguments = $NumberOfArguments + 1; # add arguments -:



# Writing in the R output file
print FILE_1 "NumberReplicates <- $number_simulations\n\n";
print FILE_1 "recombination_rate <- numeric(NumberReplicates)\n";
print FILE_1 "substitution_rate <- numeric(NumberReplicates)\n";
print FILE_1 "RHO <- numeric(NumberReplicates)\n";
print FILE_1 "THETA <- numeric(NumberReplicates)\n";



print FILE_1 "\n\n";
print FILE_1 "##### PARAMETERS FOR THE ENTIRE PROTEIN are directly specified in the print file #####\n\n";

print FILE_1 "### COALESCENT SIMULATIONS ###\n";
print FILE_1 "# ./ProteinEvolverProtABC1.2.0 -n1 -s8 765 -e1000 2 -r2.0e-06 -u5.6e-05 -\@JTT -c1 1 0 -y1 -:1\n";

print FILE_1 "# ./ProteinEvolverProtABC1.2.0 -n2 -s8 150 -e1000 2 -=4 1995 1 1 2003 4 6 1997 2 3 2001 7 8 -/1200 -g1 3 1000 1250 1000 1300 1550 2000 1560 1000 3000 -q1 4 2 2 3 1 -t3 100 800 0.002 0.001 0.003 -\%1 1 2 10000 -r2.3e-6 -o0.1 -u4.1e-5 -f20 0.04 0.06 0.05 0.05 0.08 0.02 0.05 0.05 0.03 0.07 0.04 0.06 0.05 0.05 0.05 0.05 0.05 0.05 0.04 0.06 -\@WAG -a0.7 -i0.52 -bsequences -c1 1 0 -y1 -:1\n";
print FILE_1 "\n";


print FILE_1 "##### PARAMETERS FOR EACH SIMULATION #####\n";
print FILE_1 "\nThisReplicate<-0\n";
print FILE_1 "while (ThisReplicate < NumberReplicates)\n";
print FILE_1 "	{\n";
print FILE_1 "	ThisReplicate<-ThisReplicate+1\n";
print FILE_1 "\n";

print FILE_1 "	SampleSize_SequenceLength_print <- paste (\" -s$samplesize $aa_number\",sep=\"\")\n";

if ($datedtips eq "default")
	{
	print FILE_1 "	DatedTips_print <- paste (\"\",sep=\"\")\n";	
	}
else
	{
	print FILE_1 "	DatedTips_print <- paste (\" -=$datedtips\",sep=\"\")\n";	
	}
	
if ($demogperiods eq "default")
	{
	print FILE_1 "	DemogPeriods_print <- paste (\"\",sep=\"\")\n";	
	}
else
	{
	print FILE_1 "	DemogPeriods_print <- paste (\" -g$demogperiods\",sep=\"\")\n";	
	#print "\ndemogperiods = $demogperiods\n";	
	}

if ($demogperiodsDemesNcte eq "default")
	{
	print FILE_1 "	DemogPeriodsNcte_print <- paste (\"\",sep=\"\")\n";	
	}
else
	{
	print FILE_1 "	DemogPeriodsNcte_print <- paste (\" -g$demogperiodsDemesNcte\",sep=\"\")\n";
	if ($convdemesYes == 0)
		{
		print "\n\nERROR!. Populations/species tree is required for this demographic option!\n\n";
		print "Type CTRL+C to abort the execution \n";
		my $HereError = <STDIN>;
		chop($HereError);
		exit;
		}	
	}
	
if ($demogperiodsDemesNvar eq "default")
	{
	print FILE_1 "	DemogPeriodsNvar_print <- paste (\"\",sep=\"\")\n";	
	}
else
	{
	print FILE_1 "	DemogPeriodsNvar_print <- paste (\" -g$demogperiodsDemesNvar\",sep=\"\")\n";
	if ($convdemesYes == 0)
		{
		print "\n\nERROR!. Populations/species tree is required for this demographic option!\n\n";
		print "Type CTRL+C to abort the execution \n";
		my $HereError = <STDIN>;
		chop($HereError);
		exit;
		}	
	}

if ($Migrationmodel eq "default")
	{
	print FILE_1 "	Migrationmodel_print <- paste (\"\",sep=\"\")\n";	
	}
else
	{
	print FILE_1 "	Migrationmodel_print <- paste (\" -q$Migrationmodel\",sep=\"\")\n";
	if ($migrationrateYes == 0)
		{
		print "\n\nERROR!. A migration model requires a migration rate!\n\n";
		print "Type CTRL+C to abort the execution \n";
		my $HereError = <STDIN>;
		chop($HereError);
		exit;
		}		
	}
	
if ($migrationrate eq "default")
	{
	print FILE_1 "	MigrationRate_print <- paste (\"\",sep=\"\")\n";	
	}
else
	{
	print FILE_1 "	MigrationRate_print <- paste (\" -t$migrationrate\",sep=\"\")\n";	
	}
	
if ($convdemes eq "default")
	{
	print FILE_1 "	ConvDemes_print <- paste (\"\",sep=\"\")\n";	
	}
else
	{
	print FILE_1 "	ConvDemes_print <- paste (\" -%$convdemes\",sep=\"\")\n";	
	}
	

print FILE_1 "	Hap_Dip <- paste (\" $hap_dip\",sep=\"\")\n";	
	
if ($distribution_popsize eq "fix")
	{
	print FILE_1 "	PopSize <- $value1popsize		# unif or fix\n";
	}
else
	{
	print FILE_1 "	PopSize <- sample($value1popsize:$value2popsize,1,replace=T)		# unif or fix\n";
	}

print FILE_1 "	PopSize_Hap_Dip_print <- paste (\" -e\",PopSize,Hap_Dip,sep=\"\")\n";	


if ($distribution_gtime eq "fix")
	{
	print FILE_1 "	GenTimeValue <- $value1gtime		# unif or fix\n";
	print FILE_1 "	GenTime_print <- paste (\" -\/\",GenTimeValue,sep=\"\")\n";	
	}
elsif ($distribution_gtime eq "uniform")
	{
	print FILE_1 "	GenTimeValue <- sample($value1gtime:$value2gtime,1,replace=T)		# unif or fix\n";
	print FILE_1 "	GenTime_print <- paste (\" -\/\",GenTimeValue,sep=\"\")\n";	
	}
else
	{
	print FILE_1 "	GenTime_print <- paste (\"\",sep=\"\")\n";	
	}

		

	print FILE_1 "	GrowthRate_print <- paste (\"\",sep=\"\")\n";
	if ($distribution_GrowthRate eq "fix")
		{
		print FILE_1 "	GrowthRate <- $value1GrowthRate		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
		print FILE_1 "	GrowthRate_print <- paste (\" -g0 \",GrowthRate,sep=\"\")\n";
		}
	elsif ($distribution_GrowthRate eq "uniform")
		{
		print FILE_1 "	GrowthRate <- runif(1,$value1GrowthRate,$value2GrowthRate)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
		print FILE_1 "	GrowthRate_print <- paste (\" -g0 \",GrowthRate,sep=\"\")\n";
		}
	elsif ($distribution_GrowthRate eq "normal")
		{
		if ($value3GrowthRate eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4GrowthRate \n";
			print FILE_1 "	MaxValue <- $value5GrowthRate \n";
			print FILE_1 "	GrowthRate <- 'inf' \n";
			print FILE_1 "	while (GrowthRate > MaxValue || GrowthRate < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		GrowthRate <- rnorm(1,$value1GrowthRate,$value2GrowthRate)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";		
			print FILE_1 "	GrowthRate_print <- paste (\" -g0 \",GrowthRate,sep=\"\")\n";			
			}
		else
			{
			print FILE_1 "	GrowthRate <- -1\n";
			print FILE_1 "	while (GrowthRate < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		GrowthRate <- rnorm(1,$value1GrowthRate,$value2GrowthRate)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	GrowthRate_print <- paste (\" -g0 \",GrowthRate,sep=\"\")\n";		
			}
		}
	elsif ($distribution_GrowthRate eq "exponential")
		{
		if ($value2GrowthRate eq "t") # t
			{
			print FILE_1 "	MinValue <- $value3GrowthRate \n";
			print FILE_1 "	MaxValue <- $value4GrowthRate \n";
			print FILE_1 "	GrowthRate <- 'inf' \n";
			print FILE_1 "	while (GrowthRate > MaxValue || GrowthRate < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		GrowthRate <- rexp(1,$value1GrowthRate)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			print FILE_1 "	GrowthRate_print <- paste (\" -g0 \",GrowthRate,sep=\"\")\n";
			}
		else
			{
			print FILE_1 "	GrowthRate <- -1\n";
			print FILE_1 "	while (GrowthRate < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		GrowthRate <- rexp(1,$value1GrowthRate)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	GrowthRate_print <- paste (\" -g0 \",GrowthRate,sep=\"\")\n";		
			}
		}
	elsif ($distribution_GrowthRate eq "gamma")
		{
		if ($value3GrowthRate eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4GrowthRate \n";
			print FILE_1 "	MaxValue <- $value5GrowthRate \n";
			print FILE_1 "	GrowthRate <- 'inf' \n";
			print FILE_1 "	while (GrowthRate > MaxValue || GrowthRate < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		GrowthRate <- rgamma(1,$value1GrowthRate,$value2GrowthRate)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	GrowthRate_print <- paste (\" -g0 \",GrowthRate,sep=\"\")\n";				
			}
		else
			{
			print FILE_1 "	GrowthRate <- -1\n";
			print FILE_1 "	while (GrowthRate < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		GrowthRate <- rgamma(1,$value1GrowthRate,$value2GrowthRate)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			print FILE_1 "	GrowthRate_print <- paste (\" -g0 \",GrowthRate,sep=\"\")\n";
			}
		}
	elsif ($distribution_GrowthRate eq "beta")
		{
		if ($value3GrowthRate eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4GrowthRate \n";
			print FILE_1 "	MaxValue <- $value5GrowthRate \n";
			print FILE_1 "	GrowthRate <- 'inf' \n";
			print FILE_1 "	while (GrowthRate > MaxValue || GrowthRate < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		GrowthRate <- rbeta(1,$value1GrowthRate,$value2GrowthRate)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	GrowthRate_print <- paste (\" -g0 \",GrowthRate,sep=\"\")\n";				
			}
		else
			{
			print FILE_1 "	GrowthRate <- -1\n";
			print FILE_1 "	while (GrowthRate < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		GrowthRate <- rbeta(1,$value1GrowthRate,$value2GrowthRate)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			print FILE_1 "	GrowthRate_print <- paste (\" -g0 \",GrowthRate,sep=\"\")\n";
			}
		}	


	print FILE_1 "	outgroup_print <- paste (\"\",sep=\"\")\n";
	if ($distribution_outgroup eq "fix")
		{
		print FILE_1 "	outgroup <- $value1Houtgroup		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "	outgroup_print <- paste (\" -o\",outgroup,sep=\"\")\n";
		}
	elsif ($distribution_outgroup eq "uniform")
		{
		print FILE_1 "	outgroup <- runif(1,$value1Houtgroup,$value2Houtgroup)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "	outgroup_print <- paste (\" -o\",outgroup,sep=\"\")\n";
		}
	elsif ($distribution_outgroup eq "normal")
		{
		if ($value3Houtgroup eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4Houtgroup \n";
			print FILE_1 "	MaxValue <- $value5Houtgroup \n";
			print FILE_1 "	outgroup <- 'inf' \n";
			print FILE_1 "	while (outgroup > MaxValue || outgroup < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		outgroup <- rnorm(1,$value1Houtgroup,$value2Houtgroup)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	outgroup_print <- paste (\" -o\",outgroup,sep=\"\")\n";
			}
		else
			{
			print FILE_1 "	outgroup <- -1\n";
			print FILE_1 "	while (outgroup < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		outgroup <- rnorm(1,$value1Houtgroup,$value2Houtgroup)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	outgroup_print <- paste (\" -o\",outgroup,sep=\"\")\n";
			}
		}
	elsif ($distribution_outgroup eq "exponential")
		{
		if ($value2Houtgroup eq "t") # t
			{
			print FILE_1 "	MinValue <- $value3Houtgroup \n";
			print FILE_1 "	MaxValue <- $value4Houtgroup \n";
			print FILE_1 "	outgroup <- 'inf' \n";
			print FILE_1 "	while (outgroup > MaxValue || outgroup < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		outgroup <- rexp(1,$value1Houtgroup)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			print FILE_1 "	outgroup_print <- paste (\" -o\",outgroup,sep=\"\")\n";
			}
		else
			{
			print FILE_1 "	outgroup <- -1\n";
			print FILE_1 "	while (outgroup < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		outgroup <- rexp(1,$value1Houtgroup)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	outgroup_print <- paste (\" -o\",outgroup,sep=\"\")\n";
			}
		}
	elsif ($distribution_outgroup eq "gamma")
		{
		if ($value3Houtgroup eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4Houtgroup \n";
			print FILE_1 "	MaxValue <- $value5Houtgroup \n";
			print FILE_1 "	outgroup <- 'inf' \n";
			print FILE_1 "	while (outgroup > MaxValue || outgroup < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		outgroup <- rgamma(1,$value1Houtgroup,$value2Houtgroup)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	outgroup_print <- paste (\" -o\",outgroup,sep=\"\")\n";
			}
		else
			{
			print FILE_1 "	outgroup <- -1\n";
			print FILE_1 "	while (outgroup < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		outgroup <- rgamma(1,$value1Houtgroup,$value2Houtgroup)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			print FILE_1 "	outgroup_print <- paste (\" -o\",outgroup,sep=\"\")\n";
			}
		}
	elsif ($distribution_outgroup eq "beta")
		{
		if ($value3Houtgroup eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4Houtgroup \n";
			print FILE_1 "	MaxValue <- $value5Houtgroup \n";
			print FILE_1 "	outgroup <- 'inf' \n";
			print FILE_1 "	while (outgroup > MaxValue || outgroup < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		outgroup <- rbeta(1,$value1Houtgroup,$value2Houtgroup)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";		
			print FILE_1 "	outgroup_print <- paste (\" -o\",outgroup,sep=\"\")\n";
			}
		else
			{
			print FILE_1 "	outgroup <- -1\n";
			print FILE_1 "	while (outgroup < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		outgroup <- rbeta(1,$value1Houtgroup,$value2Houtgroup)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			print FILE_1 "	outgroup_print <- paste (\" -o\",outgroup,sep=\"\")\n";
			}
		}	



	if ($distribution_homogRec eq "default")
		{
		print "\n\nERROR!!!. Recombination rate is incorrect (distribution $distribution_homogRec?)\n";
		print "Type CTRL+C to abort the execution \n";
		my $HereError = <STDIN>;
		chop($HereError);
		exit;
		}
	elsif ($distribution_homogRec eq "fix")
		{
		print FILE_1 "	HomogRec <- $value1Hrec		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
		}
	elsif ($distribution_homogRec eq "uniform")
		{
		print FILE_1 "	HomogRec <- runif(1,$value1Hrec,$value2Hrec)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
		}
	elsif ($distribution_homogRec eq "normal")
		{
		if ($value3Hrec eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4Hrec \n";
			print FILE_1 "	MaxValue <- $value5Hrec \n";
			print FILE_1 "	HomogRec <- 'inf' \n";
			print FILE_1 "	while (HomogRec > MaxValue || HomogRec < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		HomogRec <- rnorm(1,$value1Hrec,$value2Hrec)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";					
			}
		else
			{
			print FILE_1 "	HomogRec <- -1\n";
			print FILE_1 "	while (HomogRec < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		HomogRec <- rnorm(1,$value1Hrec,$value2Hrec)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";			
			}
		}
	elsif ($distribution_homogRec eq "exponential")
		{
		if ($value2Hrec eq "t") # t
			{
			print FILE_1 "	MinValue <- $value3Hrec \n";
			print FILE_1 "	MaxValue <- $value4Hrec \n";
			print FILE_1 "	HomogRec <- 'inf' \n";
			print FILE_1 "	while (HomogRec > MaxValue || HomogRec < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		HomogRec <- rexp(1,$value1Hrec)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			}
		else
			{
			print FILE_1 "	HomogRec <- -1\n";
			print FILE_1 "	while (HomogRec < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		HomogRec <- rexp(1,$value1Hrec)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";			
			}
		}
	elsif ($distribution_homogRec eq "gamma")
		{
		if ($value3Hrec eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4Hrec \n";
			print FILE_1 "	MaxValue <- $value5Hrec \n";
			print FILE_1 "	HomogRec <- 'inf' \n";
			print FILE_1 "	while (HomogRec > MaxValue || HomogRec < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		HomogRec <- rgamma(1,$value1Hrec,$value2Hrec)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";					
			}
		else
			{
			print FILE_1 "	HomogRec <- -1\n";
			print FILE_1 "	while (HomogRec < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		HomogRec <- rgamma(1,$value1Hrec,$value2Hrec)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			}
		}
	elsif ($distribution_homogRec eq "beta")
		{
		if ($value3Hrec eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4Hrec \n";
			print FILE_1 "	MaxValue <- $value5Hrec \n";
			print FILE_1 "	HomogRec <- 'inf' \n";
			print FILE_1 "	while (HomogRec > MaxValue || HomogRec < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		HomogRec <- rbeta(1,$value1Hrec,$value2Hrec)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";					
			}
		else
			{
			print FILE_1 "	HomogRec <- -1\n";
			print FILE_1 "	while (HomogRec < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		HomogRec <- rbeta(1,$value1Hrec,$value2Hrec)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			}
		}
	else
		{
		print "\n\nERROR!!!. Recombination rate is incorrect (distribution $distribution_homogRec?)\n";
		print "Type CTRL+C to abort the execution \n";
		my $HereError = <STDIN>;
		chop($HereError);
		exit;	
		}	
	print FILE_1 "	HomogRec_print <- paste (\" -r\",HomogRec,sep=\"\")\n";
	print FILE_1 "	recombination_rate[ThisReplicate]<-HomogRec\n";



	if ($distribution_subsRate eq "default")
		{
		print "\n\nERROR!!!. Substitution rate is incorrect (distribution $distribution_subsRate?)\n";
		print "Type CTRL+C to abort the execution \n";
		my $HereError = <STDIN>;
		chop($HereError);
		exit;
		}
	elsif ($distribution_subsRate eq "fix")
		{
		print FILE_1 "	SubsRate <- $value1subsR		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
		}
	elsif ($distribution_subsRate eq "uniform")
		{
		print FILE_1 "	SubsRate <- runif(1,$value1subsR,$value2subsR)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
		}
	elsif ($distribution_subsRate eq "normal")
		{
		if ($value3subsR eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4subsR \n";
			print FILE_1 "	MaxValue <- $value5subsR \n";
			print FILE_1 "	SubsRate <- 'inf' \n";
			print FILE_1 "	while (SubsRate > MaxValue || SubsRate < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		SubsRate <- rnorm(1,$value1subsR,$value2subsR)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";					
			}
		else
			{
			print FILE_1 "	SubsRate <- -1\n";
			print FILE_1 "	while (SubsRate < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		SubsRate <- rnorm(1,$value1subsR,$value2subsR)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";			
			}
		}
	elsif ($distribution_subsRate eq "exponential")
		{
		if ($value2subsR eq "t") # t
			{
			print FILE_1 "	MinValue <- $value3subsR \n";
			print FILE_1 "	MaxValue <- $value4subsR \n";
			print FILE_1 "	SubsRate <- 'inf' \n";
			print FILE_1 "	while (SubsRate > MaxValue || SubsRate < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		SubsRate <- rexp(1,$value1subsR)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			}
		else
			{
			print FILE_1 "	SubsRate <- -1\n";
			print FILE_1 "	while (SubsRate < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		SubsRate <- rexp(1,$value1subsR)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";			
			}
		}
	elsif ($distribution_subsRate eq "gamma")
		{
		if ($value3subsR eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4subsR \n";
			print FILE_1 "	MaxValue <- $value5subsR \n";
			print FILE_1 "	SubsRate <- 'inf' \n";
			print FILE_1 "	while (SubsRate > MaxValue || SubsRate < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		SubsRate <- rgamma(1,$value1subsR,$value2subsR)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";					
			}
		else
			{
			print FILE_1 "	SubsRate <- -1\n";
			print FILE_1 "	while (SubsRate < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		SubsRate <- rgamma(1,$value1subsR,$value2subsR)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			}
		}
	elsif ($distribution_subsRate eq "beta")
		{
		if ($value3subsR eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4subsR \n";
			print FILE_1 "	MaxValue <- $value5subsR \n";
			print FILE_1 "	SubsRate <- 'inf' \n";
			print FILE_1 "	while (SubsRate > MaxValue || SubsRate < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		SubsRate <- rbeta(1,$value1subsR,$value2subsR)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";					
			}
		else
			{
			print FILE_1 "	SubsRate <- -1\n";
			print FILE_1 "	while (SubsRate < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		SubsRate <- rbeta(1,$value1subsR,$value2subsR)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			}
		}
	else
		{
		print "\n\nERROR!!!. Substitution rate is incorrect (distribution $distribution_subsRate?)\n";
		print "Type CTRL+C to abort the execution \n";
		my $HereError = <STDIN>;
		chop($HereError);
		exit;	
		}	
	print FILE_1 "	SubsRate_print <- paste (\" -u\",SubsRate,sep=\"\")\n";
	print FILE_1 "	substitution_rate[ThisReplicate]<-SubsRate\n";

	

	# AA frequencies
	print FILE_1 "\n"; # by default
	if ($distribution_FreqsAA1 eq "fix")
		{
		print FILE_1 "	AAFreqs1 <- c($value1FreqsAA,$value2FreqsAA,$value3FreqsAA,$value4FreqsAA,$value5FreqsAA,$value6FreqsAA,$value7FreqsAA,$value8FreqsAA,$value9FreqsAA,$value10FreqsAA,$value11FreqsAA,$value12FreqsAA,$value13FreqsAA,$value14FreqsAA,$value15FreqsAA,$value16FreqsAA,$value17FreqsAA,$value18FreqsAA,$value19FreqsAA,$value20FreqsAA)		\n";
		}
	elsif ($distribution_FreqsAA1 eq "dirichlet")
		{
		print FILE_1 "	AAFreqs1 <- rdirichlet(1, c($value1FreqsAA,$value2FreqsAA,$value3FreqsAA,$value4FreqsAA,$value5FreqsAA,$value6FreqsAA,$value7FreqsAA,$value8FreqsAA,$value9FreqsAA,$value10FreqsAA,$value11FreqsAA,$value12FreqsAA,$value13FreqsAA,$value14FreqsAA,$value15FreqsAA,$value16FreqsAA,$value17FreqsAA,$value18FreqsAA,$value19FreqsAA,$value20FreqsAA))		\n";
		}
	else
		{
		print FILE_1 "	AAFreqs1 <- numeric(20)		\n";
		print FILE_1 "	AAFreqs1[1] <- 0.05 \n"; # by default
		print FILE_1 "	AAFreqs1[2] <- 0.05 \n"; # by default
		print FILE_1 "	AAFreqs1[3] <- 0.05 \n"; # by default
		print FILE_1 "	AAFreqs1[4] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[5] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[6] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[7] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[8] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[9] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[10] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[11] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[12] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[13] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[14] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[15] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[16] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[17] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[18] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[19] <- 0.05 \n"; # by default
        print FILE_1 "	AAFreqs1[20] <- 0.05 \n"; # by default
		}
    print FILE_1 "\n"; # by default
	print FILE_1 "	AAFreqs <- numeric(20)		\n";
	print FILE_1 "	AAFreqs[1] <- AAFreqs1[1]		\n";
	print FILE_1 "	AAFreqs[2] <- AAFreqs1[2]		\n";
	print FILE_1 "	AAFreqs[3] <- AAFreqs1[3]		\n";
	print FILE_1 "	AAFreqs[4] <- AAFreqs1[4]		\n";
	print FILE_1 "	AAFreqs[5] <- AAFreqs1[5]		\n";
	print FILE_1 "	AAFreqs[6] <- AAFreqs1[6]		\n";
	print FILE_1 "	AAFreqs[7] <- AAFreqs1[7]		\n";
	print FILE_1 "	AAFreqs[8] <- AAFreqs1[8]		\n";
	print FILE_1 "	AAFreqs[9] <- AAFreqs1[9]		\n";
	print FILE_1 "	AAFreqs[10] <- AAFreqs1[10]		\n";
	print FILE_1 "	AAFreqs[11] <- AAFreqs1[11]		\n";
	print FILE_1 "	AAFreqs[12] <- AAFreqs1[12]		\n";
    print FILE_1 "	AAFreqs[13] <- AAFreqs1[13]		\n";
    print FILE_1 "	AAFreqs[14] <- AAFreqs1[14]		\n";
    print FILE_1 "	AAFreqs[15] <- AAFreqs1[15]		\n";
    print FILE_1 "	AAFreqs[16] <- AAFreqs1[16]		\n";
    print FILE_1 "	AAFreqs[17] <- AAFreqs1[17]		\n";
    print FILE_1 "	AAFreqs[18] <- AAFreqs1[18]		\n";
    print FILE_1 "	AAFreqs[19] <- AAFreqs1[19]		\n";
    print FILE_1 "	AAFreqs[20] <- AAFreqs1[20]		\n";

	print FILE_1 "	AminoacidFrequencies_print <- paste (\" -f20 \",AAFreqs[1],\" \",AAFreqs[2],\" \",AAFreqs[3],\" \",AAFreqs[4],\" \",AAFreqs[5],\" \",AAFreqs[6],\" \",AAFreqs[7],\" \",AAFreqs[8],\" \",AAFreqs[9],\" \",AAFreqs[10],\" \",AAFreqs[11],\" \",AAFreqs[12],\" \",AAFreqs[13],\" \",AAFreqs[14],\" \",AAFreqs[15],\" \",AAFreqs[16],\" \",AAFreqs[17],\" \",AAFreqs[18],\" \",AAFreqs[19],\" \",AAFreqs[20],sep=\"\")\n";

    print FILE_1 "\n	SubsModel_print <- paste (\" -@\",\"$Subs_dip\",sep=\"\")\n";




	# +G
	print FILE_1 "	AArateHetSites_print <- paste (\"\",sep=\"\")\n";
	if ($distribution_AArateHetSites eq "fix")
		{
		print FILE_1 "	AArateHetSites <- $value1AArateHetSites		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
		print FILE_1 "	AArateHetSites_print <- paste (\" -a\",AArateHetSites,sep=\"\")\n";
		}
	elsif ($distribution_AArateHetSites eq "uniform")
		{
		print FILE_1 "	AArateHetSites <- runif(1,$value1AArateHetSites,$value2AArateHetSites)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
		print FILE_1 "	AArateHetSites_print <- paste (\" -a\",AArateHetSites,sep=\"\")\n";
		}
	elsif ($distribution_AArateHetSites eq "normal")
		{
		if ($value3AArateHetSites eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4AArateHetSites \n";
			print FILE_1 "	MaxValue <- $value5AArateHetSites \n";
			print FILE_1 "	AArateHetSites <- 'inf' \n";
			print FILE_1 "	while (AArateHetSites > MaxValue || AArateHetSites < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AArateHetSites <- rnorm(1,$value1AArateHetSites,$value2AArateHetSites)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	AArateHetSites_print <- paste (\" -a\",AArateHetSites,sep=\"\")\n";
			}
		else
			{
			print FILE_1 "  AArateHetSites <- -1\n";
			print FILE_1 "	while (AArateHetSites < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AArateHetSites <- rnorm(1,$value1AArateHetSites,$value2AArateHetSites)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	AArateHetSites_print <- paste (\" -a\",AArateHetSites,sep=\"\")\n";
			}
		}
	elsif ($distribution_AArateHetSites eq "exponential")
		{
		if ($value2AArateHetSites eq "t") # t
			{
			print FILE_1 "	MinValue <- $value3AArateHetSites \n";
			print FILE_1 "	MaxValue <- $value4AArateHetSites \n";
			print FILE_1 "	AArateHetSites <- 'inf' \n";
			print FILE_1 "	while (AArateHetSites > MaxValue || AArateHetSites < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AArateHetSites <- rexp(1,$value1AArateHetSites)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			print FILE_1 "	AArateHetSites_print <- paste (\" -a\",AArateHetSites,sep=\"\")\n";
			}
		else
			{
			print FILE_1 "  AArateHetSites <- -1\n";
			print FILE_1 "	while (AArateHetSites < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AArateHetSites <- rexp(1,$value1AArateHetSites)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	AArateHetSites_print <- paste (\" -a\",AArateHetSites,sep=\"\")\n";
			}
		}
	elsif ($distribution_AArateHetSites eq "gamma")
		{
		if ($value3AArateHetSites eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4AArateHetSites \n";
			print FILE_1 "	MaxValue <- $value5AArateHetSites \n";
			print FILE_1 "	AArateHetSites <- 'inf' \n";
			print FILE_1 "	while (AArateHetSites > MaxValue || AArateHetSites < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AArateHetSites <- rgamma(1,$value1AArateHetSites,$value2AArateHetSites)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	AArateHetSites_print <- paste (\" -a\",AArateHetSites,sep=\"\")\n";
			}
		else
			{
			print FILE_1 "  AArateHetSites <- -1\n";
			print FILE_1 "	while (AArateHetSites < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AArateHetSites <- rgamma(1,$value1AArateHetSites,$value2AArateHetSites)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			print FILE_1 "	AArateHetSites_print <- paste (\" -a\",AArateHetSites,sep=\"\")\n";
			}
		}
	elsif ($distribution_AArateHetSites eq "beta")
		{
		if ($value3AArateHetSites eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4AArateHetSites \n";
			print FILE_1 "	MaxValue <- $value5AArateHetSites \n";
			print FILE_1 "	AArateHetSites <- 'inf' \n";
			print FILE_1 "	while (AArateHetSites > MaxValue || AArateHetSites < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AArateHetSites <- rbeta(1,$value1AArateHetSites,$value2AArateHetSites)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	AArateHetSites_print <- paste (\" -a\",AArateHetSites,sep=\"\")\n";
			}
		else
			{
			print FILE_1 "  AArateHetSites <- -1\n";
			print FILE_1 "	while (AArateHetSites < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AArateHetSites <- rbeta(1,$value1AArateHetSites,$value2AArateHetSites)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			print FILE_1 "	AArateHetSites_print <- paste (\" -a\",AArateHetSites,sep=\"\")\n";
			}
		}
	else
		{
		;	
		}


	# +I
	print FILE_1 "	AAPinv_print <- paste (\"\",sep=\"\")\n";
	if ($distribution_AApinv eq "fix")
		{
		print FILE_1 "	AAPinv <- $value1AAPinv		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
		print FILE_1 "	AAPinv_print <- paste (\" -i\",AAPinv,sep=\"\")\n";
		}
	elsif ($distribution_AApinv eq "uniform")
		{
		print FILE_1 "	AAPinv <- runif(1,$value1AAPinv,$value2AAPinv)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
		print FILE_1 "	AAPinv_print <- paste (\" -i\",AAPinv,sep=\"\")\n";
		}
	elsif ($distribution_AApinv eq "normal")
		{
		if ($value3GrowthRate eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4AAPinv \n";
			print FILE_1 "	MaxValue <- $value5AAPinv \n";
			print FILE_1 "	AAPinv <- 'inf' \n";
			print FILE_1 "	while (AAPinv > MaxValue || AAPinv < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AAPinv <- rnorm(1,$value1AAPinv,$value2AAPinv)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	AAPinv_print <- paste (\" -i\",AAPinv,sep=\"\")\n";				
			}
		else
			{
			print FILE_1 "  AAPinv <- -1\n";
			print FILE_1 "	while (AAPinv < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AAPinv <- rnorm(1,$value1AAPinv,$value2AAPinv)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	AAPinv_print <- paste (\" -i\",AAPinv,sep=\"\")\n";		
			}
		}
	elsif ($distribution_AApinv eq "exponential")
		{
		if ($value2GrowthRate eq "t") # t
			{
			print FILE_1 "	MinValue <- $value3AAPinv \n";
			print FILE_1 "	MaxValue <- $value4AAPinv \n";
			print FILE_1 "	AAPinv <- 'inf' \n";
			print FILE_1 "	while (AAPinv > MaxValue || AAPinv < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AAPinv <- rexp(1,$value1AAPinv)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			print FILE_1 "	AAPinv_print <- paste (\" -i\",AAPinv,sep=\"\")\n";
			}
		else
			{
			print FILE_1 "  AAPinv <- -1\n";
			print FILE_1 "	while (AAPinv < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AAPinv <- rexp(1,$value1AAPinv)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";	
			print FILE_1 "	AAPinv_print <- paste (\" -i\",AAPinv,sep=\"\")\n";		
			}
		}
	elsif ($distribution_AApinv eq "gamma")
		{
		if ($value3GrowthRate eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4AAPinv \n";
			print FILE_1 "	MaxValue <- $value5AAPinv \n";
			print FILE_1 "	AAPinv <- 'inf' \n";
			print FILE_1 "	while (AAPinv > MaxValue || AAPinv < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AAPinv <- rgamma(1,$value1AAPinv,$value2AAPinv)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";		
			print FILE_1 "	AAPinv_print <- paste (\" -i\",AAPinv,sep=\"\")\n";			
			}
		else
			{
			print FILE_1 "  AAPinv <- -1\n";
			print FILE_1 "	while (AAPinv < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AAPinv <- rgamma(1,$value1AAPinv,$value2AAPinv)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			print FILE_1 "	AAPinv_print <- paste (\" -i\",AAPinv,sep=\"\")\n";
			}
		}
	elsif ($distribution_AApinv eq "beta")
		{
		if ($value3AAPinv eq "t") # t
			{
			print FILE_1 "	MinValue <- $value4AAPinv \n";
			print FILE_1 "	MaxValue <- $value5AAPinv \n";
			print FILE_1 "	AAPinv <- 'inf' \n";
			print FILE_1 "	while (AAPinv > MaxValue || AAPinv < MinValue) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AAPinv <- rbeta(1,$value1AAPinv,$value2AAPinv)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";		
			print FILE_1 "	AAPinv_print <- paste (\" -i\",AAPinv,sep=\"\")\n";			
			}
		else
			{
			print FILE_1 "  AAPinv <- -1\n";
			print FILE_1 "	while (AAPinv < 0) \n";
			print FILE_1 "		{ \n";
			print FILE_1 "		AAPinv <- rbeta(1,$value1AAPinv,$value2AAPinv)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)\n";
			print FILE_1 "		} \n";
			print FILE_1 "	AAPinv_print <- paste (\" -i\",AAPinv,sep=\"\")\n";
			}
		}
	else
		{
		;	
		}


	print FILE_1 "\n";
    print FILE_1 "\n";




##############################################
##############################################
# Writing Exectuable file

$counterLines = 0;
print "\n\n> Writing executable output file for simulations \"$setsdir/ProteinEvolverABC_Phase2.sh\", \"$setsdir/ProteinEvolverABC_arguments.txt\" and \"$setsdir/ProteinEvolverABC_Pop_arguments.txt\" ... ";
#print "\n> Writing output file with the true values \"TrueSimulatedValues.txt\" ... ";


# print FILE_1 "# ./ProteinEvolverProtABC1.2.0 -n2 -s8 150 -e1000 2 -=4 1995 1 1 2003 4 6 1997 2 3 2001 7 8 -/1200 -g1 3 1000 1250 1000 1300 1550 2000 1560 1000 3000 -q1 4 2 2 3 1 -t3 100 800 0.002 0.001 0.003 -%1 1 2 10000 -r2.3e-6 -o0.1 -u4.1e-5 -f20 0.04 0.06 0.05 0.05 0.08 0.02 0.05 0.05 0.03 0.07 0.04 0.06 0.05 0.05 0.05 0.05 0.05 0.05 0.04 0.06 -@WAG -a0.7 -i0.52 -bsequences -c1 1 0 -y1 -:1\n";

######################


# print ProteinEvolverABC_arguments.txt
if ($show_info == 1)
    {
    print FILE_1 "	ExecutionHeader <- paste (\"-n1\",SampleSize_SequenceLength_print,PopSize_Hap_Dip_print,DatedTips_print,GenTime_print,GrowthRate_print,DemogPeriods_print,DemogPeriodsNcte_print,DemogPeriodsNvar_print,Migrationmodel_print,MigrationRate_print,ConvDemes_print,HomogRec_print,SubsRate_print,outgroup_print,AminoacidFrequencies_print,SubsModel_print,AArateHetSites_print,AAPinv_print,\" -bsequences -c1 1 0 -y1 -:\",ThisReplicate,sep=\"\")\n";
    }
else
    {
    print FILE_1 "	ExecutionHeader <- paste (\"-n1\",SampleSize_SequenceLength_print,PopSize_Hap_Dip_print,DatedTips_print,GenTime_print,GrowthRate_print,DemogPeriods_print,DemogPeriodsNcte_print,DemogPeriodsNvar_print,Migrationmodel_print,MigrationRate_print,ConvDemes_print,HomogRec_print,SubsRate_print,outgroup_print,AminoacidFrequencies_print,SubsModel_print,AArateHetSites_print,AAPinv_print,\" -bsequences -c1 1 0 -y0 -:\",ThisReplicate,sep=\"\")\n";
    }

print FILE_1 "\n	write(paste(\"\",ExecutionHeader,\"\",sep=\"\"),\"$setsdir/ProteinEvolverABC_arguments.txt\",append=T)";

print FILE_1 "\n	Rho_sim <- 2 * $hap_dip * PopSize * HomogRec * $aa_number \n";
print FILE_1 "	HomogRec_print_Pop <- paste (\" -r\",Rho_sim,sep=\"\")\n";
print FILE_1 "\n	Theta_sim <- 2 * $hap_dip * PopSize * SubsRate * $aa_number \n";
print FILE_1 "	SubsRate_print_pop <- paste (\" -u\",Theta_sim,sep=\"\")\n";
print FILE_1 "	RHO[ThisReplicate]<-Rho_sim\n";
print FILE_1 "	THETA[ThisReplicate]<-Theta_sim\n";

if ($show_info == 1)
    {
    print FILE_1 "	ExecutionHeader_Pop <- paste (\"-n1\",SampleSize_SequenceLength_print,PopSize_Hap_Dip_print,DatedTips_print,GenTime_print,GrowthRate_print,DemogPeriods_print,DemogPeriodsNcte_print,DemogPeriodsNvar_print,Migrationmodel_print,MigrationRate_print,ConvDemes_print,HomogRec_print_Pop,SubsRate_print_pop,outgroup_print,AminoacidFrequencies_print,SubsModel_print,AArateHetSites_print,AAPinv_print,\" -bsequences -c1 1 0 -y1 -:\",ThisReplicate,sep=\"\")\n";
    }
else
    {
    print FILE_1 "	ExecutionHeader_Pop <- paste (\"-n1\",SampleSize_SequenceLength_print,PopSize_Hap_Dip_print,DatedTips_print,GenTime_print,GrowthRate_print,DemogPeriods_print,DemogPeriodsNcte_print,DemogPeriodsNvar_print,Migrationmodel_print,MigrationRate_print,ConvDemes_print,HomogRec_print_Pop,SubsRate_print_pop,outgroup_print,AminoacidFrequencies_print,SubsModel_print,AArateHetSites_print,AAPinv_print,\" -bsequences -c1 1 0 -y0 -:\",ThisReplicate,sep=\"\")\n";
    }
print FILE_1 "\n	write(paste(\"\",ExecutionHeader_Pop,\"\",sep=\"\"),\"$setsdir/ProteinEvolverABC_Pop_arguments.txt\",append=T)";

# print TrueSimulatedValues.txt
#print FILE_1 "\n	if (ThisReplicate < 2)		{\n";
#print FILE_1 "		headhere<-paste(\"Simulation number	Recombination rate	Substitution rate	Omega (dN/dS)\",sep=\"\")\n";
#print FILE_1 "		write(paste(\"\",headhere,\"\",sep=\"\"),\"$setsdir/TrueSimulatedValues.txt\",append=T)\n";
#print FILE_1 "	}\n";
#print FILE_1 "	bodyhere<-paste(\"\",ThisReplicate,\"	\",HomogRec,\"	\",SubsRate,\"	\",GY94M0dnds,sep=\"\")\n";
#print FILE_1 "	write(paste(\"\",bodyhere,\"\",sep=\"\"),\"$setsdir/TrueSimulatedValues.txt\",append=T)\n";

print FILE_1 "\n\n} # end replicates\n"; # end of replicates ThisReplicate



##############################################
##############################################
# print ProteinEvolverABC_Phase2.sh
print FILE_1 "	NumberOfProcessors <- $numberOfProcessors\n";
print FILE_1 "	NumberOfArguments <- $NumberOfArguments\n";
#if ($OperativeSystem eq "linux" && $numberOfProcessors > 1)
if (($OperativeSystem eq "linux" || $OperativeSystem eq "darwin") && $numberOfProcessors > 1) # Allowed for Linux and Mac
	{
    #print "\n\n\n PARALLEL \n\n\n\n";
	print FILE_1 "	bodyhereNew<-paste(\"more \\\"$setsdir/ProteinEvolverABC_arguments.txt\\\" | xargs -n \",NumberOfArguments,\" -P \",NumberOfProcessors,\" \\\"$maindir/bin/ProteinEvolverProtABC1.2.0\\\"\",sep=\"\")\n";
	print FILE_1 "	write(paste(\"\",bodyhereNew,\"\",sep=\"\"),\"$setsdir/ProteinEvolverABC_Phase2.sh\",append=T)\n";	
	}	
else
	{
	print FILE_1 "	bodyhereNew<-paste(\"more \\\"$setsdir/ProteinEvolverABC_arguments.txt\\\" | xargs -n \",NumberOfArguments,\" \\\"$maindir/bin/ProteinEvolverProtABC1.2.0\\\"\",sep=\"\")\n";
	print FILE_1 "	write(paste(\"\",bodyhereNew,\"\",sep=\"\"),\"$setsdir/ProteinEvolverABC_Phase2.sh\",append=T)\n";	
	}


##############################################
##############################################
# Printing prior distributions for Recombination and Substitution rates
print FILE_1 "\n\n# Printing prior distributions for Recombination and Substitution rates \n";

print FILE_1 "\nfigureName1<-paste(\"$setsdir/Histogram_PriorRecombination.pdf\",sep=\"\")\n";
print FILE_1 "pdf(figureName1)\n";
print FILE_1 "par(mfrow = c(1,1))\n";
print FILE_1 "\nif (NumberReplicates < 1000)	{\n";
print FILE_1 "	hist (recombination_rate, breaks=NumberReplicates)\n";
print FILE_1 "	} else {\n";
print FILE_1 "	hist (recombination_rate, breaks=1000)\n";
print FILE_1 "	}\n";
print FILE_1 "dev.off()\n";

print FILE_1 "\nfigureName2<-paste(\"$setsdir/Histogram_PriorSubstitution.pdf\",sep=\"\")\n";
print FILE_1 "pdf(figureName2)\n";
print FILE_1 "par(mfrow = c(1,1))\n";
print FILE_1 "\nif (NumberReplicates < 1000)	{\n";
print FILE_1 "	hist (substitution_rate, breaks=NumberReplicates)\n";
print FILE_1 "	} else {\n";
print FILE_1 "	hist (substitution_rate, breaks=1000)\n";
print FILE_1 "	}\n";
print FILE_1 "dev.off()\n";

print FILE_1 "\nfigureName4<-paste(\"$setsdir/Histogram_PriorRho.pdf\",sep=\"\")\n";
print FILE_1 "pdf(figureName4)\n";
print FILE_1 "par(mfrow = c(1,1))\n";
print FILE_1 "\nif (NumberReplicates < 1000)	{\n";
print FILE_1 "	hist (RHO, breaks=NumberReplicates)\n";
print FILE_1 "	} else {\n";
print FILE_1 "	hist (RHO, breaks=1000)\n";
print FILE_1 "	}\n";
print FILE_1 "dev.off()\n";

print FILE_1 "\nfigureName5<-paste(\"$setsdir/Histogram_PriorTheta.pdf\",sep=\"\")\n";
print FILE_1 "pdf(figureName5)\n";
print FILE_1 "par(mfrow = c(1,1))\n";
print FILE_1 "\nif (NumberReplicates < 1000)	{\n";
print FILE_1 "	hist (THETA, breaks=NumberReplicates)\n";
print FILE_1 "	} else {\n";
print FILE_1 "	hist (THETA, breaks=1000)\n";
print FILE_1 "	}\n";
print FILE_1 "dev.off()\n";

#print FILE_1 "\nfigureName1<-paste(\"$setsdir/Histogram_PriorRecombination.jpg\",sep=\"\")\n";
#print FILE_1 "jpeg(figureName1, w=930, h=420)\n";
#print FILE_1 "par(mfrow = c(1,1))\n";
#print FILE_1 "\nif (NumberReplicates < 1000)	{\n";
#print FILE_1 "	hist (recombination_rate, breaks=NumberReplicates)\n";
#print FILE_1 "	} else {\n";
#print FILE_1 "	hist (recombination_rate, breaks=1000)\n";
#print FILE_1 "	}\n";
#print FILE_1 "dev.off()\n";

#print FILE_1 "\nfigureName2<-paste(\"$setsdir/Histogram_PriorSubstitution.jpg\",sep=\"\")\n";
#print FILE_1 "jpeg(figureName2, w=930, h=420)\n";
#print FILE_1 "par(mfrow = c(1,1))\n";
#print FILE_1 "\nif (NumberReplicates < 1000)	{\n";
#print FILE_1 "	hist (substitution_rate, breaks=NumberReplicates)\n";
#print FILE_1 "	} else {\n";
#print FILE_1 "	hist (substitution_rate, breaks=1000)\n";
#print FILE_1 "	}\n";
#print FILE_1 "dev.off()\n";


close (FILE_1);



##############################################
##############################################
# Compute prior distributions
print "\n\n> Computing prior distributions ... \n";
print "This phase requires the programming language R with the installed libraries: lattice, MCMCpack, ape, graphics\n";
print "An output file is created with the R outputs ...\n";
my $returnSystem = system ("R --vanilla < \"$setsdir/MakePriorsProteinEvolverABC.r\" > \"$setsdir/Output_MakePriorsProteinEvolverABC.txt\"");




##############################################
##############################################
#### END ####

if ( $returnSystem == 0)  {
	if (-e "$setsdir/OtherOutputs")
		{
		}
	else
		{
		system ("mkdir \"$setsdir/OtherOutputs\"");
		}
	system ("mv \"$setsdir\"/Histogram_Prior* \"$setsdir/OtherOutputs\"");
	system ("mv \"$setsdir\"/Output* \"$setsdir/OtherOutputs\"");
	system ("mv \"$setsdir\"/*.r \"$setsdir/OtherOutputs\"");
	print "\n> Successful!\n\n";
}
else  {
	print "\n> Error when running R script \"$setsdir/MakePriorsProteinEvolverABC.r\"\n\n";
}
exit;

