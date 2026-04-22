###
### ProteinEvolverABC_Phase2.pl, vs 2.1
###
### Script for calculating the summary statistics of the real dataset and running the SSCPE and ProtASR2 models (if specified)
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
print "\nProteinEvolverABC_Phase2.pl";
print "\n- Calculation of the summary statistics of the real data";
print "\n- Running the SSCPE and ProtASR2 models (if specified)";
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
# Reading information


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
my $Subs_dip_MF_WT = -1;
my $SSCPE = 0; # If 0, SSCPE is not run. If 1, SSCPE runs
my $ProtASR2 = 0; # If 0, ProtASR2 is not run. If 1, ProtASR2 runs
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

# ELENA #
my $NameOfPhylipFile_PDB = "";
#######

my $NameOfPhylipFile_init = "";
my $NameOfFastaFile = "RealData.fas";

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
            my $NameOfPhylipFile_PDB = "${NameOfPhylipFile}_PDB.phy";
            print "\n    $setsdir/$NameOfPhylipFile_PDB was detected... Ok!"
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
		
## ELENA
    my $NameOfPhylipFile_PDB = "$setsdir/${NameOfPhylipFile}_PDB.phy";
    #print "ELENA 1 $NameOfPhylipFile_PDB";

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

    
            
    # Substituion model of amino acid evolution (i.e., Blosum62, CpRev, Dayhoff, DayhoffDCMUT, HIVb, HIVw, JTT, JonesDCMUT, LG, Mtart, Mtmam, Mtrev24, RtRev, VT, WAG, UserEAAM), and the SCS models (SSCPE, ProtASR2_WT, ProtASR2_MF)
    if ($_ =~ /^*SubstitutionModel=/)
        {
        ($Subs_init, $Subs_dip) = split(/=/, $_);
        $Subs_dip =~ s/\n//g; # remove the end of line
        if ($Subs_dip ne "Blosum62" && $Subs_dip ne "CpRev" && $Subs_dip ne "Dayhoff" && $Subs_dip ne "DayhoffDCMUT" && $Subs_dip ne "HIVb" && $Subs_dip ne "HIVw" && $Subs_dip ne "JTT" && $Subs_dip ne "JonesDCMUT" && $Subs_dip ne "LG" && $Subs_dip ne "Mtart" && $Subs_dip ne "Mtmam" && $Subs_dip ne "Mtrev24" && $Subs_dip ne "RtRev" && $Subs_dip ne "VT" && $Subs_dip ne "WAG" && $Subs_dip ne "UserEAAM" && $Subs_dip ne "UserEAAMsites" && $Subs_dip ne "SSCPE" && $Subs_dip ne "ProtASR2_WT" && $Subs_dip ne "ProtASR2_MF")
            {
            print "\n\nERROR!. Substitution model must be one of the following: Blosum62, CpRev, Dayhoff, DayhoffDCMUT, HIVb, HIVw, JTT, JonesDCMUT, LG, Mtart, Mtmam, Mtrev24, RtRev, VT, WAG, UserEAAM, UserEAAMsites, SSCPE ($Subs_dip), ProtASR2_WT ($Subs_dip), ProtASR2_MF ($Subs_dip) \n\n";
            print "Type CTRL+C to abort the execution \n";
            my $HereError = <STDIN>;
            chop($HereError);
            exit;
            }
        print "  Substitution model: $Subs_dip\n";
        $NumberOfArguments = $NumberOfArguments + 1; # add arguments -@

        if ($Subs_dip eq "SSCPE")
            {
            $SSCPE = 1;
            $Subs_dip = "UserEAAMsites";
            print "  Substitution model: $Subs_dip (SCS model)\n";
            }
            
        if ($Subs_dip eq "ProtASR2_WT")
            {
            $Subs_dip_MF_WT = 0;
            }
        if ($Subs_dip eq "ProtASR2_MF")
            {
            $Subs_dip_MF_WT = 1;
            }
            
        if (($Subs_dip eq "ProtASR2_WT") || ($Subs_dip eq "ProtASR2_MF"))
            {
            $ProtASR2 = 1;
            $Subs_dip = "UserEAAMsites";
            print "  Substitution model: $Subs_dip (SCS model)\n";
            }
                
        }
            
            
            
    # PDB file, required for SCS models (SSCPE, ProtASR2_WT, ProtASR2_MF)
        if ($_ =~ /^pdb=/)
        {
            ($PDB_init, $PDB_dip) = split(/=/, $_);
            $PDB_dip =~ s/\n//g; # remove the end of line
                    
            print "  PDB file: $PDB_dip\n";
            ## $NumberOfArguments = $NumberOfArguments + 1; # add arguments -@
                    
            if (($SSCPE == 1) || ($ProtASR2 == 1))
                {
                if (-e "$setsdir/$PDB_dip")
                    {
                    print "    $setsdir/$PDB_dip was detected... Ok!";
                    }
                else
                    {
                    print "\nERROR!. There is not a file \"$setsdir/$PDB_dip\" in this directory \n";
                    print "Type CTRL+C to abort the execution \n";
                    my $error = <STDIN>;
                    chop($error);
                    exit;
                    }
                }
                    
            if (($SSCPE == 1) || ($ProtASR2 == 1))
                {
                if ($PDB_dip eq "default")
                    {
                    print "\n\nERROR!. SCS model ($Subs_dip) requires a PDB file\n\n";
                    print "Type CTRL+C to abort the execution \n";
                    my $HereError = <STDIN>;
                    chop($HereError);
                    exit;
                    }
                }
            }
        
            
        # Chain of the PDB file, required for SCS models
        if ($_ =~ /^chain=/)
            {
                ($Chain_init, $Chain_dip) = split(/=/, $_);
                $Chain_dip =~ s/\n//g; # remove the end of line
                
                print "  Chain: $Chain_dip\n";
                ## $NumberOfArguments = $NumberOfArguments + 1; # add arguments -@
                
                if (($SSCPE == 1) || ($ProtASR2 == 1))
                {
                    if ($Chain_dip eq "default")
                    {
                        print "\n\nERROR!. SCS model ($Subs_dip) requires a PDB file with a Chain\n\n";
                        print "Type CTRL+C to abort the execution \n";
                        my $HereError = <STDIN>;
                        chop($HereError);
                        exit;
                    }
                }
            }
            
        
            
        # Sequence of the PDB file, required for SCS models
        if ($_ =~ /^pdbsequence=/)
            {
                ($SeqPDB_init, $SeqPDB_dip) = split(/=/, $_);
                $SeqPDB_dip =~ s/\n//g; # remove the end of line
                    
                print "  File with the PDB sequence: $SeqPDB_dip\n";
                ## $NumberOfArguments = $NumberOfArguments + 1; # add arguments -@
                    
                if (($SSCPE == 1) || ($ProtASR2 == 1))
                {
                    if ($SeqPDB_dip eq "default")
                    {
                        print "\n\nERROR!. SCS model ($Subs_dip) requires a file with the PDB sequence\n\n";
                        print "Type CTRL+C to abort the execution \n";
                        my $HereError = <STDIN>;
                        chop($HereError);
                        exit;
                    }
                }
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
                     
            # Sequence of the PDB file, required for SSCPE and ProtASR2, for SS
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
# Computes model-specific file with SSCPE and ProtASR2, using a PDB file ($PDB_dip)
# perl SSCPE.pl -pdb PDB/1zio.pdb -ali ALI/1ZIO_A_mafft_reduced.fasta -raxml
# Concatenate site-specific files
if (($SSCPE == 1) || ($ProtASR2 == 1))
    {
    if ($SSCPE == 1)
        {
        print "\n\n> Computing site-specific matrices and frequencies at the equilibrium with the SSCPE model .. \n";
        # copy the required material to the working directory
        system ("cp -r -f bin/SSCPE/* .");
        system ("cp -r -f bin/SSCPE/DIR_PROT_EVOL .");
        system ("cp -r -f bin/SSCPE/DIR_TNM .");
        system ("cp -r -f bin/SSCPE .");
        }
    if ($ProtASR2 == 1)
        {
        print "\n\n> Computing site-specific matrices and frequencies at the equilibrium with a ProtASR2 model .. \n";
        # copy the required material to the working directory
        system ("cp -r -f bin/Prot_evol_ProtASR2_exe .");
        system ("cp -r -f source/Prot_evol_ProtASR2_src/structures.in .");
        }
    # provide permissions to scripts, etc
    system ("chmod a+x *");
        
    #### ELENA #####
    system ("perl \"$maindir/scripts/PDBseq_real_data.pl\" \"$setsdir/$NameOfPhylipFile\" \"$setsdir/$SSSeqPDB_dip\" ");
    #print " ELENA 8: $setsdir/$NameOfPhylipFile_PDB $setsdir/$NameOfFastaFile";
    my $NameOfPhylipFile_PDB = "${NameOfPhylipFile}_PDB.phy";
    print "\n    $setsdir/$NameOfPhylipFile_PDB was detected... Ok!";
    # create the input file in fasta format
    system ("perl \"$maindir/scripts/Phylip2Fasta.pl\" \"$setsdir/$NameOfPhylipFile_PDB\" \"$setsdir/$NameOfFastaFile\" ");
    # For running SSCPE.pl: SSCPE should print the site-specific matrices and do not run RAxML-NG (not required), so the code of SSCPE.pl should be modified
    
    #print "\n INPUTS (SSCPE):    $setsdir/$PDB_dip";
    #print "\n INPUTS (SSCPE):    $setsdir/$NameOfFastaFile \n";
    #print "\n INPUTS (ProtASR2):    $setsdir/$PDB_dip";
    #print "\n INPUTS (ProtASR2):    $setsdir/$NameOfFastaFile";
        
        
    #### ELENA ####
    #print "ELENA $NameOfPhylipFile_PDB ----- $setsdir/$NameOfFastaFile";
    if ($SSCPE == 1)
        {
        system ("perl \"$maindir/SSCPE.pl\" -pdb \"$setsdir/$PDB_dip\" -ali \"$setsdir/$NameOfFastaFile\" -print_exch");
        }
        
    if ($ProtASR2 == 1)
        {
        # Make the input file
        my $out_M_file = "Prot_Evol_ProtASR2.in";
        if (-e "$out_M_file")
            {
            system ("rm $out_M_file");
            print "$out_M_file detected and replaced..\n";
            }
        open (OUT_M_FILE,">$out_M_file");
            
        # Here the parameters of the ProtASR2 model are specified according to the user-specifications, otherwise by default
            
            print OUT_M_FILE "#####################################################################\n";
            print OUT_M_FILE "#                                                                   #\n";
            print OUT_M_FILE "#         Parameters for the program Prot_evol                      #\n";
            print OUT_M_FILE "#                                                                   #\n";
            print OUT_M_FILE "#####################################################################\n";
            print OUT_M_FILE "#===================================================================\n";
            print OUT_M_FILE "# A) Input files defining the protein\n";
            print OUT_M_FILE "PDB=$PDB_dip # file_pdb\n";
            print OUT_M_FILE "CHAIN=  $Chain_dip\n";
            print OUT_M_FILE "#SEQ=    tpis.dna    # file_seq (optional)\n";
            print OUT_M_FILE "ALI=$NameOfFastaFile # Related proteins in FASTA format\n";
            print OUT_M_FILE "FILE_STR=structures.in\n";
            print OUT_M_FILE "## B) Substitution models\n";
            print OUT_M_FILE "MEANFIELD= 1        # Generate substitution models?\n";
            print OUT_M_FILE "# Subst. models are generated based on folding stability, struct.\n";
            print OUT_M_FILE "# conservation and combination of both. Parameters are amino acid\n";
            print OUT_M_FILE "# frequencies and selection parameter Lambda, which is determined\n";
            print OUT_M_FILE "# minimizing the KL divergence between model and regularized\n";
            print OUT_M_FILE "# distribution from PDB seq and input MSA.\n";
            print OUT_M_FILE "OPT_REG=    ! Automatically determine the regularization param.\n";
            print OUT_M_FILE "SCORE_CV=1   ! Optimize REG with Cv (1) or |KL_mod-KL_reg| (0)\n";
            print OUT_M_FILE "REG=     ! regularization param. if OPT_REG=0, starting value if OPT_REG=1\n";
            print OUT_M_FILE "MF_COMP=$Subs_dip_MF_WT    ! Perform (1) or omit (0) mean-field computations of stability constrained model (slow), otherwise only wild-type computation is performed (faster and often better performing).\n";
            print OUT_M_FILE "#===================================================================\n";
            print OUT_M_FILE "# C)  Amino acid frequencies\n";
            print OUT_M_FILE "REMUT=0        # Determine a.a. freq twice, the first time\n";
            print OUT_M_FILE "# by fitting observed frequencies with a.a. frequencies alone,\n";
            print OUT_M_FILE "# the second time fitting observed frequencies with full model\n";
            print OUT_M_FILE "# that includes selection.\n";
            print OUT_M_FILE "GET_FREQ=2        # Allowed: 0,1,2,3\n";
            print OUT_M_FILE "# 0= Use input mutation parameters\n";
            print OUT_M_FILE "# 1= Fit mutation parameters from prot sequences\n";
            print OUT_M_FILE "# 2= Combine fitted mutation model and a.a. frequencies\n";
            print OUT_M_FILE "# 3= Get background distribution from amino acid frequencies\n";
            print OUT_M_FILE "# Parameters of the mutation model if GET_FREQ=0:\n";
            print OUT_M_FILE "FREQ A 0.25        # f(A) (if GET_FREQ=0)\n";
            print OUT_M_FILE "FREQ T 0.25        # f(T) (if GET_FREQ=0)\n";
            print OUT_M_FILE "FREQ C 0.25        # f(C) (if GET_FREQ=0)\n";
            print OUT_M_FILE "FREQ G 0.25        # f(G) (if GET_FREQ=0)\n";
            print OUT_M_FILE "kCpG=2        # # Enhanced rate at CpG dinucleotides\n";
            print OUT_M_FILE "TT_RATIO=1.3        # transition-transversion ratio (>1)\n";
            print OUT_M_FILE "TWONUCMUT=0.25        # Ratio between 1-nuc and 2-nuc mutations\n";
            print OUT_M_FILE "#===================================================================\n";
            print OUT_M_FILE "# D) Exchangeability matrix\n";
            print OUT_M_FILE "EXCHANGE=FLUX    # Allowed: FLUX (default), RATE, EXCH, MUT\n";
            print OUT_M_FILE "MATRIX=JTT        # Empirical exchange matrix (JTT, WAG)\n";
            print OUT_M_FILE "#===================================================================\n";
            print OUT_M_FILE "# E) Thermodynamic model\n";
            print OUT_M_FILE "TEMP=    0.5        # Temperature\n";
            print OUT_M_FILE "SU1=    0.065        # configurational entropy per residue (unfolded)\n";
            print OUT_M_FILE "SC1=    0.065        # configurational entropy per residue (misfolded)\n";
            print OUT_M_FILE "SC0=    0.0        # configurational entropy offset (misfolded)\n";
            print OUT_M_FILE "REM=    2 (1,2,3)    # Use up to 1,2,3 moments of misfolding energy?\n";
            print OUT_M_FILE "A_LOC=    0        # Use secondary structure propensities?\n";
            print OUT_M_FILE "#===================================================================\n";
            print OUT_M_FILE "# F) Computations and output\n";
            print OUT_M_FILE "PRINT_E=1        # Print exchangeability matrix at all sites?\n";
            print OUT_M_FILE "FORMAT=PAML        # Use PAML format for exchangeability matrix\n";
            print OUT_M_FILE "ALL_MUTS=0            # Predict the effect of all nucl. mut?\n";
            print OUT_M_FILE "#===================================================================\n";
            print OUT_M_FILE "# G) Simulations of evolution\n";
            print OUT_M_FILE "TMAX=   000        # ITMAX: # of substitutions\n";
            print OUT_M_FILE "Samples= 5        # Independent trajectories simulated\n";
            print OUT_M_FILE "NEUTRAL= 0        # 1:Neutral fitness landscape 0: Fitness=1/(1+exp(DG/T))\n";
            print OUT_M_FILE "NPOP=    10        # effective population size (if MEANFIELD=0, NEUTRAL=0)\n";
            print OUT_M_FILE "#===================================================================\n";
                
        close (OUT_M_FILE);
        
            
        # Run to obtain the site-specific matrices
        system ("./\"$maindir/Prot_evol_ProtASR2_exe\" -file \"$setsdir/$out_M_file\" ");
        
            
        # Prepare the output of site-specific matrices
        # check the output files
        foreach my $file (glob("*"))
                {
                # it is a file (no directory)
                next unless -f $file;

                # pattern of the file
                if ($Subs_dip_MF_WT == 1)
                    {
                    if ($file =~ /_exchangeability_sites_/)
                        {
                        my $new_name = "SitesMatrix.txt";
                        copy($file, $new_name) or warn "The file could not be copied $file: $!";
                        print "Copy: $file -> $new_name\n";
                        }
                    }
            
                if ($Subs_dip_MF_WT == 0)
                    {
                    if ($file =~ /_exchangeability_sites_/)
                        {
                        my $new_name = "SitesMatrix.txt";
                        copy($file, $new_name) or warn "The file could not be copied $file: $!";
                        print "Copy: $file -> $new_name\n";
                        }
                    }
                }
                        
        # check format of the output file (matrices)
        system ("perl \"$maindir/scripts/makeEAAMsites_ProtASR2.pl\"");
        
        
        # cleaning
        system ("mkdir OutputsProtASR2");
        #system ("mv *_AA_profile_global.txt OutputsProtASR2");
        #system ("mv *_AArates.dat OutputsProtASR2");
        system ("mv -f *_rate_profile.dat OutputsProtASR2");
        system ("mv -f *_exchangeability_* OutputsProtASR2");
        system ("mv -f *_summary.dat OutputsProtASR2");
        #system ("mv *_likelihood.dat OutputsProtASR2");
        #system ("mv *_fitness.dat OutputsProtASR2");
        system ("mv -f *_Threading.dat OutputsProtASR2");
        system ("mv -f *_DeltaG.dat OutputsProtASR2");
        system ("mv -f Local_interactions.dat OutputsProtASR2");
        #system ("mv -f Contact_statistics_* OutputsProtASR2");
        system ("mv -f E_loc* OutputsProtASR2");
        system ("mv -f *_mutation.dat OutputsProtASR2");
        #system ("mv -f *_Contact_matrix.cm OutputsProtASR2");
        #system ("mv -f *_rate_mut.dat OutputsProtASR2");
        if ($Subs_dip_MF_WT == 0)
            {
            system ("mv -f *_WT_AA_profiles.txt OutputsProtASR2");
            }
        if ($Subs_dip_MF_WT == 1)
            {
            system ("mv -f *_MF_AA_profiles.txt OutputsProtASR2");
            }
        system ("mv -f SitesMatrix.txt OutputsProtASR2");
        system ("mv -f *.Prot_evol.log OutputsProtASR2");
        system ("mv -f *_entropy.dat OutputsProtASR2");
        system ("mv -f *_RealData.map OutputsProtASR2");
        system ("mv -f PDBseq.txt OutputsProtASR2");
        system ("mv -f *.DeltaG OutputsProtASR2");
        system ("mv -f *_Lambda.dat OutputsProtASR2");
        system ("mv -f *AA_profiles.txt OutputsProtASR2");
        system ("mv -f Prot_Evol_ProtASR2.in OutputsProtASR2");
        system ("mv -f OutputsProtASR2 OtherOutputs");
        system ("rm -f Prot_evol_ProtASR2_exe");
        }
                
    if ($SSCPE == 1)
        {
        # concatenate the site-specific matrices into a single file
        system ("perl \"$maindir/scripts/makeEAAMsites.pl\"");
        # cleaning
        system ("mkdir OutputsSSCPE");
        system ("mv *.site_*.txt OutputsSSCPE");
        system ("mv RealData.*.log OutputsSSCPE");
        system ("mv RealData.*.fasta OutputsSSCPE");
        system ("mv *.AA_profiles.txt OutputsSSCPE");
        system ("mv *.SSCPE.summary.dat OutputsSSCPE");
        system ("mv TNM_DATA OutputsSSCPE");
        system ("mv *.partitionsSite.txt OutputsSSCPE");
        system ("mv OutputsSSCPE OtherOutputs");
        system ("rm -f wag.txt Alignments.zip DE.zip Prot_evol.zip RMSD.zip tnm.zip");
        system ("rm -f SSCPE.pl script_w* script_r* script_i* script_d* README.md README_SSCPE qsubmit.pl Mutation_para.in lg.txt k2_mat.pl jtt.txt Input_TNM.in Input_Prot_evol.in");
        system ("rm -r -f DIR_TNM DIR_PROT_EVOL");
        system ("rm -f Prot_evol");
        system ("rm -f raxml-ng");
        system ("rm -r ALI");
        system ("rm -r PDB");
        }
                
    }


##############################################
##############################################
# Compute SS from input data ($file = Settings.txt)
# New version includes new SS (energies, etc)
print "\n\n> Computing summary statistics from input protein sequences file .. \n";

if ($SSPDB_dip eq "default") # SS from the protein structure are not calculated
    {
    if ($indels == 1)
        {
        system ("perl \"$maindir/scripts/ComputeProtSSfromRealData.pl\" \"$file\" \"$show_info\"");
        }
    else
        {
        system ("perl \"$maindir/scripts/ComputeProtSSfromRealData_NoIndels.pl\" \"$file\" \"$show_info\"");
        }
    }
else # SS from the protein structure are calculated
    {
    system ("cp -f bin/Prot_evol3 .");
        
    if (($SSCPE != 1) && ($ProtASR2 != 1))
        {
        # provide permissions to scripts, etc
        system ("chmod a+x *");
                
        #### ELENA #####
        system ("perl \"$maindir/scripts/PDBseq_real_data.pl\" \"$setsdir/$NameOfPhylipFile\" \"$setsdir/$SSSeqPDB_dip\" ");
        #print " ELENA 8: $setsdir/$NameOfPhylipFile_PDB $setsdir/$NameOfFastaFile";
        my $NameOfPhylipFile_PDB = "${NameOfPhylipFile}_PDB.phy";
        print "\n    $setsdir/$NameOfPhylipFile_PDB was detected... Ok!";
        # create the input file in fasta format
        system ("perl \"$maindir/scripts/Phylip2Fasta.pl\" \"$setsdir/$NameOfPhylipFile_PDB\" \"$setsdir/$NameOfFastaFile\" ");
        # For running SSCPE.pl: SSCPE should print the site-specific matrices and do not run RAxML-NG (not required), so the code of SSCPE.pl should be modified
        #print "\n INPUTS :    $setsdir/$SSSeqPDB_dip";
        #print "\n INPUTS :    $setsdir/$NameOfFastaFile \n";
        }

        
    if ($indels == 1)
        {
        system ("perl \"$maindir/scripts/ComputeProtSSfromRealData_STR.pl\" \"$file\" \"$show_info\" \"$SSPDB_dip\" \"$SSChain_dip\" \"$SSSeqPDB_dip\" ");
        }
    else
        {
        system ("perl \"$maindir/scripts/ComputeProtSSfromRealData_NoIndels_STR.pl\" \"$file\" \"$show_info\" \"$SSPDB_dip\" \"$SSChain_dip\" \"$SSSeqPDB_dip\" ");
        }
    }

# cleaning
if ($SSCPE == 1)
    {
    system ("rm *.mut_DE.dat");
    system ("rm *.mut_prof_RMSD.dat");
    system ("rm *.mut_RMSD.dat");
    system ("rm *.SSCPE.summary.dat");
    system ("rm *.tnm.summary.dat");
    }
if ($ProtASR2 == 1)
    {
    unlink $_ or warn "It did not remove $_: $!" for glob("*.Modes_*");
    #system ("rm *.Modes_*");
    system ("rm *.SSCPE.summary.dat");
    system ("rm *.mut_DE.dat");
    system ("rm *.mut_prof_RMSD.dat");
    system ("rm *.mut_RMSD.dat");
    system ("rm *.tnm.summary.dat");
    }

exit;
