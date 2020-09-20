### Miguel Arenas. 2020.
###
### Script for computing the SS from the real data
###


use strict;
#use warnings;
use File::Copy;
use File::Basename;
use Cwd;


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

my $show = $ARGV[1];


##############################################
##############################################
# Directories
my $scrtdir = dirname($0);
#print ">  AAABC directory detected: $maindir \n\n";
my $setsdir = dirname($file);
#print "> Directory of Settings.txt detected: $setsdir \n\n";


##############################################
##############################################
# variables
my $TotalNumberLines = 0;
my $counterLines = 0;

my $NameOfPhylipFile = "";
my $NameOfPhylipFile_init = "";
my $aa_number = -1;
my $samplesize = 0;


##############################################
##############################################
# Reading Settings from input file
#print "> Reading Settings from input file ... \n";
while (<FROM>) # for each line
	{
	$TotalNumberLines++;
	}
close (FROM);

open(FROM,$file);
while (<FROM>) 
	{
	$counterLines++;
	
	# Detect Target alignment input file
	if ($_ =~ /^*NameOfPhylipFile/) 
		{
		($NameOfPhylipFile_init, $NameOfPhylipFile) = split(/=/, $_);
		$NameOfPhylipFile =~ s/\n//g; # remove the end of line
		print "  Target alignment file: $NameOfPhylipFile\n";	
			
		if (-e "$setsdir/$NameOfPhylipFile")
			{
			print "  \"$setsdir/$NameOfPhylipFile\" was detected... Ok!";			
			}
		else
			{
			print "\nERROR!. There is not a file \"$setsdir/$NameOfPhylipFile\" \n";
			print "Type CTRL+C to abort the execution \n";
			my $error = <STDIN>;
			chop($error);
			exit;
			}	
		}
		
		
	} # end of lines FROM

close (FROM);



##############################################
##############################################
# Reading Settings from data file
#print "\n> Reading Settings from data file ... \n";
#print " Reading input file $NameOfPhylipFile... \n";
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
		$aa_number = $num4;
		print "\n  $NameOfPhylipFile with sample size $samplesize and length $aa_number amino acids\n\n";
						
		}	
	}
close (FROM_targetAln);		



##############################################
##############################################
# Compute SS
print "> Calculating summary statistics from the target data set ..\n";
system ("perl \"$scrtdir/SummS_from_AA_Data_NoIndels.2.pl\" \"$setsdir/$NameOfPhylipFile\" \"$show\"");
system ("mv \"$setsdir/SS$NameOfPhylipFile.ss\" \"$setsdir/SSRealData.ss\"");


print "Summary statistics for the target data were successfully computed!\n";
exit;
	
	
	
