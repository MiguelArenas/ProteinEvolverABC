###
### ComputeSSfromSimulatedData.pl, vs 1.5
### Miguel Arenas. 2013.
###
### Script for reading the input sequences file and creates a .ss file
###


use strict;
#use warnings;
use File::Copy;
use File::Basename;
use Cwd;


##############################################
##############################################
# Loading input file, "Sequences"
my $file = $ARGV[0];
unless (open(FROM,$file))
	{
	print STDERR "Cannot open file \"$file\"\n\n";
	exit;
	}
#print "> Input file uploaded: $file \n\n";	

my $save = $ARGV[1];

my $show = $ARGV[2];

#print "\n\n\n\n>>>>>> save: $save // show: $show  ..\n\n\n\n";


##############################################
##############################################
# Directories
my $scrtdir = dirname($0);
#print ">  ProtABC directory detected: $maindir \n\n";
my $setsdir = dirname($file);
#print "> Directory of Settings.txt detected: $setsdir \n\n";


##############################################
##############################################
# Reading file
#print " Reading input file $file... \n";
# variables
#my $nuc_number = -1;
#my $aa_number = -1;
#my $samplesize = 0;

#my $TotalNumberLines_target = 0;
#while (<FROM>) # for each line
#	{
#	$TotalNumberLines_target++;	
#	if ($TotalNumberLines_target == 1)
#		{
#		my $num3 = -1;
#		my $num4 = -1;	
#		my @seq1 = split(/\s+/, $_);
#		#($num3, $num4) = split(/ /, $_);	
#		$num3 = @seq1[0];
#		$num4 = @seq1[1];
#		$num4 =~ s/\n//g; # remove the end of line
#						
#		$samplesize = $num3;
#		$nuc_number = $num4;
#		print "\n  $file with sample size $samplesize and length $nuc_number nucleotides\n\n";
#		
#		$aa_number = $nuc_number/3;
#		
#		if ($nuc_number % 3 == 0) # coding data?
#			{
#			; # ok
#			}
#		else
#			{
#			print "\nERROR!. Is it coding data?. The number of nucleotides ($nuc_number) is not divisible by 3: $aa_number\n";
#			print "Type CTRL+C to abort the execution \n";
#			my $error = <STDIN>;
#			chop($error);
#			exit;	
#			}
#
#						
#		}	
#	}
#close (FROM);		


# Name of the file
my $newfileX = $file;
my $noThis = "";
my $numberThis = "";
($noThis, $numberThis) = split(/sequences/, $file);
$numberThis =~ s/\n//g; # remove the end of line	
my $newfileX = sprintf ("sequences%s", $numberThis);


##############################################
##############################################
# Compute SS
#print "> Calculating summary statistics from the target data set ..\n";
system ("perl \"$scrtdir/SummS_from_AA_Data_NoIndels.2.pl\" \"$file\" \"$show\"");


# Remove sequence?
if ($save == 0)
	{
	system ("rm \"$setsdir/$newfileX\"");
	}

#print "Summary statistics were successfully computed!\n";
exit;
	
	
	
