### Miguel Arenas. 2025.
###
### Paths Prot_evol
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



##############################################
##############################################
# Reading Settings from input file
#print "> Reading Settings from input file ... \n";
while (<FROM>) # for each line
	{
	$TotalNumberLines++;
	}
close (FROM);


my $newfileX = $file;
$newfileX = sprintf ("%s_TEMP.c", $file);
open (FILE_OUTPUT, ">$newfileX");


open(FROM,$file);
while (<FROM>) 
	{
	$counterLines++;
	
	# Detect Target alignment input file
	if ($_ =~ /^char DIR_TNM\[/)
		{
        print FILE_OUTPUT "char DIR_TNM[100]=\"./\";\n";
		}
    elsif ($_ =~ /^char FILE_STR\[/)
        {
        print FILE_OUTPUT "char FILE_STR[200]=\"structures.in\";\n";
        }
    elsif ($_ =~ /^char tnm\[/)
        {
        print FILE_OUTPUT "char tnm[100]=\"bin/SSCPE/tnm\";\n";
        }

    elsif ($_ =~ /^char tnm_mut_para\[/)
        {
        print FILE_OUTPUT "char tnm_mut_para[100]=\"Mut_para.in\";\n";
        }
    elsif ($_ =~ /mv \*.mut_DE.dat /)
        {
        print FILE_OUTPUT "// $_";
        }
    elsif ($_ =~ /mv \*.mut_RMSD.dat/)
        {
        print FILE_OUTPUT "// $_";
        }
    elsif ($_ =~ /mv \*prof_DE.dat /)
        {
        print FILE_OUTPUT "// $_";
        }
    elsif ($_ =~ /mv \*prof_RMSD.dat /)
        {
        print FILE_OUTPUT "// $_";
        }
    else
        {
        print FILE_OUTPUT $_;
        }
        
	} # end of lines FROM

close (FROM);


system ("mv $newfileX $file");


print "\nPaths Prot_evol3 done!\n";
exit;
	
	
	
