# pdb_simulated_data.pl is a script aimed to align the pdb chain that is given by the user to the sequences alignments created by the model in the simulation step

# WARNING: IF STRUCTURE AND STABILITY SUMMARY STATISTICS ARE TO BE CONSIDERED, THE PDB SEQUENCE MUST BE ALIGNED 
# AND INSERTED IN THE FIRST LINE OF THE ALIGNMENT GIVEN 


use strict;
#use warnings;
use File::Copy;
use File::Basename;
use Cwd;

my $settings_file = $ARGV[0];
unless (open(FROM,$settings_file))
    {
    print STDERR "Cannot open file \"$settings_file\"\n\n";
    exit;
    }
#print "> Input file uploaded: $settings_filesettings_file \n\n";

#my $settings_file = 'Settings.txt';


# Directories
my $scrtdir = dirname($0);
#print ">  AAABC directory detected: $maindir \n\n";
my $setsdir = dirname($settings_file);
#print "> Directory of Settings.txt detected: $setsdir \n\n";




# Verify that Settings.txt exists
die "There wasn't found any $settings_file in the working directory :(.\n" unless -e $settings_file;

# Open Settings.txt
open(my $fh, '<', $settings_file) or die "The $settings_file couldn't be open: $!\n";

my $pdb_sequence;
while (my $line = <$fh>)
    {
    chomp $line;
    # if ($line =~ /^pdbsequence=(.+)$/)
    if ($line =~ /^SSpdbsequence=(.+)$/) # Modified by miguel as this is to calculate SS, thus to allow calculation of STR SSs is sequences simulated with empirical models
        {
        $pdb_sequence = $1;
        last;
        }
    }
close $fh;

die "There wasn't found any line headed as 'SSpdbsequence=' en $settings_file.\n" unless defined $pdb_sequence;

#print "\nDone!\n";
#print "PDB sequence file extracted!: '$pdb_sequence'\n";

# Verify that the pdb seq file exists in the current directory
if (-e $pdb_sequence) {
    #print "'$pdb_sequence' exists in the working directory\n";
} else {
    print "ERROR: '$pdb_sequence' was NOT found in the working directory!\n";
    exit;
}

open(my $pfh, '<', $pdb_sequence) or die "Could not open $pdb_sequence: $!\n";

# Read all the document, which contains the aligned pdb sequence that we need
my $lines = <$pfh>;
chomp $lines;

close $pfh;

print "PDB from $pdb_sequence succesfully extracted: '$lines'\n";

# Get the simulated data alignments
my @seq_files = glob("sequences*");
die "No 'sequences*' files found in current directory.\n" unless @seq_files;

#print "Found sequences files: @seq_files\n";

# Process each simulated data alignment
foreach my $seq_file (@seq_files) {
    #print "Processing file: $seq_file\n";

    open(my $in_fh, '<', $seq_file) or die "Cannot open $seq_file: $!\n";

    # Read the first line, containing the number of sequences and its length
    my $line1 = <$in_fh>;
    chomp $line1;

    # Add one unit to the number of sequences
    # Split by tabulation marks
    my @parts = split(/\s+/, $line1);
    if (@parts < 1) {
        warn "Warning: First line of $seq_file does not seem to have numbers.\n";
        close $in_fh;
        next;
    }
    # Add one to the first number
    if ($parts[0] =~ /^\d+$/) {
        $parts[0] = $parts[0] + 1;
    } else {
        warn "Warning: First element in first line of $seq_file is not a number.\n";
        close $in_fh;
        next;
    }
    my $modified_line1 = join(" ", @parts);

    # Read the rest of the file
    my @rest_lines = <$in_fh>;
    chomp @rest_lines;
    close $in_fh;

    # Create output file, renamed as '/name of the input file/_pdb'
    my $out_file = $seq_file . '_pdb';
    open(my $out_fh, '>', $out_file) or die "Cannot create $out_file: $!\n";

    # 
    print $out_fh $modified_line1, "\n";

    # Write the PDB sequence extracted between the first and the second line
    print $out_fh "PDB       $lines\n";

    # Write the rest of the file as given
    foreach my $line (@rest_lines) {
        print $out_fh $line, "\n";
    }

    close $out_fh;
    #print "PDB simulated data alignment succesfully created: $out_file\n";
}

print "All done!\n";
#print "Simulated data structure and stability summary statistics are ready to be estimated! :)\n\n";

exit ;
