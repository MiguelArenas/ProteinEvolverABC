# pdb_real_data.pl is a script aimed to add the pdb chain that is given by the user to the real data alignment given by the user

# WARNING: IF STRUCTURE AND STABILITY SUMMARY STATISTICS ARE TO BE CONSIDERED, THE PDB SEQUENCE MUST BE ALIGNED AND HAVE NO GAPS -INDELS
# IT WILLL BE INSERTED IN THE FIRST LINE OF THE ALIGNMENT GIVEN 

use strict;
use warnings;
use File::Copy;
use File::Basename;
use Cwd;

# Get the real data alignment and the PDB sequence

my $phylip_file = $ARGV[0];
open(my $in_fh, '<', $phylip_file) or die "Cannot open file \"$phylip_file\": $!\n";

my $seq_pdb = $ARGV[1];
open(my $pdb_fh, '<', $seq_pdb) or die "Cannot open $seq_pdb: $!\n";
my $pdb_line = <$pdb_fh>;
chomp $pdb_line;
close $pdb_fh;

# Read the first line, containing the number of sequences and its length
my $header = <$in_fh>;
chomp $header;
my @header_parts = split(/\s+/, $header);

die "Invalid header line in $phylip_file\n" unless @header_parts >= 2;
die "First item in header is not a number.\n" unless $header_parts[0] =~ /^\d+$/;

$header_parts[0]++;
my $new_header = join(" ", @header_parts);

    # Read the rest of the file
    my @rest_lines = <$in_fh>;
    chomp @rest_lines;
    close $in_fh;

    # Create output file, renamed as '/name of the input file/_PDB'
    my $out_file = $phylip_file . '_PDB.phy';
    open(my $out_fh, '>', $out_file) or die "Cannot create $out_file: $!\n";
    print $out_fh $new_header, "\n";

    # Write the PDB sequence extracted from the phylip file (real data alignment) between the first and the second line
    print $out_fh "SeqPDB    $pdb_line\n";

    # Write the rest of the file as given
    foreach my $line (@rest_lines) 
    	{
   		print $out_fh $line, "\n";
   		}

close $out_fh;
print "\nDone!\n";
print "Modified alignment with PDB sequence written to: $out_file\n";


exit;
