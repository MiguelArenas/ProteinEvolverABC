#!/usr/bin/env perl 

use strict ;  # Robust syntax
use warnings ;   # Display warnings if any trouble in the code

# Print presentation

print "Preparing the site-specific exchangeable matrices from ProtASR2 ..\n";

my $input  = "SitesMatrix.txt";
my $output = "UserEAAMsites";

# Remove a possible existing file with that name (i.e., from a previous execution)
if (-e $output) {
    unlink $output or die "It could not delete $output: $!";
    print "Deleted: $output\n";
}


open(my $in,  "<", $input)  or die "It could not open $input: $!";
open(my $out, ">", $output) or die "It could not create $output: $!";

# Reading the lines
my @lines = <$in>;
close($in);

# New file
print $out "# Exchangeability matrices\n";
print $out "A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";

# Exclude first two and last three lines of the original file

for (my $i = 2; $i < @lines - 3; $i++) {

    my $line = $lines[$i];

    if ($line =~ /^SITE\s+(\d+)\b/) {
        print $out "site $1\n";
    } else {
        print $out $line;
    }
}

close($out);

print "File saved: $output\n";
print "Done! \n";
print "An output file with all the matrices per site has been successfully created \n";

exit;



