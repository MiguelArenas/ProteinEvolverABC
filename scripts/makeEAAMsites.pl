#!/usr/bin/env perl 

use strict ;  # Robust syntax
use warnings ;   # Display warnings if any trouble in the code

# Print presentation

print "Concatenating site-specific exchangeable matrices into a single file ..\n";

# Open an output file in writing mode. Write headings

my $output_file = 'UserEAAMsites' ; 		# Declare the variable of the output file
open ( my $out_fh, '>', $output_file ) 	
	or die "Error opening the output file: $!"; 	# In case there is an error, $! prints information about the causes
	
print $out_fh "# Exchangeability matrices\n";
print $out_fh "A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";

# Obtain the files beginning for *.site_NUMBER.txt, being the number a variable field
# Use glob to find all the files finishing by .txt 

my @todos = glob("*.txt");

# Filter that dataset by using grep for the common part and d+ for the numerical part. $ indicates the end of the chain

my @matrix_per_site = grep { /\.site_\d+\.txt$/ } @todos; 

# Extract the number from every file name and compare it numerically, by using <=>

my @ordenados = sort {
    ($a =~ /\.site_(\d+)\.txt$/)[0] <=> ($b =~ /\.site_(\d+)\.txt$/)[0]
} @matrix_per_site;

# For each input that matches the grid, open it and read the information contained

foreach my $input_file (@ordenados) {
    if (-e $input_file) {
        open(my $in_fh, '<', $input_file) or die "Error opening the input file $input_file: $!";
        
# Extract the number of the file (biological meaning: protein site) and print it in the output file 

        if ($input_file =~ /\.site_(\d+)\.txt$/) {
            my $site_number = $1; 
            print $out_fh "site $site_number\n";  
        }
         # Read all lines from the input file
        my @lines = <$in_fh>;
        
        # Remove the last line 
        pop @lines if @lines > 1;
        
        # Concatenate the lines (without the last one) in the output file
        foreach my $line (@lines) {
            print $out_fh $line;  
        }

        close($in_fh);
    } else {
        warn "Warning: the file '$input_file' could not be found \n";
    }
}

# Close the output file
# Inform that the process has been completed successfully

close($out_fh);
print "Done! \n";
print "An output file with all the matrices per site has been successfully created \n";

exit;



