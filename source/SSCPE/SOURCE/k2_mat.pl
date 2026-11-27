#!/usr/bin/perl
#
# Ktreedist
# Calculation of the minimum branch length distance (K tree score) between phylogenetic trees
#
# Copyright (C) 2007 Victor Soria-Carrasco & Jose Castresana
# Institute of Molecular Biology of Barcelona (IBMB), CSIC, Jordi Girona 18, 08034 Barcelona, Spain
# vscagr@ibmb.csic.es (VSC) & jcvagr@ibmb.csic.es (JC)
#
# The present version has been modified by Ugo Bastolla ubastolla@cbm.csic.es
# Centro de Biologia Molecular Severo Ochoa (CSIC-UAM), Madrid Spain
# The program determines two types of branches and fits two different scales between the trees,
# inspired by branches where the protein structure is conserved or positively selected to vary
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# Requirements:
#	Perl v.5.8.x or greater (previous version have not been tested)
#
# Usage:
#  K2_mat.pl <options> -rt <reference tree file> -ct <comparison tree/s file>
#     Options:
#          -s [<filename>]  Output file for comparison tree/s after scaling
#          -t [<filename>]  Output file for table of results
#          -p [<filename>]  Output file for partitions table
#          -r  Output symmetric difference (Robinson-Foulds)
#          -n  Output number of partitions in the comparison tree/s
#          -a  Equivalent to all options
#          -write_graph Write reference tree in graph format
#          -write_fit   Write br. len fit for every node of comparison trees
#          -write_class Write classes joining nodes for which k2 is small
#
# For example:
# Ktreedist.pl -rt example_ref.tree -ct example_comp.tree -a
# or, if you have not made  Ktreedist.pl executable in a UNIX system:
# perl Ktreedist.pl -rt example_ref.tree -ct example_comp.tree -a

use warnings;
use strict;


# Ktreedist  version 1.0
my $normalize=2; # Normalize K and K2 score by average branch length?
my $normalize_RF=1; # Normalize RF distance by maximum value
my $version = "2.0";
my $monthyear = "May 2024";
my $write_fit=0;
my $write_graph=0;
my $write_class=0;

# //////////////////////////////////////////////////////////////////////////////
# //////////////////////////////  HIDDEN OPTIONS  //////////////////////////////
# //////////////////////////////////////////////////////////////////////////////

# The variable hard_polytomies is "off" by default. By default, zero branch length partitions count in the calculation of the symmetric difference.
# With hard_polytomies "on", zero branch length partitions are collapsed and treated as polytomies for this calculation.
my $hard_polytomies = "off";

# Debug mode is "off" by default. With debug mode "on" you do not need to enter -rt and -ct for the reference tree and
# comparison tree/s files, respectively, but combinations of parameters to be entered may be less flexible.
my $debugmode = "off";


# Safe mode is "off" by default. With safe mode "on" you are warned in case of overwriting files. When safe mode is "on"
# a new visible option appears to control overwriting:
#          -f  Overwrite output files
my $safemode = "off";


# Number of decimal places for branch lenghts in scaled trees
my $precision = 10;

my $supporting=0;

# //////////////////////////////////////////////////////////////////////////////
# ///////////////////////////////  PRELIMINARS  ////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////

#print "\nKtreedist version $version -  $monthyear";
#print "\nCalculation of the minimum branch length distance (K tree score) between trees\n\n";
print "\nK2treedist version 1.0 -  Oct 2023";
print "\nCalculation of the minimum branch length distance between trees after 2 branch types fitting\n\n";

# Get command line arguments
my ($reftree, $comptrees, $scaled, $table, $partitions, $robinson, $npartitions, $overwrite) = get_cmd_line(@ARGV);
for(my $i=0; $i<scalar(@ARGV); $i++){
    if($ARGV[$i] eq "-write_graph"){$write_graph=1;}
    elsif($ARGV[$i] eq "-write_fit"){$write_fit=1;}
    elsif($ARGV[$i] eq "-write_class"){$write_class=1;}
}

# DEFAULT DEFINITIONS, PATHS & FILENAMES
#-------------------------------------------------------------------------------
# Get paths & filenames
my ($filename_reftree, $filename_reftree_woext, $path_reftree) = get_path_filename($reftree);
my ($filename_comptrees, $filename_comptrees_woext, $path_comptrees) = get_path_filename($comptrees);

#This is really necessary in case of no full path
$reftree = "$path_reftree$filename_reftree" if ($reftree ne "");
$comptrees = "$path_comptrees$filename_comptrees" if ($comptrees ne "");

my $default_scaled = "$path_comptrees$filename_comptrees.scaled";
my $default_table = "$path_comptrees$filename_comptrees.tab";
my $default_partitions = "$path_comptrees$filename_comptrees.part";

# EXTRAS
$scaled = $default_scaled if ($scaled eq "yes");
$table = $default_table if ($table eq "yes");
$partitions = $default_partitions if ($partitions eq "yes");

# Do initial checkings and get trees if there are no errors
print "Checkings:\n";
my ($reftr, $comptr, $nametr, $txbl, $root, $bl, $t0,  $tr1) = initial_check ($reftree, $comptrees, $scaled, $table, $partitions);
$reftree = ($$reftr);
$comptrees = ($$comptr);
my %nametrees = (%$nametr);
my $taxblock = ($$txbl);
my $rooted = ($$root);
my $bk_line = ($$bl);
my $tree0 = ($$t0);
my @trees1 = (@$tr1);


# //////////////////////////////////////////////////////////////////////////////
# ///////////////////////////////////  MAIN  ///////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////


if (scalar(@trees1) == 1){
    print "\nK tree score calculation:\n'$filename_reftree' vs '$filename_comptrees'\n";
}
else{
    print "\nK tree score calculation:\n'$filename_reftree' vs '$filename_comptrees' (".scalar(@trees1)." trees)\n";
}

my ($claA, $claB, $brl, $sp) = get_partitions_tree($tree0, "yes");
my @clados0A = (@$claA);
my @clados0B = (@$claB);
my %brlen0 = (%$brl);
my @species0 = (@$sp);
my $nsp0=scalar(@species0);
my @clados0 = (@clados0A, @clados0B);
my @nodes0=get_nodes(\@clados0A, \%brlen0, $filename_reftree, 0, 1,
		     $write_graph);
#print "Nodes_ref:";
#for(my $i=0; $i<scalar(@nodes0); $i++){print " ", $nodes0[$i][1];}
#print "\n";
Put_sons(\@nodes0);

my $ntree = 0;
my @tab = ();
my @part = ();
my @scaled_trees = ();

my $nameref = "";
if (exists $nametrees{$reftree}{0}){ $nameref = $nametrees{$reftree}{0}; }
else{ $nameref = $filename_reftree_woext; }
my $rtd = "";
if ($rooted eq "yes"){ $rtd = "rooted";}
else{ $rtd = "unrooted";}

my $sum_br_ref=0;
my $min=0;
foreach my $br (values %brlen0){
    $sum_br_ref+=$br;
    if ($br == 0){$min++;}
}
$sum_br_ref/=2;

my $n_clados0 = scalar(@clados0)/2;

if ($hard_polytomies eq "on"){
    print "Zero branch length partitions, if present, are collapsed.\n";
    $n_clados0 = $n_clados0 - $min/2;
    $min = 0;
}

print "Reference tree ($nameref) is $rtd and has ".scalar(@species0).
" tips and ".$n_clados0." partitions.\n";
print "Number of partitions on comparison tree/s are listed below.\n"
if ($npartitions eq "yes");

my $initial_lines = 
"Comparison_tree  0p-fit(BS) 1p-fit(KS) Scale 2p-fit(K2S) Scale_0 Scale_1 ";
$initial_lines = $initial_lines."  RF" if ($robinson eq "yes");
$initial_lines = $initial_lines."  N_br" if ($npartitions eq "yes");
$initial_lines = $initial_lines."  Support(K2)";
$initial_lines = $initial_lines."  Sum_br";
$initial_lines = "\n  $initial_lines  \n-";

my $len = (length($initial_lines)-3);
for (my $i=0; $i < $len; $i++){ $initial_lines = "-".$initial_lines."-";}
$initial_lines =~ s/-$/\n/g;
print "$initial_lines";

my $sum_BS=0; my $sum_Ks1=0; my $sum_Ks2=0; my $sum_RF=0;
my $sum_K1=0; my @sum_K2; $sum_K2[0]=0; $sum_K2[1]=0; 
my $sum_npart=0; my $sum_support=0; my $sum_n=0;
my $num_tree=0;

foreach my $tree1 (@trees1){
	my ($cla, $clb, $brl) = get_partitions_tree($tree1, "yes");
	my @clados1A = (@$cla);
	my @clados1B = (@$clb);
	my @clados1 = (@clados1A, @clados1B);
	my %brlen1 = (%$brl);
	my $n_clados1 = scalar (@clados1)/2;

	my @nodes1=get_nodes(\@clados1A, \%brlen1, $filename_comptrees,
			     $ntree, scalar(@trees1), $write_graph);

	if ($hard_polytomies eq "on"){
	    my $min=0;
	    foreach my $br (values %brlen1){
		if ($br == 0){$min++;}
	    }
	    $n_clados1 = $n_clados1 - $min/2;		
	}
	
	# Get matrix for each pair of trees
	my @matrix = get_matrix (\@clados0A,\@clados0B,
				 \@clados1A,\@clados1B, 
				 \%brlen0, \%brlen1,
				 \@nodes0, \@nodes1);
	# Matrix array structure
	# First index: Index of partition
	# Second index:
	# 0			      -	partition (tips)
	# 1			      -	symmetric partition (tips)
	# 2			      -	branch length reference tree
	# 3			      -	branch length comparison tree
	# 4			      - branch length reference tree "NONE"
	# 5			      - branch length comparison tree "NONE"
	# 6                           - node name
	# 7                           - parent node name
	# After scale factor
	# 8			      - branch length comparison tree after scaling 1-par
	# 9			      - branch length comparison tree after scaling "NONE"

	# Get K-score
	my ($K_value, $sum_br) = get_K_value(@matrix);
	my ($ks, $rf, $bs, $ma) = get_K_score (\@matrix, \$K_value, $nsp0);
	my $K_score = ($$ks);  
	my $drf = ($$rf);      # RF score
	my $BS = ($$bs); 	#Brach Lenght Score.
	@matrix = (@$ma);

	# K2 score
	my ($kk, $zz) = 
	    get_K_value_EM(\@matrix, $filename_comptrees, $ntree,
			   scalar(@trees1));
	my @z2=(@$zz);
	my @k2=(@$kk);
	if($k2[0]==0 && $k2[1]==0){next;}
	$num_tree++;

	if($write_class){
	    # Join if z2==0
	    Write_classification(\@nodes0, \@matrix, \@z2,
				 $filename_reftree, $num_tree);
	}

	my $k2_score = get_k2_score(\@matrix, \@k2, \@z2);

	# Scale comparison tree & output if -s was invoked
	# my ($tr, $nwbrl, $sct) = modify_brlen_tree(\$scaled, \$tree1, \%brlen1, \$K_value);
	my ($tr, $sct) = modify_brlen_tree_2(\$scaled, \$tree1, \@k2, \@z2);
#$nwbrl, \%brlen1, 
	$tree1 = ($$tr);
	#%brlen1 = (%$nwbrl);
	push (@scaled_trees, ($$sct));

	# //// OUTPUT DISPLAY ////
	my $output = "  ".output_display(\$comptrees, \$ntree, \%nametrees, 
					 \$BS, \$K_score, \$K_value,
					 \$k2_score, \@k2,
					 \$robinson, \$drf, \$n_clados1,
					 \$sum_br);
	
	print "$output";
	$sum_n++; $sum_RF+=$drf;
	$sum_BS+=$BS; $sum_Ks1+=$K_score; $sum_K1+=$K_value;
	$sum_Ks2+=$k2_score; $sum_K2[0]+=$k2[0]; $sum_K2[1]+=$k2[1];
	$sum_npart+=$n_clados1; $sum_support+=$supporting;
	#$sum_br_ref+=$sum_br;

	# //// GENERATE ARRAYS FOR OUTPUT FILES ////
	
	push (@tab, get_table ( \$reftree, \$comptrees, \%nametrees, 
				\$filename_comptrees_woext, \$ntree,
				\$npartitions, 
				\$BS, \$K_score, \$K_value, \$k2_score, \@k2,
				\$robinson, \$drf,
				\$n_clados0, \$n_clados1, \$sum_br))
	    if ($table ne "no");
	
	push (@part,
	      get_parts_table ( \$comptrees, \%nametrees, \$ntree, \@matrix)) 
	    if ($partitions ne "no");

	$ntree++;
}


$sum_RF/=$sum_n;
$sum_BS/=$sum_n; $sum_Ks1/=$sum_n; $sum_K1/=$sum_n;
$sum_Ks2/=$sum_n; $sum_K2[0]/=$sum_n; $sum_K2[1]/=$sum_n;
$sum_npart/=$sum_n; $sum_support/=$sum_n;

if ($table ne "no"){
    $nametrees{$comptrees}{$ntree}="(Mean)";
    $supporting=$sum_support;
    push (@tab,
	  get_table ( \$reftree, \$comptrees, \%nametrees, 
		      \$filename_comptrees_woext, \$ntree,
		      \$npartitions, 
		      \$sum_BS, \$sum_Ks1, \$sum_K1, \$sum_Ks2, \@sum_K2,
		      \$robinson, \$sum_RF,
		      \$n_clados0, \$sum_npart, \$sum_br_ref));
}

my $final_line = "";
for (my $i=0; $i < $len; $i++){ $final_line = $final_line."-"; }
$final_line = $final_line."\n";
print $final_line;


# //////////////////////////////////////////////////////////////////////////////
# ///////////////////////////////  OUTPUT  FILES ///////////////////////////////
# //////////////////////////////////////////////////////////////////////////////


print "\nFile saving:" if (($table ne "no") || ($partitions ne "no") || $scaled ne "no");
print "\n";

# Output scaled trees
if ($scaled ne "no"){
	print "Generating output file for comparison tree/s after scaling... ";
	#Avoid writing to this file, will not be used.
	#---------------------------------------------
	#open (FILE, ">$scaled") or die ("Cannot write to $scaled");
	#	foreach my $st (@scaled_trees){
	#		print FILE "$st\n";
	#	}
	#close (FILE);
	newick2nexus (\$scaled, \$taxblock, \%nametrees, \$comptrees) if (exists ($nametrees{$comptrees}{0}));
	unix2macdos($bk_line, $scaled) if ($bk_line =~ m/^MAC|DOS$/);
	print "Done\n";
}

# 	Output results in a tabulated text file
if ($table ne "no"){
	print "Generating table of results... ";
	open (FILE, ">$table") or die ("Cannot write to $table");
	my $line = "Trees\t0-fit(BS)\t1-fit(KS)\tScale\t2-fit(K2S)\tScale_0\tScale_1";
	$line = $line."\tRF" if ($robinson eq "yes");
	$line = $line."\tN_partitions_ref_tree" if ($npartitions eq "yes"); #\tN_partitions_comp_tree
	$line = $line."\tSupport(K2)";
	$line = $line."\tsum_br";
	print FILE "$line\n";
	#my $i=0;
	foreach my $l (@tab){
	    print FILE "$l\n";
	    #$i++; if($i==$ntree){print FILE "# Average results:\n";}
	}
	if($normalize){
	    print FILE "# K and K2 score normalized by mean branch length\n";
	}
	if($normalize_RF){
	    print FILE "# RF distance normalized by maximum value\n";
	}
	close (FILE);
	unix2macdos($bk_line, $table) if ($bk_line =~ m/^MAC|DOS$/);
	print "Done\n";
}


#Avoid writing to this file, will not be used.
#---------------------------------------------
# 	Output partitions table
#if ($partitions ne "no"){
	#print "Generating table of partitions... ";
	#open (FILE, ">$partitions") or die ("Cannot write to $partitions");
	#	print FILE "File generated by Ktreedist version $version\n\n";
	#	print FILE "Reference file: '$reftree'\n";
	#	print FILE "Comparison file: '$comptrees'\n";
	#	print FILE "Zero branch-length partitions, if present, are collapsed.\n" if ($hard_polytomies eq "on");
	#	print FILE "\n";
	#	print FILE "Meaning of each column:\n";
	#	print FILE "\tBrl_ref_tree: Branch length of this partition on the reference tree.\n";
	#	print FILE "\tBrl_comp_tree: Branch length of this partition on the original comparison tree.\n";
	#	print FILE "\tBrl_comp_tree_K: Branch length of this partition on the comparison tree after scaling.\n";
	#	print FILE "\tPartition: Tip nodes that constitute this partition.\n\n";
	#	foreach my $l (@part){
	#		print FILE "$l";
	#	}
	#close(FILE);
	#unix2macdos($bk_line, $partitions) if ($bk_line =~ m/^MAC|DOS$/);
	#print "Done\n";
#}
print "\n";

exit();


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------- SUBROUTINES ----------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Get paths & filenames (with and without extension) ///
sub get_path_filename{
	my $filename = shift(@_);
	my $path = $filename;
	$filename =~ s/.+[\/|\\]//g;
	my $filename_woext = $filename;
	$filename_woext =~ s/\..+//g;
	$path =~ s/\/*$filename$//g;
	$path = "./".$path if ($path !~ m/^[A-Z]*\:*[\/|\.\/|\\]/);
	$path = $path."/" if ($path !~ m/[\/|\\]$/);

	return($filename, $filename_woext, $path);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Print usage howto ///
sub usage{
    my $wrong_option = shift(@_);
    print "WARNING: Option '$wrong_option' not recognized.\n\n" if (defined ($wrong_option));
    print "Usage:\n";
    if ($debugmode eq "on"){
	print " Ktreedist.pl <reference tree file> <comparison tree/s file> [<options>]\n";
    }
    else{
	print " Ktreedist.pl -rt <reference tree file> -ct <comparison tree/s file> [<options>]\n";
    }
    print "    Options:\n";
    print "         -t [<filename>]  Output file for table of results\n";
    print "         -p [<filename>]  Output file for table of partitions\n";
    print "         -s [<filename>]  Output file for comparison tree/s after scaling\n";
    print "         -r  Output symmetric difference (Robinson-Foulds)\n";
    print "         -n  Output number of partitions in the comparison tree/s\n";
    print "         -f  Overwrite output files\n" if ($safemode eq "on");
    print "         -write_graph Write ref tree in graph format\n";
    print "         -write_fit   Write fit for each branch of comp trees\n";
    print "         -write_class Write classes joining nodes with small k2 (z2[i]=0)\n";
    print "         -a  Equivalent to all options\n";
    print "         -norm  Not normalize K score with branch length\n";
    print "         \n";
    exit();
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Get command line parameters ///
sub get_cmd_line{
	my $get = "";
	my $reference_tree = "";
	my $comp_trees = "";
	my $scaled = "no";
	my $table = "no";
	my $partitions = "no";
	my $npartitions = "no";
	my $robinson = "no";
	my $overwrite = "yes";
	$overwrite = "no" if ($safemode eq "on");

	my @cmd = @_;

	if ($debugmode eq "on"){
# 		print "\nCaution: Debug mode is on.\n\n";
	    $reference_tree = $ARGV[0] if (exists ($ARGV[0]));
	    $comp_trees = $ARGV[1] if (exists ($ARGV[1]));
	    
	    if (scalar(@ARGV) == 2){
		@cmd = ();
		push (@cmd, "-r");
		push (@cmd, "-n");
	    }
	}

	for (my $i=0; $i < scalar (@cmd); $i++){
	    if ($cmd[$i] =~ m/^-[a|A]$/){
		$scaled = "yes";
		$table = "yes";
		$partitions = "yes";
		$robinson = "yes";
		$npartitions = "yes";
		$overwrite = "yes";
		next;
	    }
	    if ($cmd[$i] =~ m/^-[h|H]$/){
		usage();
	    }
	    if ($cmd[$i] =~ m/^-[s|S]$/){
		$get = "scaled";
		$scaled = "yes";
		next;
	    }
	    elsif ($cmd[$i] =~ m/^-[t|T]$/){
		$get = "table";
		$table = "yes";
		next;
	    }
	    elsif ($cmd[$i] =~ m/^-[p|P]$/){
		$get = "partitions";
		$partitions = "yes";
		next;
	    }
	    elsif ($cmd[$i] =~ m/^-[r|R]$/){
		$robinson = "yes";
		next;
	    }
	    elsif ($cmd[$i] eq "-norm"){
		$normalize=2; next;
	    }
	    elsif ($cmd[$i] =~ m/^-[n|N]$/){
		$npartitions = "yes";
		next;
	    }
	    elsif ($cmd[$i] =~ m/^-[r|R][t|T]$/){
		$get = "reference";
		next;
	    }
	    elsif ($cmd[$i] =~ m/^-[c|C][t|T]$/){
		$get = "comp";
		next;
	    }
	    elsif ($cmd[$i] =~ m/^-[f|F]$/){
		$overwrite = "yes";
		next;
	    }
	    
	    #*****
	    
	    if ($get eq "scaled" && $scaled eq "yes"){
		$scaled = $cmd[$i];
		next;
	    }
	    if ($get eq "table" && $table eq "yes"){
		$table = $cmd[$i];
		next;
	    }
	    if ($get eq "partitions" && $partitions eq "yes"){
		$partitions = $cmd[$i];
		next;
	    }
	    elsif ($i<2 && $debugmode eq "on") {next;}
	    if ($get eq "reference" && $reference_tree eq ""){
		$reference_tree = $cmd[$i];
		next;
	    }
	    elsif ($i<2 && $debugmode eq "on") {next;}
	    if ($get eq "comp" && $comp_trees eq ""){
		$comp_trees = $cmd[$i] ;
		next;
	    }
	    
	    if ($cmd[$i] !~ m/^-[a|h|s|t|c|f]$/i){
		if($cmd[$i] ne "-write_graph" &&
		   $cmd[$i] ne "-write_fit" &&
		   $cmd[$i] ne "-write_class"){
		    usage($cmd[$i]);
		}
	    }
	}
	usage () if (($reference_tree eq "") || ($comp_trees eq ""));
	return ($reference_tree, $comp_trees, $scaled, $table, $partitions, $robinson, $npartitions, $overwrite);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# This subroutine checks several parameters before starting tree comparison:
#	- if input and output files exist (and can be writen)
#	- kind of line break
#	- input files format (nexus needs internal conversion to newick)
#	- if all trees are rooted or unrooted
#	- if trees have branch length information
#	- if trees have same tip names
sub initial_check{
    my $reftree = shift(@_);
    my $comptrees = shift(@_);
    my $scaled = shift(@_);
    my $table = shift(@_);
    my $partitions = shift(@_);
    my %nametrees = ();
    my $nexus = "";
    my $taxblock = "no";
    my $out = "no";

    my $rooted = "no";
    
    my $ERR = 0;
    my $ERROR = "";
    my $WARNING = "";
    
    # =====================================
    # Check existence of files
    # =====================================
    if (($comptrees eq "") || ($reftree eq "")){
	usage();
    }
    elsif ((!-e $reftree) || (!-e $comptrees)){
	if ((!-e $reftree) && (!-e $comptrees)){
	    $ERROR = "     Files '$reftree' & '$comptrees' do not exist\n";
	}
	elsif (!-e $reftree){
	    $ERROR = "     File '$reftree' does not exist\n";
	}
	else{
	    $ERROR = "     File '$comptrees' does not exist\n";
	}
	die("ERROR:\n$ERROR\n");
    }
    
    # =====================================
    # =====================================
    
    # =====================================
    #Check line break
    # =====================================
    print "Line break...  ";
    my @lines_ref = break_line($reftree);
    my @lines_comp = break_line($comptrees);
    my $breakline_ref = shift(@lines_ref);
    my $breakline_comp = shift(@lines_comp);
    my $bk_line_line = "OK ($breakline_ref & $breakline_comp)";
    
    print "$bk_line_line\n";
    my $bk_line = $breakline_comp;
    
    # =====================================
    # =====================================
    
    # =====================================
    # Check output files
    # =====================================
    if (($scaled ne "no") || ($table ne "no") || ($partitions ne "no")){
	if ($overwrite eq "no"){
	    my $printout = "";
	    if (-e $scaled){
		$printout = "WARNING: File '$scaled' exists. Do you want to overwrite it? (y/n) ";
	    }
	    if (-e $table){
		if ($printout ne ""){
		    $printout = "WARNING: Files '$scaled' & '$table' exist. Do you want to overwrite them? (y/n) ";
		}
		else{
		    $printout = "WARNING: File '$table' exists. Do you want to overwrite it? (y/n) ";
		}
	    }
	    if (-e $partitions){
		if ($printout ne ""){
		    if (grep(/them/,$printout)){
			$printout = "WARNING: Files '$scaled', '$table' & '$partitions' exist. Do you want to overwrite them? (y/n) ";
		    }
		    else{
			if (grep (/$table/, $printout)){
			    $printout = "WARNING: Files '$table' & '$partitions' exist. Do you want to overwrite them? (y/n) ";
			}
			else{
			    $printout = "WARNING: Files '$scaled' & '$partitions' exist. Do you want to overwrite them? (y/n) ";
			}
		    }
		}
		else{
		    $printout = "WARNING: File '$partitions' exists. Do you want to overwrite it? (y/n) ";
		}
	    }
	    if ($printout ne ""){
		my $answer = "";
		while ($answer !~ /y|Y|n|N/){
		    print $printout;
		    $answer = <STDIN>;
		}
		chomp($answer);
		if ($answer !~ m/^y|Y$/){
		    die("ERROR:\n     Output files overwrite not allowed.\n\n");
		}
		else{
		    unlink ("$scaled");
		    unlink ("$table");
		    unlink ("$partitions");
		}
	    }
	}
	else{
	    unlink ("$scaled");
	    unlink ("$table");
	    unlink ("$partitions");
	}
    }
    # =====================================
    # =====================================
    
    # =====================================
    # CHECK if files are nexus or newick
    # =====================================
    print "File format...  ";
    
    my @lines = ([@lines_ref], [@lines_comp]);
    my $fileformat_line = "";
    my $formaterr = "";
    my @reftrees = ();
    my @comptrees = ();
    for (my $l=0; $l < scalar(@lines); $l++){
	my $newick = "yes";
	my @aux = @{$lines[$l]};
	
	# Check newick
	# -------------------------------------
	my $save = "";
	foreach my $a (@aux){
	    next if ($a =~ m/^[\n|\/|\#|\>]/);
	    $a =~ s/^\[.+\]//g if ($a =~ m/^\[.+\]/);
	    $save = $save.$a;
	    chomp($save);
	}
	my @trs = split (/;/, $save);
	foreach (@trs){
	    $_ = $_.";";
	    my @aux = split(/\,/, $_);
	    my $n_intercom = scalar(@aux);
	    @aux = split(/ *\: *[0-9]+\.*[0-9]*e*E*\-*[0-9]* *| *\, *| *\( *| *\) *[0-9]*\.*[0-9]* *| *\;/, $_);
	    my $n_spp = 0;
	    foreach my $a (@aux){
		$n_spp++ if ($a ne "");
	    }
	    @aux = split (/\(/,$_);
	    my $o_par = scalar(@aux);
	    @aux = split (/\)/,$_);
	    my $c_par = scalar(@aux);
	    if ((!/^\(/ && !/\;$/) || ($n_intercom != $n_spp) || ($o_par != $c_par)){
		$newick = "no";
		last;
	    }
	}
	if ($newick eq "yes"){
	    if ($fileformat_line eq ""){
		$fileformat_line = "newick";
	    }
	    else{
		$fileformat_line = $fileformat_line." & newick";
	    }
	    if ($l == 0){ @reftrees = @trs; }
	    else { @comptrees = @trs; }
	    
	}
	
	
	# Check nexus
	# -------------------------------------
	if ($newick eq "no"){
	    my $nex = "no";
	    my $confirmnex = "no";
	    foreach (@aux){
		next if (/^[\n|\/|\>]|^[ |\t]+$/);
		if (/^[ |\t]*\#NEXUS/){
		    $nex = "yes";
		    next;
		}
		if ($nex eq "yes"){
		    if (/TREE/i){
			if ($fileformat_line eq ""){
			    $fileformat_line = "nexus";
			}
			else{
			    $fileformat_line = $fileformat_line." & nexus";
			}
			$confirmnex = "yes";
			last;
		    }
		}
	    }
	    
	    if (($nex eq "no") || ($nex eq "yes" && $confirmnex eq "no")){
		$nexus = $nexus."no";
	    }
	    else{
		my ($nametr, $txbl, $trs) = nexus2newick(\@aux);
		my %auxname = (%$nametr);
		$taxblock = ($$txbl);
		if ($l == 0){ @reftrees = (@$trs); }
		else{ @comptrees = (@$trs); }
		$fileformat_line = $fileformat_line." with taxa block" if ($taxblock eq "yes");
		if ($l == 0){ $nametrees{$reftree} = \%auxname; }
		else { $nametrees{$comptrees} = \%auxname; }
		$nexus = $nexus."yes";
	    }
	}
	else{ $nexus = $nexus."no"; }
	
	if (($newick eq "no") && ($nexus =~ m/no$/)){
	    if ($formaterr ne ""){ $formaterr = $formaterr."\n          &\n          format of file "; }
	    else{ $formaterr = "Format of file "; }
	    if ($l == 0) { $formaterr = $formaterr." '$reftree' not recognized"; }
	    else { $formaterr = $formaterr." '$comptrees' not recognized"; }
	}
    }
    
    if ($formaterr ne ""){
	print "\n";
	die ("     ERROR:\n          $formaterr.\n\n");
    }
    else{
	print "OK ($fileformat_line)\n";
    }
    # =====================================
    # =====================================
    
    # =====================================
    # Reading trees
    # =====================================
    my @rtrs = ();
    my @ctrs = ();
    my @auxtrees = (@reftrees, @comptrees);
    
    for (my $t=0; $t < scalar(@auxtrees); $t++){
	my @trees = ();
	my $line = "";
	chomp($auxtrees[$t]);
	next if ($auxtrees[$t] =~ m/^[\n|\/|\#|\>]/);
	$auxtrees[$t] =~ s /^\[.+\]//g;
	if ($auxtrees[$t] =~ m/\;$/){
	    $line = $line.$auxtrees[$t];
	    if ($t < scalar (@reftrees)){ push (@rtrs, $line); }
	    else { push (@ctrs, $line); }
	    $line = "";
	}
	else{
	    $line = $line.$auxtrees[$t];
	}
    }
    # =====================================
    # =====================================
    
    # =====================================
    # Check if there is more than one reference tree
    # =====================================
    print "Number of reference trees...  ";
    
    if (scalar(@rtrs) > 1){
	$WARNING = "\n     WARNING:  There is more than one tree in the reference file. Only the first one will be used.";
	print "$WARNING\n";
    }
    else{
	print "OK\n";
    }
    my $t0 = shift(@rtrs); #ÊThe first reference tree will be used
    
    # =====================================
    # =====================================
    
    # =====================================
    # Check branch lengths
    # =====================================
    print "Checking branch lengths...  ";
    my @trs = ($t0, @ctrs);
    my $nobrlen = "no";
    my $wrongbrlen = "no";	
    my $indexes = "";
    my $indexes_wrbrlen = "";
    my $n_wo_brlen = 0;
    my $n_wr_wo_brlen = 0;	
    my $refbr = "yes";
    my $wrrefbr = "yes";
    for (my $i=0; $i < scalar(@trs); $i++){
	if ($i == 0){ # It's reftree
	    if (!grep (/\:/, $trs[$i])){
		$refbr = "no";
		$nobrlen = "yes";
	    }
	}
	else{ # It's comptrees
	    if (!grep (/\:/, $trs[$i])){
		if ($nexus =~ m /yes$/){
		    $indexes = $indexes.$nametrees{$comptrees}{$i-1}.", ";
		}
		else{
		    $indexes = $indexes.($i).", ";
		}
		$nobrlen = "yes";
		$n_wo_brlen++;
	    }
	}
	
	#Checks branch lengths are right (not negative)
	my @aux = split (/\:|\)|\,/, $trs[$i]);
	foreach my $a (@aux){
	    if ($a=~ m/^\-[0-9]+/){
		$wrongbrlen = "yes";
		if ($i == 0){					
		    $wrrefbr = "no";
		}
		else{
		    if ($nexus =~ m /yes$/){
			$indexes_wrbrlen = $indexes_wrbrlen.$nametrees{$comptrees}{$i-1}.", ";
		    }
		    else{
			$indexes_wrbrlen = $indexes_wrbrlen.($i).", ";
		    }
		    $n_wr_wo_brlen++;
		}
	    }
	}		
    }
    $indexes =~ s/, $//g;
    $indexes_wrbrlen =~ s/, $//g;
    
    if ($nobrlen eq "yes"){
	my $brlenerr = "";
	if ($refbr eq "no") {
	    if ($indexes eq ""){
		$brlenerr = "          Tree '$reftree' has no branch lengths.\n";
	    }
	    else{
		if ($n_wo_brlen > 1){
		    $brlenerr = "          Tree '$reftree' has no branch lengths\n          &\n          Trees '$indexes' in file '$comptrees' have no branch lengths.\n";
		}
		else{
		    $brlenerr = "          Tree '$reftree' has no branch lengths\n          &\n          Tree '$indexes' in file '$comptrees' has no branch lengths.\n";
		}
	    }
	}
	else{
	    if ($n_wo_brlen > 1){
		$brlenerr = "          Trees '$indexes' in file '$comptrees' have no branch lengths.\n";
	    }
	    else{
		$brlenerr = "          Tree '$indexes' in file '$comptrees' has no branch lengths.\n";
	    }
	}
	
	$ERROR = $ERROR."\n     ERROR (branch lengths):\n$brlenerr\n";
	print "$ERROR";
	$ERR = 1;
    }
    else{
	
	$ERROR = "$ERROR\n\n" if ($ERROR ne "");
	my $brlenerr = "";	
	if ($wrongbrlen eq "yes"){
	    if ($wrrefbr eq "no"){
		if ($indexes_wrbrlen eq ""){
		    $brlenerr = "          Branch lengths of tree '$reftree' are incorrect (negative values).\n";
		}
		else{
		    if ($n_wo_brlen > 1){
			$brlenerr = "          Branch lengths of trees '$reftree' are incorrect (negative values)\n          &\n          Branch lengths of trees '$indexes_wrbrlen' in file '$comptrees' are incorrect (negative values).\n";
		    }
		    else{
			$brlenerr = "          Branch lengths of tree '$reftree' are incorrect (negative values)\n          &\n          Branch lengths of tree '$indexes_wrbrlen' in file '$comptrees' are incorrect (negative values).\n";
		    }
		}
	    }
	    else{
		if ($n_wr_wo_brlen > 1){
		    $brlenerr = "          Branch lengths of trees '$indexes_wrbrlen' in file '$comptrees' are incorrect (negative values).\n";
		}
		else{
		    $brlenerr = "          Branch lengths of tree '$indexes_wrbrlen' in file '$comptrees' are incorrect (negative values).\n";
		}			
	    }
	    if ($nobrlen eq "yes"){
		$ERROR = $ERROR."\n\n$brlenerr\n";
	    }
	    else{
		$ERROR = "\n     ERROR (branch lengths):\n$brlenerr\n";			
	    }
	    print "$ERROR";
	    $ERR = 1;			
	}
	else{
	    print "OK\n";
	}
    }	
    $ERROR = "";
    # =====================================
    # =====================================
    
    # =====================================
    # Check species - tips
    # =====================================
    print "Tips in trees...  ";
    my @species = ();
    foreach my $t (@trs){
	my @sps = ();
	my @auxsps = ();
	if (grep (/\:/, $t)){
	    @auxsps = split (/ *\: *\-*[0-9]+\.*[0-9]*e*E*-*[0-9]* *| *\, *| *\( *| *\) *[0-9]*\.*[0-9]* *| *\; */, $t);
	}
	else{
	    @auxsps = split (/ *\, *| *\( *| *\) *| *\; */, $t);
	}
	foreach my $s (@auxsps){
	    my $auxs = $s;
	    $auxs =~ s/\'|\"//g;
	    $auxs =~ s/ /_/g;
	    $t =~ s/$s/$auxs/g;
	    push (@sps, $auxs) if ($auxs ne "");
	}
	push (@species, [@sps]);
	$t =~ s/ //g;
    }
    $t0 = $trs[0]; # Update reference tree
    @ctrs = @trs; # Update comparison trees
    shift(@ctrs);
    my @dif_n_sps = ();
    my @dif_sps = ();
    for (my $i=1; $i < scalar(@species); $i++){
	my $hits = 0;
	if (scalar(@{$species[0]}) != scalar (@{$species[$i]})){
	    push (@dif_n_sps, $i);
	}
	else{
	    for (my $j=0; $j < scalar(@{$species[0]}); $j++){
		if (grep (/^$species[0][$j]$/, @{$species[$i]})){
		    $hits++;
		}
	    }
	    if ($hits != scalar(@{$species[0]})){
		push (@dif_sps, $i);
	    }
	}
    }
    
    my $spperr = "";
    $indexes = "";
    my $indexes_spp = "";
    
    if (scalar (@dif_n_sps) > 0){
	foreach my $dns (@dif_n_sps){
	    if (exists ($nametrees{$comptrees}{($dns-1)})){
		$indexes = $indexes.$nametrees{$comptrees}{($dns-1)}.", ";
	    }
	    else{
		$indexes = $indexes.($dns).", ";
	    }
	    $indexes_spp = $indexes_spp.scalar(@{$species[$dns]}).", ";
	}
	$indexes =~ s/, $//g;
	$indexes_spp =~ s/, $//g;
	my @aux_ind = split (/, /, $indexes_spp);
	$indexes_spp = "";
	for (my $i=0; $i < $#aux_ind; $i++){
	    $indexes_spp = $indexes_spp.$aux_ind[$i].", ";
	}
	$indexes_spp =~ s/, $/ /g;
	$indexes_spp = $indexes_spp."and " if ($indexes_spp ne "");
	$indexes_spp = $indexes_spp."$aux_ind[$#aux_ind]";
	my $n_ind = scalar(@aux_ind);
	if ($n_ind > 1){
	    $spperr = $spperr."          Reference tree has ".scalar(@{$species[0]})." tips\n          &\n          trees '$indexes' in '$comptrees' have ".$indexes_spp." tips, respectively.";
	}
	else{
	    $spperr = $spperr."          Reference tree has ".scalar(@{$species[0]})." tips\n          &\n          tree '$indexes' in '$comptrees' has ".$indexes_spp." tips.";
	}
    }
    
    $indexes = "";
    if (scalar (@dif_sps) > 0){
	foreach my $ds (@dif_sps){
	    if (exists ($nametrees{$comptrees}{($ds-1)})){
		$indexes = $indexes.$nametrees{$comptrees}{($ds-1)}.", ";
	    }
	    else{
		$indexes = $indexes.($ds).", ";
	    }
	}
	$indexes =~ s/, $//g;
	my @aux_ind = split (/, /, $indexes);
	my $n_ind = scalar(@aux_ind);
	if ($n_ind > 1){
	    $spperr = $spperr."\n          &\n" if ($spperr ne "");
	    $spperr = $spperr."          The tips of trees '".$indexes."' in '$comptrees' are different to those of reference tree.";
	}
	else{
	    $spperr = $spperr."\n          &\n" if ($spperr ne "");
	    $spperr = $spperr."          The tips of tree '".$indexes."' in '$comptrees' are different to those of reference tree.";
	}
    }
    if ($spperr ne ""){
	$spperr = "     ERROR (species - tips):\n".$spperr;
	$ERROR = "$ERROR\n$spperr\n\n";
	print "$ERROR";
	$ERR = 1;
    }
    else{
	print "OK\n";
    }
    
    $ERROR = "";
    
    # =====================================
    # =====================================
    
    # =====================================
    # Check root
    # =====================================
    print "All trees rooted or all unrooted...  ";
    my $n_rooted = 0;
    my $n_trees = 0;
    my @treenumber = ();
    foreach my $t (@trs){
	my $aux = $t;
	$aux =~ s/\)[0-9]+\:[0-9]+\.*[0-9]*/\)/g;# for bootstrap
	$aux =~ s/^\(|\)\;$|\:[0-9]+\.*[0-9]*//g;
	my @aux2 = split (/(\(|\))/,$aux);
	my $opened = 0;
	my $closed = 0;
	my $chg = "";
	foreach my $a (@aux2){
	    next if ($a eq "\n" || $a eq "");
	    $opened ++ if ($a eq "\(");
	    
	    $closed ++ if ($a eq "\)");
	    
	    if ($opened != 0){
		$a = "\\(" if ($a eq "\(");
		$a = "\\)" if ($a eq "\)");
		$chg = $chg.$a;
	    }
	    if (($opened == $closed) && ($opened != 0) && ($chg ne $a)){
		$aux =~ s/$chg/X/g;
		$chg = "";
		$opened = 0;
		$closed = 0;
	    }
	}
	my @cm = split (/\,/, $aux);
	if (scalar (@cm) == 2) {
	    push (@treenumber, $n_trees."R");
	    $n_rooted++;
	}
	else {
	    push (@treenumber, $n_trees."U");
	}
	$n_trees++;
    }
    
    if (($n_rooted != scalar (@trs)) && ($n_rooted != 0)){
	my $about_trees = "";
	my $indexes = "";
	
	if ($treenumber[0] =~ /R$/){ # Reference tree is rooted
	    $rooted = "yes";
	    my $n_unrooted = 0;
	    foreach my $tn (@treenumber){
		if (grep (/U/,$tn)){
		    $tn =~ s/U//g;
		    if ($nexus =~ m/yes$/){
			$indexes = $indexes.$nametrees{$comptrees}{$tn-1}.", ";
		    }
		    else{
			$indexes = $indexes.$tn.", ";
		    }
		    $n_unrooted++;
		}
	    }
	    shift (@treenumber);
	    if ($n_unrooted == scalar(@treenumber)){
		if (scalar(@treenumber) == 1){
		    $about_trees = "          Reference tree '$reftree' is rooted\n          &\n          '$comptrees' is unrooted";
		}
		else{
		    $about_trees = "          Reference tree '$reftree' is rooted\n          &\n          all trees in '$comptrees' are unrooted";
		}
	    }
	    else{
		$indexes =~ s/\, $//g;
		my @aux = split(/,/,$indexes);
		my $n_ind = scalar(@aux);
		if ($n_ind > 1){
		    $about_trees = "          Reference tree '$reftree' is rooted\n          &\n          trees '$indexes' in '$comptrees' are unrooted";
		}
		else{
		    $about_trees = "          Reference tree '$reftree' is rooted\n          &\n          tree '$indexes' in '$comptrees' is unrooted";
		}
	    }
	}
	else {
	    $n_rooted = 0;
	    foreach my $tn (@treenumber){
		if (grep (/R/,$tn)){
		    $tn =~ s/R//g;
		    if ($nexus =~ m/yes$/){
			$indexes = $indexes.$nametrees{$comptrees}{$tn-1}.", ";
		    }
		    else{
			$indexes = $indexes.$tn.", ";
		    }
		    $n_rooted++;
		}
	    }
	    shift (@treenumber);
	    if ($n_rooted == scalar(@treenumber)){
		if (scalar(@treenumber) == 1){
		    $about_trees = "          Reference tree '$reftree' is unrooted tree \n          &\n          '$comptrees' is rooted";
		}
		else{
		    $about_trees = "          Reference tree '$reftree' is unrooted \n          &\n          all trees in '$comptrees' are rooted" ;
		}
	    }
	    else{
		$indexes =~ s/\, $//g;
		my @aux = split(/,/,$indexes);
		my $n_ind = scalar(@aux);
		if ($n_ind > 1){
		    $about_trees = "          Reference tree '$reftree' is unrooted \n          &\n          trees '$indexes' in '$comptrees' are rooted";
		}
		else{
		    $about_trees = "          Reference tree '$reftree' is unrooted \n          &\n          tree '$indexes' in '$comptrees' is rooted";
		}
	    }
	}
	if ($ERROR ne ""){
	    $ERROR = "$ERROR\n\nERROR (rooted - unrooted):\n$about_trees.\n\n     Trees should be all rooted or all unrooted.\n\n";
	}
	else{
	    $ERROR = "\n     ERROR (rooted - unrooted):\n$about_trees.\n\n          Trees should be all rooted or all unrooted.\n\n";
	}
	print "$ERROR";
	$ERR = 1;
    }
    else{
	if ($treenumber[0] =~ /R$/){
	    $rooted = "yes";
	}
	print "OK\n";
    }
    $ERROR = "";
    # =====================================
    # =====================================
    
    
    my $errormsg = "\n".'*********************************************'."\n".'  There were one or more ERRORS. See above.  '."\n".'*********************************************'."\n\n";
    
    die ("$errormsg") if ($ERR == 1);

    return (\$reftree, \$comptrees, \%nametrees, \$taxblock, \$rooted, \$bk_line, \$t0, \@ctrs);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Check kind of line break (and change it to UNIX if necessary) ///
sub break_line{
    my $file = shift(@_);
    my $break = "";
    my $br = "\n";
    open (FILE, "$file") or die ("Cannot open file $file");
    my @file = (<FILE>);
    close(FILE);
    foreach (@file){
	if (/\r\n/){ $break = 'DOS'; $br = "\r\n"; last; }
	if (/\r/)  { $break = 'MAC'; $br = "\r"; last; }
	if (/\n/)  { $break = 'UNIX'; last; }
    }
    my @lines = ();
    if ($break ne 'UNIX'){
	my $all_lines = "";
	foreach (@file){
	    s/$br/\n/g;
	    $all_lines = $all_lines.$_;
	}
	@lines = split (/\n/, $all_lines);
    }
    else{
	@lines = @file;
    }
    return($break, @lines);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Converts nexus tree to newick ///
sub nexus2newick{
	my ($a, $b) = @_;
	my @aux = (@$a);
	my @newicktrees = ();
	my %nametrees = ();
	my $taxblock = "no";
	my $treeblock = "no";
	if (grep (/translate/i, @aux)){ # Nexus with taxa block
		$taxblock = "yes";
		my %taxa = ();
		my $n = 0;
		my $continue = "no";
		foreach (@aux){

			if (/^Begin +TREES/i){
				$treeblock = "yes";
			}
			next if ($treeblock eq "no");
			last if (/End\;/ && $treeblock eq "yes");

			$continue = "yes" if (/Translate/i);
			next if ($continue eq "no");
			
			if (/^( *|\t*)[0-9]+/){
				s/^( +|\t+)//g;
				my @aux2 = split(/^([0-9]+)/, $_);
				my $number = $aux2[1];
				my $taxon = $aux2[2];
				$taxon =~ s/\t+| +|\,|\;|\n//g;
				$taxa{"$number"} = $taxon;
			}
			
			if (/^ *\t*[U|R]*TREE/i){
				my $tree = $_;
				my $nametree = $tree;
				$nametree =~ s/^ *\t*[U|R]*TREE *\t*|\[\&r*u*\]| *\t*(\[.+\])* *\t*\=.*|\'|\n//ig;
				$tree =~ s/.+=[ |\t|\[|\]|\&|a-z]*\(/\(/ig;
				if (grep (/\:/,$tree)){
					my @auxtree = split (/([0-9]+\:)/, $tree);
					$tree = "";
					foreach my $a (@auxtree){
						if ($a =~ m/^[0-9]+\:/){
							$a =~ s/\://g;
							$a = $taxa{$a};
							$a = $a.":";
						}
						$tree = $tree.$a;
					}
				}
				else{# No branch lengths
					my @auxtree = split (/([0-9]+)/, $tree);
					$tree = "";
					foreach my $a (@auxtree){
						$a = $taxa{$a} if ($a =~ m/^[0-9]+/);
						$tree = $tree.$a;
					}
				}
				push (@newicktrees, $tree);
				$nametrees{$n} = $nametree;
				$n++;
			}
		}
	}
	else{ # Nexus without taxa block
		my $n = 0;
		foreach (@aux){
			if (/^Begin +TREES/i){
				$treeblock = "yes";
			}
			next if ($treeblock eq "no");
			last if (/End\;/ && $treeblock eq "yes");

			if (/^ *\t*[U|R]*TREE/i){
				my $tree = $_;
				my $nametree = $tree;
				$nametree =~ s/^ *\t*[U|R]*TREE *\t*|\[\&r*u*\]| *\t*(\[.+\])* *\t*\=.*|\'|\n//ig;
				$tree =~ s/.+=[ |\t|\[|\]|\&|a-z]*\(/\(/ig;
				push (@newicktrees, $tree);
				$nametrees{$n} = $nametree;
				$n++;
			}
		}
	}

	return (\%nametrees, \$taxblock, \@newicktrees);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Converts a newick file to nexus ///
sub newick2nexus{
	my ($sc, $txbl, $nt, $cpt) = @_;
	my $scaled = ($$sc);
	my $taxblock = ($$txbl);
	my %nametrees = (%$nt);
	my $comptrees = ($$cpt);
	open (FILE, "$scaled") or die ("Cannot open $scaled");
		my @auxfile = (<FILE>);
	close (FILE);
	open (FILE, ">$scaled") or die ("Cannot write to $scaled");
		print FILE "#NEXUS\n\n";
		print FILE "\[\!\n";
		print FILE "File generated by Ktreedist version $version\n\n";
		print FILE "Reference file: '$reftree'\n";
		print FILE "Comparison file: '$comptrees'\n";
		print FILE "]\n\n";
		print FILE "Begin TREES;\n\n";
		my $ntips = 1;
		my %nametips = ();
		if ($taxblock eq "yes"){
			print FILE "\tTranslate\n";
			my @aux = ();
			foreach my $af (@auxfile){
				if ($af =~ m/^\(/){
					@aux = split (/\,|\(|\)|\;/, $af);
					last;
				}
			}
			my $ind = 1;
			my $taxblock_output = "";
			foreach my $a (@aux){
				if ($a =~ m/^.+\:/){
					$a =~ s/\:.+//g;
					$taxblock_output = $taxblock_output."\t\t$ind\t'$a',\n";
					$nametips{$ind} = $a;
					$ind++;
				}
			}
			$ntips = $ind;
			$taxblock_output =~ s/,\n$/\n\t;\n\n/g;
			print FILE $taxblock_output;
		}
		for (my $i=0; $i < scalar(@auxfile); $i++){
			if ($auxfile[$i] =~ m/^\(/){
				if ($taxblock eq "yes"){
					for (my $j=1; $j < $ntips; $j++){
						$auxfile[$i] =~ s/$nametips{$j}/$j/g;
					}
				}
				print FILE "TREE '$nametrees{$comptrees}{$i}' = $auxfile[$i]";
			}
			else{
				print FILE $auxfile[$i];
			}
		}
		print FILE "\nEnd;\n";
	close (FILE);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# /// Get partitions ///
sub get_partitions_tree{
	my $tree = shift(@_);
	my $separated = shift(@_);
	$tree =~ s/\) *[0-9]*\.*[0-9]*/\)/g;
	my @aux = split(/(:|,|\)|\;|\()/, $tree);
	my @treesplit = ();
	foreach my $a (@aux){
		push (@treesplit, $a) if ($a ne "");
	}

	my %brlen = ();
	my %aux_brlen = ();
	my @cladosA = ();
	my @cladosB = ();
	my @species = ();
	my @opened = ();
	my $n_comas = 0;
	my $cont = -1;
	my $rooted = "no";

	# Get internal partitions
	for (my $i=0; $i < scalar(@treesplit); $i++) {
	    $n_comas ++ if (($treesplit[$i] eq "\,") && ($#opened == 0));
	    if ($treesplit[$i] eq "\("){
		$cont++;
		push (@opened, $cont);
	    }
	    if ($treesplit[$i] eq "\:"){
		if ($treesplit[$i-1] =~ m/\w+/){
		    if ($treesplit[$i-2] ne "\)"){
                        # it's a tip and not a bootstrap value
			push (@species, $treesplit[$i-1]);
			$aux_brlen{$species[$#species]} = $treesplit[$i+1];
			foreach my $o (@opened){
			    push (@{$cladosA[$o]}, $treesplit[$i-1]);
			}
		    }
		}
	    }
	    if ($treesplit[$i] eq "\)"){
		if ($treesplit[$i+1] ne ";"){
		    $brlen{$cladosA[$opened[$#opened]]} = $treesplit[$i+2];
		}
		pop(@opened);
	    }
	}
	$rooted = "yes" if ($n_comas == 1);

	my $index = $#cladosA;
	shift(@cladosA); # First partition is all the tree

	#Get symmetric ones
	for (my $i=0; $i < scalar(@cladosA); $i++){
		foreach my $sp (@species){
			if (!grep(/^$sp$/, @{$cladosA[$i]})){
				push (@{$cladosB[$i]}, $sp);
			}
		}
		push (@{$cladosB[$i]}, "root") if ($rooted eq "yes");
		$brlen{$cladosB[$i]} = $brlen{$cladosA[$i]};
	}
	# Get tip partitions
	foreach my $sp (@species){
		$cladosA[$index][0] = $sp;
		$brlen{$cladosA[$#cladosA]} = $aux_brlen{$sp};
		foreach my $sp2 (@species){
			push (@{$cladosB[$index]}, $sp2) if ($sp2 ne $sp);
		}
		push (@{$cladosB[$index]}, "root") if ($rooted eq "yes");
		$brlen{$cladosB[$#cladosB]} = $aux_brlen{$sp};
		$index++;
	}

	my @clados = ();
	my %brlen_aux = ();

	# Join partitions when requested
	if ($separated eq "no"){
	    for (my $i=0; $i < scalar (@cladosA); $i++){
		push (@clados, [@{$cladosA[$i]}]);
		$brlen_aux{$clados[$#clados]} = $brlen{$cladosA[$i]};
		push (@clados, [@{$cladosB[$i]}]);
		$brlen_aux{$clados[$#clados]} = $brlen{$cladosB[$i]};
	    }
	    %brlen = %brlen_aux;
	    return (\@clados, \%brlen, \@species, \$rooted);
	}
	else{
	    return (\@cladosA, \@cladosB, \%brlen, \@species, \$rooted);
	}
}

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

sub Set_root{
    my $node= $_[0];
    my @Root;
    my $n_nodes=scalar(@$node);

    $Root[1]="Root";
    for(my $a=2; $a<=9; $a++){$Root[$a]=-1;}
    my $nson=0;
    my $all_sons=0;
    my @sons;
    for(my $i=0; $i<$n_nodes; $i++){
	if($node->[$i][2]<0){
	    $nson++;
	    for(my $k=0; $k<$node->[$i][0]; $k++){
		$sons[$k+$all_sons]=$node->[$i][4]->[$k];
	    }
	    $all_sons+=$node->[$i][0];
	    if($Root[7]<0){$Root[7]=$i;}
	    elsif($Root[8]<0){$Root[8]=$i;}
	    elsif($Root[9]<0){$Root[9]=$i;}
	    else{
		print "WARNING, >3 sons for Root: ",
		$Root[7]," ",$Root[8]," ",$Root[9]," $i\n";
	    }
	}

    }
    $Root[4]=\@sons;
    $Root[0]=$all_sons;
    print "Root has $nson son nodes and $all_sons descendants\n";
    return(@Root);
}

sub Put_sons{
    my $node= $_[0];
    my $n_nodes=scalar(@$node);
    for(my $i=0; $i<$n_nodes; $i++){
	$node->[$i][7]=-1;
	$node->[$i][8]=-1;
	$node->[$i][9]=-1;
    }
    my $triple=0; my $assigned=0; my $no_par=0;
    for(my $i=$n_nodes-1; $i>=0; $i--){
	my $p=$node->[$i][2];
	if($p>=0){
	    $assigned++;
	    if($node->[$p][7]<0){$node->[$p][7]=$i;}
	    elsif($node->[$p][8]<0){$node->[$p][8]=$i;}
	    elsif($node->[$p][9]<0){$node->[$p][9]=$i; $triple++;}
	    else{
		print "WARNING, >3 sons for node ${p}: ",
		$node->[$p][7]," ",$node->[$p][8]," ",$node->[$p][9]," $i\n";
		$triple++;
		#die;
	    }
	}else{
	    $no_par++;
	}
    }
    print "$assigned sons assigned out of $n_nodes nodes ($no_par without parent), ",
    "$triple nodes with 3 or more sons\n";
}

sub Write_classification{
    my ($node, $matrix, $z2, $file, $num_tree)= @_;

    my @nodes=@{$node};
    my $n_node=scalar(@nodes);
    my $nseq=($n_node+3)/2;
    my $n_diff=0;
    for (my $i=0; $i < $n_node; $i++){if($z2->[$i]>0.5){$n_diff++;}}
    print "$n_diff branches with different function and ",$n_node-$n_diff,
    " with same function\n";

    my @Root=Set_root($node);
    for(my $a=0; $a<=9; $a++){$nodes[$n_node][$a]=$Root[$a];}
    $z2->[$n_node]=1;
    $n_node++;

    my $nameout=Remove_ext($file, 1);
    if($num_tree>1){$nameout=sprintf("%s_%d", $nameout, $num_tree);}
    $nameout=sprintf("%s.func", $nameout);
    print "Printing functional classification in $nameout\n";
    open(my $fo, '>', $nameout);
    
    print $fo
	"# Functional classification of $nseq proteins in structure tree ",
	Remove_ext($file, 1), " sequence divergence from tree $num_tree\n";
    
    my $cluster=0;
    my $nprot;

    # Set names
    my @names; $nseq=0;
    for (my $i=0; $i < $n_node; $i++){
	if($nodes[$i][0]==1){$names[$nseq]=$nodes[$i][1]; $nseq++;}
    }
    print "$nseq sequences exist among $n_node nodes\n";
    my @printed; for(my $i=0; $i<$nseq; $i++){$printed[$i]=0;}

    if(0){
	print $fo "# z2: ";
	for (my $i=0; $i < $n_node; $i++){
	    print $fo " $nodes[$i][1] ",$z2->[$i];
	}
	print $fo "\n";
    }

    for(my $sons=1; $sons<=$nseq; $sons++){ # number of sons
	for (my $i=0; $i < $n_node; $i++){
	    if($nodes[$i][0]!=$sons){next;}
	    if($z2->[$i]<=0.5){next;} # no split
	    my $split=0; my @tmp;
	    for(my $k=0; $k<$sons; $k++){
		my $seq=$nodes[$i][4]->[$k];
		my $ik=label($seq, \@names, $nseq);
		if($ik<0){
		    print "ERROR, sequence $seq not found in tree\n"; die;
		}
		if($printed[$ik]==0){
		    $tmp[$split]=$seq; $split++;
		    $printed[$ik]=1;
		}
	    }
	    if($split==0){next;} # No change of function
	    #print "Node $nodes[$i][1] $nodes[$i][0] has different function\n";
	    print $fo "# cluster $cluster $nodes[$i][1] $split proteins\n";
	    for(my $k=0; $k<$split; $k++){
		print $fo $tmp[$k], "\n";
	    }
	    $cluster++;
	    $nprot+=$split;
	} # end print cluster
    } # end nodes
    my $out="# $nprot proteins over $nseq have been classified\n".
    	"# $cluster clusters printed in file $nameout\n";
    print $fo $out; print $out;
    close($fo);

    my $tmp="tmp.txt";
    open($fo, '>', $tmp);
    print $fo "# num_clusters= $cluster num_protein= $nprot\n" ;
    close($fo);
    `cat $nameout >> $tmp`;
    `mv -f $tmp $nameout`;
}

sub label{
    my ($name, $names, $n)=@_;
    for(my $i=0; $i<$n; $i++){
	if($name eq $names->[$i]){return($i);}
    }
    return(-1);
}

#/// Get matrix ////



sub get_matrix{
	my ($cl0A, $cl0B, $cl1A, $cl1B, $brl0, $brl1, $nod0, $nod1) = @_;
	my @clados0A = (@$cl0A);
	my @clados0B = (@$cl0B);
	my @clados1A = (@$cl1A);
	my @clados1B = (@$cl1B);
	my @clados1 = (@clados1A, @clados1B);
	my %brlen0 = (%$brl0);
	my %brlen1 = (%$brl1);
	my @nodes0 = (@$nod0);
	my @nodes1 = (@$nod1);

	my @matrix = ();
	my $m = 0;
	# Add partitions there are in tree 0
	for (my $i=0; $i < scalar(@clados0A); $i++){
	    my $priorm = $m;
	    for (my $j=0; $j < scalar(@clados1); $j++){
		if (scalar(@{$clados0A[$i]}) == scalar(@{$clados1[$j]})){
		    my $n_equal = 0;
		    foreach my $cl0 (@{$clados0A[$i]}){
			if (grep (/^$cl0$/, @{$clados1[$j]})){
			    $n_equal++;
			}
		    }
		    if ($n_equal == scalar (@{$clados0A[$i]})){
			if (($hard_polytomies eq "off") ||
			    ($hard_polytomies eq "on" && $brlen0{$clados0A[$i]} != 0 &&
			     $brlen1{$clados1[$j]} != 0)){
			    push (@{$matrix[$m]}, [@{$clados0A[$i]}]); # index 0, partition
			    push (@{$matrix[$m]}, [@{$clados0B[$i]}]); # index 1, symmetric
			    push (@{$matrix[$m]}, $brlen0{$clados0A[$i]}); # index 2, brlen
			    push (@{$matrix[$m]}, $brlen1{$clados1[$j]});  # index 3, brlen
			    push (@{$matrix[$m]}, $brlen0{$clados0A[$i]}); # "NONE" column
			    push (@{$matrix[$m]}, $brlen1{$clados1[$j]}); # "NONE" column
			    push (@{$matrix[$m]}, $nodes0[$i][1]); # name of node in ref.tree
			    push (@{$matrix[$m]}, $nodes0[$i][5]); # name of parent node in ref.tree
			    $nod0->[$i][6]=$m;
			    $m++;
			}
			last;
		    }
		}
	    }
	    if ($priorm == $m){
		if ($hard_polytomies eq "off" ||
		    ($hard_polytomies eq "on" && $brlen0{$clados0A[$i]} != 0)){
		    push (@{$matrix[$m]}, [@{$clados0A[$i]}]); 
		    push (@{$matrix[$m]}, [@{$clados0B[$i]}]);
		    push (@{$matrix[$m]}, $brlen0{$clados0A[$i]});
		    push (@{$matrix[$m]}, 0);
		    push (@{$matrix[$m]}, $brlen0{$clados0A[$i]}); # "NONE" column
		    push (@{$matrix[$m]}, "NONE"); # "NONE" column
		    push (@{$matrix[$m]}, $nodes0[$i][1]); # name of node in ref.tree
		    push (@{$matrix[$m]}, $nodes0[$i][5]); # name of parent node in ref.tree
		    $nod0->[$i][6]=$m;
		    $m++;
		}
	    }
	}

	# Add partitions there are in tree 1 but not in tree 0
	for (my $i=0; $i < scalar(@clados1A); $i++){
	    my $n_nonequal = 0;
	    for (my $j=0; $j < scalar(@matrix); $j++){
		my $n_equal = 0;
		for (my $k=0; $k < 2; $k++){
		    $n_equal = 0;
		    if (scalar(@{$clados1A[$i]}) == scalar (@{$matrix[$j][$k]})){
			foreach my $cl1 (@{$clados1A[$i]}){
			    if (grep (/^$cl1$/, @{$matrix[$j][$k]})){
				$n_equal++;
			    }
			}
			last if ($n_equal == scalar(@{$clados1A[$i]}));
		    }
		}
		if ($n_equal != scalar(@{$clados1A[$i]})){
		    $n_nonequal ++;
		}
		else{
		    last;
		}
	    }
	    if ($n_nonequal == scalar (@matrix)){
		if (($hard_polytomies eq "off") ||
		    ($hard_polytomies eq "on" && $brlen1{$clados1A[$i]} != 0)){
		    push (@{$matrix[$m]}, [@{$clados1A[$i]}]);
		    push (@{$matrix[$m]}, [@{$clados1B[$i]}]);
		    push (@{$matrix[$m]}, 0);
		    push (@{$matrix[$m]}, $brlen1{$clados1A[$i]});
		    push (@{$matrix[$m]}, "NONE"); # "NONE" column
		    push (@{$matrix[$m]}, $brlen1{$clados1A[$i]}); # "NONE" column
		    push (@{$matrix[$m]}, $nodes1[$i][1]); # name of node in ref.tree
		    push (@{$matrix[$m]}, $nodes1[$i][5]); # name of parent node in ref.tree
		    $nod1->[$i][6]=$m;
		    $m++;
		}
	    }
	}
	return (@matrix);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#Ê//// Calculate K-value ////

sub get_K_value{
    my @matrix = @_;
    
    # K-value calculation 2=reference tree 3= comparison tree
    my $numerator = 0;
    my $denominator = 0;
    my $sum_br=0;
    for (my $i=0; $i < scalar(@matrix); $i++){
	$numerator = $numerator + $matrix[$i][2]*$matrix[$i][3];
	$denominator = $denominator+($matrix[$i][3]**2);
	$sum_br += $matrix[$i][3];
    }
    my $kval = $numerator/$denominator;
    return ($kval, $sum_br);
}

sub get_K_value_EM{ #\@matrix, $filename, $ntree, $tot_trees
    my ($ma) = @_;
    my @matrix = (@$ma);
    my $n=scalar(@matrix);

    my @x; # ratio of branch lengths (0: ref 1: comp)
    my @wb; # weights, given by b'_i^2
    #my @w;  # fake

    my @k2 = (); $k2[0]=0; $k2[1]=0;
    my @z_out= (); 

    # Ratio of branch lengths
    my $norm=0;
    for (my $i=0; $i <$n;  $i++){
	if($matrix[$i][3]>0.00001){
	   $x[$i]=$matrix[$i][2]/$matrix[$i][3];
	   $wb[$i]=$matrix[$i][3]*$matrix[$i][3];
	}else{
	   $x[$i]=9.9999; $wb[$i]=0;
	}
	$norm+=$wb[$i];
    }
    $norm/=$n;
    for (my $i=0; $i <$n;  $i++){$wb[$i]/=$norm;}

    my $x1=0; my $x2=0; $norm=0;
    #for (my $i=0; $i<$n; $i++){$w[$i]=1;}
    for (my $i=0; $i <$n;  $i++){
	$x1+=$x[$i]*$wb[$i];
	$x2+=$x[$i]*$x[$i]*$wb[$i];
	$norm+=$wb[$i];
    }
    $x1/=$norm; # norm=n
    $x2=$x2/$norm-$x1*$x1;

    # Output file
    my $nameout=Remove_ext($_[1], 1); my $fo;
    {
	my $ntree=$_[2];
	my $tot_trees=$_[3];
	if($tot_trees>1){$nameout=sprintf("%s.%d", $nameout, $ntree+1);}
	$nameout=sprintf("%s.node", $nameout);
	print "Computing expectation maximization for tree ", $ntree+1,
	" with $n branches\n";
	if($write_fit){
	    print "Writing results in file $nameout\n";
	    open($fo, '>', $nameout);
	}
    }
    
    # Mean and coefficient of variation of branch lengths
    my $diff2=0;
    my $b0_1=0; my $b0_2=0; my $n0=0;
    my $b1_1=0; my $b1_2=0; my $n1=0;
    for (my $i=0; $i<$n;  $i++){
	my $d;
	my $b=$matrix[$i][2]; # branch of ref.tree
	if($b){$b0_1+=$b; $b0_2+=$b*$b; $n0++; $d=$b;}
	$b=$matrix[$i][3];    # branch of comp.tree
	if($b){$b1_1+=$b; $b1_2+=$b*$b; $n1++; $d-=$b;}
	$diff2+=$d*$d;
    }
    if($n0){$b0_1/=$n0;} $b0_2=($b0_2-$n0*$b0_1*$b0_1);
    if($n0>1){$b0_2=sqrt($b0_2/($n0-1));}
    my $out=sprintf("# Ref. tree: %d branches, mean: %.4g CoV: %.3g\n",
		    $n0, $b0_1, $b0_2/$b0_1);
    if($write_fit){print $fo $out;} #print $out; 
    if($n0){
	$diff2=sqrt($diff2/$n0)/$b0_1;
	if($diff2< 0.0001){
	    $out="# Ref. tree and comp. tree have RMSD= $diff2 <b>, exiting\n";
	    if($write_fit){print $fo $out;} print $out;
	    return(\@k2, \@z_out);
	}
    }

    if($n1){$b1_1/=$n1;} $b1_2=($b1_2-$n1*$b1_1*$b1_1);
    if($n1>1){$b1_2=sqrt($b1_2/($n1-1));}
    $out=sprintf("# Comp tree: %d branches, mean: %.4g CoV: %.3g\n",
		 $n1, $b1_1, $b1_2/$b1_1);
    if($write_fit){print $fo $out;} #print $out; 
    
    # Expectation Maximization, single branch type
    my $par1=2;
    my ($log_lik_1, $zz1, $t1, $m1, $s1)=EM_d1(\@x, \@wb, $n, $x1, $x2, $fo, 1);
    # t1: a-priori probability of each type of branch
    # m1: mean value of x for branches of type j
    # s1: variance of x for branches of type j
    my $AIC_1=2*($par1-$log_lik_1); # 2: AIC
    my $AICc_1=$AIC_1+(2*$par1*$par1+2*$par1)/($n-$par1-1); # 3: AICc
    my $BIC_1=log($n)*$par1-2*$log_lik_1; # 4: BIC
    my @mu1 = (@$m1);
    my @sig1 = (@$s1);
    my @z1 = (@$zz1);
    $out="# One type of branch $par1 par:";
    $out=sprintf("%s mean b/b'= %.3g", $out, $mu1[0]);
    $out=sprintf("%s s.dev.= %.3g\n", $out, $sig1[0]);
    $out=sprintf("%s# Log_lik= %.4g AIC= %.4g AICc= %.4g BIC= %.4g",
		 $out, $log_lik_1, $AIC_1, $AICc_1, $BIC_1);
    if($write_fit){print $fo $out, "\n";}
    
    # Expectation Maximization, two branch types
    my $par2=5;
    my ($log_lik_2, $zz2, $t2, $m2, $s2)=EM_d1(\@x, \@wb, $n, $x1, $x2, $fo, 2);
    my $AIC_2=2*($par2-$log_lik_2); # 2: AIC
    my $AICc_2=$AIC_2+(2*$par2*$par2+2*$par2)/($n-$par2-1); # 3: AICc
    my $BIC_2=log($n)*$par2-2*$log_lik_2; # 4: BIC
    my @mu2 = (@$m2);
    my @sig2 = (@$s2);
    my @tau2 = (@$t2);
    my @z2 = (@$zz2);
    $out="# Two types of branch $par2 par:";
    $out=sprintf("%s mean b/b'= %.3g %.3g", $out, $mu2[0], $mu2[1]);
    $out=sprintf("%s s.dev.= %.3g %.3g ", $out, $sig2[0], $sig2[1]);
    $out=sprintf("%s tau0= %.3g\n", $out, $$t2[0]);
    $out=sprintf("%s# Log_lik= %.4g AIC= %.4g AICc= %.4g BIC= %.4g",
		 $out, $log_lik_2, $AIC_2, $AICc_2, $BIC_2);
    if($write_fit){print $fo $out, "\n";}
    
    my @criteria; my $n_not=0;
    if($AIC_2  > $AIC_1){push (@criteria, "AIC"); $n_not++;}
    if($AICc_2 > $AICc_1){push (@criteria, "AICc"); $n_not++;}
    if($BIC_2  > $BIC_1){push (@criteria, "BIC"); $n_not++;}
    $supporting=3-$n_not;
    print "2 types of branches are supported by $supporting criteria\n";
    if($write_fit){
	print $fo "# Two types of branches ";
	if($n_not){
	    print $fo "not supported by";
	    for(my $i=0; $i<$n_not; $i++){print $fo " ",$criteria[$i];}
	    print $fo "\n";
	}else{
	    print $fo "supported by all criteria (AIC, AICc, BIC)\n";
	}
    }

    my @norm2 = (); $norm2[0]=0; $norm2[1]=0; my $j;
    
    for(my $i=0; $i <$n; $i++){
	if($wb[$i]==0){$z_out[$i]=-1; next;}
	if($z2[$i]>0.5){$j=0;}
	else{$j=1;}
	$z_out[$i]=$j; $k2[$j]+= $wb[$i]*$x[$i]; $norm2[$j]+=$wb[$i];
    }
    if($norm2[0]){$k2[0]/=$norm2[0]};
    if($norm2[1]){$k2[1]/=$norm2[1]}; 
    $out="# Scale_factors: ".sprintf(" %.4f %.4f\n", $k2[0], $k2[1]);

    if($write_fit){
	print $fo $out;
	print $fo "# The branch presents structural acceleration if ".
	    "P(type 0)<0.5\n";
	print $fo
	    "# node\tparent\tb_ref\tb_fit\tb_comp\tb_r/b_c\tProb(type 0)\n";
    }
    for (my $i=0; $i <$n; $i++){
	my $b_fit=$matrix[$i][3];
	my $j=$z_out[$i];
        if($j>=0){$b_fit*=$k2[$j];}
	if($write_fit){
	    my $out=sprintf("%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.3f\n",
			    $matrix[$i][6], $matrix[$i][7], $matrix[$i][2],
			    $b_fit, $matrix[$i][3], $x[$i], $z2[$i]);
	    print $fo $out;
	}
    }
    # print "Results written in $nameout\n";
    if($write_fit){close($fo);}

    return (\@k2, \@z_out);
}


#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#//// Calculate K-score ////
sub get_K_score{
	my ($ma, $kv, $nsp0) = @_;
	my @matrix = (@$ma);
	my $K_value = ($$kv);
	my $BLD = 0;
	my $BS	= 0;
	my $drf = 0;

	my $b_tot=0; my $nbr=0;
	for (my $i=0; $i < scalar (@matrix); $i++){
	    $matrix[$i][8] = $matrix[$i][3]*$K_value;
	    $BLD += ($matrix[$i][2] - $matrix[$i][8])**2;
	    $BS  += ($matrix[$i][2] - $matrix[$i][3])**2;
	    if($matrix[$i][2]){$nbr++;}
	    if($normalize==2){
		$b_tot += $matrix[$i][2]*$matrix[$i][2];
	    }elsif($normalize==1){
		$b_tot += $matrix[$i][2];
	    }

	    if ($matrix[$i][4] eq "NONE" || $matrix[$i][5] eq "NONE"){
		$drf++;
		$matrix[$i][9] = "NONE" if ($matrix[$i][5] eq "NONE");
	    }
	}
	if($normalize_RF){
	    $drf /= (2*$nsp0 -6);
	    $drf = sprintf("%.4f", $drf);
	}

	$BLD = sqrt($BLD/$nbr);
	$BS = sqrt($BS/$nbr);
	if($normalize){
	    $b_tot/=$nbr;
	    if($normalize==2){$b_tot=sqrt($b_tot);}
	    $BLD /= $b_tot;
	    $BS/= $b_tot;
	}
	return (\$BLD, \$drf, \$BS, \@matrix);
}

sub get_k2_score{
	my ($ma, $kk, $zz) = @_;
	my @matrix = (@$ma);
	my @z2 = (@$zz);
	my @k2 = (@$kk);

	my $K2s=0;
	my $b_tot=0;
	my $nbr=0;
	for (my $i=0; $i < scalar (@matrix); $i++){
	    my $j=$z2[$i]; my $k;
	    if($j >= 0){$k=$k2[$j];}
            else{$k=0;}
	    $K2s += ($matrix[$i][2] - $k*$matrix[$i][3])**2;
	    if($normalize==2){$b_tot += $matrix[$i][2]*$matrix[$i][2];}
	    elsif($normalize==1){$b_tot += $matrix[$i][2];}

	    if($matrix[$i][2]){$nbr++;}
	}
	$K2s = sqrt($K2s/$nbr);
	if($normalize){
	    $b_tot/=$nbr;
	    if($normalize==2){$b_tot=sqrt($b_tot);}
	    $K2s /= $b_tot;
	}
	return ($K2s);
}

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#Ê//// Scale the tree by K-value ////

sub modify_brlen_tree{
	my ($fi, $tr, $brl, $Ka) = @_;
	my $file = ($$fi);
	my $tree = ($$tr);
	my %brlen = (%$brl);
	my $K = ($$Ka);

	my $scaled_tree = "";
	for (values %brlen){
		$_ = sprintf ("%.".$precision."f",($_*$K));
	}

	my @aux = split(/(:|,|\)|\;|\()/, $tree);
	my @treesplit = ();
	foreach my $a (@aux){
		push (@treesplit, $a) if ($a ne "");
	}
	$tree = "";
	for (my $i=0; $i < scalar(@treesplit); $i++) {
		$treesplit[$i+1] = sprintf ("%.".$precision."f",($treesplit[$i+1]*$K)) if ($treesplit[$i] eq "\:");
		$tree = $tree.$treesplit[$i];
	}

	$scaled_tree = $tree if ($file ne "no");  # output scaled tree

	return (\$tree, \%brlen, \$scaled_tree);
}

sub modify_brlen_tree_2{
	my ($fi, $tr, $kk, $zz) = @_; #$brl, 
	my $file = ($$fi);
	my $tree = ($$tr);
	#my %brlen = (%$brl);
	my @k2 = (@$kk);
	my @z2 = (@$zz);

	my $scaled_tree = ""; my $l=0;
	#for (values %brlen){
	#    my $K; if($z2[$l]==0){$K=$k2[0];}else{$K=$k2[1];} $l++;
	#    $_ = sprintf ("%.".$precision."f",($_*$K));
	#}

	my @aux = split(/(:|,|\)|\;|\()/, $tree);
	my @treesplit = ();
	foreach my $a (@aux){
		push (@treesplit, $a) if ($a ne "");
	}
	$tree = ""; $l=0;
	for (my $i=0; $i < scalar(@treesplit); $i++) {
	    if ($treesplit[$i] eq "\:"){
		my $K; if($z2[$l]==0){$K=$k2[0];}else{$K=$k2[1];} $l++;
		$treesplit[$i+1] = sprintf ("%.".$precision."f",($treesplit[$i+1]*$K));
	    }
	    $tree = $tree.$treesplit[$i];
	}
	$scaled_tree = $tree if ($file ne "no");  # output scaled tree
	return (\$tree, \$scaled_tree); #\%brlen, 
}


#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# Change line break from UNIX to DOS or MAC
sub unix2macdos{
	my $kind = shift(@_);
	my $file = shift(@_);
	my $bl = "\n";
	if ($kind eq "MAC"){ $bl = "\r"; }
	elsif ($kind eq "DOS"){	$bl = "\r\n"; }
	open (FILE, "$file") or die ("Cannot open $file");
		my @aux=(<FILE>);
	close(FILE);
	open (FILE, ">$file") or die ("Cannot write to $file");
	foreach (@aux){
		s/\n/$bl/g;
		print FILE $_;
	}
	close(FILE);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#Ê/// Generate lines for display output ///
sub output_display{
	my ($ct, $nt, $namt, $bs, $ks, $kv, $k2s, $k2v, $r, $dr,
	    $nc1, $sum) = @_;
	my $comptrees = ($$ct);
	my $ntree = ($$nt);
	my %nametrees = (%$namt);
	my $BS = ($$bs);
	my $K_score = ($$ks);
	my $K_value = ($$kv);		
	my $K2_score = ($$k2s);
	my @K2 = (@$k2v);		
	my $robinson = ($$r);		
	my $drf = ($$dr);		
	my $n_clados1 = ($$nc1);
	my $sum_br = ($$sum);

	# Comparison tree name	
	my $name = "";
	my $output = "";
	if (exists($nametrees{$comptrees}{$ntree})){ $name = $nametrees{$comptrees}{$ntree}; }
	else{ $name = "Tree ".($ntree+1); }
	$name = substr($name, 0, 15) if (length($name) > 15);
	my $pre_name_bl = "";
	my $post_name_bl = "";
	for (my $i=0; $i < int((15-length($name))/2); $i++){ $post_name_bl = $post_name_bl." "; }
	my $check = ((15-length($name))%2);
	if ($check != 0){
		for (my $i=0; $i < int(((15-length($name))/2)+1); $i++){ $pre_name_bl = $pre_name_bl." "; }
	}
	else{
		$pre_name_bl = $post_name_bl;
	}
	$output = $output."$name$pre_name_bl$post_name_bl";

	# K-score & K-value columns
	$output = $output.
	    sprintf("   %.4f     %.4f    %.4f ", $BS, $K_score, $K_value);
	$output = $output.
	    sprintf("   %.4f     %.4f    %.4f ", $K2_score, $K2[0], $K2[1]);

	# Symmetric difference score column
	if ($robinson eq "yes"){
		my $pre_rf_bl = "";
		my $post_rf_bl = "";
		for (my $i=0; $i < int((15-length($drf))/2); $i++){
		    $post_rf_bl = $post_rf_bl." ";
		}
		my $check = ((15-length($drf))%2);
		if ($check != 0){
			for (my $i=0; $i < int(((15-length($drf))/2)+1); $i++){
			    $pre_rf_bl = $pre_rf_bl." ";
			}
		}
		else{
		    $pre_rf_bl = $post_rf_bl;
		}
		$output = $output.sprintf("  $pre_rf_bl%.4g$post_rf_bl", $drf);
	}

	# Number of partitions column
	if ($npartitions eq "yes"){
		my $parts = "$n_clados1";
		my $pre_parts_bl = "";
		my $post_parts_bl = "";
		for (my $i=0; $i < int((12-length($parts))/2); $i++){
		    $post_parts_bl = $post_parts_bl." ";
		}
		$check = ((12-length($parts))%2);
		if ($check != 0){
			for (my $i=0; $i<int(((12-length($parts))/2)+1); $i++){
			    $pre_parts_bl = $pre_parts_bl." "; 
			}
		}
		else{
		    $pre_parts_bl = $post_parts_bl;
		}
		$output = $output."  $pre_parts_bl$parts$post_parts_bl";
	}

        # Supporting criteria
	$output = $output.sprintf("  %.3g",$supporting);
	$output = $output.sprintf("  %.4g",$sum_br);
	$output = $output."  \n";		
	
	return ($output);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#Ê/// Generate lines for table file output ///
sub get_table {
	my ($rt, $ct, $namt, $fn, $nt, $np, $bs, $ks, $kv, $k2s, $k2v,
	    $r, $dr, $nc0, $nc1, $sum) = @_;
	my $reftree = ($$rt);
	my $comptrees = ($$ct);
	my %nametrees = (%$namt);
	my $ntree = ($$nt);		
	my $filename_comptrees_woext = ($$fn);	
	my $npartitions = ($$np);
	my $BS = ($$bs);		
	my $K_score = ($$ks);
	my $K_value = ($$kv);
	my $K2_score =($$k2s);
	my @K2 = (@$k2v);
	my $robinson = ($$r);
	my $drf = ($$dr);
	my $n_clados0 = ($$nc0);
	my $n_clados1 = ($$nc1);
	my $sum_br = ($$sum);

	my $linea = "";
	if (exists ($nametrees{$reftree}{0})){
		$linea = "$nametrees{$reftree}{0}";
# 		$linea = $filename_reftree_woext."."."$nametrees{$reftree}{0} (ref)"; # For nexus, in case you want the filename to appear in the table
	}
	else{
	    $linea = "$filename_reftree_woext";
	}
	if (exists ($nametrees{$comptrees}{$ntree})){
		$linea = $linea." - ".$nametrees{$comptrees}{$ntree};
#		$linea = $linea." - ".$filename_comptrees_woext.".".$nametrees{$comptrees}{$ntree}; # For nexus, in case you want the filename to appear in the table
	}
	else{
		$linea = $linea." - ".$filename_comptrees_woext;
#." ("."Tree ".($ntree+1).")";
	}
	$linea = $linea.sprintf("\t%.3g\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
				$BS, $K_score, $K_value,
				$K2_score, $K2[0], $K2[1]);
	$linea = $linea.sprintf("\t%.4f",$drf) if ($robinson eq "yes");
	#$linea = $linea."\t$n_clados0\t$n_clados1" if ($npartitions eq "yes");
	$linea = $linea."\t$n_clados0" if ($npartitions eq "yes");
	$linea = $linea.sprintf("\t%.3g", $supporting);
	$linea = $linea.sprintf("\t%.4g", $sum_br);
	return ($linea);
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
#Ê/// Generate lines for partitions file output ///
sub get_parts_table {
	my ($ct, $namt, $nt, $ma) = @_;

	my $comptrees = ($$ct);
	my %nametrees = (%$namt);
	my $ntree = ($$nt);		
	my @matrix = (@$ma);				
	
	my $whatree = "";
	if (exists ($nametrees{$comptrees}{$ntree})){
		$whatree = $nametrees{$comptrees}{$ntree};
	}
	else{
		$whatree = " Tree ".($ntree+1)." ";
	}
	my $chars = length($whatree);
	my $upanddown = "";
	for (my $i=0; $i < ($chars+2); $i++){ $upanddown = $upanddown."-"; }
	my $linea = "\n$upanddown\n\|$whatree\|\n$upanddown\n";
	$linea = $linea."Brl_ref_tree\tBrl_comp_tree\tBrl_comp_tree_K\tPartition\n";
	for (my $i=0; $i < scalar(@matrix); $i++){
		if ($matrix[$i][4] eq "NONE"){
			$linea = $linea."$matrix[$i][4]\t".
			    sprintf ("%.5f\t%.5f\t", $matrix[$i][5],$matrix[$i][8]);
		}
		elsif ($matrix[$i][5] eq "NONE"){
			$linea = $linea.sprintf ("%.5f\t", $matrix[$i][4]).
			    "$matrix[$i][5]\t$matrix[$i][9]\t";
		}
		else{
			$linea = $linea.sprintf ("%.5f\t%.5f\t%.5f\t", $matrix[$i][2],$matrix[$i][3],$matrix[$i][8]);
		}

		$linea = $linea."(";
		my @part = ($matrix[$i][0]);
		for (my $h=0; $h < scalar (@part); $h++){
			for (my $q=0; $q < scalar (@{$part[$h]}); $q++){
 					$linea = $linea. "$part[$h][$q], ";
			}
		}
		$linea =~ s/, $//g;
		$linea = $linea.") / (";
		@part = ($matrix[$i][1]);
		for (my $h=0; $h < scalar (@part); $h++){
			for (my $q=0; $q < scalar (@{$part[$h]}); $q++){
 					$linea = $linea. "$part[$h][$q], ";
			}
		}
		$linea =~ s/, $/)\n/g;

	}

	return ($linea);		
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
	
sub EM_d1 
#Input: (\@x, \@w, $n, $x1, $x2, $fo, $nstates) Output: $log_lik, $z, $tau, $mu, $sig2
# 1 or 2 states, 1dim variable x
{ 
    my ($xx, $ww, $n, $x1, $x2, $fo, $nstates) = @_;
    my @x=(@$xx);
    my @w=(@$ww);
    my $i; # sample
    my $j; # state

    my @tau = ();
    my @mu = ();
    my @sig2 = ();

    my $l2pi=1.8378; # log(2pi)
    if($x2<0.00001){$x2=0.00001;}
    for($j=0; $j<$nstates; $j++){
	$tau[$j] = 1/$nstates;
	$mu[$j] = $x1;
	$sig2[$j]=$x2;
    }
    my $del2=$x2/10;
    my $del=sqrt($del2); # bin width of the distribution
    my $ld2=log($del2);
    my $log_lik=-0.5*$n*(1+$l2pi+log($x2)-$ld2);

    my @z=();
    my @tt = (); # membership
    if($nstates==1){
	for ($i=0; $i <$n;  $i++){$z[$i]=1;}
	return($log_lik, \@z, \@tau, \@mu, \@sig2);
    }

    # nstates =2
    $mu[0]=$x1-sqrt($x2);
    $mu[1]=$x1+sqrt($x2);

    # EM algorithm
    my $it_max=30; my $iter=0; my $eps=0.0001;
    my $log_lik_old=0;
    my $out=sprintf("# log_lik= %.2f (1 type)\n", $log_lik);
    #if($write_fit){print $fo $out;}
    for ($i=0; $i <$n;  $i++){$tt[0][$i]=0;}
    for($iter=0; $iter<$it_max; $iter++){
	# Expectation (membership probabilities)
	my @fact;
	for($j=0; $j<$nstates; $j++){
	    if($sig2[$j]>$del2){
		$fact[$j]=$del/sqrt($sig2[$j]);
	    }else{$fact[$j]=1;} # delta distribution
	}
	for ($i=0; $i <$n;  $i++){
	    if($w[$i]==0){next;}
	    my $sum=0;
	    for($j=0; $j<$nstates; $j++){
		my $d2=$x[$i]-$mu[$j]; $d2*=$d2;
		if($sig2[$j]>0){$d2/=$sig2[$j];}
		elsif($d2){$d2=100;}
		$tt[$j][$i]=$tau[$j]*$fact[$j]*exp(-0.5*$d2);
		$sum+=$tt[$j][$i];
	    }
	    for($j=0; $j<$nstates; $j++){if($sum){$tt[$j][$i]/=$sum;}}
	}
        # Maximization (new parameters)
	my $sum=0;
	for($j=0; $j<$nstates; $j++){
	    $tau[$j]=0; $mu[$j]=0; $sig2[$j]=0; 
	    for ($i=0; $i <$n;  $i++){
                if($w[$i]==0){next;}
		my $wt=$w[$i]*$tt[$j][$i];
		$tau[$j]+=$wt;
		$mu[$j]+=$wt*$x[$i];
		$sig2[$j]+=$wt*$x[$i]*$x[$i];
	    }
	    $mu[$j]/=$tau[$j];
	    $sig2[$j]=$sig2[$j]/$tau[$j]-$mu[$j]*$mu[$j];
	    $sum+=$tau[$j];
	}
	$log_lik=-0.5*$sum*(1+$l2pi-$ld2);
	for($j=0; $j<$nstates; $j++){
	    if($tau[$j]){$log_lik+=$tau[$j]*log($tau[$j]/$sum);}
	    if($sig2[$j]>$x2){$log_lik+=$tau[$j]*(-0.5*log($sig2[$j]));}
	    else{$log_lik+=$tau[$j]*(-0.5*$ld2);}
	    $tau[$j]/=$sum;
	}
	if(0 && $write_fit){
	    $out=sprintf("# log_lik= %.2f\n", $log_lik);
	    print $fo $out;
	}
	if($iter && ($log_lik-$log_lik_old)<$eps){last;}
	$log_lik_old=$log_lik;
    }
    if($iter==$it_max){
	print "WARNING, EM algorithm did not converge in $it_max iterations\n";
    }
    for ($i=0; $i <$n;  $i++){$z[$i]=$tt[0][$i];}
    return($log_lik_old, \@z, \@tau, \@mu, \@sig2);
}

sub get_nodes{
    my  ($cl, $brl, $file, $ntree, $tot_trees, $write) = @_;
    my @clados = (@$cl);
    my %brlen = (%$brl);
    my $n_nodes=scalar(@clados);
    #print "file= $file tree $ntree out of $tot_trees write= $write\n";

    my @node = ();
    # structure of array nodes:
    # index 0 = n.species
    # index 1: name
    # index 2: parent
    # index 3: branch length
    # index 4: list of species
    # index 5: name of parent
    # index 6: matrix index m (assigned in get_matrix)
    # index 7: son 1
    # index 8: son 2
    # index 9: son 3 (only for unresolved nodes, ex. root)

    my $max=1; my $nsp=0;
    for (my $i=0; $i < $n_nodes; $i++){
	$node[$i][0]=scalar(@{$clados[$i]});
	if($node[$i][0]==1){$node[$i][1]=$clados[$i][0]; $nsp++;}
	elsif($node[$i][0]>$max){$max=$node[$i][0];}
	$node[$i][3]=$brlen{$clados[$i]};  
	$node[$i][4]=$clados[$i];
    }

   # Name internal nodes
    my $nint=1;
    for(my $n=$max; $n>1; $n--){
	for (my $i=0; $i < $n_nodes; $i++){
	    if($node[$i][0]==$n){
		$node[$i][1]="Node$nint"; $nint++;
	    }
	}
    }
    $nint--;

    # Look for parent
    my $npar=0; my $nroot=0;
    for (my $i=0; $i < $n_nodes; $i++){
	$node[$i][2]=-1; # parent
	my $min=$n_nodes;
	for(my $j=0; $j<$n_nodes; $j++){
	    if($node[$j][0]<=$node[$i][0]){next;}
	    my $n_equal = 0;
	    foreach my $cl0 (@{$clados[$i]}){
		if (grep (/^$cl0$/, @{$clados[$j]})){
		    $n_equal++;
		}
	    }
	    if($n_equal==$node[$i][0] && $node[$j][0]<$min){
		$node[$i][2]=$j; $min=$node[$j][0];
	    }
	} 
	if($node[$i][2] >=0){$node[$i][5]=$node[$node[$i][2]][1]; $npar++;}
	else{ $node[$i][5]="Node0"; $nroot++;}
    }

    my $out="# $nsp species, $nint internal nodes, $n_nodes total nodes.";
    if($n_nodes==(2*$nsp-3)){$out=sprintf("%s Tree is unrooted", $out);}
    elsif($n_nodes==(2*$nsp-2)){$out=sprintf("%s Tree is rooted", $out);}
    else{print "$out\nI do not know if tree is rooted or unrooted\n"; die;}
    $out=sprintf("%s\n# %d parent nodes have been determined + %d with root\n",
		 $out, $npar, $nroot);
    #print $out;

    if($write){
	my $nameout=Remove_ext($file, 1);
	if($tot_trees>1){$nameout=sprintf("%s.%d", $nameout, $ntree+1);}
	$nameout=sprintf("%s.graph", $nameout);
	print "Printing tree in graph format in $nameout\n";
	open(my $fo, '>', $nameout);
	print $fo $out;
	my $b1=0; my $b2=0;
	for (my $i=0; $i < $n_nodes; $i++){
	    my $b=$node[$i][3]; $b1+=$b; $b2+=$b*$b;
	}
	$b1/=$n_nodes; $b2=($b2-$n_nodes*$b1*$b1);
	if($n_nodes>1){$b2=sqrt($b2/($n_nodes-1));}
	my $out=sprintf("# %d branches, mean: %.4g CoV: %.3g\n",
			$n_nodes, $b1, $b2/$b1);
	print $fo $out;
	print $fo "#node\tparent\tbranch\tsons\n";
	for (my $i=0; $i < $n_nodes; $i++){
	    $out=$node[$i][1]."\t".$node[$i][5].
		sprintf("\t%.4g",$node[$i][3])."\t".$node[$i][0];
	    #for(my $j=0; $j<$node[$i][0]; $j++){$out=$out."\t".$node[$i][4]->[$j];}
	    print $fo $out, "\n";
	}
	close($fo);
    }

    return(@node);
}

sub Remove_ext{
    my $skip=$_[1];
    my @word=split(/\./, $_[0]);
    my $out=$word[0];
    for(my $i=1; $i<(scalar(@word)-$skip); $i++){
	$out=sprintf("%s.%s", $out, $word[$i]);
    }
    return($out);
}
