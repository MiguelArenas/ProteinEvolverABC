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
#  script_read_tree.pl <options> -tree <Newick file> 
#

use warnings;
use strict;

# Ktreedist	version 1.0
# Last modified 25 June 2007
my $version = "2.0";
my $monthyear = "November 2023";


# //////////////////////////////////////////////////////////////////////////////
# //////////////////////////////  HIDDEN OPTIONS  //////////////////////////////
# //////////////////////////////////////////////////////////////////////////////

# The variable hard_polytomies is "off" by default.
# By default, zero branch length partitions count in the calculation of the symmetric difference.
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

# //////////////////////////////////////////////////////////////////////////////
# ///////////////////////////////  PRELIMINARS  ////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////

print "\n$0 version 1.0 -  Nov 2023";
print "\nFor all nodes in the tree, prints parent branch length ingroup\n";
print "Based on Ktreedist by Castresana's group\n";

# Get command line arguments
my $tree="";
my $print_sons=1;
chomp($ARGV[scalar(@ARGV)-1]);
for(my $i=0; $i<scalar(@ARGV); $i++){
    if($ARGV[$i] eq "-tree"){
	$i++; $tree=$ARGV[$i];
    }elsif($ARGV[$i] eq "-nosons"){
	$print_sons=0;
    }else{
	print "WARNING, unknown option ",$ARGV[$i],"\n";
    }
}
if($tree eq ""){
    print "ERROR, option -tree is mandatory\n";
    usage();
}


# DEFAULT DEFINITIONS, PATHS & FILENAMES
#-------------------------------------------------------------------------------
# Get paths & filenames
my ($filename_tree, $filename_tree_woext, $path_tree) = 
    get_path_filename($tree);
#This is really necessary in case of no full path
$tree = "$path_tree$filename_tree" if ($tree ne "");

# Do initial checkings and get trees if there are no errors
print "Checkings:\n";
my ($reftr, $nametr, $txbl, $root, $bl, $t0) = initial_check ($tree);
my $reftree = ($$reftr);
my %nametrees = (%$nametr);
my $taxblock = ($$txbl);
my $rooted = ($$root);
my $bk_line = ($$bl);
my $tree0 = ($$t0);

# //////////////////////////////////////////////////////////////////////////////
# ///////////////////////////////////  MAIN  ///////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////


my ($claA, $claB, $brl, $sp) = get_partitions_tree($tree0, "yes");
my @clados0A = (@$claA);
my @clados0B = (@$claB);
my %brlen0 = (%$brl);
my @species0 = (@$sp);
my @clados0 = (@clados0A, @clados0B);

my @nodes0=get_nodes(\@clados0A, \%brlen0, $filename_tree, 0, 1, 1, $rooted);


exit();


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------- SUBROUTINES ----------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

sub help{

    print "Program $0 For each node of the tree prints: parent, branch length, number of sons, sons\n";
    print
	"USAGE: ",$0,"\n",
	" -tree <Newick file> (mandatory)\n\n";
}

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
    print " $0 -tree <reference tree file>/s file> [<options>]\n";
    print "    Options:\n";
    print "         -nosons ! Do not print list of sons\n";
    print "         \n";
    exit();
}
#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

#-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

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
    if ($reftree eq ""){usage();}
    elsif (!-e $reftree){
	$ERROR = "     File '$reftree' does not exist\n";
	die("ERROR:\n$ERROR\n");
    }
    
    # =====================================
    # =====================================
    
    # =====================================
    #Check line break
    # =====================================
    print "Line break...  ";
    my @lines_ref = break_line($reftree);
    my $breakline_ref = shift(@lines_ref);
    my $bk_line_line = "OK ($breakline_ref)";
    print "$bk_line_line\n";
    
    # =====================================
    # =====================================
    
    # =====================================
    # Check output files
    # =====================================

    # =====================================
    # =====================================
    
    # =====================================
    # CHECK if files are nexus or newick
    # =====================================
    print "File format...  ";
    
    my @lines = ([@lines_ref]);
    my $fileformat_line = "";
    my $formaterr = "";
    my @reftrees = ();
    if(scalar(@lines)>1){
	print "WARNING, ",scalar(@lines)," trees read, only the first one will be used\n";
    }

    for (my $l=0; $l < 1; $l++){ #scalar(@lines)
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
	    }else{
		$fileformat_line = $fileformat_line." & newick";
	    }
	    if ($l == 0){ @reftrees = @trs; }
	}
	
	if ($newick eq "no"){
	    # Check nexus
	    # -------------------------------------
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
	    }else{
		my ($nametr, $txbl, $trs) = nexus2newick(\@aux);
		my %auxname = (%$nametr);
		$taxblock = ($$txbl);
		if ($l > 0){last;}
		@reftrees = (@$trs);
		$fileformat_line = $fileformat_line." with taxa block" if ($taxblock eq "yes");
		$nametrees{$reftree} = \%auxname;
		$nexus = $nexus."yes";
	    }
	}else{ 
	    $nexus = $nexus."no";
	}
	
	if (($newick eq "no") && ($nexus =~ m/no$/)){
	    if ($formaterr ne ""){
		$formaterr = $formaterr."\n          &\n          format of file "; 
	    }else{
		$formaterr = "Format of file "; 
	    }
	    if ($l == 0) {
		$formaterr = $formaterr." '$reftree' not recognized";
	    }
	}
    }
    if ($formaterr ne ""){
	print "\n";
	die ("     ERROR:\n          $formaterr.\n\n");
    }else{
	print "OK ($fileformat_line)\n";
    }
    # =====================================
    # =====================================
    
    # =====================================
    # Reading trees
    # =====================================
    my @rtrs = ();
    my @auxtrees = (@reftrees);
    
    for (my $t=0; $t < scalar(@auxtrees); $t++){
	my @trees = ();
	my $line = "";
	chomp($auxtrees[$t]);
	next if ($auxtrees[$t] =~ m/^[\n|\/|\#|\>]/);
	$auxtrees[$t] =~ s /^\[.+\]//g;
	if ($auxtrees[$t] =~ m/\;$/){
	    $line = $line.$auxtrees[$t];
	    if ($t < scalar (@reftrees)){ push (@rtrs, $line); }
	    $line = "";
	}else{
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
    }else{
	print "OK\n";
    }
    my $t0 = shift(@rtrs); #The first reference tree will be used
    
    # =====================================
    # =====================================
    
    # =====================================
    # Check branch lengths
    # =====================================
    print "Checking branch lengths...  ";
    my @trs = ($t0);
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
	}else{ # It's comptrees
	    last;
	}
	
	#Checks branch lengths are right (not negative)
	my @aux = split (/\:|\)|\,/, $trs[$i]);
	foreach my $a (@aux){
	    if ($a=~ m/^\-[0-9]+/){
		$wrongbrlen = "yes";
		if ($i == 0){					
		    $wrrefbr = "no";
		}else{
		    last;
		}
	    }
	}		
    }
    $indexes =~ s/, $//g;
    $indexes_wrbrlen =~ s/, $//g;
    
    if ($nobrlen eq "yes"){
	my $brlenerr = "";
	if ($refbr eq "no") {
	    $brlenerr = "          Tree '$reftree' has no branch lengths.\n";
	}
	$ERROR = $ERROR."\n     ERROR (branch lengths):\n$brlenerr\n";
	print "$ERROR";
	$ERR = 1;
    }else{
	
	$ERROR = "$ERROR\n\n" if ($ERROR ne "");
	my $brlenerr = "";	
	if ($wrongbrlen eq "yes"){
	    if ($wrrefbr eq "no"){
		$brlenerr = "          Branch lengths of tree '$reftree' are incorrect (negative values).\n";			
	    }
	    if ($nobrlen eq "yes"){
		$ERROR = $ERROR."\n\n$brlenerr\n";
	    }else{
		$ERROR = "\n     ERROR (branch lengths):\n$brlenerr\n";			
	    }
	    print "$ERROR";
	    $ERR = 1;			
	}else{
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
	}else{
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
        
    # =====================================
    # =====================================
    
    # =====================================
    # Check root
    # =====================================
    print "All trees rooted or all unrooted...  \n";
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
	}else {
	    push (@treenumber, $n_trees."U");
	}
	$n_trees++;
    }
    

        
    my $errormsg = "\n*********************************************\n".
	'  There were one or more ERRORS. See above.  '.
	"\n*********************************************\n\n";
    
    die ("$errormsg") if ($ERR == 1);

    return (\$reftree, \%nametrees, \$taxblock, \$rooted, \$bk_line, \$t0);
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
				if ($treesplit[$i-2] ne "\)"){# it's a tip and not a bootstrap value
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

	#ÊGet symmetric ones
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

sub get_nodes{
    my  ($cl, $brl, $file, $ntree, $tot_trees, $write, $rooted) = @_;
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

    my $max=1; my $nsp=0;
    for (my $i=0; $i < $n_nodes; $i++){
	$node[$i][0]=scalar(@{$clados[$i]});
	if($node[$i][0]==1){$node[$i][1]=$clados[$i][0]; $nsp++;}
	elsif($node[$i][0]>$max){$max=$node[$i][0];}
	$node[$i][3]=$brlen{$clados[$i]};  
	$node[$i][4]=$clados[$i];
    }

   # Name internal nodes
    my $nint;
    if($rooted eq "yes"){$nint=0;}
    else{$nint=1;}
    for(my $n=$max; $n>1; $n--){
	for (my $i=0; $i < $n_nodes; $i++){
	    if($node[$i][0]==$n){
		$node[$i][1]="Node$nint"; $nint++;
	    }
	}
    }
    if($rooted ne "yes"){$nint--;}

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
	if($node[$i][2] >=0){
	    $node[$i][5]=$node[$node[$i][2]][1]; $npar++;
	}else{
	    if($rooted ne "yes"){$node[$i][5]="Node0";} $nroot++;
	}
    }

    my $out="# $nsp species, $nint internal nodes, $n_nodes total nodes.";
    if($n_nodes==(2*$nsp-3)){$out=sprintf("%s Tree is unrooted", $out);}
    elsif($n_nodes==(2*$nsp-2)){$out=sprintf("%s Tree is rooted", $out);}
    else{print "$out\nI do not know if tree is rooted or unrooted\n"; die;}
    $out=sprintf("%s\n# %d parent nodes have been determined + %d with root\n",
		 $out, $npar, $nroot);
    #print $out;

    if($write){
	my $nameout=Remove_ext($file);
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
	    $out=sprintf("%s\t%s\t%.4g\t%d",
			 $node[$i][1], $node[$i][5], $node[$i][3], $node[$i][0]);
	    if($print_sons){
		for(my $j=0; $j<$node[$i][0]; $j++){$out=$out."\t".$node[$i][4]->[$j];}
	    }
	    print $fo $out,"\n";
	}
	close($fo);
    }

    return(@node);
}

sub Remove_ext{
    my @word=split(/\./, $_[0]);
    my $out=$word[0];
    for(my $i=1; $i<(scalar(@word)-1); $i++){
	$out=sprintf("%s.%s", $out, $word[$i]);
    }
    return($out);
}
