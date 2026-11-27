#!/usr/bin/env perl 
use strict;
use warnings;
my $input="Prot_evol.c";
my $output="tmp";
my $nchange=4;

my $dir=`pwd`; chomp ($dir);
print "Working directory: ",$dir,"\n";

chdir("..");
my $maindir=`pwd`; chomp ($maindir);
print "Installation directory: ",$maindir,"\n";
chdir($dir);

open(my $fo, '>', $output);
open(my $fh, '<', $input)
or die "Could not open file '$input' $!";
my $written=0;
while (my $row = <$fh>){
    if($written < $nchange){
	if(substr($row, 0, 13) eq "char FILE_STR"){
	    print $fo "char FILE_STR[]=\"$maindir/structures.in\";\n";
	    $written++; next;
	}elsif(substr($row, 0, 9) eq "char tnm["){
	    print $fo "char tnm[]=\"$maindir/tnm\";\n";
	    $written++; next;
	}elsif(substr($row, 0, 18) eq "char tnm_mut_para["){
	    print $fo "char tnm_mut_para[]=\"$maindir/Mut_para.in\";\n";
	    $written++; next;
	}elsif(substr($row, 0, 12) eq "char DIR_TNM"){
	    my $TNM_DATA=$maindir."/TNM_DATA";
	    print $fo "char DIR_TNM[]=\"$TNM_DATA\";\n";
	    `mkdir $TNM_DATA`;
	    print "mkdir $TNM_DATA\n";
	    $written++; next;
	}
    }
    print $fo $row;
}
close $fh;
close $fo;


if($written<$nchange){
    print "WARNING, $input was modified at $written out of $nchange lines\n";
}else{
    print "$input was modified at $nchange lines\n";
}
if($written){`mv -f $output $input`;}
