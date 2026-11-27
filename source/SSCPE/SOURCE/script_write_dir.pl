#!/usr/bin/env perl 
use strict;
use warnings;
my $SSCPE="SSCPE.pl";
my $run_SSCPE="script_run_SSCPE_models.pl";
my $Input_TNM="Input_TNM.in";
my $Input_Prot_evol="Input_Prot_evol.in";

my $dir=`pwd`; chomp ($dir);
print "Working directory: ",$dir,"\n";

# Write directory in SSCPE.pl
Write_dir("my \$dir_SSCPE", $dir, $SSCPE);
Write_dir("my \$dir_SSCPE", $dir, $run_SSCPE);
`chmod u+x $SSCPE`;
`chmod u+x $run_SSCPE`;

my $file=$dir."/Mutation_para.in";
Write_dir("MUT_PARA", $file, $Input_TNM);
Write_dir("FILE_STR", $file, $Input_Prot_evol);
print "All done\n";

sub Write_dir{
    my ($string, $dir, $file)=@_;

    my $modified=0;
    my $tmp="tmp.out";
    open(my $fo, '>', $tmp);
    open(my $fh, '<:encoding(UTF-8)', $file)
	or die "Could not open file '$file' $!";
    my $l=length($string);
    while (my $row = <$fh>){
	if(substr($row, 0, $l) eq $string){
	print $fo "$string=\"",$dir,"\";\n";
	$modified++;
    }else{
	print $fo $row;
    }
    }
    close $fh;
    close $fo;
    if($modified){
	`mv -f $tmp $file`;
	print "$file was modified\n";
    }else{
	print "WARNING, $file, was not modified\n";
    }
}
