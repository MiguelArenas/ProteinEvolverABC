#!/usr/bin/env perl 
use strict;
use warnings;

my $nrun=50;
my $List="";

# Constant parameters (change when installing on a different machine)
my $script_SSCPE="/data/ubastolla/BIN/SSCPE/SSCPE.pl";
my $script_raxml="/data/ubastolla/BIN/SSCPE/raxml-ng";
my $script_Prot_evol="/data/ubastolla/BIN/SSCPE/Prot_evol";
my $submit="qsubmit.pl -q x86_64 -n 8 -s";
my $FLUX=0;

# Modifiable input:
my $pdbdir="/ngs/databases/pdb";
my $ali_dir=".";
my @msa;
my @pdb;
my @chain="A";
my $setrate=1;
my @models=("MF","WT","DE","RMSD","DEMF","DEWT","RMSDMF","RMSDWT");
#my @models=("WT");

my @mod_emp=("JTT","WAG","LG");
#my @mod_emp=();
my @mod2=("","+F","+G","+F+G");
my $reuse=1;
my $gamma=0;

print "Script $0 Author: ",
    "Ugo Bastolla <ubastolla\@cbm.csic.es>\n",
    "It needs as input a list of proteins (-list) with three columns:\n",
    "pdb.code, MSA containing the sequence of the PDB, chain (optional)\n",
    "The first line (optional) is ALI= <folder where MSA are stored>\n",
    "it runs the programs tnm and Prot_evol (or reuses previous results ",
    "if they exist, unless -noreuse is set), it generates site-specific ",
    "stability and structure constrained (SSCPE) substitution matrices, ",
    "it transforms them in format readable by the program RAxML-NG ",
    "and it runs RAxML-NG\n\n";

if(scalar(@ARGV)<2){help();}
for(my $i=0; $i<scalar(@ARGV); $i++){
    if($ARGV[$i] eq "-list"){
	$i++; $List=$ARGV[$i];
    }elsif($ARGV[$i] eq "-n"){
       $i++; $nrun=$ARGV[$i];
    }elsif($ARGV[$i] eq "-gamma"){
	$gamma=1;
    }elsif($ARGV[$i] eq "-pdbdir"){
       $i++; $pdbdir=$ARGV[$i];
    }elsif($ARGV[$i] eq "-alidir"){
       $i++; $ali_dir=$ARGV[$i];
    }elsif($ARGV[$i] eq "-rate"){
# 0: BU=1; 1: BU=1./site-spec subst.rate; 2: BU=rate",
	$i++; $setrate=$ARGV[$i];
	print "Site-specific substitution rates: $setrate\n";
    }elsif($ARGV[$i] eq "-h"){
	help();
    }else{
	print "WARNING, unknown option ",$ARGV[$i],"\n";
    }
}

if($List eq ""){
    print "ERROR, the option -list is mandatory\n";
    help();
}

open(my $fh, '<', $List)
or die "Could not open file '$List' $!";
my $nprot=0;
while (my $row = <$fh>){
    if(substr($row, 0, 1) eq "#"){next;}
    chomp($row); my @word=split(/\s+/, $row);
    if(scalar(@word)<2){next;}
    if(substr($row, 0, 3) eq "ALI"){
	$ali_dir=$word[1]; next;
    }
    $msa[$nprot]=$word[0];
    $pdb[$nprot]=$word[1];
    if(0 && scalar(@word)>2){$chain[$nprot]=$word[2];}
    else{$chain[$nprot]="A";}
    $nprot++;
}
print "$nprot alignments read in file $List\n";
close $fh;
my $nmod=scalar(@models)+scalar(@mod_emp)*scalar(@mod2);
my $nmod2=$nmod+8;
my $nmin=$nprot; if($nrun < $nmin){$nmin=$nrun;}
print "Inferring trees for $nmod models and $nmin alignments\n";

my $n=0;
for(my $i=0; $i<$nprot; $i++){
    my $name=Remove_extension($msa[$i]);
    my @tree=`ls -1 ${name}*.bestTree`;
    my $done=scalar(@tree);
    if($done==$nmod2 && $tree[0] ne ""){
	    print "$done trees already inferred for alignment $msa[$i], ",
	    "skipping\n";
	    next;
    }
    my $pdbtype="pdb";
    my $pdbfile="";
    if($pdbdir ne ""){
	$pdbfile=$pdbdir;
	if(substr($pdbdir, length($pdbdir)-1, 1) ne "/"){
	    $pdbfile=$pdbfile."/";
	}
    }
    $pdbfile=$pdbfile."${pdb[$i]}.pdb";
    if(!-e $pdbfile){
	print "WARNING, $pdbfile does not exist\n";
	if(-e $pdb[$i] && $pdbdir ne ""){
	    print "$pdb[$i] exists, interpreting it as a list of PDB\n";
	    $pdbfile=$pdb[$i]; $pdbtype="pdblist";
	}else{
	    print "ERROR, neither $pdbfile nor $pdb[$i] exist\n"; die;
	}
    }

    my $script="tmp.$name";
    print "Writing script $script for alignment $msa[$i]\n";
    open(my $fo, '>', $script);
    my $out;
# Check if Prot_evol results exist, otherwise run it
    my $summ_file="";
    if($pdbtype eq "pdb"){$summ_file="${pdb[$i]}";}
    else{$summ_file="${name}.";}
    my @summary=`ls -1 ${summ_file}*summary.dat`;
    if($reuse==0 || scalar(@summary)==0 || $summary[0] eq ""){
	$out="$script_Prot_evol -ali $msa[$i] -$pdbtype $pdbfile\n";
	print $fo $out; print $out;
    }
    foreach my $model (@models){
	my $ali="${name}.${model}.HB1.FL${FLUX}.rate${setrate}.fasta";
	$out="cp $msa[$i] $ali\n";
	$out=$out.
	    "$script_SSCPE -ali $ali -pdbdir ${pdbdir} -$pdbtype $pdbfile".
	    " -model $model -rate $setrate";
	if($FLUX){$out=$out." -flux";}
	$out=$out." -raxml\n";
	$out=$out."rm -f $ali\n";
	print $out; print $fo $out; 
    }
    foreach my $model (@mod_emp){
    foreach my $variant (@mod2){
	my $ali="${name}.${model}";
	if($variant eq "+F+G"){$ali=$ali."_F_G";}
	elsif($variant eq "+G"){$ali=$ali."_G";}
	elsif($variant eq "+F"){$ali=$ali."_F";}
	$ali=$ali.".fasta";
	$out="cp $msa[$i] $ali\n";
	$out=$out."$script_raxml --msa $ali -model $model$variant\n";
	$out=$out."rm -f $ali\n";
	print $out; print $fo $out; 
    }
    }
	
    close $fo;
    `chmod u+x $script`;
    `$submit $script`;
    $n++;
    if($n==$nrun){
       print "$n computations sent to the cluster, exiting\n"; last;
    }
}
if($n!=$nrun){
    print "$n computations sent to the cluster, exiting\n";
}

sub help{
    print
	"USAGE: ",$0,
	" -list <list of proteins> (1:MSA 2:pdb or pdb_list 3: chain)\n",
	" -n number of MSA to run, default $nrun\n",
	" -pdbdir <path to PDB files>\n",
	" -alidir <path to MSA files>\n",
	" -rate <0,1,2> 1: BU=site-spec subst.rate; 0: BU=1; 2: BU=1/rate\n",
	" -gamma ! Use +G option in RAxML-NG\n",
	" -noreuse (Do not reuse TNM or Prot_evol results ",
	"even if they already exist)\n",
	"FORMAT of PDB list: 1: pdb code 2: chain\n",
	"\n";
    die;
}

sub Remove_extension{
    my @word=split(/\./, $_[0]);
    my $out=$word[0];
    for(my $i=1; $i<(scalar(@word)-1); $i++){
	$out=sprintf("%s.%s",$out,$word[$i]);
    }
    return $out;
}
