#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(getcwd);
use File::Copy;
use File::Basename;
use Config;
use FindBin;

# Authors: Ugo Bastolla <ubastolla@cbm.csic.es> 

# path to programs and databases, modify as needed:
my $dir_SSCPE="/data/ubastolla/BIN/SSCPE";

# install programs and data in the above directory:
my $Mut_para="$dir_SSCPE/Mutation_para.in";
my $Prot_evol="$dir_SSCPE/Prot_evol";
my $raxmlng="$dir_SSCPE/raxml-ng";
my $tnm="$dir_SSCPE/tnm";
my $script_get_br="$dir_SSCPE/script_get_br_len.pl";
my $structures="$dir_SSCPE/structures.in";
my $mat_LG="$dir_SSCPE/lg.txt";
my $mat_WAG="$dir_SSCPE/wag.txt";
my $mat_JTT="$dir_SSCPE/jtt.txt";
my $pdbdir="";
my $pdblist="";
my $Join_nomod=1;  # Join all not modelled sites in one partition
my $Emp_nomod=1; # Sites not modelled with PDB are
# described by empirical model instead of P_mut
my $F_nomod=0;  # Sites not modelled with PDB run with option +F
my $nostr=0;
my $last_model=-1;
my $tree="";
my $print_exch=0;

my @MODELS=("LG","WAG","JTT",   # Empirical
	    "MF","WT",    # STAB "DDG",
	    "DE","RMSD",  # STRUCT
	    "DEWT","RMSDWT", # SSCPE
	    "DEMF","RMSDMF", # SSCPE
	    "ML");
	    #"DE2","RMSD2","DE2DDG","RMSD2DDG" # $allmut==1

my $nmod1=scalar(@MODELS);

print "Script $0 Author: ",
    "Ugo Bastolla <ubastolla\@cbm.csic.es> ",
    "with the collaboration of Miguel Arenas <marenas\@uvigo.es>\n",
    "It needs as input an MSA file (-ali) and a list of pdb files (-pdblist) ",
    "and their path (-pdbdir) or a single pdb file (-pdb)\n",
    "It runs the programs Prot_evol (or reuses previous results ",
    "if they exist, unless -noreuse or -noreuse2 is set), which generates ",
    "site-specific structure and folding-stability constrained (SSCPE) ",
    "substitution models, transforms them in format readable by the ",
    "programs RAxML-NG and PAML and infers ML trees for the model -model <> ",
    "with the program raxml-ng upon request (-raxml) ",
    "or infers the branches of a given tree (-evaluate -tree <>) ",
    "if raxml-ng is installed in the directory $dir_SSCPE.\n",
    "If -print_exch is used, it prints the exchangeability matrices of the ",
    "specified model",
    "\n\n";

my @symbols=("A","C","D","E","F","G","H","I","K","L",
	     "M","N","P","Q","R","S","T","V","W","Y");

# Input arguments
my $msa="";
my $pdb="";
my $chain="";
my $reuse_TNM=1; # Reuse previous TNM results if they exist
my $reuse_PE=1;  # Reuse previous Prot_evol results if they exist
my $runraxml=0;  # Run RAxML-NG for phylogenetic analysis
my $thread=4;    # Number of threads of raxmlng
my $clean=1;     # Clean site-specific matrices after running raxmlng
my $setrate=1;
my $model="ML";
my $matrix="OPT";
my $index="";
my $option="";
my $gamma=0;
my $allmut=0;
my $iwt=0;
my $kmod=-1;
my $compute_exch=1;
my $flux=1;
my $HB=1;
my $evaluate=0;

if(scalar(@ARGV)<2){help();}
#chomp($ARGV[scalar(@ARGV)-1]);
for(my $i=0; $i<scalar(@ARGV); $i++){
    if($ARGV[$i] eq "-ali"){
	$i++; $msa=$ARGV[$i];
    }elsif($ARGV[$i] eq "-pdb"){
	$i++; $pdb=$ARGV[$i];
    }elsif($ARGV[$i] eq "-chain"){
	$i++; $chain=$ARGV[$i];
    }elsif($ARGV[$i] eq "-pdblist"){
	$i++; $pdblist=$ARGV[$i];
    }elsif($ARGV[$i] eq "-pdbdir"){
   	$i++; $pdbdir=$ARGV[$i];
    }elsif($ARGV[$i] eq "-raxml"){
	$runraxml=1;
    }elsif($ARGV[$i] eq "-print_exch"){
	$print_exch=1;
    }elsif($ARGV[$i] eq "-evaluate"){
	$evaluate=1;
    }elsif($ARGV[$i] eq "-tree"){
	$i++; $tree=$ARGV[$i];

    }elsif($ARGV[$i] eq "-model"){
	$i++; $kmod=-1; my $k=0;
	foreach my $mod (@MODELS){
	    if($ARGV[$i] eq $mod){$kmod=$k; last;} $k++;
	}
	if($kmod<0){
	    print "ERROR, $ARGV[$i] is not an allowed model\nAllowed:";
	    foreach my $mod (@MODELS){print " ",$mod;}
	    die;
	    #print "\nUsing default $model\n"; $k=0;
	    foreach my $mod (@MODELS){
		if($model eq $mod){$kmod=$k; last;} $k++;
	    }
	}else{
	    $model=$ARGV[$i];
	    print "Using model $model of Prot_evol\n";
	}
    }elsif($ARGV[$i] eq "-gamma"){
	$gamma=1;
    }elsif($ARGV[$i] eq "-allmut"){
	$allmut=1;
    }elsif($ARGV[$i] eq "-iwt"){
	$i++;
	if($ARGV[$i]<0){
	    print "ERROR, iwt=$ARGV[$i] not allowed\n"; die;
	    $i--;
	}else{
	    $iwt=$ARGV[$i]; # $reuse_PE=0;
	}
    }elsif($ARGV[$i] eq "-matrix"){
	$i++; 
	if($ARGV[$i] ne "OPT" && 
	   $ARGV[$i] ne "JTT" &&
	   $ARGV[$i] ne "WAG" &&
	   $ARGV[$i] ne "LG"){
	    print "ERROR, $ARGV[$i] is not an allowed substitution matrix\n",
	    "Only allowed: OPT (default) LG WAG JTT\n";
	    #"Using default $matrix\n";
	    die;
	}else{
	    $matrix=$ARGV[$i];
	    print "Using substitution matrix $matrix in Prot_evol\n";
	}
    }elsif($ARGV[$i] eq "-nohb"){
	$HB=0;
    }elsif($ARGV[$i] eq "-noflux"){
	$flux=0;
    }elsif($ARGV[$i] eq "-index"){
	$i++; $index=$ARGV[$i];
    }elsif($ARGV[$i] eq "-noclean"){
	$clean=0;
    }elsif($ARGV[$i] eq "-option"){
	$i++; $option=$ARGV[$i];
    }elsif($ARGV[$i] eq "-thread"){
    	$i++; $thread=$ARGV[$i];
    }elsif($ARGV[$i] eq "-noreuse"){
	$reuse_PE=0; $reuse_TNM=1;
   }elsif($ARGV[$i] eq "-nostr"){
	$nostr=1;
    }elsif($ARGV[$i] eq "-noreuse2"){
	$reuse_PE=0; $reuse_TNM=0;
    }elsif($ARGV[$i] eq "-h"){
	help();
    }elsif($ARGV[$i] eq "-rate"){
# 1: BU=site-spec subst.rate; 0: BU=1; -1: BU=1/rate",
	$i++; $setrate=$ARGV[$i];
	print "Site-specific substitution rates: $setrate\n";
    }else{
	print "WARNING, unknown option ",$ARGV[$i],"\n";
    }
}
if($pdb eq "" && $pdblist eq ""){
    print "ERROR, PDB file must be provided either as list (-pdblist -pdbdir) ",
    "or as single file (-pdb -chain)\n";
    help();
}
if($msa eq ""){
    print "ERROR, alignment file must be provided with option -ali <file>\n";
    help();
}
if($kmod>= $nmod1 && $allmut==0){
    print "ERROR, the required model $model requires the option -allmut\n";
    die;
}

if($evaluate){
    if($tree eq ""){
	print "ERROR, with option -evaluate you have to specify a tree -tree\n";
	die;
    }elsif($runraxml==0){
	print "WARNING, with option -evaluate you need -raxml\n";
	$runraxml=1;
    }
}

# Prepare directories
getcwd();
my $maindir=Print_dir(1);

my $dirstore=".";
#$dirstore="${maindir}/$pdbcode";

unless(-d $dirstore){
    `mkdir $dirstore`;
    my $out="Creating folder ${dirstore}\n"; print $out;
}


####################
#   MSA
#####################
my $msa_file=File_name($msa);
if(substr($msa,0,1) ne "/"){$msa="$maindir/$msa";}
# copy MSA to directory dirstore
my $tmp="${dirstore}/$msa_file";
unless(-e $tmp){`cp $msa $dirstore`;}
my $out="MSA file: $msa\n";
print $out;

# RaxML options
if($option && -e $option){
    my $doption="${dirstore}/$option";
    unless(-e $doption){`cp $option $dirstore`;}
}

chdir $dirstore;
my $workdir=Print_dir(0);

#################################
#      Define the model 
#################################

# Model name
#my @msa_tmp=split(/\./, $msa_file);
#my $MyModel=$msa_tmp[0];
my $MyModel=Remove_extension($msa_file);
my $SSCPE=1;
if($model eq "LG" || $model eq "WAG" || $model eq "JTT"){
    $MyModel="${MyModel}.$model";
    if($gamma){$MyModel="${MyModel}_G";}
    if($index){$MyModel="${MyModel}.$index";}
    $out="Empirical model, TNM and Prot_evol computations are not needed\n";
    print $out; 
    $SSCPE=0; goto end_SSCPE;
}

######################
# SSCPE model 
#######################
my $pdbcode="";
my $prot="";
my $dir_tmp="";
my $name_SSCPE="";
my $name_prot="";

if($pdbdir && substr($pdbdir,0,1) ne "/"){$pdbdir="$maindir/$pdbdir";}

if($pdb){
    $pdbcode=File_name($pdb);
    $pdbcode=substr($pdbcode,0,4);
    $prot="${pdbcode}${chain}";
    $name_SSCPE=$prot;
    $name_prot=$prot;
    if(substr($pdb,0,1) ne "/"){
       if($pdbdir){$pdb="${pdbdir}/$pdb";}
       else{$pdb="${maindir}/$pdb";}
    }
    $out="PDB file: $pdb\n";
    print $out;
    $dir_tmp="tmp_$pdbcode";
}elsif($pdblist){
    print "Reading list of PDB files in $pdblist path= $pdbdir\n";
    $name_SSCPE=Remove_extension($pdblist);
    $name_prot="pdblist";
}else{
    print "ERROR, pdb file must be provided either in list (-pdblist) ",
    "or as unique file (-pdb)\n";
    die;
}

unless(-e $Prot_evol){
    print "ERROR, wrong path to the program Prot_evol, $Prot_evol\n";
    print "Probably the folder name $dir_SSCPE in SSCPE.pl is wrong\n";
    print "Please rerun the installation script script_install_SSCPE.sh ",
    "in the folder where you downloaded SSCPE.zip\n";
    die;
}

# Input file names
my $file_RMSD="";
my $file_DE=""; 
my $file_RMSD_all="";
my $file_DE_all=""; 
my $file_summ_tnm="";
my $file_summ_Prot_evol="";
my $file_exch_SSCPE="";
my $file_freq_SSCPE="";
my $file_rate_SSCPE="";
my $buffer="";

my $start=$workdir."/";
if($name_SSCPE){$start=$start.$name_SSCPE;}
if($reuse_PE){
# Check if Prot_evol results exist, otherwise run Prot_evol
    $out="Checking if Prot_evol results exist in ${workdir}\n";
    print $out; $buffer=$buffer.$out;
    $file_summ_Prot_evol=Get_file($start, "SSCPE.summary.dat");
}else{
    $out="Prot_evol results not reused even if they exist\n";
    print $out; $buffer=$buffer.$out;
}

if(-e $file_summ_Prot_evol){
    $out="Prot_evol results $file_summ_Prot_evol exist in $workdir\n";
    print $out; $buffer=$buffer.$out;
}else{
    # Prot-evol results do not exist, run the program
    $out="Prot_evol results recomputed, running Prot_evol\n";
    print $out; $buffer=$buffer.$out;
    Run_Prot_evol($pdb, $chain, $pdblist, $pdbdir, $msa);
    $file_summ_Prot_evol=Get_file($start, "SSCPE.summary.dat");
}

if($model eq "ML"){
    ($model, $HB, $flux, $setrate)=Read_model($file_summ_Prot_evol);
}
my $mod_name=$model;

Get_files_Prot_evol(".", $name_SSCPE, $model);

$out="Prot_evol results: $file_summ_Prot_evol\n$file_freq_SSCPE\n";
if($file_exch_SSCPE){$out=$out."${file_exch_SSCPE}\n";}
if($file_rate_SSCPE){$out=$out."${file_rate_SSCPE}\n";}

print $out; $buffer=$buffer.$out;
$out="Using Prot_evol results $file_summ_Prot_evol";
if($file_freq_SSCPE){$out="$out ${file_freq_SSCPE}\n";}
else{$out="$out ${file_exch_SSCPE} ${file_rate_SSCPE}\n";}
print $out; $buffer=$buffer.$out;

#my $file_map=""; my $subst_mod="";
#my $npdb=Read_summary($file_map, $subst_mod, $file_summ_Prot_evol);

my $subst_mod=Read_subst_mod($file_summ_Prot_evol);

#$MyModel="${MyModel}.${name_prot}";
$MyModel="${MyModel}.$mod_name";
$MyModel="${MyModel}.HB$HB";
#$MyModel="${MyModel}.$subst_mod";
#if($flux){$MyModel="${MyModel}_flux";}
#else{$MyModel="${MyModel}_exch";}
$MyModel="${MyModel}.FL$flux";
$MyModel="${MyModel}.rate${setrate}";
if($index){$MyModel="${MyModel}.$index";}

# log file
my $nameout="${MyModel}.log";
open(my $fo, '>', $nameout);
print "Writing log in $nameout\n";
print $fo $buffer;
$buffer="";

# Alignment
my $msa_new="${MyModel}.fasta";
FastaManyLines_FastaOneLine($msa_file, $msa_new);
my ($ntaxa, $lali, $sq, $naa)=Remove_all_gaps($msa_new);
my @seq=(@$sq);
my @num_aa=(@$naa);

$out="Results are stored in folder ${dirstore}\n";
print $out; print $fo $out;
$out="Alignment: $msa l=$lali\n";
print $out; print $fo $out;
if($prot){$out="Reference structure: $prot\n";}
else{$out="PDB list: $pdblist\n";}
print $out; print $fo $out;

if($file_freq_SSCPE eq "" &&
   ($file_exch_SSCPE eq "" || $file_rate_SSCPE eq "")){
    $out="Model $mod_name was not computed by Prot_evol, please change it\n";
    print $out; print $fo $out;
    die;
}

# Clean Prot_evol and TNM results
if(-d $dir_tmp){
    `rm -rf $dir_tmp`;
    $out="Deleting Prot_evol and TNM results in ${dir_tmp}\n";
    print $out; print $fo $out;
}

# Align and count sites modelled by Prot_evol
my $len_model=0;
my @ali_pdb;
for(my $i=0; $i<$lali; $i++){$ali_pdb[$i]=-1;}
my $file_in=$file_freq_SSCPE; #$file_exch_SSCPE
open(my $fh, '<:encoding(UTF-8)', $file_in) 
    or die "Could not open file '$file_in' $!";
while (my $row = <$fh>){
    if(substr($row, 0, 1) eq "#"){next;}
    my @word=split(/\s+/, $row); my $i;
    if(substr($row, 0, 1) ne " "){$i=$word[0];}
    else{$i=$word[1];}
    $ali_pdb[$i]=$len_model;
    $len_model++;
}
close $fh;
$out="$len_model amino-acid sites found in file ${file_in}\n";
print $out; print $fo $out;
my $lpdb=$len_model;

#my $lpdb=0;
#my @ali_pdb=Align_to_pdb($lpdb, $ntaxa, $lali, \@seq, $npdb, $file_map);
#if($lpdb != $len_model){
#    $out="ERROR, different number of sites in model $file_freq_SSCPE ".
#	"(${len_model}) and in mapping of PDB seq of MSA $msa_new (${lpdb})\n";
#    print $out; print $fo $out;
#    die;
#}


my $partition_file="${MyModel}.partitionsSite.txt";
my @rate; 
if($runraxml || $print_exch){
    if($compute_exch){
	$out="Reading site-specific frequencies in $file_freq_SSCPE\n";
	print $out; print $fo $out;
	my ($fr, $fm)=Read_frequencies($file_freq_SSCPE, $len_model);
	print $out; print $fo $out;
	my @freq=(@$fr); my @fmut=(@$fm); 
	$out="Computing site-specific rates and exchangeability\n";
	print $out; print $fo $out;
	@rate=Compute_print_exch($MyModel, $len_model, \@freq, \@fmut,
				 $subst_mod, $HB, $flux);
    }else{
	$out="Reading site-specific rates in $file_rate_SSCPE\n";
	print $out; print $fo $out;
	@rate=Read_rates($file_rate_SSCPE, $len_model);
	$out="Reading site-specific exchangeability in $file_exch_SSCPE\n";
	print $out; print $fo $out;
	Matrices_FromProtEvol_to_PAML($MyModel, $lali, $file_exch_SSCPE);
	$rate[$len_model]=1;
    }
    $last_model=$len_model;

    if($SSCPE==0){
	Identical_partitions($model, $lali, $gamma, $partition_file);
    }else{
	if($setrate==0){$out="Uniform rate at each site\n";}
	elsif($setrate==2){$out="BU_i=1/rate_i\n";}
	else{
	    $out="";
	    if($setrate!=1){
		$out="WARNING, rate $setrate not allowed, using default 1\n";
	    }
	    $setrate=1;
	    $out=$out."BU_i=rate_i\n";
	}
	print $out; print $fo $out;
	Make_partition_new($MyModel, $lali, \@ali_pdb, \@rate, $subst_mod,
			   $partition_file, \@num_aa, $last_model);
	
    }
    print "End of exchangeability computations\n";
} 

 end_SSCPE:
 

# Run RAxML-NG
if($runraxml){
    $out="\nRunning RAxML-NG Input name: ${msa_new}\n"; 
    print $out; print $fo $out; 
    die;

    my $branches=2*$ntaxa-3;
    $out="$ntaxa taxa, $branches branches, $lali positions\n";
    print $out; print $fo $out;
    
    my $raxml_option="";
    if(-e $option){
	open(my $fh, '<:encoding(UTF-8)', $option);
	$raxml_option = <$fh>;
	chomp($raxml_option);
	close $fh;
	$out="Reading other RAxML options in ${option}, found:\n".
	    "${raxml_option}\n\n";
	print $out; print $fo $out;
    }
    
    # Launch get_br_len for tracking last tree
    my $command;
    if(0){
	$command="$script_get_br -ali $msa_new &";
	print $command,"\n"; print $fo $command,"\n";
	system($command);
    }

    # Launch raxml
    if($evaluate==0){
	$command="$raxmlng --search";
    }else{
	$command="$raxmlng --evaluate --tree $tree";
    }
    if($thread){$command=$command." --thread $thread";}
    $command=$command." --msa $msa_new --model $partition_file";
    if($gamma){$command="$command --force";}

    # Rescale the branches of partitions
    if($setrate==0){
        # All part. have equal branch scale
    	$command="$command --brlen linked";
    }else{
	# Use provided branch scales +BU Do not optimize branch scales
	$command="$command --brlen scaled";
	$command="$command --opt-model off";
    }
    if(0){$command="$command --blopt nr_safe";}
    if($raxml_option){$command="$command $raxml_option";}
    print $command,"\n"; print $fo $command,"\n";
    system($command);

    #`$raxml --search --msa $msa_new --model JTT +FC +G --opt-model off`;
    # For JTT

    my $raxmlout="${msa_new}.raxml.log";
    #my $raxmlout=`ls -1 ${msa_new}.*.log`; chomp($raxmlout);
    if(-e $raxmlout){
	$out="\nReading RAxML-NG results in ${raxmlout}\n";
	print $out; print $fo $out;
	$fh=open(my $fh, '<', $raxmlout);
	my $LogLik=0;
	while (my $row = <$fh>){
	    if(substr($row, 0, 12) eq "Final LogLik"){
		my @word=split(/\s+/, $row);
		$LogLik=$word[2]; last;
	    }
	}
	close $fh;
	$out=sprintf("LogLikelihood= %.3f /(n*branches): %.5g\n",
		     $LogLik, $LogLik/($lali*$branches));
    }else{
	$out="WARNING RAxML-NG results not found in ${raxmlout}\n";
    }
    print $out; print $fo $out;

    if($clean){
	if($print_exch==0){
	    `rm -f $MyModel*.site_*.txt`;
	    `rm -f $partition_file`;
	}
	my $tmp="${msa_new}.raxml.rba"; `rm -f $tmp`;
	$tmp="${msa_new}.raxml.reduced.*"; `rm -f $tmp`;
	$tmp="${msa_new}.raxml.startTree"; `rm -f $tmp`;
	$tmp="${msa_new}.raxml.bestModel"; `rm -f $tmp`;
	#$tmp="${msa_new}.raxml.mlTrees"; `rm -f $tmp`;
	$tmp="${msa_new}.raxml.log"; `rm -f $tmp`;
	$tmp="${msa_new}.fasta"; `rm -f $tmp`;
    }

}
print $fo $buffer;
print "log information written in $nameout\n";
close $fo;
exit;

sub help{

    print "ERROR, PDB file and alignment file must be specified\n";
    print
	"USAGE: ",$0,"\n",
	" -ali <MSA file> (FASTA)\n",
	" -pdblist <list of PDB files> line: pddbcode chain\n",
	" -pdb <single PDB file>\n",
	" -chain <single PDB chain> (default: first chain)\n",
	" -pdbdir <path to PDB> (optional)\n",
	" -model <MOD> (allowed: MF WT DE RMSD DEWT RMSDWT DEMF RMSDMF default: $model)\n",
	" -raxml ! infer ML tree running RAxML-NG\n",
	" -evaluate ! Compute with RAxML-NG branch lengths for given tree\n",
	" -tree  <tree to be evaluated>\n",
	" -print_exch ! Prints and does not delete exchangeability matrices\n",
	" -index <index of output files (optional)>\n",
#	" -iwt   ! RMSD of a_wt: mean over lowest iwt ones (def $iwt)\n",
	" -nostr  (Do not run TNM if its results do not exist)\n",
	" -noreuse  (Do not reuse Prot_evol results even if exist)\n",
	" -noreuse2 (Do not reuse Prot_evol & TNM results even if exist)\n",
	" -option <one line file with other RAxML options> (optional)\n",
	" -thread (number of processors used by RAxML-NG, default $thread)\n",
	" -noclean (Do not clean RAxML partitions after use)\n",
	"Options fot the exchangeability matrix:\n",
	" -matrix <global subst. matrix: LG WAG JTT OPT (default $matrix)>\n",
	" -rate <1,0> 1: BU=site-spec subst.rate; 0: BU=1;\n",
	" -nohb   (Do not include fixation prob. with Halpern-Bruno formula)\n",
	" -noflux (Do not use flux of emp. model instead of exchangeability)\n",
#	" -gamma ! use +G option in RAxML-NG\n",
#	" -allmut ! Comp. deformation of all possible mutations (DE2,RMSD2)\n",
	"\n";
    print "Allowed models (-model):";
    foreach my $mod (@MODELS){print " ",$mod;}
    print "\n";

    die;
}

sub Get_file{ #($dirstore, ".tnm.summary.")
    my ($start, $mid) = @_;
    my $tmp=`ls -1 ${start}*${mid}*`; chomp($tmp);
    if(substr($tmp,0,3) eq "ls:"){return("");}
    my @word=split(/\//, $tmp); my $w=scalar(@word);
    #print "File: ",$tmp,"\n"; 
    #print $w; if($w>1){print " ", $word[$w-1];} print "\n";
    if($w>1){return($word[$w-1]);}
    else{return("");}
}

sub Identical_partitions{ # $model $length $gamma, $partition_file

    my $model=$_[0]; my $lali=$_[1];
    my $gamma=$_[2]; my $output=$_[3]; 
    print "Making partitions for columns 1-",$lali,
    ", output: ",$dirstore,"/",$output,"\n";
    open (my $fo,'>',$output);
    for (my $s=1; $s<=$lali; $s++){
	my $out=$model;
	if($gamma){$out="${out}+G";}
	$out="${out}, p${s} = ${s}-${s}\n";
	print $fo $out;
    }
    close $fo;
    return;
}

sub Make_partition{ # $MyModel $length, \@ali_pdb, $partition_file

    my ($model, $lali, $ali, $rate, $Emp_matr, $output) = @_;
    print "Making partitions for columns 1-",$lali,
    ", output: ",$dirstore,"/",$output,"\n";
    open (my $fo,'>',$output);
    for (my $i=0; $i<$lali; $i++){
	my $out; my $s=$i+1; my $b; my $spdb=$ali->[$i];
	if($spdb>=0){
	    if($setrate==0){$b=1;}
	    elsif($setrate==1){$b=$rate->[$spdb];}
	    elsif($setrate==2){$b=1/$rate->[$spdb];}
	    $spdb++;
	    $out="PROTGTR{${model}.site_${spdb}.txt}";
	    if($gamma){$out="${out}+G";}
	    if($setrate){$out=sprintf("%s+BU{%.3g}", $out, $b);}
	}else{
	    $out="$Emp_matr";
	}
	$out="${out}, p${s} = ${s}-${s}\n";
	print $fo $out;
    }
    close $fo;
    return;
}

sub Make_partition_new{ 
    # The nomodel partition contains all sites without structural model
    my ($model, $lali, $ali, $rate, $Emp_matr, $output, $num_aa, $last_model)
	= @_;
    my $threshold=100; # Minimum number of aa for the option +F
    print "Making partitions for columns 1-",$lali,
    ", output: ",$dirstore,"/",$output,"\n";
    my $nomodel="nomodel="; my $ini_nomod=-1;
    my $num_aa_nomod=0; my $print_nomod=0;
    open (my $fo,'>',$output);
    for (my $i=0; $i<$lali; $i++){
	my $out; my $s=$i+1; my $spdb=$ali->[$i]; my $b=1;
	if($spdb>=0){
	    if($setrate==1){$b=$rate->[$spdb];}
	    elsif($setrate==2){$b=1/$rate->[$spdb];}
	    $spdb++;
	    $out="PROTGTR{${model}.site_${spdb}.txt}";
	    if($gamma){$out="${out}+G";}
	    if($setrate){$out=sprintf("%s+BU{%.4g}", $out, $b);}
	    $out="${out}, p${s} = ${s}-${s}\n";
	    print $fo $out;
	    if($ini_nomod>=0){
		if($print_nomod){$nomodel=$nomodel.",";}
		else{$print_nomod=1;}
		$nomodel=$nomodel."${ini_nomod}-$i";
		$ini_nomod=-1;
	    }
	}else{
	    if($Join_nomod){
		$num_aa_nomod+=$num_aa->[$i];
		if($ini_nomod<0){$ini_nomod=$s;}
	    }else{
		if($Emp_nomod){
		    $out="$Emp_matr";
		    if($setrate){$out=sprintf("%s+BU{1.0}", $out);}
		}else{
		    my $spdb=$last_model; my $b;
		    if($setrate==0){$b=1;}
		    elsif($setrate==1){$b=$rate->[$spdb];}
		    elsif($setrate==2){$b=1/$rate->[$spdb];}
		    $spdb++;
		    $out="PROTGTR{${model}.site_${spdb}.txt}";
		    if($setrate){$out=sprintf("%s+BU{%.4g}", $out, $b);}
		}
		$out="${out}, p${s} = ${s}-${s}\n";
		print $fo $out;
	    }
	}

    }
    if($Join_nomod && $num_aa_nomod){
	if($Emp_nomod){
	    $out="$Emp_matr";
	    if($F_nomod && $num_aa_nomod>$threshold){
		$out=$out."+F";
	    }
	    if($setrate){$out=sprintf("%s+BU{1.0}", $out);}
	}else{
	    my $spdb=$last_model; my $b;
	    if($setrate==0){$b=1;}
	    elsif($setrate==1){$b=$rate->[$spdb];}
	    elsif($setrate==2){$b=1/$rate->[$spdb];}
	    $spdb++;
	    $out="PROTGTR{${model}.site_${spdb}.txt}";
	    if($gamma){$out="${out}+G";}
	    if($F_nomod && $num_aa_nomod>$threshold){$out="${out}+F";}
	    if($setrate){$out=sprintf("%s+BU{%.4g}", $out, $b);}
	}
	if($ini_nomod>=0){
	    if($print_nomod){$nomodel=$nomodel.",";}
	    $nomodel=$nomodel."${ini_nomod}-$lali";
	}
	$out=$out.", $nomodel\n";
	print $fo $out;
    }
    close $fo;
    return;
}

sub FastaManyLines_FastaOneLine{ #$_[0]=Input $_[1]=Output

    print "Converting MSA to one line format\n";

# opening lines of input files
    my $Input = $_[0];
    unless (open(FROM,$Input)){
	print STDERR "Cannot open file \"$Input\"\n\n";
    }

# Reading input file 
    print "Loading $Input\n";

    my $newfile_tree = $_[1];
    open (FILE_OUTPUT,">$newfile_tree");

    my $NumberSites = 0;
    my @AALong_SequenceHere1;
    my @AALong_HeadSeqHere1;

    my $SeqNumber = 0;    
    while (<FROM>){
	if ($_ =~ />/){
	    $SeqNumber++;
	    if ($SeqNumber > 1){print FILE_OUTPUT "\n";}
        
	    my $TaxaName = $_;
	    $TaxaName =~ s/\n//g; # remove the end of line
	    $TaxaName =~ s/\s//g; # remove blank site
	    $AALong_HeadSeqHere1[$SeqNumber] = $TaxaName;
            
	    print FILE_OUTPUT $_;   
        }else{
	    my $sequenceHere = $_;
	    $sequenceHere =~ s/\n//g; # remove the end of line
	    $sequenceHere =~ s/\s//g; # remove blank site
	    $NumberSites = length($sequenceHere);
	    $AALong_SequenceHere1[$SeqNumber] = $sequenceHere;
	    print FILE_OUTPUT "$sequenceHere";
        }		
    }
    print FILE_OUTPUT "\n";
    close (FROM);
    close (FILE_OUTPUT);
    print "Finished\n\n";
    return;
}
	
sub RemakeAlnsWithoutGapsFirstSeq_Ugo{

# It does not assume that the fasta file contains only one line per sequence
    my $AAFileName1=$_[0];
    print "Remaking alignment $_[0] without columns with gaps in first seq.\n";
    
# opening lines of input files
    unless (open(FROM,$AAFileName1)){
	print STDERR "Cannot open file \"$AAFileName1\"\n\n"; die;
    }

    my $newfile_ALN = $_[1];
    open (FILE_OUTPUT,">$newfile_ALN");

    my @gapPos;
    my $NumberGapsPos = 0;
    my $SeqNumber = 0;
    my $num_aa=0;
    open(FROM,$AAFileName1);
    while (<FROM>){ # for each line
	if ($_ =~ />/){
	    if($SeqNumber){print FILE_OUTPUT "\n";}
	    $SeqNumber++;
	    $num_aa=0;
	    print FILE_OUTPUT $_;
	}else{
# identify sites with gaps and print only the sequence without gaps
	    my $ThisSeq_1 = $_;
	    $ThisSeq_1 =~ s/\n//g; # remove the end of line
	    $ThisSeq_1 =~ s/\s//g; # remove blank site
	    my @AASeq = split ('', $ThisSeq_1);
	    if($SeqNumber==1){
		# Identify columns with gaps
		my $l=length($ThisSeq_1); $num_aa+=$l;
		for (my $AAn = 0; $AAn < $l; $AAn++){
		    if ($AASeq[$AAn] ne "-"){$gapPos[$AAn]=0;}
		    else{$gapPos[$AAn]=1; $NumberGapsPos++;}
		}
	    }
	    for (my $AAn = 0; $AAn < length($ThisSeq_1); $AAn++){
		if ($gapPos[$AAn]==0){		    
		    print FILE_OUTPUT "$AASeq[$AAn]";
                }
            }
        }

    }
    print FILE_OUTPUT "\n";

    # Print columns with gap
    print "Columns with gaps: ";
    for (my $n = 0; $n < $num_aa; $n++){
	if($gapPos[$n]){print "$n ";}
    }
    print "\n";
    close (FROM);
    close (FILE_OUTPUT);
    return $num_aa;
}

sub RemakeAlnsWithoutGapsFirstSeq{
# It assumes that the fasta file contains only one line per sequence

    my $AAFileName1=$_[0];
    print "Remaking alignment $_[0] without columns with gaps in first seq.\n";

    my $newfile_ALN = $_[1];
    open (FILE_OUTPUT,">$newfile_ALN");
    
# opening lines of input files
    unless (open(FROM,$AAFileName1)){
	print STDERR "Cannot open file \"$AAFileName1\"\n\n"; die;
    }

    my $col_number = -1;
    my $aa_number1=0;
    my $SeqNumber = 0;
    
    my @gapPos;
    my $NumberGapsPos = 0;
    my $counterLines=0;
    while (<FROM>){
	$counterLines++;
	if($_ =~ />/){
	    print FILE_OUTPUT $_; $SeqNumber++;
	}else{
	    my $ThisSeq_1 = $_;
	    $ThisSeq_1 =~ s/\n//g; # remove the end of line
	    $ThisSeq_1 =~ s/\s//g; # remove blank site
	    my @AASeq = split ('', $ThisSeq_1);

	    if($SeqNumber == 1){
		# identify columns with gaps and print only col without gaps
		$col_number = length($ThisSeq_1);
		print "Total number of columns: $col_number \n";
		
		for (my $AAn = 0; $AAn < $col_number; $AAn++){
		    if ($AASeq[$AAn] eq "-"){
			$NumberGapsPos++;
			$gapPos[$AAn] = 1;
		    }else{
			$gapPos[$AAn] = 0;
		    }
		}
		$aa_number1 = $col_number-$NumberGapsPos;
		print "Excluding gaps in first sequence: ", $aa_number1,"\n";
	    }
            
	    for (my $AAn = 0; $AAn < $col_number; $AAn++){
                if ($gapPos[$AAn] == 0){
		    print FILE_OUTPUT "$AASeq[$AAn]";
                }
            }
	    print FILE_OUTPUT "\n";
        }
    }

    close (FROM);
    close (FILE_OUTPUT);
    
    print "Columns with gaps: ";
    for (my $n = 0; $n < $col_number; $n++){
	if($gapPos[$n]){print "$n ";}
    }
    print "\n\n";
    return $aa_number1;
}


sub ModelsFromDATtoRAxMLng2{
### Programmer: Miguel Arenas Busto. March, 2019.
### Version 2.0 
###
    print "Converting to RAxML format\n";
    my @filesDat = <*.dat>;
    foreach my $file (@filesDat){ # for each .dat file
	my $newfile_rates = $file;
	$newfile_rates = "${file}.txt";
	open (FILE_OUTPUT,">$newfile_rates");
	my @EntireFile = ();
	my $countLinesEntireFile = 0;
	my $line_position = 0;
	unless (open(FROM,$file)){
	    print STDERR "Cannot open file \"$file\"\n\n";
	}
	while (<FROM>){
	    $countLinesEntireFile++;
	    $EntireFile[$countLinesEntireFile] = $_;
	} 
	# now print the files
	for (my $n = 1; $n <= $countLinesEntireFile; $n++){
	    my $sequenceLineHere = $EntireFile[$n];
	    if ($n < 22){
		print FILE_OUTPUT $sequenceLineHere;
	    }
	}      
	close (FILE_OUTPUT);
	close (FROM);
    }
    #system ("rm *.dat");
    print "Finished\n\n";
    return;
}

sub Get_files_TNM{ #0=$dirstore
    my ($workdir, $pdbcode) = @_;
    my $start=$workdir."/";
    if($pdbcode){$start=$start.$pdbcode;}
    $file_summ_tnm=Get_file($start, ".summary.");
    $file_RMSD=Get_file($start, "mut_RMSD.dat");
    $file_DE=  Get_file($start, "mut_DE.dat");
    if($allmut){
	$file_RMSD_all=Get_file($start, "mut_RMSD_all.dat");
	$file_DE_all=  Get_file($start, "mut_DE_all.dat");
    }
}

sub Get_files_Prot_evol{ #0=$dirstore 1=model
    my ($workdir, $pdbcode, $model) = @_;
    my $start=$workdir."/";
    if($pdbcode){$start=$start.$pdbcode;}
    $file_summ_Prot_evol=Get_file($start, "SSCPE.summary.dat");
    $file_freq_SSCPE=Get_file($start, ".${model}.AA_profiles");
    #$file_exch_SSCPE=Get_file($start, ".${model}.exchangeability_sites");
    #$file_rate_SSCPE=Get_file($start, ".${model}.rate_profile");
    return;
}

sub Run_Prot_evol{
    my ($pdb, $chain, $pdblist, $pdbdir, $msa) = @_;

    my $dir_tmp;
    if($pdb){$dir_tmp="tmp_$pdbcode";}
    else{$dir_tmp="tmp_list";}
    
    unless(-d $dir_tmp){
	`mkdir $dir_tmp`;
	$out="Creating folder ${dir_tmp}\n";
	print $out; $buffer=$buffer.$out;
    }
    if($pdblist ne ""){`cp $pdblist $dir_tmp`;}

    chdir $dir_tmp;
    my $subdir=Print_dir(1);
    
    # Prepare configuration file for Prot_evol
    my $command="$Prot_evol -ali $msa";
    if($iwt>0){$command=$command." -iwt $iwt";}
    if($pdbdir){$command=$command." -pdbdir $pdbdir";}
    if($nostr){$command=$command." -nostr";}
    my $input="";
    if($pdb){
	#$input=Write_input_Prot_evol($pdb, $msa, $chain, $pdbcode);
	#$command=$command." -file $input";
	$command=$command." -pdb $pdb";
	if($chain){$command=$command." -chain $chain";}
	$command=$command." > ${pdbcode}.SSCPE.log\n";
    }else{
	$command=$command." -pdblist $pdblist";
	$command=$command." > ${pdblist}.SSCPE.log\n";
    }

    $out="Running Prot_evol: $command";
    print $out; $buffer=$buffer.$out;
    `$command`;
    print "Finished\n\n";

    # Store results of Prot_evol
    Move_files_Prot_evol("../");
    if($pdblist){
	`mv *.mut_DE* ../`;
	`mv *.mut_RMSD* ../`;
	`rm -f $pdblist`;
    }
    if($input){`mv $input ../`;}
    chdir "..";
    $subdir=Print_dir(0);

    # Clean
    #`rm -rf $dir_tmp`; print "Deleting dir ",$dir_tmp,"\n";
}

sub Move_files_Prot_evol
{
    my $dir=$_[0];
    `mv *.AA_profile* $dir`;
    `mv *SSCPE.summary* $dir`;
    #`mv *.map $dir`;
    #`mv *.target $dir`;
    #`mv *.exchangeability_sites* $dir`;
    #`mv *.rate_profile* $dir`;    
}

sub Run_TNM{ 

    my $pdb=$_[0];
    my $chain=$_[1];
    my $Mut_para=$_[2];

    unless(-e $tnm){
	print 
	    "ERROR, wrong path to the program tnm $tnm\n",
	    "Probably the folder name $dir_SSCPE in SSCPE.pl is wrong\n",
	    "Please try to rerun the installation script ",
	    "script_install_SSCPE.sh in the folder where you downloaded ",
	    "SSCPE.zip\n";
	die;
    }
    unless(-e $Mut_para){
	print
	    "ERROR, wrong path to the mutation parameters $Mut_para\n",
	    "Probably the folder name $dir_SSCPE in SSCPE.pl is wrong\n",
	    "Please try to rerun the installation script ",
	    "script_install_SSCPE.sh in the folder where you downloaded ",
	    "SSCPE.zip\n";
	die;
    }

    unless(-d $dir_tmp){
	`mkdir $dir_tmp`;  print "Creating folder ",$dir_tmp,"\n";
    }
    chdir $dir_tmp;
    Print_dir(1);

    # Prepare configuration file for TNM
    my $input=Write_input_TNM($pdb, $chain, $Mut_para);
    my $command="$tnm $input > ${pdbcode}.tnm.log\n";
    $out="Running TNM: $command";
    print $out; $buffer=$buffer.$out;
    `$command`; # run the tnm program
    $out="Finished\n";
    print $out; $buffer=$buffer.$out;

    # Store TNM results
    my $name=substr(`ls -1 *summary.dat`,0,5);
    if(substr($name,0,3) eq "ls:"){
	print "ERROR, TNM program could not be run\n"; die;
    }else{
	print "TNM executed for protein chain ",$name,"\n";
    }
    $file_RMSD="${name}.mut_RMSD.dat";
    $file_DE="${name}.mut_DE.dat";
    $file_RMSD_all="${name}.mut_RMSD_all.dat";
    $file_DE_all="${name}.mut_DE_all.dat";
    my $file_tmp="${name}.summary.dat";
    $file_summ_tnm="${name}.tnm.summary.dat";
    `mv $file_tmp ../$file_summ_tnm`;
    `mv $file_RMSD ../`;
    `mv $file_DE ../`;
    if($allmut){
	`mv $file_RMSD_all ../`;
	`mv $file_DE_all ../`;
    }
    `mv $input ../`;

    chdir "..";
    Print_dir(0);
}

sub Remove_extension{
    my @word=split(/\./, $_[0]);
    my $out=$word[0];
    for(my $i=1; $i<(scalar(@word)-1); $i++){
	$out="${out}.$word[$i]";
    }
    return $out;
}

sub File_name{
    my $path=$_[0];
    my @word=split(/\//, $path);
    return($word[scalar(@word)-1]);
}


sub Write_input_Prot_evol{ #($pdb, $msa, $chain, $pdbcode)  
    my ($pdb, $msa, $chain, $pdbcode) = @_;
    my $input="Input_Prot_evol_${pdbcode}.in";

    open(my $fo, '>', $input);
    print $fo "ALI= $msa\n";
    print $fo "PDB= ",$pdb,"\n";
    if($chain){print $fo "CHAIN= ",$chain,"\n";}
    if($file_RMSD){print $fo "STR_MUT= ../",$file_RMSD,"\n";}
    if($file_DE){print $fo "STR_MUT= ../",$file_DE,"\n";}
    if($file_RMSD_all){print $fo "STR_MUT= ../",$file_RMSD_all,"\n";}
    if($file_DE_all){print $fo "STR_MUT= ../",$file_DE_all,"\n";}
    if($iwt>0){print $fo "IWT= $iwt\n";}
    print $fo "MF_COMP=1\n";
    print $fo "MATRIX=$matrix\n";
    print $fo "EXCHANGE=FLUX\n";
    print $fo "MEANFIELD= 1\n";
    print $fo "PRINT_E=1\n";
    print $fo "OPT_REG= 0\n";
    print $fo "SCORE_CV= 0\n";
    print $fo "REG=0.01\n";
    print $fo "REMUT=1\n";
    print $fo "GET_FREQ=3\n";
    print $fo "FILE_STR= ",$dir_SSCPE,"/structures.in\n";
    print $fo "A_LOC= 1            # Use secondary structure propensities?\n";
    print $fo "TEMP=  1.0	   # Temperature\n";
    print $fo "SU1=   0.13	   # config. entropy per res (unfold)\n";
    print $fo "SC1=  0.065	   # config. entropy per res (misfold)\n";
    print $fo "SC0=  0.0	   # config. entropy offset (misfold)\n";
    print $fo "REM=   2	           # Use 0,1,2 moments of misfolding energy\n";
    close $fo;
    return($input);
}

sub Write_input_TNM{ #($pdb, $chain, $Mut_para)

    my $pdb=$_[0];
    my $chain=$_[1];
    my $Mut_para=$_[2];

    my $input="Input_TNM_${pdbcode}.in";
    open(my $fo, '>', $input);
    print $fo "PDB1= ",$pdb,"\n";
    if($chain){print $fo "CH1= ",$chain,"\n";}
    print $fo "FIT_B=1           ! Fit force constant from B factors\n";
    if($allmut) {
	print $fo "PRED_MUT=2    ! Predict RMSD of all possible mutations\n";
    }else{
	print $fo "PRED_MUT=1    ! Predict RMSD of mutations from WT\n";
    }
    #if($iwt>0){print $fo "IWT= $iwt\n";}
    print $fo "MUT_PARA=",$Mut_para,"   ! File with mutation parameters\n";
    close $fo;
    return($input);
}

sub Print_dir{
    my $dir=`pwd`; chomp($dir);
    if($_[0]){
	my $out="Working directory: ${dir}\n";
	print $out; #print $fo $out;
    }
    return($dir);
}

sub Read_summary{
    my $npdb=-1; my $mat=""; my $file="";
    my $input=$_[2];
    open(my $fh, '<', $input)
	or die "Could not open file '$input' $!";
    while (my $row = <$fh>){
	if(substr($row,0,17) eq "# PDB seq in MSA:"){
	    chomp($row); $npdb=substr($row,18)-1;
	    my $out="PDB sequence is present in MSA at line $npdb\n";
	    print $out; $buffer=$buffer.$out;
	}elsif(substr($row,0,9) eq "# Mapping"){ # Mapping written in 
	    chomp($row); $file=substr($row, 21);
	    my $out="Mapping between MSA and PDB read in $file\n";
	    print $out; $buffer=$buffer.$out;
	}elsif(substr($row,0,9) eq "# Adopted"){
	    chomp($row); my @word=split(/\s+/, $row);
	    $mat=$word[4];
	    my $out="Substitution matrix $mat read in $input\n";
	    print $out; $buffer=$buffer.$out;
	    last;
	}
    }
    close $fh;
    $_[0]=$file;
    if($mat eq "LG" || $mat eq "WAG" || $mat eq "JTT"){
	$_[1]=$mat;
    }elsif($matrix eq "LG" || $matrix eq "WAG" || $matrix eq "JTT"){
	$_[1]=$matrix;
    }else{
	print "ERROR, could not find substitution matrix\n"; die;
    }
    if($npdb<0){
	my $out="WARNING, index of PDB sequence in MSA not found in $_[0]\n";
	print $out; $buffer=$buffer.$out;
    }
    return($npdb);
}

sub Read_subst_mod{
    my $input=$_[0];
    my $mat="";
    open(my $fh, '<', $input)
	or die "Could not open file '$input' $!";
    while (my $row = <$fh>){
	if(substr($row,0,9) eq "# Adopted"){
	    chomp($row); my @word=split(/\s+/, $row);
	    $mat=$word[4];
	    my $out="Substitution matrix $mat read in $input\n";
	    print $out; $buffer=$buffer.$out;
	    last;
	}
    }
    close $fh;
    if($mat eq "LG" || $mat eq "WAG" || $mat eq "JTT"){
	return($mat);
    }elsif($matrix eq "LG" || $matrix eq "WAG" || $matrix eq "JTT"){
	return($matrix);
    }else{
	print "ERROR, could not find substitution matrix\n"; #die;
	return("WAG");
    }
}

sub Read_model{
    my $input=$_[0];
    my $model=""; my $HB=-1; my $FL=1; my $rate=1;
    open(my $fh, '<', $input)
	or die "Could not open file '$input' $!";
    print "Reading model in $input\n";
    while (my $row = <$fh>){
	if(substr($row,0,12) eq "# Best model"){
	    chomp($row); my @word=split(/\s+/, $row);
	    $model=$word[3]; $HB=$word[5]; $FL=$word[7]; $rate=$word[9];
	    my $out="Best model $model HB= $HB flux= $FL rate= $rate ".
		"read in $input\n";
	    print $out; $buffer=$buffer.$out;
	    last;
	}
    }
    close $fh;
    if($HB<0){
	print "ERROR, Best model not found in $input\n"; die;
    }
    return($model, $HB, $FL, $rate);
}

#my @seq=
sub Remove_all_gaps{ #$ntaxa, $lali, $msa_new);
    my $file=$_[0];

    my @aa; my @seq; my @names;
    open(my $fh, '<', $file)
	or die "Could not open file '$file' $!";
    my $lali=0; my $n=-1;
    while (my $row = <$fh>){
	chomp($row);
	if(substr($row,0,1) eq "#"){
	    next;
	}elsif(substr($row,0,1) eq ">"){
	    $n++; $seq[$n]=""; $names[$n]=substr($row,1);
	}elsif($n>=0){
	    $seq[$n]=$seq[$n].$row;
	}
    }
    close $fh;
    my $ntaxa=$n+1;

    for(my $n=0; $n<$ntaxa; $n++){
	my $ln=length($seq[$n]); 
	if($lali==0){
	    $lali=$ln;
	    for(my $i=0; $i<$lali; $i++){$aa[$i]=0;}
	}elsif($ln != $lali){
	    print "ERROR in MSA ",$file," different number of columns ",
	    $ln, " at line ",$n+1," != ",$lali,"\n"; die;
	}
	for(my $i=0; $i<$ln; $i++){
	    my $a=substr($seq[$n],$i,1); my $isaa=0;
	    foreach my $s (@symbols){
		if($a eq $s){$isaa=1; last;}
	    }
	    if($isaa){$aa[$i]++;}
	    elsif($a eq "X"){substr($seq[$n],$i,1)="A"; $aa[$i]++;}
	    else{substr($seq[$n],$i,1)="-";}
	}
    }
    my $allgaps=0;
    for(my $i=0; $i<$lali; $i++){if($aa[$i]==0){$allgaps++;}}
    my $out="Found $ntaxa sequences of length $lali and $allgaps columns".
	" with all gaps in $file\n";
    print $out; $buffer=$buffer.$out;

    if($allgaps){
	$out="Removing $allgaps columns with all gaps in $file\n";
	print $out; $buffer=$buffer.$out;
	open(my $fo, '>', $file);
	for(my $n=0; $n<$ntaxa; $n++){
	    print $fo ">",$names[$n],"\n";
	    my $k=0;
	    for(my $i=0; $i<$lali; $i++){
		if($aa[$i]){
		    my $a=substr($seq[$n],$i,1);
		    print $fo $a;
		    if($k!=$i){substr($seq[$n],$k,1)=$a;}
		    $k++;
		}
	    }
	    print $fo "\n";
	}
	close $fo;
	my $newfile=Remove_extension($file)."_nogap.fas";
	`cp $file $newfile`;
	my $k=-1;
	for(my $i=0; $i<$lali; $i++){
	    if($aa[$i]){
		$k++; if($k!=$i){$aa[$k]=$aa[$i];}
	    }
	}
    }

    $lali-=$allgaps;
    return($ntaxa, $lali, \@seq, \@aa);
}

sub Get_matrix{
    my $default="LG";
    my @word=split(/_/, $_[0]);
    for(my $i=0; $i<scalar(@word); $i++){
	if($word[$i] eq "LG" ||
	   $word[$i] eq "JTT"||
	   $word[$i] eq "WAG"){
	    return($word[$i]);
	}
    }
    return($default);
}


sub Read_subst_model{ #($subst_mat)
    my $input=$_[0];
    open(my $fh, '<', $input) or die "Could not open file '$input' $!";
    print "Reading substitution model in $input\n";
    my @f; my @exch;
    for(my $a=0; $a<20; $a++){
	$f[$a]=0; for(my $b=0; $b<20; $b++){$exch[$a][$b]=0;}
    }
    my $a=1;
    while (my $row = <$fh>){
	if(substr($row,0,1) eq "#"){next;}
	chomp($row); my @word=split(/\s+/, $row);
	if($a<20){
	    for(my $b=0; $b<$a; $b++){
		$exch[$a][$b]=$word[$b];
		$exch[$b][$a]=$word[$b];
		print $exch[$a][$b]," ";
	    }
	    print "\n";
	}elsif($a==21){
	    for(my $b=0; $b<20; $b++){$f[$b]=$word[$b]; print $f[$b]," ";}
	    print "\n";
	}
	$a++;
    }
    return(\@f, \@exch);
}

sub Compute_print_exch {
    my ($model, $len_model, $freq, $fmut, $subst_mod, $HB, $flux) = @_;
    print "Computing substitution models emp= $subst_mod HB= $HB flux= $flux\n";
    my @rate;

    my $subst_mat=""; 
    if($subst_mod eq "LG"){$subst_mat=$mat_LG;}
    elsif($subst_mod eq "WAG"){$subst_mat=$mat_WAG;}
    elsif($subst_mod eq "JTT"){$subst_mat=$mat_JTT;}
    else{
	print "ERROR, empirical model $subst_mod not defined\n"; die;
    }
    my ($f, $exch) = Read_subst_model($subst_mat);

    my @f_emp = (@$f); my @exch_emp = (@$exch);
    my @flux_emp=Normalize_model(\@f_emp, \@exch_emp);

    # Not modelled sites have f_i=f_mut
    for(my $a=0; $a<20; $a++){$freq->[$len_model][$a]=$fmut->[$a];}
    my $L_model=$len_model+1;

    my @exch_site;
    if($HB == 0){
	if($flux){
	    @exch_emp = Flux_matrix(\@flux_emp, $freq, $L_model);
	}
    }else{
	my @P_fix = P_fix($freq, \@f_emp, $L_model); #$fmut
	if($flux){
	    my @exch_flux=Flux_matrix_HB(\@flux_emp,$freq,\@P_fix,$L_model);
	    @exch_site=Multiply_matrix(\@exch_flux, \@P_fix, $L_model);
	}else{
	    @exch_site=Multiply_matrix(\@exch_emp, \@P_fix, $L_model);
	}
    }

    for(my $i=0; $i<$len_model; $i++){
	# name of model: site_i+1.txt Corresponding to column: ali_pdb[c]=i
	my $output = sprintf("%s.site_%d.txt", $model, $i+1);
	$rate[$i]=Print_model($output,$freq,\@exch_emp,\@exch_site,$i,$HB);
    }

    return(@rate);
}

sub Compute_rate{ #(\@f_i, \@exch)
    my ($f_i, $exch) = @_;
    my $rate=0;
    for(my $a=1; $a<20; $a++){
	for(my $b=0; $b<$a; $b++){
	    $rate += $f_i->[$a]*$f_i->[$b]*$exch->[$a][$b];
	}
    }
    return(2*$rate);
}

sub Normalize_model{ #(\@f_emp, \@exch_emp)
    my ($f, $exch) = @_;
    my @flux_mat;
    my $rate=0;
    for(my $a=1; $a<20; $a++){
	for(my $b=0; $b<$a; $b++){
	    $flux_mat[$a][$b] = $f->[$a]*$f->[$b]*$exch->[$a][$b];
	    $rate +=$flux_mat[$a][$b];
	}
    }
    $rate*=2;
    for(my $a=0; $a<20; $a++){
	$flux_mat[$a][$a]=0;
	for(my $b=0; $b<$a; $b++){
	    $exch->[$a][$b]/=$rate;
	    $exch->[$b][$a]=$exch->[$a][$b];
	    $flux_mat[$a][$b]/=$rate;
	    $flux_mat[$b][$a]=$flux_mat[$a][$b];
	}
    }
    return(@flux_mat);
}

sub Flux_matrix{ #(\@flux_emp, $freq, $len_model)

    my ($flux, $freq, $len_model) = @_;
    my @exch_flux;
    for(my $a=0; $a<20; $a++){
	$exch_flux[$a][$a]=0;
	for(my $b=0; $b<$a; $b++){
	    my $sum=0;
	    for(my $i=0; $i<$len_model; $i++){
		$sum += $freq->[$i][$a]*$freq->[$i][$b];
	    }
	    $sum/=$len_model;
	    $exch_flux[$a][$b]=$flux->[$a][$b]/$sum;
	    $exch_flux[$b][$a]=$exch_flux[$a][$b];
	}
    }
    return(@exch_flux);
}

sub Flux_matrix_HB{ #(\@flux_emp, $freq, \@exch_site, $len_model)

    my ($flux, $freq, $exch_site, $len_model) = @_;
    my @exch_flux;
    for(my $a=0; $a<20; $a++){
	$exch_flux[$a][$a]=0;
	for(my $b=0; $b<$a; $b++){
	    my $sum=0;
	    for(my $i=0; $i<$len_model; $i++){
		$sum += $freq->[$i][$a]*$freq->[$i][$b]
		    *$exch_site->[$i][$a][$b];
	    }
	    $sum/=$len_model;
	    $exch_flux[$a][$b]=$flux->[$a][$b]/$sum;
	    $exch_flux[$b][$a]=$exch_flux[$a][$b];
	}
    }
    return(@exch_flux);
}
 
sub Multiply_matrix{
    my ($exch_glob, $P_fix, $len_model) = @_;
    my @exch_i;
    for(my $i=0; $i<$len_model; $i++){
	for(my $a=0; $a<20; $a++){
	    $exch_i[$i][$a][$a]=0;
	    for(my $b=0; $b<$a; $b++){
		$exch_i[$i][$a][$b]=
		    $exch_glob->[$a][$b]*$P_fix->[$i][$a][$b];
		$exch_i[$i][$b][$a]=$exch_i[$i][$a][$b];
	    }
	}
    }
    return(@exch_i);
}

sub P_fix{ #($freq, $fmut, $len_model)
    my ($freq, $P_mut, $len_model) = @_;
    my $P_MIN=0.001; my $log_MIN=log($P_MIN);
    my @P_fix;

    for(my $i=0; $i<$len_model; $i++){
	my @P_sel; my @log_P_sel;
	for(my $a=0; $a<20; $a++){
	    $P_sel[$a]=$freq->[$i][$a]/$P_mut->[$a];
	    if($P_sel[$a]>$P_MIN){$log_P_sel[$a]=log($P_sel[$a]);}
	    else{$log_P_sel[$a]=$log_MIN;}
	}
	for(my $a=1; $a<20; $a++){
	    # Fixation probability from a to b
	    my $Fix;
	    $P_fix[$i][$a][$a]=1;
	    for($b=0; $b<$a; $b++){
		if($P_sel[$a]!=$P_sel[$b]){
		    $Fix=($log_P_sel[$a]-$log_P_sel[$b])/
			($P_sel[$a]-$P_sel[$b]);
		}else{
		    $Fix=1;
		}
		$P_fix[$i][$a][$b]=$Fix;
	    }
	}
    }
    return(@P_fix);
}

sub Print_model{ #$output, \@freq, \@exch_emp, $i
    my ($output, $freq, $exch, $exch_site, $i, $HB) = @_;
    open (my $fo,'>',$output);
    # print "Writing $output\n";
    my $rate=0; my $e;
    for(my $a=1; $a<20; $a++){
	my $out="";
	my $ra=0;
	for(my $b=0; $b<$a; $b++){
	    if($b){$out="${out} ";}
	    if($HB==0){$e=$exch->[$a][$b];}
	    else{$e=$exch_site->[$i][$a][$b];}

	    $out=sprintf("%s%.4f", $out, $e);
	    $ra += $e*$freq->[$i][$b];
	}
	$rate += $freq->[$i][$a]*$ra;
	print $fo $out,"\n";
    }
    my $out="\n";
    for(my $a=0; $a<20; $a++){
	if($a){$out="${out} ";}
	$out=sprintf("%s%.4f", $out, $freq->[$i][$a]);
    }
    print $fo $out,"\n\n";
    close $fo;
    return(2*$rate);
}

sub Read_likelihood_script{
    my ($msa_new, $fo) = @_;
    my $namelog=Remove_extension($msa_new).".log";
    my $raxmllog=$msa_new.".raxml.log";
    print $fo "grep \"Loaded alignment\" $raxmllog >> $namelog\n";
    print $fo "grep LogLik $raxmllog >> $namelog\n";
}

sub Read_frequencies{
    my ($input, $len_model) = @_;
    print "Reading $input\n";
    open(my $fh, '<:encoding(UTF-8)', $input)
	or die "Could not open file '$input' $!";

    my @freq; my @fmut; my $s=0;
    while (my $row = <$fh>){
	if(substr($row, 0, 6) eq "#P_MUT"){
	    my @word=split(/\s+/, $row);
	    for(my $a=0; $a<20; $a++){$fmut[$a]=$word[$a+2];}
	    next;
	}elsif(substr($row, 0, 1) eq "#"){next;}
	my @word=split(/\s+/, $row);
	my $shift; my $i;
	if(substr($row, 0, 1) eq " "){$i=$word[1]; $shift=3;}
	else{$i=$word[0]; $shift=2;}
	#for(my $a=0; $a<20; $a++){$freq[$i][$a]=$word[$a+$shift];}
	for(my $a=0; $a<20; $a++){$freq[$s][$a]=$word[$a+$shift];}
	$s++; 
    }
    print "$s freqs read\n";
    if($s != $len_model){
	print "ERROR, $s sites read but $len_model sites modelled\n"; die;
    }
    close $fh;
    return(\@freq, \@fmut);
}

##############################
# OBSOLETE
############################## 

#my @ali_pdb=
sub Align_to_pdb{ #$lpdb, $ntaxa, $lali, \@seq, $npdb $file_map
    my $ntaxa=$_[1]; my $lali=$_[2];
    my $seq=$_[3]; my $npdb=$_[4];
    my $file_map=$_[5];
    my @ali_pdb;

    my $out="$ntaxa sequences of length $lali in MSA\n";
    print $out; $buffer=$buffer.$out;
    my $read=0; my $ipdb=0;
    
    if($file_map && open(my $fh, $file_map)){
	my $i=0;
	while(my $row = <$fh>){
	    if(substr($row,0,1) eq "#"){next;}
	    chomp($row); if($row>$ipdb){$ipdb=$row;}
	    $ali_pdb[$i]=$row; $i++;
	}
	close $fh;
	if($i!=$lali){
	    $out="ERROR in file ${file_map}, ".
		"read $i columns instead of $lali\n";
	}else{
	    $out="Read mapping between MSA and PDB in $file_map\n";
	    $read=1;
	}
	print $out; $buffer=$buffer.$out;
    }
    if($read){$_[0]=$ipdb+1; return(@ali_pdb);}

    $ipdb=0;
    if($npdb >=0 && $npdb < $ntaxa){
	$out="PDB sequence: $seq->[$npdb]\n";
	print $out; $buffer=$buffer.$out;
	for(my $i=0; $i<$lali; $i++){
	    if(substr($seq->[$npdb],$i,1) ne "-"){
		$ali_pdb[$i]=$ipdb; $ipdb++;
	    }else{
		$ali_pdb[$i]=-1;
	    }
	    #print substr($seq->[$npdb],$i,1);
	}
	$out="contains $ipdb residues\n";
	print $out; $buffer=$buffer.$out;
    }else{
	$out="WARNING, PDB sequence $npdb not found in MSA. ".
	    "I assume it has no gaps\n";
	print $out; $buffer=$buffer.$out;
	$ipdb=$lali;
	for(my $i=0; $i<$lali; $i++){$ali_pdb[$i]=$i;}
    }
    $_[0]=$ipdb;
    return(@ali_pdb);
}

sub Read_rates{ 
    my ($input, $len_model) = @_;
    print "Reading $input\n";
    open(my $fh, '<:encoding(UTF-8)', $input)
	or die "Could not open file '$input' $!";

    my @rate; my $s=0;
    while (my $row = <$fh>){
	if(substr($row, 0, 1) eq "#"){next;}
	my @word=split(/\s+/, $row);
	$rate[$s]=$word[scalar(@word)-1];
	$s++; 
    }
    print "$s rates read\n";
    if($s != $len_model){
	print "ERROR, $s sites read but $len_model sites modelled\n"; die;
    }
    close $fh;
    return(@rate);
}

sub Matrices_FromProtEvol_to_PAML{ #$_[0]=$file_exch_SSCPE $_[1]=$length
    print "Converting Matrix files From ProtEvol to PAML format\n";
    # Author: Miguel Arenas (2014, CBMSO, CSIC-UAM Madrid)"

# Loading input file with matrices from Prot_evol
    my $Model= $_[0]; my $length= $_[1]; my $file = $_[2]; 
    unless (open(FROM,$file)){
	print STDERR "Cannot open file \"$file\"\n\n"; die;
    }
    print "Matrices file uploaded: $file $length sites\n";	

    print "Writing one exchangeability matrix per site: ";
    my $WorkingSite=0;
    my $Site=0;
    my $Line=0;
    open(FROM, $file);
    while (<FROM>){
	#$Lines++;
	# Detect site
	if ($_ =~ /^SITE/){
	    my @seq1 = split(/\s+/, $_);
	    my $PDBSite = $seq1[1];
	    if ($Site){close (FILE_2);}
	    $Site++; print "$Site ";
	    my $createdFile = "${Model}.site_${Site}.txt";
	    open (FILE_2,'>',$createdFile);
	    $Line=0;
	}else{
	    $Line++; 
	    if($Site && $Line<=22){print FILE_2 $_;}
	}	
    } # end of lines FROM
    print "\nModels have been written in $dirstore/$Model\_site_<n>.txt\n";
    print "Finished\n\n";

    close (FILE_2);
    close (FROM);
}
