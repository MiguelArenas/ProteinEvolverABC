#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(getcwd);
use File::Copy;
use File::Basename;
use Config;

# Send to cluster: Modify as needed for your cluster
my $run="perl qsubmit.pl -q x86_64 -n $nthread -s"; # 

# /data/ubastolla/BIN/SSCPE/SSCPE.pl -pdbdir /ngs/databases/pdb/
# -ali $msa -pdblist $pdblist -model $model
# -nohb -norun -rate 1 
#>perl qsubmit.pl -q x86_64 -n 8 -s

my $runmodel=1;
my $reuse=1;
my $examine=0;
my $coeff_slope=0.667; # Ratio between Cv fit and OLS fit
my $exp_br=1; # Mean (1) or mean square (2) branch?
my $nthread=4;

my @iwts=(0); #(5,10,20); #my @iwts=(0,3,5,7,10,20);
my @fluxes=(1,0);
my @Halpern=(1,0);
my @rates=(1,0); # 1,0,2
my @MODELS=("JTT","JTT+F","JTT+G","JTT+F+G",
	    "LG","LG+F","LG+G","LG+F+G",  
	    "WAG","WAG+F","WAG+G","WAG+F+G", # Empirical
	    "MF", "WT",    # STAB "DDG",
	    "DE","RMSD",  # STRUCT
	    "DEWT","RMSDWT", # SSCPE
	    "DEMF","RMSDMF"
	   ); # $allmut==1

my $wmodel="DE"; # First model for whicih iwt is relevant
my $smodel="MF";

my $dir_SSCPE="/data/ubastolla/BIN/SSCPE";
my $Prog_SSCPE="${dir_SSCPE}/SSCPE.pl";
my $Prog_K2="${dir_SSCPE}/k2_mat.pl";
my $Prog_raxml="${dir_SSCPE}/raxml-ng";
#my $script_get_br="$dir_SSCPE/script_get_br_len.pl";

my $pdbdir="";
my $pdblist="";

my $ref_tree="";
my $msa="";
my $pdb="";
my $chain="";
my $ntaxa=0;
my $nbranch=0;
my $nsite=0;

# Results
my @KS0; my @KS1; my @KS2;
my @Kv1; my @Kv21;my @Kv22;
my @RF;  my @support;
my @ll;
my @br_len;
my @models;

print "Script $0 author ",
    "Ugo Bastolla <ubastolla\@cbm.csic.es>\n",
    "It needs an input MSA (-ali), a pdb file either list (-pdblist -pdbdir) ",
    "or single (-pdb -chain) and reference tree (-tree).\n",
    "It runs the program SSCPE.pl, which generates site-specific ",
    "folding-stability and structure constrained substitution models ",
    "and adopts them for inferring a phylogenetic tree with RAxML-NG ",
    "and compares the inferred tree with the reference tree with k2_mat.pl ",
    "for all models and different options\n\n";

# Input arguments
if(scalar(@ARGV)<2){help();}
for(my $i=0; $i<scalar(@ARGV); $i++){
    if($ARGV[$i] eq "-ali"){
	$i++; $msa=$ARGV[$i];
    }elsif($ARGV[$i] eq "-pdbdir"){
	$i++; $pdbdir=$ARGV[$i];
   }elsif($ARGV[$i] eq "-pdblist"){
	$i++; $pdblist=$ARGV[$i];
   }elsif($ARGV[$i] eq "-pdb"){
	$i++; $pdb=$ARGV[$i];
    }elsif($ARGV[$i] eq "-chain"){
	$i++; $chain=$ARGV[$i];
    }elsif($ARGV[$i] eq "-tree"){
	$i++; $ref_tree=$ARGV[$i];
    }elsif($ARGV[$i] eq "-norun"){
	$runmodel=0;
    }elsif($ARGV[$i] eq "-examine"){
	$examine=1;
    }else{
	print "WARNING, unknown option ",$ARGV[$i],"\n";
    }
}
if($pdb eq "" && $pdblist eq ""){
    print "ERROR, PDB file(s) must be provided with option -pdb or -pdblist\n";
    help();
}
if($msa eq ""){
    print "ERROR, alignment file must be provided with option -ali <file>\n";
    help();
}
my $msa_file=File_name($msa);
my @word=split(/\./, $msa_file);
my $ali_name=$word[0];
print "Superfamily $ali_name\n";

my $prot;
if($pdb){
    if($chain eq ""){$chain="A";}
    my $pdb_file=File_name($pdb);
    $prot=Remove_extension($pdb_file)."${chain}";
    if($pdbdir && substr($pdb,0,1) ne "/"){$pdb="${pdbdir}/$pdb";}
}else{
    $prot="pdblist";
    #$name=File_name($pdblist);
    #$name=Remove_extension($name)
}

my @ext=("bestTree","minBrTree");
my $type_max=1; # type_max=2 for computing minBrTree

my $nmod=0;  my $iw=0; 
foreach my $iwt (@iwts){
    $iw++; my $imod=0; my $wmod=0;
    foreach my $model (@MODELS){
	$imod++; my $kmod=0;
	if($model eq $wmodel){$wmod=1;}
	if($wmod==0 && $iw > 1){next;}
	print "Model $model $imod\n";
	$kmod=0;
	my $SSCPE=1; 
	my $mod2=substr($model, 0, 2);
	if($mod2 eq "LG" || $mod2 eq "WA" || $mod2 eq "JT"){
	    $SSCPE=0;
	}
	my $label="${model}";
	if($SSCPE==0){
	    print "kmod= $kmod nmod= $nmod\n";
	    for(my $j=2; $j<=5; $j++){
		if(length($label)<=$j){last;}
		if(substr($label,$j,1) eq "+"){substr($label,$j,1)="_";}
	    }
	    my $msa_new=Remove_extension($msa_file).".${label}.fasta";
	    my $command="--model $model";
	    if($examine){
		Examine_tree($command, $label, $msa_new, $nmod);
	    }else{
		Run_model($SSCPE, $command, $label, $msa_new, $runmodel);
	    }
	    $kmod++; $nmod++;
	    next;
	} 
	# SSCPE model
	foreach my $HB (@Halpern){
	    foreach my $flux (@fluxes){
		#if($HB==0 && $flux==0){next;}
	        foreach my $rate(@rates){
		    if($HB==0 && $flux==0 && $rate==0){next;}
		    #if($HB==0 && $flux==1 && $rate==0){next;}
		    if($HB==0 && $flux==0 && $rate==1){next;}

		    print "kmod= $kmod nmod= $nmod\n";

		    my $command="-raxml -thread $nthread";
		    if($chain){$command="$command -chain $chain";}

		    #my $label="${prot}.${model}";
		    my $label=$model;
		    $command="$command -model $model";
		    if($wmod && $iwt){
			$label="${label}.iwt$iwt";
			$command="${command} -iwt $iwt";
		    }
		    $label="${label}.HB${HB}.FL${flux}.rate${rate}";
		    if($HB==0){
			$command="${command} -nohb";
		    }
		    if($flux==0){
			$command="${command} -noflux";
		    }
		    $command="${command} -rate $rate";
		    my $msa_new=
			Remove_extension($msa_file).".${label}.fasta";

		    my $run_Prot_evol=0;
		    if($reuse==0 && $kmod==0 && 
		       ($model eq $smodel || ($model eq $wmodel && $iw>1))){
			$command="${command} -noreuse";
			$run_Prot_evol=1;
		    }
		    if($examine){
			Examine_tree($command, $label, $msa_new, $nmod);
		    }else{
			Run_model($SSCPE,$command,$label,$msa_new,$runmodel);
			if($runmodel){
			    if($run_Prot_evol){sleep 500; $run_Prot_evol=0;}
			    else{sleep 0.1;}
			}
		    }
		    $kmod++; $nmod++;
		}
	    }
	}
    }
}
if($examine){
    for(my $i=0; $i<$type_max; $i++){
	Print_Results($i, $nmod);
    }
}

sub Examine_tree{
    my ($command, $label, $msa_new, $kmod) = @_;

    my $model=$label;
    my @word=split(/\./, $model);
    if(scalar(@word)>=5){
	$model=$word[1];
	for(my $i=2; $i<scalar(@word); $i++){$model=$model.".$word[$i]";}
    }

    $models[$kmod]=$model;
    $ll[$kmod]=0; $br_len[$kmod]=0;
    $KS0[$kmod]=0; $KS1[$kmod]=0;  $KS2[$kmod]=0;
    $Kv1[$kmod]=0; $Kv21[$kmod]=0; $Kv22[$kmod]=0;
    $RF[$kmod]=0; $support[$kmod]=0;


    my $Diff="Diff.${ali_name}.${label}";
    my $out;

    my $file_log=Remove_extension($msa_new).".log";
    my $raxmllog=$msa_new.".raxml.log";
    if(-e $raxmllog && !-e $file_log){
	$out="grep \"Loaded alignment\" $raxmllog >> $file_log\n".
	    "grep LogLik $raxmllog >> $file_log\n";
	print $out;
	`$out`;
    }

    my $LogLik=0;
    my $file=$file_log;
    my $out_taxa="";
    print "Reading $file\n";
    if(!-e $file){
	print "WARNING, file $file does not exist\n";
    }else{
	open(my $fh, '<', $file);
	while (my $row = <$fh>){
	    my @word=split(/\s+/, $row);
	    my $nw=scalar(@word);
	    if($nw<=1){
		next;
	    }elsif(substr($word[0],0,6) eq "LogLik"){
		$LogLik=$word[1]; last;
	    }elsif($nw>2 && substr($word[1],0,6) eq "LogLik"){
		$LogLik=$word[2]; last;
	    }elsif($nw<5){
		next;
	    }elsif($ntaxa==0){
		if($word[1] eq "taxa,"){
		    $ntaxa=$word[0]; $nsite=$word[4];
		}elsif($nw > 7 && $word[5] eq "taxa"){
		    $ntaxa=$word[4];  $nsite=$word[7];
		}else{
		    next;
		}
	    }
	}
	close $fh;
    }
    if($LogLik >= 0){
	my $tmp_log="tmp.${label}.log";
	print "grep LogLik $tmp_log >> $file\n";
	`grep LogLik $tmp_log >> $file`;
	if(-e $file){
	    open(my $fh, '<', $file);
	    while (my $row = <$fh>){
		my @word=split(/\s+/, $row);
		my $nw=scalar(@word);
		if($nw>1 && substr($word[0],0,6) eq "LogLik"){
		    $LogLik=$word[1]; if($LogLik<0){last;}
		}elsif($nw>2 && substr($word[1],0,6) eq "LogLik"){
		    $LogLik=$word[2]; if($LogLik<0){last;}
		}
	    }
	    close $fh;
	}
    }
    $ll[$kmod]=$LogLik;
 
    for(my $i=0; $i<$type_max; $i++){ # $type_max=2 for computing minBrTree
	my $tree="${msa_new}.raxml.$ext[$i]";
	if(!-e $tree){
	    print "WARNING, tree $tree does not exist\n"; next;
	}
	$br_len[$kmod]=Read_br_len($tree, $exp_br);

	my $Diff_i="${Diff}.${ext[$i]}.txt";
	if(!-e $Diff_i){
	    $out="perl $Prog_K2 -rt $ref_tree -ct $tree ".
		"-r -t $Diff_i\n";
	    print $out;
	    `$out`;
	}
	$file=$Diff_i;
	if(!-e $file){
	    print "WARNING, file $file was not written\n";
	    return;
	}else{
	    print "Reading $file\n";
	    my $l=0;
	    open(my $fh, '<', $file);
	    while (my $row = <$fh>){
		$l++;
		if($l!=2){next;}
		chomp($row);
		my @word=split(/\t/, $row);
		$KS0[$kmod]=$word[1];
		$KS1[$kmod]=$word[2]; $KS2[$kmod]=$word[4];
		$Kv1[$kmod]=$word[3];
		$Kv21[$kmod]=$word[5]; $Kv22[$kmod]=$word[6];
		$RF[$kmod]=$word[7]; $support[$kmod]=$word[8];
#$br_len=$word[9];
		#print "$KS0 $KS1 $KS2 $Kv1 $Kv21 $Kv22 $RF $support\n";
#Ploop_C2 - Ploop_C2	1.91	1.0768	0.5282	0.8445	0.4580	1.0983	2
	    }
	    close $fh;
	}
    }
    print "$kmod\t$models[$kmod]",
    sprintf("\t%.4g",$ll[$kmod]),"\t$KS1[$kmod]\t$KS2[$kmod]\t$RF[$kmod]",
    "\t$Kv1[$kmod]\t$Kv21[$kmod]\t$Kv22[$kmod]\t$support[$kmod]",
    sprintf("\t%.4g\n",$br_len[$kmod]);

}

sub Print_Results{
    my ($type, $nmod)=@_;
    my $name="Results.${ali_name}.${ext[$type]}.txt";
    print "Writing all results in ${name}\n";
    open(my $f_all, '>', $name);

    my $nbr=2*$ntaxa-3; my $norm=$nsite*$nbr;
    my $out_taxa="# ntaxa= $ntaxa nbranch= $nbr sites=$nsite $nmod models\n";
    print $f_all $out_taxa;

    my @x; my @y; my $nn=0;
    for(my $m=0; $m<$nmod; $m++){
	if($ll[$m]==0 || $br_len[$m]==0){next;}
	if($norm){$ll[$m]/=$norm;}
	$x[$nn]=$br_len[$m];
	$y[$nn]=$ll[$m];
	$nn++;
    }
    my ($r, $slope, $offset)=Regression(\@x, \@y, $nn);
    my $out="# Correlation between br.len/nbr and log.lik/(ns*nbr):";
    $out=$out.sprintf(" %.3f", $r)." slope=".sprintf(" %.3g", $slope).
	" offset=".sprintf(" %.3g", $offset)." n= $nn\n";
    print $out;
    print $f_all $out;

    $slope*=$coeff_slope;
    print $f_all "# reg=-log.lik/(ns*nbr)+",
    "${coeff_slope}*fit_slope*Mean_br_len\n";
    print $f_all "# ${coeff_slope}*fit_slope= ",sprintf("%.4g\n", $slope);

    print $f_all
	"#1=model 2=LogLik 3=KS1 4=KS2 5=RF 6=Kv1 7=Kv21 8=Kv22 ",
	"9=support 10=br_len 11=reg\n";
    my ($n_all, $ave_br_len, $n_LL, )=(0, 0, 0);
    my ($ave_LL, $ave_K1, $ave_K2, $ave_RF, $ave_reg)=(0, 0, 0, 0, 0);
    my ($ave_Kv1, $ave_Kv21, $ave_Kv22, $ave_supp)=(0, 0, 0, 0);

    for(my $k=0; $k<$nmod; $k++){
	if($ll[$k]==0 || $br_len[$k]==0){
	    print $f_all "#${models[$k]}\n"; next;
	}
	my $lc=sprintf("%.4g", $ll[$k]);
	my $reg=-$ll[$k]+$slope*$br_len[$k];

	print $f_all "${models[$k]}\t$lc\t$KS1[$k]\t$KS2[$k]\t$RF[$k]",
	"\t$Kv1[$k]\t$Kv21[$k]\t$Kv22[$k]\t$support[$k]",
	sprintf("\t%.4g", $br_len[$k]),
	sprintf("\t%.4g\n", $reg);
	$ave_LL+=$ll[$k]; $n_LL++;
	$n_all++; 
	$ave_K1+=$KS1[$k]; $ave_K2+=$KS2[$k]; $ave_RF+=$RF[$k];
	$ave_Kv1+=$Kv1[$k]; $ave_Kv21+=$Kv21[$k]; $ave_Kv22+=$Kv22[$k];
	$ave_supp+=$support[$k]; $ave_br_len+=$br_len[$k];
	$ave_reg+=$reg;
    }

    if($n_all==0){
	print "WARNING, no results found\n";
    }else{
	my $out="#Mean:";
	if($n_LL){$ave_LL/=$n_LL;}
	$out=$out.sprintf("\t%.4g", $ave_LL);
	$out=$out.sprintf("\t%.4g", $ave_K1/$n_all);
	$out=$out.sprintf("\t%.4g", $ave_K2/$n_all);
	$out=$out.sprintf("\t%.4g", $ave_RF/$n_all);
	$out=$out.sprintf("\t%.4g", $ave_Kv1/$n_all);
	$out=$out.sprintf("\t%.4g", $ave_Kv21/$n_all);
	$out=$out.sprintf("\t%.4g", $ave_Kv22/$n_all);
	$out=$out.sprintf("\t%.4g", $ave_supp/$n_all);
	$out=$out.sprintf("\t%.4g", $ave_br_len/$n_all);
	$out=$out.sprintf("\t%.4g", $ave_reg/$n_all);
	print $f_all $out, "\n";
    }
    close $f_all;
    print "Results for all models printed in file $name\n";
    return;
}

sub Run_model{
    my ($SSCPE, $command, $label, $msa_new, $runmodel) = @_;

    my $Diff="Diff.${ali_name}.${label}";

    print "label= $label\n";
    print "command= $command\n";
    print "msa= $msa_new\n";

    my $file_log=Remove_extension($msa_new).".log";
    my $run_this=1;
    if(-e $file_log){
	print "WARNING, log file $file_log exists, skipping\n";
	$run_this=0;
	#return;
    }

    my $script="tmp.${label}";
    open(my $fo, '>', $script);
    if($SSCPE==0){
	print $fo "cp -f $msa $msa_new\n";
	# Launch program that reads minimum branch length tree
	#print $fo "$script_get_br -ali $msa_new &\n"; # Not needed anymore
	print $fo "$Prog_raxml --msa $msa_new $command\n";
	my $raxmllog=$msa_new.".raxml.log";
	my $out="grep \"Loaded alignment\" $raxmllog >> $file_log\n".
	    "grep LogLik $raxmllog >> $file_log\n";
	#print $out;
	print $fo $out;
	print $fo "perl $Prog_SSCPE -ali $msa";
    }elsif($pdblist){
	 print $fo " -pdblist $pdblist -pdbdir $pdbdir";
    }else{
	print $fo " -pdb $pdb";
    }
    if($SSCPE){print $fo " $command\n";}

    for(my $i=0; $i<$type_max; $i++){
	my $Diff_i="${Diff}.${ext[$i]}.txt";
	my $out="perl $Prog_K2 -rt $ref_tree -ct ${msa_new}.raxml.$ext[$i] ".
	    "-r -t $Diff_i\n";
	#print $fo "echo $out";
	print $fo $out;
    }

    print $fo "rm -f ${msa_new}.raxml.rba\n";
    print $fo "rm -f ${msa_new}.raxml.bestModel\n";
    print $fo "rm -f ${msa_new}.raxml.startTree\n";
    print $fo "rm -f ${msa_new}.raxml.reduced.phy\n";
    #print $fo "rm -f ${msa_new}.raxml.mlTrees\n";

    close $fo;
    `chmod u+x $script`;
    if($run_this && $runmodel){`$run $script`;}
}

sub Remove_extension{
    my @word=split(/\./, $_[0]);
    my $out=$word[0];
    for(my $i=1; $i<(scalar(@word)-1); $i++){
	$out="${out}.$word[$i]";
    }
    return $out;
}

sub help{

    print "ERROR, PDB file and alignment file must be specified\n";
    print
	"USAGE: $0 \n",
	" -ali <MSA file> (FASTA)\n",
	" -pdbdir <Path to PDB>\n",
	" -pdblist <list of PDB files> line: pddbcode chain\n",
	" -pdb   <single PDB file>\n",
	" -chain <single PDB chain> (default: first chain)\n",
	" -tree <Reference tree>\n",
	" -examine <Examine previous results without running SSCPE>\n",
	" -norun <Write the scripts without running them>\n",
	"\n",
	"Examples:\n",
	"perl script_run_SSCPE_models.pl -ali Ploop.PCAli.fas -tree Ploop.PCAli.unroot.tree -pdbdir /ngs/databases/pdb -pdblist Ploop.PCAli.pdb.list\n",
	"perl script_run_SSCPE_models.pl -ali Ploop.PCAli.fas -tree Ploop.PCAli.unroot.tree -pdb /ngs/databases/pdb/2bwj.pdb -chain A\n",
	"perl script_run_SSCPE_models.pl -ali Ploop.PCAli.fas -tree Ploop.PCAli.unroot.tree -pdb /ngs/databases/pdb/2bwj.pdb -chain A -examine\n";
    die;
}

sub File_name{
    my $path=$_[0];
    my @word=split(/\//, $path);
    return($word[scalar(@word)-1]);
}

sub Read_br_len{ # Newick format

    my ($file, $exp) = @_;
    chomp($file);
    if(!-e $file){
	print "WARNING, tree $file does not exist\n";
	return(0);
    }

    print "Reading $file\n";
    my $pr=0;

    open(my $fh, '<', $file);
    my $sum=0; my $nbr=0;
    while (my $row = <$fh>){
	if(substr($row, 0, 1) eq "#"){next;}
	$nbr=0; $sum=0; my $lbr=0; my $ini=-1;
	print "Branches of tree $file: ";
	for(my $l=0; $l<length($row); $l++){
	    my $c=substr($row, $l, 1);
	    if($c eq ":"){
		$ini=$l+1;
	    }elsif($ini>=0 &&
		   ($c eq "," || $c eq "(" || $c eq ")" ||
		    $c eq ";" || $c eq "\n")){
		$lbr=substr($row, $ini, $l-$ini);
		if($exp==2){$lbr*=$lbr;} # squared
		$sum+=$lbr; $nbr++;
		$ini=-1;
		if($pr){print sprintf(" %.5f", $lbr);}
	    }
	}
	$sum/=$nbr; if($exp==2){$sum=sqrt($sum);}
	print " ($nbr branches), ";
	if($exp==2){print "RMS(lbr)= ";}
	else{print "Mean(lbr)= ";}
	print sprintf("%.4g\n",$sum);
    }

    return($sum);
}

sub Regression{
    my ($x, $y, $n)=@_;
    my $x1=0; my $x2=0; my $y1=0; my $y2=0; my $xy=0;
    for(my $i=0; $i<$n; $i++){
	my $xx=$x->[$i]; my $yy=$y->[$i];
	$x1+=$xx; $x2+=$xx*$xx; 
	$y1+=$yy; $y2+=$yy*$yy;
	$xy+=$xx*$yy;
    }
    $x1/=$n; $x2=$x2/$n-$x1*$x1;
    $y1/=$n; $y2=$y2/$n-$y1*$y1;
    $xy=$xy/$n-$x1*$y1;
    my $r=$xy/sqrt($x2*$y2);
    my $slope=$xy/$x2;
    my $offset=$y1-$slope*$x1;
    return($r, $slope, $offset);
}
