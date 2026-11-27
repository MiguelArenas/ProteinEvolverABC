#!/usr/bin/env perl 
use warnings; 
use Cwd 'abs_path';
use File::Basename;

sub help{
	die "Usage: qsubmit.pl [-q queue]  [-g number_of_gpus] [-n numberCores] [--mem GBsMemory] [-m mailaddres] [-S] [-s scriptFile] [-I]\n\t -S: Split the job in 8 cores per node\n\t Default memory is 1.8 GBs per core\n\t-I: Indicates interative job\n";
}

$queue = 'x86_64';
$script = '';
$mail = '';
$cores = 1;
$memory=0;
$def_mem=1845;
$gpu = 0;
$max_np = 1;
$interactive = 0;

%qcores = ();
open (QLIST, "qstat -Q|") or die "qstat not available in this machine";
while (<QLIST>){
    if ($_=~/^----/){
        last;
    }
}
while (<QLIST>){
    if ($_!~/^\s/){
        @parms = split(/\s+/,$_);
        $qcores{$parms[0]}=$parms[1];
    }
}
if ( keys%qcores eq '' ){
    die "No queues available\n";
}
if (@ARGV) {
	while (@ARGV) {
		if ($ARGV[0] eq '-q') {
	 		if ($ARGV[1]) {
	 			shift@ARGV;
	 			$queue = shift@ARGV;
				if ( !defined($qcores{$queue}) ){
					print "Available queues:\n";
					foreach $queue (keys%qcores){
						print "\t$queue\n";
					}
					die " Specified queue '$queue' doesn't exit\n";
				}
	 		} else {
	 			die "A queue name must be provided\n";
	 		}
	 	}elsif ($ARGV[0] eq '-s') {
	 		if ($ARGV[1]) {
				 shift@ARGV;
				 $script=shift@ARGV;
				if ($script =~ /(.*?)\s+/){
					$exe = $1;
				}else{
					$exe = $script;
				}
				if ( ! -r $exe ) {
					die "Script file '$exe' not found\n";
				}
				if ( ! -x $exe ){
					die "Script file '$exe' is not executable, use chmod!\n";
				}
			 } else {
				 die "A script file must be provided\n";
			 }
		}elsif ($ARGV[0] eq '-n') {
			$cores=$ARGV[1];
			if (defined($cores)) { 
				shift@ARGV; 
				$cores=shift@ARGV;
				if ($cores =~ /\D/ ){
					die "Number of cores must be indicated by an integer\n";
				}
				if ($cores == 0){
					die "Number of cores must be >=1\n";
				}
			}else{
				print STDERR "The number of required cores must be specified\n";
				help();
			}
		}elsif ($ARGV[0] eq '--mem') {
			if ($ARGV[1]) { 
				shift@ARGV; 
				$memory=shift@ARGV;
				if ($memory =~ /\D/ ){
					die "GBs of memory for the job must be indicated by an integer\n";
				}
			}
		}elsif ($ARGV[0] eq '-S') {
				shift@ARGV; 
				$max_np=0;
		}elsif ($ARGV[0] eq '-I') {
				shift@ARGV; 
				$interactive=1;
		}elsif ($ARGV[0] eq '-g') {
			$gpu=$ARGV[1];
			if (defined($gpu)) { 
				shift@ARGV; 
				$gpu=shift@ARGV;
				if ($gpu =~ /\D/ ){
					die "Number of gpus must be indicated by an integer\n";
				}
			}else{
				print STDERR "The number of required GPUS must be specified\n";
				help();
			}
	        }elsif ($ARGV[0] eq '--mail') {
			if ($ARGV[1]) { 
				shift@ARGV; 
				$mail=shift@ARGV;
				if ($mail !~ /(.*)\@cbm.uam.es/ and $mail !~ /(.*)\@cbm.csic.es/ ){
					die "Only CBM addresses are allowed\n";
				}
			}
		}else {
		 	print "Unrecognized argument $ARGV[0]\n";
			 help();
		 }
	}
}else {
	help();
}
if ($queue eq '' ){
	print STDERR "Queue not specified. Available queues: ", join(' ',keys%qcores),"\n";
	help();
}
if ( $script eq '' and $interactive == 0){
	print STDERR "Missing script to send to queue\n";
	help();
}
if ( $cores > $qcores{$queue} ){
    die "Queue $queue has a limit of $qcores{$queue} cores\n";
}

$pwd = $ENV{'PWD'};

if ($script){
	# Getting script full path
	@script = split(' ',$script);
	$binary = abs_path(shift(@script));
	$path = dirname($binary);
	unshift(@script, basename($binary));
	$script = join(' ',@script);

	# Getting script parameters
	$scriptparm = $script;
	$scriptparm =~ s/[\s\/\'\"\:,\+\!=]/_/g;
	$scriptparm =~ s/-+/-/g;
	$scriptparm =~ s/_+/_/g;
	$scriptparm =~ s/_-/_/g;
	#$scriptparm =~ s/\//_/g;
	#$scriptparm =~ s/\'/_/g;
	#$scriptparm =~ s/\"/_/g;
	#$scriptparm =~ s/\:/_/g;
	#$scriptparm =~ s/,/_/g;
	$scriptparm =~ s/_+/_/g;
	if (length($scriptparm) > 100){
	    $head = substr($scriptparm,0,75);
	    $tail = substr($scriptparm,-25,25);
	    $scriptparm = "${head}-${tail}";
	}
}
#Calculate number of nodes:
#open (QPARMS,"echo p q $queue|qmgr|") or die "There's no access to qmgr\n";
#$ppn = 0;
#while (<QPARMS>){
#    if ($_=~/resources_default.nodes = (\d+):/){
#        $ppn = $1;
#        last;
#    }
#}

#Obtain type of nodes:
open (QPARMS,"echo p q $queue|qmgr|") or die "There's no access to qmgr\n";
while (<QPARMS>){
	#if ($_=~/resources_default.nodes = (\d+):(.*)/){
	if ($_=~/default_chunk.queue_name = (.*)/){
        	#$type = $1;
        	last;
	}
}

# >>>>> Some x86_64 machines got more than 8 processors

if ($queue eq 'i686'){
	$ppn = 2;
}else{
	if ($max_np) {
		if ( $gpu > 0 ) {
			$ppn = ($cores > 16 )?16:$cores;
		}else{
			$ppn = ($cores > 120)?120:$cores;
		}
	}else {
		$ppn = 8
	}
}
	
#print "Cores per node:$ppn\n";
if ($ppn == 0) {
    die "It hasn't been possible to determine the number of cores per node on queue $queue\n";
}else{
	if ($gpu > 0){
		if ($cores < $gpu){
			die "At least 1 core per GPU must be selected\n";
		}
		$nnodes = $gpu;
		$ppn =int($cores/$nnodes);
		if ($ppn > 16){
			die "Unable to allocate more than 16 cores per GPU\n";
		}
	}else{
	    $nnodes  = int ($cores / $ppn);
	}
	$extracores = $cores - ($nnodes*$ppn);
}
#$nodesstring = 'select=';
#print "ppn: $ppn, nodes: $nnodes, extracores: $extracores\n";
#die;
$nodesstring = "select=${nnodes}:ncpus=$ppn";
#if ($nnodes>1){
	#$nodesstring .=  "${nnodes}:$type:ppn=$ppn";
#	$nodesstring .=  "node=${nnodes}:ppn=$ppn";
#}else{
#	$nodesstring .= "${nnodes}:ncpus=$ppn";
#}
if ($extracores > 0){
	#if ($nnodes){
	    #$nodesstring .= "+1:$type:ppn=$extracores";
	$nodesstring .= "+1:ncpus=$extracores";
	#}else{
	    #$nodesstring .= "1:$type:ppn=$extracores";
	    #$nodesstring .= "1:ncpus=$extracores";
	    #}
	$nnodes++;
}
if ($memory == 0) {
	if ($cores >= 1 and $cores < 12){
		$memory = $def_mem*$cores;
		$nodesstring .= ":mem=${memory}mb";
	}
}else{
	$nodesstring .= ":mem=${memory}gb";
}

if ( $gpu > 0 ) {
    $gpu = 1; #There's only one GPU per host
    $nodesstring .= ":ngpus=$gpu";
    $nnodes=$gpu;
} 
if ($nnodes > 1 or $extracores > 0){
	print STDERR "\nWarning: This job can't be assigned to a single node, splitting into $nnodes nodes\n\n";
}
$host=$ENV{'HOSTNAME'};
if ($host!~/eth$/){
	$host=~s/chiron./chiron1eth/;
}
#$host=~/(.*?)\./;
#$host=$1;
#die "Host: $host\n";
if (!$scriptparm){
	$scriptparm='interactiveJob';
	$script='';
	$path='.';
}
open (PBS,">${scriptparm}.pbs") or die "Can't write PBS script\n";
print PBS 
"#!/bin/bash
# queue name, one of { ", join(',',keys%qcores), " }
#PBS -q $queue\n";

if ( $gpu > 0 ) {
    print PBS "#PBS -l $nodesstring \n";
}
if ($interactive == 1){
	print PBS "#PBS -I\n";
	print PBS "#PBS -l walltime=04:00:00\n";
	print STDERR "\nWarning: Interactive jobs have a time limit of 4h\n\n";
	$interactive='-I -X -V'
}else{
	print PBS "# path/filename for standard error.
#PBS -e $host:$pwd/$scriptparm.error
# path/filename for standard output.
#PBS -o $host:$pwd/$scriptparm.log
	\n";
	$interactive='';
}
if ($mail ne ''){
    print PBS
    "# send me e-mail when job begins (b), when job ends (e) and when job aborts (a)
#PBS -m bea
#PBS -M $mail\n";
    $qsub_mail_flag='';
}else{
    print PBS "#PBS -m n\n";
    $qsub_mail_flag='-m n';
}
print PBS
"# This jobs's working directory            
echo '============= PBS Job on UBCBMSO Cluster ==========='
echo 'Working directory is :' $pwd
cd $pwd
echo 'Job Name             :' \$PBS_JOBNAME
echo 'JobID                :' \$PBS_JOBID 
echo 'Node                 :' `hostname`
echo 'Number of cores      :' \$NCPUS
echo 'Time is              :' `date`
echo 'Actual directory is  :' `pwd`
echo 'Command              :' \"${path}/$script\"
echo '===================================================='

echo \"\`date\`: Started job from \$USER at queue $queue\" >> /data/guest/pbs_stats.txt
START=\$(date +%s)
${path}/$script
END=\$(date +%s)
DIFF=\$(( \$END - \$START ))
echo \"\`date\`: Finished job from \$USER at queue $queue in \$DIFF seconds\" >> /data/guest/pbs_stats.txt\n";

close (PBS);

#$command = "/usr/bin/qsub -l nodes=$nodesstring  $scriptparm.pbs\n";
$command = "/opt/pbs/bin/qsub -l $nodesstring $qsub_mail_flag $interactive $scriptparm.pbs\n";
print "$command";
$ret=system($command);
print $ret;
