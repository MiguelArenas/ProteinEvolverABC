### Prepare_Files_for_ABC.pl
### Prepare files to be uploaded in R in order to perform ABC
###



use strict;
#use warnings;
use File::Copy;
use File::Basename;

my $settings;    #handle to Settings file
my $params_in;   #handle to Params input file
my $sstats_in;   #handle to SStats input file
my $params_out;  #handle to Params output file
my $sstats_out;  #handle to SStats output file



##############################################
##############################################
# Loading input file, "Settings"
my $file = $ARGV[0];
unless (open($settings,$file))  {
	print STDERR "Cannot open file \"$file\"\n\n";
	exit;
}
my $setsdir = dirname($file);
close($settings);


##############################################
##############################################
# Declare fixed variables
my $numberSummStats = 16;
my $numberParams = 2;
my $fileInSummStats = "$setsdir/SSsimulations.ss";
my $fileOutSummStats = "$setsdir/SSsimulations.csv";
my $fileInParams = "$setsdir/ProteinEvolverABC_Pop_arguments.txt"; # let's use pop level
my $fileOutParams = "$setsdir/PSimulations.csv";



##############################################
##############################################
# Check for incomplete Summary Statistics files and then reduce file size

unless (open($params_in, '<', $fileInParams))  {
	print STDERR "Cannot open Summary Statistics file \"$fileInParams\"\n\n";
}
unless (open($sstats_in, '<', $fileInSummStats))  {
	print STDERR "Cannot open Parameters file \"$fileInSummStats\"\n\n";
}
unless (open($params_out,'>',$fileOutParams))  {
	print STDERR "Cannot create Summary Statistics output file \"$fileOutParams\"\n\n";
}
unless (open($sstats_out,'>',$fileOutSummStats))  {
	print STDERR "Cannot create Parameters output file \"$fileOutSummStats\"\n\n";
}
print "Working on files \"$fileInParams\" and \"$fileInSummStats\"\n";

#Print headers to files
if ($numberParams == 2)  {
	print $params_out "Rho,Theta\n";
}
else  {
	print "\n\nERROR!!!. The number of parameters must be 2 ($numberParams)\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
}

if ($numberSummStats == 16)  {
	print $sstats_out "Phi,NSS,ChiSq,p_av,p_sd,p_sk,p_ku,H_av,H_sd,H_sk,H_ku,S,si_av,si_sd,si_sk,si_ku\n";
}
else  {
	print "\n\nERROR!!!. The number of parameters must be 16 ($numberSummStats)\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
}



#Check position of values in files
my @sstatsAux = (0);  #position of the values in the Summary Statistics file
for (my $i = 0; $i < 16; $i++)  {
    $sstatsAux[$i] = $i*2+1;
}
my @paramsAux = (0);  #position of the values in the Summary Statistics file
my $pline = <$params_in>;
my @plist = split(" ", $pline);
my $i = 0;
for my $w (@plist)  {
    if ( grep /^-r/, $w ){
		$paramsAux[0] = $i;
    }
    elsif( grep /^-u/, $w){
    	$paramsAux[1] = $i ;
    }
    #elsif( grep /^-u/, $w){
    #	$paramsAux[2] = $i;
    #}
    $i++
}
seek $params_in, 0, 0;  #   Move pointer back to start of file

print "Simulations with incomplete calculation of Summary Statistics will be skipped..";
#Going through files
my $count = 0;
my $flag = 1;
while(!eof($params_in) and !eof($sstats_in))  {
	$count++;
    my $pline = <$params_in>;
    my $sline = <$sstats_in>;
	if (scalar(split(" ", $sline)) == $numberSummStats*2)  {
	    #print only values to Parameters file
	    my @plist = split(" ", $pline);
		@plist = @plist[@paramsAux];
		$plist[0] = substr($plist[0],2);
		$plist[1] = substr($plist[1],2);
	    print $params_out join(",",@plist),"\n";
	    #print only values to Summary Statistics file
	    my @slist = split(" ", $sline);
	    @slist = @slist[@sstatsAux];
	    print $sstats_out join(",",@slist),"\n";
	}
	else  {
		if ($flag == 1)  {
		    print "\nskipped line(s): ";
		    $flag = 0;
		}
		print "$count ";
	}
}
print $sstats_out "\n";

close ($params_in);
close ($sstats_in);
close ($params_out);
close ($sstats_out);
#### END ####

print "\nDone!\n\n";
exit;

