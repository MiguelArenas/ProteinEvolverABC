### SS_from_AA_Data.2.pl
### SS From MSA of AA sequences
### Programmer: Miguel Arenas Busto. 2020.
### Version 1.0 
###



use strict;
#use warnings;
use File::Copy;
use File::Basename;
use Cwd;



#print "*** SS_from_AA_Data.2.pl ***";
#print "\n Miguel Arenas Busto. 2020. \n";
#print "Input file in sequential phylip format: HEAD = 10 characters, then the sequence \n";




##############################################
##############################################
# Loading input file, "Settings"
my $file = $ARGV[0];
unless (open(FROM,$file))
	{
	print STDERR "Cannot open file \"$file\"\n\n";
	exit;
	}

my $show = $ARGV[1];
#print "\n\n\n\n>>>>>> YES in SS /// show: $show  ..\n\n\n\n";
if ($show == 1)
{
print "> Input file uploaded: $file \n\n";
}


##############################################
##############################################
# Directories
my $scrtdir = dirname($0);
#print ">  CodABC directory detected: $scrtdir \n\n";
my $setsdir = dirname($file);
#print "> Directory of Settings.txt detected: $setsdir \n\n";

my $setsbase = basename($file);



######################################################################################################
######################################### LOADING FILE #########################################
######################################################################################################

#my $outFileName = sprintf ("XX%s", $file);	
#open FILE_OUTPUT, sprintf ('>%s', $outFileName);

# opening lines of input files
unless (open(FROM,$file))
	{
	print STDERR "Cannot open file \"$file\"\n\n";
	}
if ($show == 1)
    {
    print "Working on \"$file\"\n";
    print "Gaps (-) or unknown states (?) will be considered as a new state\n";
    }
# EVOLVING
print "\rSSs dataset: \"$file\"";
#

my $aa_number = -1;
my $TotalNumberSequences = -1;
my $line_position = 0;

my @Long_SequenceHere;
my @Long_HeadSeqHere;
my $CurrentNumberSeq = 0;

while (<FROM>) # FOR EACH LINE OF THE FILE
	{
	$line_position++;
	
	if ($line_position == 1)
		{
		my $num3 = -1;
		my $num4 = -1;
		($num3, $num4) = split(/ /, $_);	
		$TotalNumberSequences = $num3;
		$aa_number = $num4;
		$aa_number =~ s/\n//g; # remove the end of line
        if ($show == 1)
            {
            print "$TotalNumberSequences sequences, $aa_number amino acids";
            }
		}
		
	if ($line_position > 1 && $line_position <= $TotalNumberSequences+1)
		{
		#print "\n $line_position \n";
		$CurrentNumberSeq++;
		
		my $line_sequence = '';
		my @SequenceHereX;
		my @HeadSeqHereX;
		my $SequenceHere = '';
		my $HeadSeqHere = '';
		my $LenghtChLine = -1;
		$LenghtChLine = $aa_number + 10;
		
		$line_sequence = $_;
		$line_sequence =~ s/\n//g; # remove the end of line
			
		#my @EndSequence = ();   # No with this
		my @InitialSequence;
		@InitialSequence = split ('', $line_sequence );
		
		my $cccc = 0;
		for ($cccc = 0; $cccc <= $LenghtChLine-1; $cccc++)
			{
			if ($cccc <= 9) # Head is the 10 first positions
				{
				$HeadSeqHereX[$cccc]	= $InitialSequence[$cccc]; # HEAD in array
				}
			else	# The rest is the sequence
				{
				$SequenceHereX[$cccc]	= $InitialSequence[$cccc];	# SEQUENCES in array
				
				if ($cccc >= 10 && $cccc <= $LenghtChLine-1) # checking
					{
					if ($InitialSequence[$cccc] eq 'A' || $InitialSequence[$cccc] eq 'R' || $InitialSequence[$cccc] eq 'N' || $InitialSequence[$cccc] eq 'D' || $InitialSequence[$cccc] eq 'C' || $InitialSequence[$cccc] eq 'Q' || $InitialSequence[$cccc] eq 'E' || $InitialSequence[$cccc] eq 'G' || $InitialSequence[$cccc] eq 'H' || $InitialSequence[$cccc] eq 'I' || $InitialSequence[$cccc] eq 'L' || $InitialSequence[$cccc] eq 'K' || $InitialSequence[$cccc] eq 'M' || $InitialSequence[$cccc] eq 'F' || $InitialSequence[$cccc] eq 'P' || $InitialSequence[$cccc] eq 'S' || $InitialSequence[$cccc] eq 'T' || $InitialSequence[$cccc] eq 'W' || $InitialSequence[$cccc] eq 'Y' || $InitialSequence[$cccc] eq 'V' || $InitialSequence[$cccc] eq '-' || $InitialSequence[$cccc] eq '?')
						{
						; # Ok	 A R N D C Q E G H I L K M F P S T W Y V
						}	
					else
						{
						print "\nThe position $cccc ($InitialSequence[$cccc]) is not a A R N D C Q E G H I L K M F P S T W Y V or -.\n Check the name of the sequences, they must have 10 characteres (NO MORE NEITHER LESS) in phylip sequencial format, so the sequence starts in character #11 \n";
						exit;		
						}
					} # end of checking
					
				}
			}
		
		
		
	#	print "\nLoading sequence number: $CurrentNumberSeq ..";
		for ($cccc = 0; $cccc <= $LenghtChLine; $cccc++)
			{
				
			if ($cccc <= 9) # Head is the 10 first positions
				{
				#if ($cccc == 0)
				#	{
				#	print "\n HEAD: "
				#	}
					
				my $letter = '';
				$letter = $HeadSeqHereX[$cccc]; # HEAD in vector
				$HeadSeqHere .= $letter;
				
				my $Pos = 0;
				$Pos = ($CurrentNumberSeq * 10) + $cccc;
				$Long_HeadSeqHere[$Pos] = $HeadSeqHereX[$cccc];
				#print "$Long_HeadSeqHere[$Pos]";
				}
			else	# The rest is the sequence
				{
				#if ($cccc == 10)
				#	{
				#	print "\n SEQUENCE: "
				#	}
				my $letter = '';
				$letter = $SequenceHereX[$cccc]; # SEQUENCES in vector
				$SequenceHere .= $letter;
				
				my $Pos = 0;
				$Pos = ($CurrentNumberSeq * $aa_number) + $cccc - 10;
				$Long_SequenceHere[$Pos] = $SequenceHereX[$cccc];	
				#print "$Long_SequenceHere[$Pos]";
				}
				
			}		
		
		
					
				
		} # end of lines > 1
				
	} # end of lines

if ($show == 1)
    {
    print "\n.. File \"$file\" loaded in memory ..\n";
    }




# ************** AA ALIGNMENT: FOR VIEWING AND WORKING *********************
my $n = 0;
my $Pos = 0;

# HEADS
my $HeadNumber = 0;
$n = 0;

for ($HeadNumber = 1; $HeadNumber <= $TotalNumberSequences; $HeadNumber++) # each sequence
	{
    #print "\nHEAD $HeadNumber: ";
	for ($n = 0; $n <= 9; $n++) # each character
		{
		$Pos = 0;
		$Pos = ($HeadNumber * 10) + $n;
		
        #print "$Long_HeadSeqHere[$Pos]";
		}
	}

# SEQUENCES
my $SeqNumber = 0;
$n = 0;
for ($SeqNumber = 1; $SeqNumber <= $TotalNumberSequences; $SeqNumber++) # each sequence
	{
    #print "\nSEQUENCE $SeqNumber: ";
		
	for ($n = 0; $n <= $aa_number-1; $n++) # each character
		{
		$Pos = 0;
		$Pos = ($SeqNumber * $aa_number) + $n;
		
        #print "$Long_SequenceHere[$Pos]";
		}
		
	}
		
	
#print "\nOutput file: $outFileName\n";	
close FROM;
#close FILE_OUTPUT;	



#######################################################################################################
########################################### LOADING AA FILE ###########################################
#######################################################################################################

my @AALong_SequenceHere = @Long_SequenceHere;
my @AALong_HeadSeqHere = @Long_HeadSeqHere;
my $AATotalNumberSequences = $TotalNumberSequences;

# ************** AA ALIGNMENT: FOR VIEWING AND WORKING *********************
my $AAn = 0;
my $AAPos = 0;

# HEADS
my $AAHeadNumber = 0;
$AAn = 0;

for ($AAHeadNumber = 1; $AAHeadNumber <= $AATotalNumberSequences; $AAHeadNumber++) # each sequence
	{
    #print "\nHEAD $AAHeadNumber: ";
	for ($AAn = 0; $AAn <= 9; $AAn++) # each character
		{
		$AAPos = 0;
		$AAPos = ($AAHeadNumber * 10) + $AAn;
		
        #print "$AALong_HeadSeqHere[$AAPos]";
		}
	}

# SEQUENCES
my $AASeqNumber = 0;
$AAn = 0;
for ($AASeqNumber = 1; $AASeqNumber <= $AATotalNumberSequences; $AASeqNumber++) # each sequence
	{
    #print "\nSEQUENCE $AASeqNumber: ";
		
	for ($AAn = 0; $AAn <= $aa_number-1; $AAn++) # each character
		{
		$AAPos = 0;
		$AAPos = ($AASeqNumber * $aa_number) + $AAn;
		
        #print "$AALong_SequenceHere[$AAPos]";
		}
	}
		
	



#####################################################################################################################
########################################### extracting Summary Statistics ###########################################
#####################################################################################################################

my $newfileX = $file;
$newfileX = sprintf ("%s/SS%s.ss", $setsdir, $setsbase);
open (FILE_OUTPUT, ">$newfileX");




#####################################################################
# Different number of AA AND Heterozygosity per site

# set to 0
my $SSnumber = 0;
my @NumbDiffAA = ();
my @HeterozygosityAA = ();
my $h = 0;
for ($h = 0; $h <= $aa_number-1; $h++) # each character
{
    $NumbDiffAA[$h] = 0;
    $HeterozygosityAA[$h] = 0;
}
my $NumberOfAA_SS = 0;


$AAPos = 0;
$AASeqNumber = 0;
$AAn = 0;
# AA: A C D E F G H I K L M N P Q R S T V W Y
my $IsA = 0;
my $IsC = 0;
my $IsD = 0;
my $IsE = 0;
my $IsF = 0;
my $IsG = 0;
my $IsH = 0;
my $IsI = 0;
my $IsK = 0;
my $IsL = 0;
my $IsM = 0;
my $IsN = 0;
my $IsP = 0;
my $IsQ = 0;
my $IsR = 0;
my $IsS = 0;
my $IsT = 0;
my $IsV = 0;
my $IsW = 0;
my $IsY = 0;
my $Isgap = 0;
my $Isint = 0;
my $FreqA = 0;
my $FreqC = 0;
my $FreqD = 0;
my $FreqE = 0;
my $FreqF = 0;
my $FreqG = 0;
my $FreqH = 0;
my $FreqI = 0;
my $FreqK = 0;
my $FreqL = 0;
my $FreqM = 0;
my $FreqN = 0;
my $FreqP = 0;
my $FreqQ = 0;
my $FreqR = 0;
my $FreqS = 0;
my $FreqT = 0;
my $FreqV = 0;
my $FreqW = 0;
my $FreqY = 0;
my $Freqgap = 0;
my $Freqint = 0;
my $Freq2A = 0;
my $Freq2C = 0;
my $Freq2D = 0;
my $Freq2E = 0;
my $Freq2F = 0;
my $Freq2G = 0;
my $Freq2H = 0;
my $Freq2I = 0;
my $Freq2K = 0;
my $Freq2L = 0;
my $Freq2M = 0;
my $Freq2N = 0;
my $Freq2P = 0;
my $Freq2Q = 0;
my $Freq2R = 0;
my $Freq2S = 0;
my $Freq2T = 0;
my $Freq2V = 0;
my $Freq2W = 0;
my $Freq2Y = 0;
my $Freq2gap = 0;
my $Freq2int = 0;


# Reading in sequences
for ($AAn = 0; $AAn <= $aa_number-1; $AAn++) # each site
{
    
    # AA DIFF
    my $AANumDiffHere = 0;
    $IsA = $IsC = $IsD = $IsE = $IsF = $IsG = $IsH = $IsI = $IsK = $IsL = $IsM = $IsN = $IsP = $IsQ = $IsR = $IsS = $IsT = $IsV = $IsW = $IsY = $Isgap = $Isint = 0;
    
    # Hetergzt
    my $AATotalDiffHere = 0;
    my $AAHeterozygosityAAThisSite = 0;
    
    
    for ($AASeqNumber = 1; $AASeqNumber <= $AATotalNumberSequences; $AASeqNumber++) # each sequence
    {
        $AAPos = 0;
        $AAPos = ($AASeqNumber * $aa_number) + $AAn;
        #print "\nCharacter: $AAn  // in Sequence: $AASeqNumber // Value: ";
        #print "$AALong_SequenceHere[$AAPos]";
        
        if ($AALong_SequenceHere[$AAPos] eq 'A')
        {
            $IsA++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'C')
        {
            $IsC++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'D')
        {
            $IsD++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'E')
        {
            $IsE++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'F')
        {
            $IsF++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'G')
        {
            $IsG++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'H')
        {
            $IsH++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'I')
        {
            $IsI++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'K')
        {
            $IsK++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'L')
        {
            $IsL++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'M')
        {
            $IsM++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'N')
        {
            $IsN++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'P')
        {
            $IsP++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'Q')
        {
            $IsQ++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'R')
        {
            $IsR++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'S')
        {
            $IsS++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'T')
        {
            $IsT++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'V')
        {
            $IsV++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'W')
        {
            $IsW++;
        }
        if ($AALong_SequenceHere[$AAPos] eq 'Y')
        {
            $IsY++;
        }
        if ($AALong_SequenceHere[$AAPos] eq '-')
        {
            $Isgap++;
        }
        if ($AALong_SequenceHere[$AAPos] eq '?')
        {
            $Isint++;
        }
        
    } # end of AA for this site
    
    
    # AA DIFF - counting the number of different AA
    if ($IsA > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsC > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsD > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsE > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsF > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsG > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsH > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsI > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsK > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsL > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsM > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsN > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsP > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsQ > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsR > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsS > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsT > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsV > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsW > 0)
    {
        $AANumDiffHere++;
    }
    if ($IsY > 0)
    {
        $AANumDiffHere++;
    }
    if ($Isgap > 0)
    {
        $AANumDiffHere++;
    }
    if ($Isint > 0)
    {
        $AANumDiffHere++;
    }
    $NumbDiffAA[$AAn] =	$AANumDiffHere;
    #print "Site $AAn  // NumAAdifHere: $AANumDiffHere\n";
    
    
    # Hetergzt
    $AATotalDiffHere = $IsA + $IsC + $IsD + $IsE + $IsF + $IsG + $IsH + $IsI + $IsK + $IsL + $IsM + $IsN + $IsP + $IsQ + $IsR + $IsS + $IsT + $IsV + $IsW + $IsY + $Isgap + $Isint;
    if ($AATotalDiffHere != $AATotalNumberSequences) # checking
    {
        print "\n\n An error have been found in SS_from_AA_and_Codon_Data.1.pl, part of AA: number of total diff = $AATotalDiffHere should not be different than the number of sequences: $AATotalNumberSequences \n\n";
        exit;
    }
    $FreqA = $IsA / $AATotalDiffHere;
    $FreqC = $IsC / $AATotalDiffHere;
    $FreqD = $IsD / $AATotalDiffHere;
    $FreqE = $IsE / $AATotalDiffHere;
    $FreqF = $IsF / $AATotalDiffHere;
    $FreqG = $IsG / $AATotalDiffHere;
    $FreqH = $IsH / $AATotalDiffHere;
    $FreqI = $IsI / $AATotalDiffHere;
    $FreqK = $IsK / $AATotalDiffHere;
    $FreqL = $IsL / $AATotalDiffHere;
    $FreqM = $IsM / $AATotalDiffHere;
    $FreqN = $IsN / $AATotalDiffHere;
    $FreqP = $IsP / $AATotalDiffHere;
    $FreqQ = $IsQ / $AATotalDiffHere;
    $FreqR = $IsR / $AATotalDiffHere;
    $FreqS = $IsS / $AATotalDiffHere;
    $FreqT = $IsT / $AATotalDiffHere;
    $FreqV = $IsV / $AATotalDiffHere;
    $FreqW = $IsW / $AATotalDiffHere;
    $FreqY = $IsY / $AATotalDiffHere;
    $Freqgap = $Isgap / $AATotalDiffHere;
    $Freqint = $Isint / $AATotalDiffHere;
    $Freq2A = $FreqA*$FreqA;
    $Freq2C = $FreqC*$FreqC;
    $Freq2D = $FreqD*$FreqD;
    $Freq2E = $FreqE*$FreqE;
    $Freq2F = $FreqF*$FreqF;
    $Freq2G = $FreqG*$FreqG;
    $Freq2H = $FreqH*$FreqH;
    $Freq2I = $FreqI*$FreqI;
    $Freq2K = $FreqK*$FreqK;
    $Freq2L = $FreqL*$FreqL;
    $Freq2M = $FreqM*$FreqM;
    $Freq2N = $FreqN*$FreqN;
    $Freq2P = $FreqP*$FreqP;
    $Freq2Q = $FreqQ*$FreqQ;
    $Freq2R = $FreqR*$FreqR;
    $Freq2S = $FreqS*$FreqS;
    $Freq2T = $FreqT*$FreqT;
    $Freq2V = $FreqV*$FreqV;
    $Freq2W = $FreqW*$FreqW;
    $Freq2Y = $FreqY*$FreqY;
    $Freq2gap = $Freqgap*$Freqgap;
    $Freq2int = $Freqgap*$Freqint;
    
    my $sumXi2 = 0;
    $sumXi2 = $Freq2A + $Freq2C + $Freq2D + $Freq2E + $Freq2F + $Freq2G + $Freq2H + $Freq2I + $Freq2K + $Freq2L + $Freq2M + $Freq2N + $Freq2P + $Freq2Q + $Freq2R + $Freq2S + $Freq2T + $Freq2V + $Freq2W + $Freq2Y + $Freq2gap + $Freq2int;
    my $restH = 0;
    $restH = 1 - $sumXi2;
    my $prefixCocientH = 0;
    $prefixCocientH = $AATotalDiffHere - 1;
    my $prefixH = 0;
    if ($prefixCocientH != 0)
    {
        $prefixH = $AATotalDiffHere / $prefixCocientH;
    }
    else
    {
        $prefixH = "inf";	
    }
    $AAHeterozygosityAAThisSite = $prefixH * $restH;
    $HeterozygosityAA[$AAn] = $AAHeterozygosityAAThisSite;
    
    
    ### number of AA segregating sites
    if ($AANumDiffHere > 1)
    {
        $NumberOfAA_SS++; 
    }
    #print " NumberdiffHere: $AANumDiffHere SScum: $NumberOfAA_SS\n";
    
} # end of sites



#####################################################################
### Phi
if ($show == 1)
{
print "\n.. Calculating Phi, NSS and MaxChiSquare SS ..";
}
    
# Name of the file
my $newfileM = $file;
my $noThis = "";
my $numberThis = "";
my $windowPhi1 = 50;
my $windowPhi2 = 50;
($noThis, $numberThis) = split(/sequences/, $file);
$numberThis =~ s/\n//g; # remove the end of line	
my $newfileM = sprintf ("sequences%s", $numberThis);

#print "\n\n$newfileM\n\n";

if (-e "TempPhiResults$newfileM.txt")
	{
	system("rm TempPhiResults$newfileM.txt");
	}

if ($aa_number > 300 && $NumberOfAA_SS > 10)
    {
    system("\"$scrtdir/../bin/Phi\" -s \"$file\" -t A -p -o -v > TempPhiResults$newfileM.txt");
    if ($show == 1)
        {
        print "\n\"$scrtdir/../bin/Phi\" -s \"$file\" -t A -p -o -v > TempPhiResults$newfileM.txt\n";
        }
    }
else
    {
    $windowPhi1 = $aa_number / 3; # too short sequence for the window size specified by default in Phi, need for reducing the window size
    if ($NumberOfAA_SS < 10)
        {
        $windowPhi1 = $windowPhi1 / 2;
        }
    my $windowPhi2 = sprintf "%.0f", $windowPhi1;
    system("\"$scrtdir/../bin/Phi\" -s \"$file\" -t A -p -w \"$windowPhi2\" -o -v > TempPhiResults$newfileM.txt");
    if ($show == 1)
        {
        print "\n\"$scrtdir/../bin/Phi\" -s \"$file\" -t A -p -w \"$windowPhi2\" -o -v > TempPhiResults$newfileM.txt\n";
        }
    }
open (FILE_PHI,"TempPhiResults$newfileM.txt");
#print FILE_OUTPUT "Phi:	$PhiActual	NSS:	$NSSActual	MaxChiSquare:	$MaxChiSqActual	";
my $PhiLine = 0;
my $PhiActual = "";
my $NSSActual = "";
my $MaxChiSqActual = "";

while (<FILE_PHI>)
{
    $PhiLine++;
    
    if ($_ =~ /Mean:        /)
    {
        #print $_;
        my @partsX = split(/\s+/, $_);
        $PhiActual = $partsX[1];
    }
    if ($_ =~ /The Neighbour Similarity score is/)
    {
        #print $_;
        my @partsNSS = split(/\s+/, $_);
        $NSSActual = $partsNSS[5];
    }
    if ($_ =~ /Value of maximum breakpoint is:/)
    {
        #print $_;
        my @partsMaxChi = split(/\s+/, $_);
        $MaxChiSqActual = $partsMaxChi[5];
    }
				
}
close (FILE_PHI);
# Activar lo de abajo
system("rm Phi.inf.list Phi.inf.sites Phi.log TempPhiResults$newfileM.txt Phi.poly.unambig.sites");

print FILE_OUTPUT "Phi:	$PhiActual	NSS:	$NSSActual	MaxChiSquare:	$MaxChiSqActual	";
if ($show == 1)
{
print ".. Phi:	$PhiActual	NSS:	$NSSActual	MaxChiSquare:	$MaxChiSqActual	";
}



#####################################################################
## AA DIFF - Obtaining results for Different AA per site.. :
my $CumNumber_of_Diff_AA = 0;
my $Average_of_Diff_AA = 0;
my $Sd_of_Diff_AA = 0;
my $Skewness_of_Diff_AA = 0;
my $Kurtosis_of_Diff_AA = 0;

# Average
$SSnumber = $SSnumber + 1;
for ($h = 0; $h <= $aa_number-1; $h++) # each site
	{
	$CumNumber_of_Diff_AA = $CumNumber_of_Diff_AA + $NumbDiffAA[$h];	# $NumbDiffAA[$h] = this value is the number of AA differences for each site
	}
$Average_of_Diff_AA = $CumNumber_of_Diff_AA / $aa_number;
print FILE_OUTPUT "AADiff_Avrge:	$Average_of_Diff_AA	";
if ($show == 1)
{
print "\n.. Calculating Average of different amino acids SS ($SSnumber) ..";
print "\t$Average_of_Diff_AA	";
}

# SD
$SSnumber = $SSnumber + 1;
my $value_for_SD = 0;
my $cum_value_for_SD = 0;
my $variance = 0;
for ($h = 0; $h <= $aa_number-1; $h++) # each site
	{
	$value_for_SD = abs($NumbDiffAA[$h] - $Average_of_Diff_AA);
	$value_for_SD = ($value_for_SD * $value_for_SD);
	$cum_value_for_SD = $cum_value_for_SD + $value_for_SD;
	}
$variance = ($cum_value_for_SD / $aa_number);
$Sd_of_Diff_AA = sqrt($variance);

$h = 0;
$value_for_SD = 0;
$cum_value_for_SD = 0;
$variance = 0;
print FILE_OUTPUT "AADiff_Sd:	$Sd_of_Diff_AA	";
if ($show == 1)
{
print "\n.. Calculating SD of different amino acids SS ($SSnumber) ..";
print "\t\t$Sd_of_Diff_AA	";
}

# Skewness
$SSnumber = $SSnumber + 1;
my $value_for_SK = 0;
my $cum_value_for_SK = 0;
my $cocient_SK = 0;
for ($h = 0; $h <= $aa_number-1; $h++) # each site
	{
	$value_for_SK = abs($NumbDiffAA[$h] - $Average_of_Diff_AA);
	$value_for_SK = ($value_for_SK * $value_for_SK * $value_for_SK);
	$cum_value_for_SK = $cum_value_for_SK + $value_for_SK;
	}
$cocient_SK = $Sd_of_Diff_AA * $Sd_of_Diff_AA * $Sd_of_Diff_AA;
$cocient_SK = ($aa_number-1) * $cocient_SK;
if ($cocient_SK != 0)
	{
	$Skewness_of_Diff_AA = $cum_value_for_SK / $cocient_SK;
	}
else
	{
	$Skewness_of_Diff_AA = "inf";	
	}

$h = 0;
$value_for_SK = 0;
$cum_value_for_SK = 0;
$cocient_SK = 0;
print FILE_OUTPUT "AADiff_Skewnss:	$Skewness_of_Diff_AA	";
if ($show == 1)
{
print "\n.. Calculating Skewness of different amino acids SS ($SSnumber) ..";
print "	$Skewness_of_Diff_AA	";
}

# Kurtosis
$SSnumber = $SSnumber + 1;
my $value_for_KU = 0;
my $cum_value_for_KU = 0;
my $cocient_KU = 0;
for ($h = 0; $h <= $aa_number-1; $h++) # each site
	{
	$value_for_KU = abs($NumbDiffAA[$h] - $Average_of_Diff_AA);
	$value_for_KU = ($value_for_KU * $value_for_KU * $value_for_KU * $value_for_KU);
	$cum_value_for_KU = $cum_value_for_KU + $value_for_KU;
	}
$cocient_KU = $Sd_of_Diff_AA * $Sd_of_Diff_AA * $Sd_of_Diff_AA * $Sd_of_Diff_AA;
$cocient_KU = ($aa_number-1) * $cocient_KU;
if ($cocient_KU != 0)
	{
	$Kurtosis_of_Diff_AA = $cum_value_for_KU / $cocient_KU;
	}
else
	{
	$Kurtosis_of_Diff_AA = "inf";	
	}

$h = 0;
$value_for_KU = 0;
$cum_value_for_KU = 0;
$cocient_KU = 0;
print FILE_OUTPUT "AADiff_Kurtosis:	$Kurtosis_of_Diff_AA	";
if ($show == 1)
{
print "\n.. Calculating Kurtosis of different amino acids SS ($SSnumber) ..";
print "	$Kurtosis_of_Diff_AA	";
}


## AA Heterozygosity - Obtaining results for Different AA per site.. :
my $CumNumber_of_Het_AA = 0;
my $Average_of_Het_AA = 0;
my $Sd_of_Het_AA = 0;
my $Skewness_of_Het_AA = 0;
my $Kurtosis_of_Het_AA = 0;

# Average
$SSnumber = $SSnumber + 1;
for ($h = 0; $h <= $aa_number-1; $h++) # each site
	{
	$CumNumber_of_Het_AA = $CumNumber_of_Het_AA + $HeterozygosityAA[$h];	# $HeterozygosityAA[$h] = this value is the number of AA Heterozygosity for each site
	}
$Average_of_Het_AA = $CumNumber_of_Het_AA / $aa_number;
print FILE_OUTPUT "AAHetrzgs_Avrge:	$Average_of_Het_AA	";
if ($show == 1)
{
print "\n.. Calculating Average of Heterozygosity amino acids SS ($SSnumber) ..";
print "	$Average_of_Het_AA	";
}

# SD
$SSnumber = $SSnumber + 1;
$value_for_SD = 0;
$cum_value_for_SD = 0;
$variance = 0;
for ($h = 0; $h <= $aa_number-1; $h++) # each site
	{
	$value_for_SD = abs($HeterozygosityAA[$h] - $Average_of_Het_AA);
	$value_for_SD = ($value_for_SD * $value_for_SD);
	$cum_value_for_SD = $cum_value_for_SD + $value_for_SD;
	}
$variance = ($cum_value_for_SD / $aa_number);
$Sd_of_Het_AA = sqrt($variance);

$h = 0;
$value_for_SD = 0;
$cum_value_for_SD = 0;
$variance = 0;
print FILE_OUTPUT "AAHetrzgs_Sd:	$Sd_of_Het_AA	";
if ($show == 1)
{
print "\n.. Calculating SD of Heterozygosity amino acids SS ($SSnumber) ..";
print "	$Sd_of_Het_AA	";
}

# Skewness
$SSnumber = $SSnumber + 1;
$value_for_SK = 0;
$cum_value_for_SK = 0;
$cocient_SK = 0;
for ($h = 0; $h <= $aa_number-1; $h++) # each site
	{
	$value_for_SK = abs($HeterozygosityAA[$h] - $Average_of_Het_AA);
	$value_for_SK = ($value_for_SK * $value_for_SK * $value_for_SK);
	$cum_value_for_SK = $cum_value_for_SK + $value_for_SK;
	}
$cocient_SK = $Sd_of_Het_AA * $Sd_of_Het_AA * $Sd_of_Het_AA; # eye ..
$cocient_SK = ($aa_number-1) * $cocient_SK;
if ($cocient_SK != 0)
	{
	$Skewness_of_Het_AA = $cum_value_for_SK / $cocient_SK;
	}
else
	{
	$Skewness_of_Het_AA = "inf";	
	}

$h = 0;
$value_for_SK = 0;
$cum_value_for_SK = 0;
$cocient_SK = 0;
print FILE_OUTPUT "AAHetrzgs_Skewnss:	$Skewness_of_Het_AA	";
if ($show == 1)
{
print "\n.. Calculating Skewness of Heterozygosity amino acids SS ($SSnumber) ..";
print "	$Skewness_of_Het_AA	";
}

# Kurtosis
$SSnumber = $SSnumber + 1;
$value_for_KU = 0;
$cum_value_for_KU = 0;
$cocient_KU = 0;
for ($h = 0; $h <= $aa_number-1; $h++) # each site
	{
	$value_for_KU = abs($HeterozygosityAA[$h] - $Average_of_Het_AA);
	$value_for_KU = ($value_for_KU * $value_for_KU * $value_for_KU * $value_for_KU);
	$cum_value_for_KU = $cum_value_for_KU + $value_for_KU;
	}
$cocient_KU = $Sd_of_Het_AA * $Sd_of_Het_AA * $Sd_of_Het_AA * $Sd_of_Het_AA; # eye ..
$cocient_KU = ($aa_number-1) * $cocient_KU;
if ($cocient_KU != 0)
	{
	$Kurtosis_of_Het_AA = $cum_value_for_KU / $cocient_KU;
	}
else
	{
	$Kurtosis_of_Het_AA = "inf";	
	}

$h = 0;
$value_for_KU = 0;
$cum_value_for_KU = 0;
$cocient_KU = 0;
print FILE_OUTPUT "AAHetrzgs_Kurtosis:	$Kurtosis_of_Het_AA	";
if ($show == 1)
{
print "\n.. Calculating Kurtosis of Heterozygosity amino acids SS ($SSnumber) ..";
print "	$Kurtosis_of_Het_AA	";
}


# Number of AA segregating sites
$SSnumber = $SSnumber + 1;
if ($show == 1)
{
print "\n.. Calculating Amino Acid segregating sites ($SSnumber) ..";
print "	$NumberOfAA_SS	";
}
print FILE_OUTPUT "AA_SegSites:	$NumberOfAA_SS	";

# Round number to 3 digits after decimal point
#$rounded = sprintf("%.3f", $number);




#####################################################################
################### NEW SS (added in 2020) ##########################
#SS based on sequence identity between pairs

#print "\n\n--- NEW Summary Statistics (2020) ---\n\n";


# ************** AA ALIGNMENT *********************
$AAn = 0;
$AAPos = 0;
my @SeqIDbetweenSequences;
my $countSeqID = 0;

# SEQUENCES
my $AASeqNumber_a = 0;
my $AASeqNumber_b = 0;
my $AAPos_a = 0;
my $AAPos_b = 0;

$AAn = 0;
my $HereDist = 0;
my $nextSequence = 0;

for ($AASeqNumber_a = 1; $AASeqNumber_a <= $AATotalNumberSequences; $AASeqNumber_a++) # each sequence
    {
    #print "\n>SEQUENCE $AASeqNumber_a \n";
    $nextSequence = $AASeqNumber_a + 1;
    
    for ($AASeqNumber_b = $nextSequence; $AASeqNumber_b <= $AATotalNumberSequences; $AASeqNumber_b++) # each sequence
        {
        #print " compared with SEQUENCE $AASeqNumber_b ";
        $HereDist = 0;
        
        
        for ($AAn = 0; $AAn <= $aa_number-1; $AAn++) # each character
            {
            $AAPos_a = 0;
            $AAPos_a = ($AASeqNumber_a * $aa_number) + $AAn;
            
            $AAPos_b = 0;
            $AAPos_b = ($AASeqNumber_b * $aa_number) + $AAn;
            
            if ($AALong_SequenceHere[$AAPos_a] eq $AALong_SequenceHere[$AAPos_b])
                {
                #print "\nEqual states: $AALong_SequenceHere[$AAPos_a] $AALong_SequenceHere[$AAPos_b]\n";
                }
            else
                {
                $HereDist = $HereDist + 1;
                #print "\nDifferent states: $AALong_SequenceHere[$AAPos_a] $AALong_SequenceHere[$AAPos_b] (HereDist = $HereDist)\n";
                }
            
            #print "$AALong_SequenceHere[$AAPos]";
            }
        
        #print "\nHereDist = $HereDist\n";
        $countSeqID = $countSeqID + 1;
        $SeqIDbetweenSequences[$countSeqID] = $HereDist;
            
        }
        
    }


## AA SeqID - Obtaining results for Different AA per site.. :
my $CumNumber_of_SeqID_AA = 0;
my $Average_of_SeqID_AA = 0;
my $Sd_of_SeqID_AA = 0;
my $Skewness_of_SeqID_AA = 0;
my $Kurtosis_of_SeqID_AA = 0;

# Average
$SSnumber = $SSnumber + 1;
my $ccc = 0;
my $cumm = 0;
my $meanCC = 0;
for ($ccc = 1; $ccc <= $countSeqID; $ccc++) # each sequence
    {
    $cumm = $cumm + $SeqIDbetweenSequences[$ccc];
    #print "$SeqIDbetweenSequences[$ccc] ";
    }
$meanCC = $cumm / $countSeqID;
#print "\nItems = $countSeqID /// cum = $cumm /// mean = $meanCC\n";
print FILE_OUTPUT "AA_SeqID_Avrge:	$meanCC	";
if ($show == 1)
{
print "\n.. Calculating Average of pairwise SeqID SS ($SSnumber) ..";
print "	$meanCC	";
}

# SD
$SSnumber = $SSnumber + 1;
my $value_for_SD = 0;
my $cum_value_for_SD = 0;
my $variance = 0;
for ($ccc = 1; $ccc <= $countSeqID; $ccc++) # each pair of sequences
    {
    $value_for_SD = abs($SeqIDbetweenSequences[$ccc] - $meanCC);
    $value_for_SD = ($value_for_SD * $value_for_SD);
    $cum_value_for_SD = $cum_value_for_SD + $value_for_SD;
    }
$variance = ($cum_value_for_SD / $countSeqID);
$Sd_of_SeqID_AA = sqrt($variance);

$h = 0;
$value_for_SD = 0;
$cum_value_for_SD = 0;
$variance = 0;
print FILE_OUTPUT "AA_SeqID_Sd:	$Sd_of_SeqID_AA	";
if ($show == 1)
{
print "\n.. Calculating SD of pairwise SeqID SS ($SSnumber) ..";
print "\t\t$Sd_of_SeqID_AA	";
}


# Skewness
$SSnumber = $SSnumber + 1;
my $value_for_SK = 0;
my $cum_value_for_SK = 0;
my $cocient_SK = 0;
for ($ccc = 1; $ccc <= $countSeqID; $ccc++) # each pair of sequences
    {
    $value_for_SK = abs($SeqIDbetweenSequences[$ccc] - $meanCC);
    $value_for_SK = ($value_for_SK * $value_for_SK * $value_for_SK);
    $cum_value_for_SK = $cum_value_for_SK + $value_for_SK;
    }
$cocient_SK = $Sd_of_SeqID_AA * $Sd_of_SeqID_AA * $Sd_of_SeqID_AA; # eye ..
$cocient_SK = $countSeqID * $cocient_SK;
if ($cocient_SK != 0)
    {
    $Skewness_of_SeqID_AA = $cum_value_for_SK / $cocient_SK;
    }
else
    {
    $Skewness_of_SeqID_AA = "inf";
    }

$h = 0;
$value_for_SK = 0;
$cum_value_for_SK = 0;
$cocient_SK = 0;
print FILE_OUTPUT "AA_SeqID_Skewnss:	$Skewness_of_SeqID_AA	";
if ($show == 1)
{
print "\n.. Calculating Skewness of pairwise SeqID SS ($SSnumber) ..";
print "	$Skewness_of_SeqID_AA	";
}
    

# Kurtosis
$SSnumber = $SSnumber + 1;
my $value_for_KU = 0;
my $cum_value_for_KU = 0;
my $cocient_KU = 0;
for ($ccc = 1; $ccc <= $countSeqID; $ccc++) # each pair of sequences
    {
    $value_for_KU = abs($SeqIDbetweenSequences[$ccc] - $meanCC);
    $value_for_KU = ($value_for_KU * $value_for_KU * $value_for_KU * $value_for_KU);
    $cum_value_for_KU = $cum_value_for_KU + $value_for_KU;
    }
$cocient_KU = $Sd_of_SeqID_AA * $Sd_of_SeqID_AA * $Sd_of_SeqID_AA * $Sd_of_SeqID_AA; # eye ..
$cocient_KU = $countSeqID * $cocient_KU;
if ($cocient_KU != 0)
    {
    $Kurtosis_of_SeqID_AA = $cum_value_for_KU / $cocient_KU;
    }
else
    {
    $Kurtosis_of_SeqID_AA = "inf";
    }

$h = 0;
$value_for_KU = 0;
$cum_value_for_KU = 0;
$cocient_KU = 0;
print FILE_OUTPUT "AA_SeqID_Kurtosis:	$Kurtosis_of_SeqID_AA	";
if ($show == 1)
{
print "\n.. Calculating Kurtosis of pairwise SeqID SS ($SSnumber) ..";
print "	$Kurtosis_of_SeqID_AA	";
print "\nDone!\n";
}




###########################

close (FILE_OUTPUT);
#### END ####

exit;
	
	
	


