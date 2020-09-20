### ProteinEvolverABC_Phase4.pl
###


use strict;
use warnings;
use File::Basename;

##############################################
##############################################
# Credits
system("clear");
print "\n****************************************************************************************************************";
print "\nProteinEvolverABC_Phase4.pl";
print "\nProteinEvolverABC Phase 4 does:";
print "\n- Read Settings file";
print "\n- Perform ABC analysis using simulated data and target data (requires R with library: abc)";
print "\n****************************************************************************************************************\n\n";


##############################################
##############################################
# Loading input file, "Settings"
my $file = $ARGV[0];
unless (open(FROM,$file))  {
	print STDERR "Cannot open file \"$file\"\n\n";
	exit;
}
#print "> Input file uploaded: $file \n\n";


##############################################
##############################################
# Directories
#my $maindir = dirname($0);
#print "> ProteinEvolverABC directory detected: $maindir \n\n";
my $setsdir = dirname($file);
#print "> Settings directory detected: $setsdir \n\n";


##############################################
##############################################
# variables
my $nIter = 0;
my $abcNIter = 0;
my $abcNTol = 0;
my $abcMethod = "";
my $abcTransf = "";
my $abcHCorr = 0;
my $summaryStatisticsAux = "";
my @summaryStatistics = (0);
my $multiPage = 0;

my $dump = "";


##############################################
##############################################
# Reading Settings from input file
print "> Reading Settings from input file ... \n";
open(FROM,$file);
while (<FROM>)  {
	# Number of simulations
	if ($_ =~ /^\*NumberOfSimulations=/)  {
		($dump, $nIter) = split(/=/, $_);
		$nIter =~ s/\n//g; # remove the end of line
		if ($nIter < 1)  {
			print "\n\nERROR!!!. The number of replicates must be higher than 0 ($nIter)\n";
			print "Type CTRL+C to abort the execution \n";
			my $HereError = <STDIN>;
			chop($HereError);
			exit;
		}
		print "  Number of simulations: $nIter\n";
	}


	### ABC settings
	# Number of iterations
	if ($_ =~ /^\*ABCIterations=/)  {
		($dump, $abcNIter) = split(/=/, $_);
		$abcNIter =~ s/\n//g; # remove the end of line
		if ($abcNIter < 1)  {
			print "\n\nERROR!!!. The number of replicates to consider must be higher than 0 ($abcNIter)\n";
			print "Type CTRL+C to abort the execution \n";
			my $HereError = <STDIN>;
			chop($HereError);
			exit;
		}
		print "  Number of iterations (simulations considered for the estimation phase): $abcNIter\n";
	}

	# ABC tolerance
	if ($_ =~ /^\*ABCTolerance=/)  {
		($dump, $abcNTol) = split(/=/, $_);
		$abcNTol =~ s/\n//g; # remove the end of line
		if ($abcNTol < 1)  {
			print "\n\nERROR!!!. The tolerance number must be higher than 0 ($abcNTol)\n";
			print "Type CTRL+C to abort the execution \n";
			my $HereError = <STDIN>;
			chop($HereError);
			exit;
		}
		print "  Tolerance number: $abcNTol\n";
	}


	# ABC method
	if ($_ =~ /^\*ABCMethod=/)  {
		($dump, $abcMethod) = split(/=/, $_);
		$abcMethod =~ s/\n//g; # remove the end of line
		if ($abcMethod ne "rejection" && $abcMethod ne "loclinear")  {
			print "\n\nERROR!. ABCMethod must be \"rejection\" or \"loclinear\" ($abcMethod)\n\n";
			print "Type CTRL+C to abort the execution \n";
			my $HereError = <STDIN>;
			chop($HereError);
			exit;
		}
		print "  ABCMethod: $abcMethod\n";
	}

	# ABC transformation
	if ($_ =~ /^\*ABCTransf=/)  {
		($dump, $abcTransf) = split(/=/, $_);
		$abcTransf =~ s/\n//g; # remove the end of line
		if ($abcTransf ne "none" && $abcTransf ne "log" && $abcTransf ne "logit")  {
			print "\n\nERROR!. ABCTransf must be \"none\", \"log\" or \"logit\" ($abcTransf)\n\n";
			print "Type CTRL+C to abort the execution \n";
			my $HereError = <STDIN>;
			chop($HereError);
			exit;
		}
		print "  ABCTransf: $abcTransf\n";
	}

	# ABC correction
	if ($_ =~ /^\*ABCHCorr=/)  {
		($dump, $abcHCorr) = split(/=/, $_);
		$abcHCorr =~ s/\n//g; # remove the end of line
		if ($abcHCorr != 0 && $abcHCorr != 1)  {
			print "\n\nERROR!. ABCHCorr must be 0 or 1 ($abcHCorr)\n\n";
			print "Type CTRL+C to abort the execution \n";
			my $HereError = <STDIN>;
			chop($HereError);
			exit;
		}
		print "  ABCHCorr: $abcHCorr\n";
	}

	# Summary Statistics
	if ($_ =~ /^\*SummaryStatistics=/)  {
		($dump, $summaryStatisticsAux) = split(/=/, $_);
		$summaryStatisticsAux =~ s/\n//g; # remove the end of line
		@summaryStatistics = split(' ', $summaryStatisticsAux);
		if (scalar grep( $_ < 1 || $_ > 16, @summaryStatistics) > 0) {
			print "\n\nERROR!. SummaryStatistics must be integers between 1 and 16 ($abcHCorr)\n\n";
			print "Type CTRL+C to abort the execution \n";
			my $HereError = <STDIN>;
			chop($HereError);
			exit;
		}
		print "  SummaryStatistics: @summaryStatistics\n";
	}



	### Graphical settings
	# Multiple pages
	if ($_ =~ /^\*MultiPage=/)  {
		($dump, $multiPage) = split(/=/, $_);
		$multiPage =~ s/\n//g; # remove the end of line
		if ($multiPage != 0 && $multiPage != 1)  {
			print "\n\nERROR!. MultiPage must be 0 or 1 ($multiPage)\n\n";
			print "Type CTRL+C to abort the execution \n";
			my $HereError = <STDIN>;
			chop($HereError);
			exit;
		}
		print "  MultiPage: $multiPage\n\n";
	}

} # end of lines FROM

close (FROM);


##############################################
##############################################
# Declare fixed variables
my $numberSummStats = 16;
my $numberParams = 2;
my $fileRScript = "$setsdir/PerformABC.r";
my $fileOutSummStats = "$setsdir/SSsimulations.csv";
my $fileOutParams = "$setsdir/PSimulations.csv";
my $fileRealData = "$setsdir/SSRealData.ss";
my $fileABCSumm1 = "$setsdir/Summary_ABC_rejection.csv";
my $fileABCSumm2 = "$setsdir/Summary_ABC_loclinear.csv";
my $fileABCFig1 = "$setsdir/Posteriors_rejection";
my $fileABCFig2 = "$setsdir/Posteriors_loclinear";
my $fileABCFig3 = "$setsdir/Histograms_SStats";
my $fileABCFig4 = "$setsdir/PCA_SStats.pdf";
my $fileABCFig5 = "$setsdir/Scaterplots_SStatsVSParams";


##############################################
##############################################
# Remove previous main output
if (-e "$setsdir/MainOutputs")  {
	system ("rm -rf \"$setsdir/MainOutputs\"");
}


##############################################
##############################################
# Writing R script

# Checking values
if ($abcNIter < 1 || $abcNIter > $nIter)  {
	print "\n\nERROR!!!. The number of replicates to consider must be higher than 0 and lower or equal the total number of replicates($abcNIter)\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
}
if ($abcNTol < 1 || $abcNTol > $abcNIter)  {
	print "\n\nERROR!!!. The ABC tolerance number must be higher than 0 and lower than number of replicates ($abcNTol)\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
}
if ($abcNTol < scalar @summaryStatistics)  {
	print "\n\nERROR!!!. The ABC tolerance number must be higher than the number of used summary statistics ($abcNTol)\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
}
if ($abcMethod ne "rejection" && $abcMethod ne "loclinear")  {
	print "\n\nERROR!. ABCMethod must be \"rejection\" or \"loclinear\" ($abcMethod)\n\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
}
if ($abcTransf ne "none" && $abcTransf ne "log" && $abcTransf ne "logit")  {
	print "\n\nERROR!. ABCTransf must be \"none\", \"log\" or \"logit\" ($abcTransf)\n\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
}
if ($abcHCorr != 0 && $abcHCorr != 1)  {
	print "\n\nERROR!. ABCHCorr must be 0 or 1 ($abcHCorr)\n\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
}
if (scalar grep( $_ < 1 || $_ > 16, @summaryStatistics) > 0) {
	print "\n\nERROR!. SummaryStatistics must be integers between 1 and 16 ($abcHCorr)\n\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
} 
if ($multiPage != 0 && $multiPage != 1)  {
	print "\n\nERROR!. MultiPage must be 0 or 1 ($multiPage)\n\n";
	print "Type CTRL+C to abort the execution \n";
	my $HereError = <STDIN>;
	chop($HereError);
	exit;
} 

# Writing in the R script
open (FILE_1,">$fileRScript");
print "> Writing R output file \"$fileRScript\" ...\n\n";

print FILE_1 "# R script to perform the ABC analysis \n";
print FILE_1 "# \n\n";
print FILE_1 "# Load libraries\n";
print FILE_1 "suppressPackageStartupMessages(library(abc))\n\n";
print FILE_1 "##### GENERAL SETTINGS #####\n";
print FILE_1 "nsstats <- $numberSummStats\n";
print FILE_1 "nparams <- $numberParams\n";
print FILE_1 "sstats.set <- c(",join(",",@summaryStatistics),")\n";
print FILE_1 "sstats.fname <- \"$fileOutSummStats\"\n";
print FILE_1 "params.fname <- \"$fileOutParams\"\n";
print FILE_1 "target.fname <- \"$fileRealData\"\n";
print FILE_1 "abc.niter <- $abcNIter\n";
print FILE_1 "abc.ntol <- $abcNTol\n";
print FILE_1 "abc.method <- \"$abcMethod\"  #loclinear or rejection\n";
print FILE_1 "abc.transf <- \"$abcTransf\"  #log, logit or none\n";
print FILE_1 "abc.hcorr <- TRUE\n";
print FILE_1 "multipage <- TRUE\n";
print FILE_1 "summ1.fname <- \"$fileABCSumm1\"\n";
print FILE_1 "summ2.fname <- \"$fileABCSumm2\"\n";
print FILE_1 "fig1.fname <- \"$fileABCFig1\"\n";
print FILE_1 "fig2.fname <- \"$fileABCFig2\"\n";
print FILE_1 "fig3.fname <- \"$fileABCFig3\"\n";
print FILE_1 "fig4.fname <- \"$fileABCFig4\"\n";
print FILE_1 "fig5.fname <- \"$fileABCFig5\"\n";
print FILE_1 "############################\n\n";
print FILE_1 "#read sstats data\n";
print FILE_1 "sstats <- matrix(scan(sstats.fname,skip=1,nlines=abc.niter,sep=\",\",na.strings=\"--\",what=double(),quiet=T),ncol=nsstats,byrow=T)\n";
print FILE_1 "colnames(sstats) <- scan(sstats.fname,skip=0,nlines=1,sep=\",\",what=\"character\",quiet=T)\n";
print FILE_1 "cat(\" Summary Statistics file read\\n\")\n\n";
print FILE_1 "#read params data\n";
print FILE_1 "params <- matrix(scan(params.fname,skip=1,nlines=abc.niter,sep=\",\",what=double(),quiet=T),ncol=nparams,byrow=T)\n";
print FILE_1 "colnames(params) <- scan(params.fname,skip=0,nlines=1,sep=\",\",what=\"character\",quiet=T)\n";
print FILE_1 "cat(\"\ Parameters file read\\n\")\n\n";
print FILE_1 "#read target data\n";
print FILE_1 "target <-as.numeric(scan(target.fname,what=\"raw\",quiet=T)[-c(seq(1,2*nsstats-1,2))])\n\n";
print FILE_1 "cat(\"\ Sequences data file read\\n\")\n\n";
print FILE_1 "#select sstats set\n";
print FILE_1 "sstats<-sstats[,sstats.set]\n";
print FILE_1 "target<-target[sstats.set]\n\n";
print FILE_1 "#get rid of NA's and Inf's\n";
print FILE_1 "sstats.na <- apply(sstats,1,function(x) all(!is.na(x)) && all(is.finite(x)))\n";
print FILE_1 "if(abc.niter - sum(sstats.na) > 0){\n";
print FILE_1 "    params<-params[sstats.na,]\n";
print FILE_1 "    sstats<-sstats[sstats.na,]\n";
print FILE_1 "    cat(paste(\"\\n   Simulations whose summary statistics have NA or Inf values discarded: \",abc.niter - sum(sstats.na),\"\\n\",sep=\"\"))\n";
print FILE_1 "}\n";
print FILE_1 "if(dim(params)[1] < abc.ntol){\n";
print FILE_1 "    abc.ntol <- dim(params)[1]\n";
print FILE_1 "    cat(\"\\n   Tolerance number is larger than number of simulations!\\n\")\n";
print FILE_1 "    cat(\"\\n   Tolerance set to number of simulations\\n\")\n";
print FILE_1 "}\n\n";
print FILE_1 "#get range of params\n";
print FILE_1 "params.min <- apply(params,2,min)\n";
print FILE_1 "params.max <- apply(params,2,max)\n";
print FILE_1 "sstats.min <- apply(sstats,2,function(x) quantile(x,names=F,probs=0.001))\n";
print FILE_1 "sstats.max <- apply(sstats,2,function(x) quantile(x,names=F,probs=0.999))\n\n";
print FILE_1 "####### ABC ANALYSIS #######\n";
print FILE_1 "abc.obj <- suppressWarnings(abc(target,params,sstats,abc.ntol/dim(params)[1],abc.method,hcorr=abc.hcorr,abc.transf,logit.bounds=cbind(params.min,params.max)))\n";
print FILE_1 "if(abc.method == \"loclinear\" && any(is.na(abc.obj\$adj.values))){\n";
print FILE_1 "    cat(\"\\n  Problems while performing the local linear regression!\\n\")\n";
print FILE_1 "    cat(\"  Only presenting estimates from rejection step\\n\")\n";
print FILE_1 "    cat(\"  Try increasing number of simulation and/or tolerance\\n\")\n";
print FILE_1 "    abc.method <- \"rejection\"\n";
print FILE_1 "}\n";
print FILE_1 "cat(\"\\n> ABC method performed, printing results ...\\n\")\n\n";
print FILE_1 "#get estimates from rejection method\n";
print FILE_1 "tab.rej <- summary(abc.obj,unadj=TRUE,intvl=0.95,print=FALSE)\n";
print FILE_1 "write(\"Statistic,Rho,Theta\",file=summ1.fname,app=F)\n";
print FILE_1 "write.table(tab.rej,file=summ1.fname,sep=\",\",col.names=FALSE,app=T)\n";
print FILE_1 "cat(\"\\nEstimates from the rejection step:\\n\")\n";
print FILE_1 "print(tab.rej)\n\n";
print FILE_1 "#get estimates from loclinear method\n";
print FILE_1 "if(abc.method == \"loclinear\"){\n";
print FILE_1 "    tab.reg <- summary(abc.obj,unadj=FALSE,intvl=0.95,print=FALSE)\n";
print FILE_1 "    write(\"Statistic,Rho,Theta\",file=summ2.fname,app=F)\n";
print FILE_1 "    write.table(tab.reg,file=summ2.fname,sep=\",\",col.names=FALSE,app=T)\n";
print FILE_1 "    cat(\"\\nEstimates from the local linear regression step:\\n\")\n";
print FILE_1 "    print(tab.reg)\n"; # MA: print(tab.rej) -> print(tab.reg)
print FILE_1 "}\n\n";
print FILE_1 "#plot posteriors from rejection method\n";
print FILE_1 "niter <- ifelse(dim(params)[1]<1000,dim(params)[1],1000)\n";
print FILE_1 "if(multipage)\n";
print FILE_1 "    pdf(file=paste(fig1.fname,\".pdf\",sep=\"\"))\n";
print FILE_1 "for(i in 1:nparams){\n";
print FILE_1 "    param.name <- colnames(params)[i]\n";
print FILE_1 "    if(!multipage)\n";
print FILE_1 "        pdf(file=paste(fig1.fname,\"_\",param.name,\".pdf\",sep=\"\"))\n";
print FILE_1 "    h.plot <- hist(abc.obj\$unadj.values[,i],plot=F)\n";
print FILE_1 "    maxy <- 1.2*max(h.plot\$density)\n";
print FILE_1 "    hist(abc.obj\$unadj.values[,i],xlim=c(params.min[i],params.max[i]),ylim=c(0,maxy),freq=F,xlab=param.name,main=\"\")\n";
print FILE_1 "    posterior.d <- density(abc.obj\$unadj.values[,i],from=params.min[i],to=params.max[i])\n";
print FILE_1 "    lines(posterior.d,xlim=c(params.min[i],params.max[i]),col=\"blue\")\n";
print FILE_1 "    prior.d <- density(params[1:niter,i],from=params.min[i],to=params.max[i])\n";
print FILE_1 "    lines(prior.d,xlim=c(params.min[i],params.max[i]),col=\"black\")\n";
print FILE_1 "    if(!multipage)\n";
print FILE_1 "        invisible(dev.off())\n";
print FILE_1 "}\n";
print FILE_1 "if(multipage)\n";
print FILE_1 "    invisible(dev.off())\n";
print FILE_1 "cat(\"\\n Plots of posterior distributions from the rejection step created\\n\")\n";
print FILE_1 "#plot posteriors from loclinear method\n";
print FILE_1 "if(abc.method == \"loclinear\"){\n";
print FILE_1 "    if(multipage)\n";
print FILE_1 "        pdf(file=paste(fig2.fname,\".pdf\",sep=\"\"))\n";
print FILE_1 "    for(i in 1:nparams){\n";
print FILE_1 "        param.name <- colnames(params)[i]\n";
print FILE_1 "        if(!multipage)\n";
print FILE_1 "            pdf(file=paste(fig2.fname,\"_\",param.name,\".pdf\",sep=\"\"))\n";
print FILE_1 "        h.plot <- hist(abc.obj\$adj.values[,i],plot=F)\n";
print FILE_1 "        maxy <- 1.2*max(h.plot\$density)\n";
print FILE_1 "        hist(abc.obj\$adj.values[,i],xlim=c(params.min[i],params.max[i]),ylim=c(0,maxy),freq=F,xlab=param.name,main=\"\")\n";
print FILE_1 "        posterior.d <- density(abc.obj\$adj.values[,i],weight=abc.obj\$weights/sum(abc.obj\$weights),from=params.min[i],to=params.max[i])\n";
print FILE_1 "        lines(posterior.d,xlim=c(params.min[i],params.max[i]),col=\"blue\")\n";
print FILE_1 "        prior.d <- density(params[1:niter,i],from=params.min[i],to=params.max[i])\n";
print FILE_1 "        lines(prior.d,xlim=c(params.min[i],params.max[i]),col=\"black\")\n";
print FILE_1 "        if(!multipage)\n";
print FILE_1 "            invisible(dev.off())\n";
print FILE_1 "    }\n";
print FILE_1 "    if(multipage)\n";
print FILE_1 "        invisible(dev.off())\n";
print FILE_1 "    cat(\"\ Plots of posterior distributions from the local linear regression step created\\n\")\n";
print FILE_1 "}\n";
print FILE_1 "######### CHECK ABC ########\n";
print FILE_1 "#check target on summstats\n";
print FILE_1 "if(multipage)\n";
print FILE_1 "    pdf(file=paste(fig3.fname,\".pdf\",sep=\"\"))\n";
print FILE_1 "for(i in 1:dim(sstats)[2]){\n";
print FILE_1 "    if(!multipage)\n";
print FILE_1 "        pdf(file=paste(fig3.fname,\"_\",colnames(sstats)[i],\".pdf\",sep=\"\"))\n";
print FILE_1 "    hist(abc.obj\$ss[,i],xlim=c(sstats.min[i],sstats.max[i]),main=\"\",xlab=colnames(sstats)[i])\n";
print FILE_1 "    abline(v=target[i],lwd=2,lty=1,col=\"blue\")\n";
print FILE_1 "    if(!multipage)\n";
print FILE_1 "        invisible(dev.off())\n";
print FILE_1 "}\n";
print FILE_1 "if(multipage)\n";
print FILE_1 "    invisible(dev.off())\n\n";
print FILE_1 "#check target on summstats 2\n";
print FILE_1 "niter <- ifelse(dim(params)[1]<10000,dim(params)[1],10000)\n";
print FILE_1 "pcaparams <- prcomp(sstats[1:niter,],scale=TRUE)\n";
print FILE_1 "sstats.pca <- pcaparams\$x\n";
print FILE_1 "sstats.rej.pca <- t(apply(abc.obj\$ss,1,function(x) ((t(x)-pcaparams\$center)/pcaparams\$scale)\%*\%pcaparams\$rotation))\n";
print FILE_1 "target.pca <- ((target-pcaparams\$center)/pcaparams\$scale)\%*\%pcaparams\$rotation\n";
print FILE_1 "pdf(file=fig4.fname)\n";
print FILE_1 "plot(sstats.pca[1:niter,1],sstats.pca[1:niter,2],col=\"black\",pch=21,xlab=\"PC1\",ylab=\"PC2\")\n";
print FILE_1 "points(sstats.rej.pca[,1],sstats.rej.pca[,2],col=\"red\",pch=21)\n";
print FILE_1 "points(target.pca[1],target.pca[2],col = \"black\",pch=21)\n";
print FILE_1 "points(target.pca[1],target.pca[2],col = \"yellow\",pch=20)\n";
print FILE_1 "invisible(dev.off())\n\n";
print FILE_1 "#check target on summstats vs params\n";
print FILE_1 "niter <- ifelse(abc.ntol<1000,abc.ntol,1000)\n";
print FILE_1 "if(multipage){\n";
print FILE_1 "    pdf(file=paste(fig5.fname,\".pdf\",sep=\"\"))\n";
print FILE_1 "    oldpar<-par(mfrow=c(1,nparams),mar=c(4,2,4,1))\n";
print FILE_1 "}\n";
print FILE_1 "for(i in 1:dim(sstats)[2]){\n";
print FILE_1 "    if(!multipage){\n";
print FILE_1 "        pdf(file=paste(fig5.fname,\"_\",colnames(sstats)[i],\".pdf\",sep=\"\"))\n";
print FILE_1 "        oldpar<-par(mfrow=c(1,nparams),mar=c(4,2,4,1))\n";
print FILE_1 "    }\n";
print FILE_1 "    for(j in 1:nparams){\n";
print FILE_1 "        param.name <- colnames(params)[j]\n";
print FILE_1 "        plot(abc.obj\$unadj.values[1:niter,j],abc.obj\$ss[1:niter,i],pch=20,ylab=\"\",ylim=c(sstats.min[i],sstats.max[i]),xlab=param.name)\n";
print FILE_1 "        abline(h=target[i],lwd=1,lty=1,col=\"blue\")\n";
print FILE_1 "        mtext(colnames(sstats)[i], side=3, line=-2, outer=TRUE)\n";
print FILE_1 "    }\n";
print FILE_1 "    if(!multipage){\n";
print FILE_1 "        par(oldpar)\n";
print FILE_1 "        invisible(dev.off())\n";
print FILE_1 "    }\n";
print FILE_1 "}    \n";
print FILE_1 "if(multipage){\n";
print FILE_1 "    par(oldpar)\n";
print FILE_1 "    invisible(dev.off())\n";
print FILE_1 "}\n";
print FILE_1 "cat(\"\ Plots to assess ABC analysis integrity created\\n\")\n";
print FILE_1 "############################\n\n";
close (FILE_1);
print "> R output file created: \"$fileRScript\" ...\n\n";


##############################################
##############################################
# Run R script to perform ABC
print "> Performing ABC... \n";
print "This phase requires the programming language R with the library: abc\n";
my $returnSystem = system ("Rscript --vanilla \"$fileRScript\"");


##############################################
##############################################
#### END ####
if ( $returnSystem == 0)  {
	if (-e "$setsdir/MainOutputs")  {
	}
	else  {
		system ("mkdir \"$setsdir/MainOutputs\"");
	}
	system ("mv \"$fileABCSumm1\"* \"$setsdir/MainOutputs\"");
	if (-e  "$fileABCSumm2" )  {
		system ("mv \"$fileABCSumm2\" \"$setsdir/MainOutputs\"");
	}
	system ("mv \"$fileABCFig1\"* \"$setsdir/MainOutputs\"");

	if ( $setsdir !~ /\s/ )  { #No whitespaces
		if ( glob("$fileABCFig2*") )  { 
			system ("mv \"$fileABCFig2\"* \"$setsdir/MainOutputs\"");
		}
	}
	else  { #Has whitespaces
		if ( glob("\"$fileABCFig2\"*") )  { 
			system ("mv \"$fileABCFig2\"* \"$setsdir/MainOutputs\"");
		}
	}

	system ("mv \"$fileABCFig3\"* \"$setsdir/MainOutputs\"");
    system ("mv \"$fileABCFig4\" \"$setsdir/MainOutputs\"");
	system ("mv \"$fileABCFig5\"* \"$setsdir/MainOutputs\"");
	if (-e "$setsdir/OtherOutputs")  {
	}
	else  {
		system ("mkdir \"$setsdir/OtherOutputs\"");
	}
	system ("mv \"$fileRScript\" \"$setsdir/OtherOutputs\"");
	
#	system ("mv \"*.csv\" \"$setsdir/OtherOutputs\""); # MA (move some output files to OtherOutputs)
#	system ("mv \"*.ss\" \"$setsdir/OtherOutputs\""); # MA (move some output files to OtherOutputs)

	
	print "\n> Successful!\n\n";
}
else  {
	print "\n> Error when running R script \"$fileRScript\"!\n\n";
}

system ("cp $setsdir/SSsimulations.csv \"$setsdir/OtherOutputs\"");
system ("cp $setsdir/PSimulations.csv \"$setsdir/OtherOutputs\"");
system ("cp $setsdir/SSRealData.ss \"$setsdir/OtherOutputs\"");

exit;

