# R script to apply priors and make executable files "ProteinEvolverABC_Phase2.sh" and "ProteinEvolverABC_arguments.txt" for ProteinEvolverProtABC1.2.0
# By Miguel Arenas
# 2020
# 

## Distributions included ##
# fix			# A value is fixed	(integer or non integer)			# v <- 2 
# sample		# Sample (integer): lowest highest 						# v <- sample(0:23,1,replace=T)
# unif			# Uniform (non integer): lowest highest; t truncated	# v <- runif(1,0.1,0.9) // v <- runif(1,0.1,0.9) + t
# norm			# Normal: mean, sd; t truncated							# v <- rnorm(1,0.9,0.1) // v <- rnorm(1,0.9,0.1) + t
# exp			# Exponential: rate; t truncated						# v <- rexp(1,0.4) // v <- rexp(1,0.4) + t
# gamma			# Gamma: shape, rate (1/scale); t truncated				# v <- rgamma(1,0.3,0.4) // v <- rgamma(1,0.3,0.4) + t
# beta			# Beta: shape1, shape2 (1/scale); t truncated			# v <- rbeta(1,0.3,0.4) // v <- rbeta(1,0.3,0.4) + t
# dirichlet		# Dirichlet: alpha (vector)								# v <- rdirichlet(1, c(1,1,1,1))

# Load libraries
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(MCMCpack))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(graphics))

unlink ("./ProteinEvolverABC_arguments.txt", recursive = FALSE)
unlink ("./ProteinEvolverABC_Pop_arguments.txt", recursive = FALSE)
unlink ("./ProteinEvolverABC_Phase2.sh", recursive = FALSE)
##### GENERAL SETTINGS #####
NumberReplicates <- 250

recombination_rate <- numeric(NumberReplicates)
substitution_rate <- numeric(NumberReplicates)
RHO <- numeric(NumberReplicates)
THETA <- numeric(NumberReplicates)


##### PARAMETERS FOR THE ENTIRE PROTEIN are directly specified in the print file #####

### COALESCENT SIMULATIONS ###
# ./ProteinEvolverProtABC1.2.0 -n1 -s8 765 -e1000 2 -r2.0e-06 -u5.6e-05 -@JTT -c1 1 0 -y1 -:1
# ./ProteinEvolverProtABC1.2.0 -n2 -s8 150 -e1000 2 -=4 1995 1 1 2003 4 6 1997 2 3 2001 7 8 -/1200 -g1 3 1000 1250 1000 1300 1550 2000 1560 1000 3000 -q1 4 2 2 3 1 -t3 100 800 0.002 0.001 0.003 -%1 1 2 10000 -r2.3e-6 -o0.1 -u4.1e-5 -f20 0.04 0.06 0.05 0.05 0.08 0.02 0.05 0.05 0.03 0.07 0.04 0.06 0.05 0.05 0.05 0.05 0.05 0.05 0.04 0.06 -@WAG -a0.7 -i0.52 -bsequences -c1 1 0 -y1 -:1

##### PARAMETERS FOR EACH SIMULATION #####

ThisReplicate<-0
while (ThisReplicate < NumberReplicates)
	{
	ThisReplicate<-ThisReplicate+1

	SampleSize_SequenceLength_print <- paste (" -s35 604",sep="")
	DatedTips_print <- paste ("",sep="")
	DemogPeriods_print <- paste ("",sep="")
	DemogPeriodsNcte_print <- paste ("",sep="")
	DemogPeriodsNvar_print <- paste ("",sep="")
	Migrationmodel_print <- paste ("",sep="")
	MigrationRate_print <- paste ("",sep="")
	ConvDemes_print <- paste ("",sep="")
	Hap_Dip <- paste (" 2",sep="")
	PopSize <- 1000				# unif or fix
	PopSize_Hap_Dip_print <- paste (" -e",PopSize,Hap_Dip,sep="")
	GenTime_print <- paste ("",sep="")
	GrowthRate_print <- paste ("",sep="")
	outgroup_print <- paste ("",sep="")
	HomogRec <- runif(1,0,5.56e-6)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)
	HomogRec_print <- paste (" -r",HomogRec,sep="")
	recombination_rate[ThisReplicate]<-HomogRec
	SubsRate <- runif(1,0,1.67e-4)		# fix, unif, norm(t), exp(t), gamma(t), beta(t)
	SubsRate_print <- paste (" -u",SubsRate,sep="")
	substitution_rate[ThisReplicate]<-SubsRate

	AAFreqs1 <- c(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05)		

	AAFreqs <- numeric(20)		
	AAFreqs[1] <- AAFreqs1[1]		
	AAFreqs[2] <- AAFreqs1[2]		
	AAFreqs[3] <- AAFreqs1[3]		
	AAFreqs[4] <- AAFreqs1[4]		
	AAFreqs[5] <- AAFreqs1[5]		
	AAFreqs[6] <- AAFreqs1[6]		
	AAFreqs[7] <- AAFreqs1[7]		
	AAFreqs[8] <- AAFreqs1[8]		
	AAFreqs[9] <- AAFreqs1[9]		
	AAFreqs[10] <- AAFreqs1[10]		
	AAFreqs[11] <- AAFreqs1[11]		
	AAFreqs[12] <- AAFreqs1[12]		
	AAFreqs[13] <- AAFreqs1[13]		
	AAFreqs[14] <- AAFreqs1[14]		
	AAFreqs[15] <- AAFreqs1[15]		
	AAFreqs[16] <- AAFreqs1[16]		
	AAFreqs[17] <- AAFreqs1[17]		
	AAFreqs[18] <- AAFreqs1[18]		
	AAFreqs[19] <- AAFreqs1[19]		
	AAFreqs[20] <- AAFreqs1[20]		
	AminoacidFrequencies_print <- paste (" -f20 ",AAFreqs[1]," ",AAFreqs[2]," ",AAFreqs[3]," ",AAFreqs[4]," ",AAFreqs[5]," ",AAFreqs[6]," ",AAFreqs[7]," ",AAFreqs[8]," ",AAFreqs[9]," ",AAFreqs[10]," ",AAFreqs[11]," ",AAFreqs[12]," ",AAFreqs[13]," ",AAFreqs[14]," ",AAFreqs[15]," ",AAFreqs[16]," ",AAFreqs[17]," ",AAFreqs[18]," ",AAFreqs[19]," ",AAFreqs[20],sep="")

	SubsModel_print <- paste (" -@","JTT",sep="")
	AArateHetSites_print <- paste ("",sep="")
	AAPinv_print <- paste ("",sep="")


	ExecutionHeader <- paste ("-n1",SampleSize_SequenceLength_print,PopSize_Hap_Dip_print,DatedTips_print,GenTime_print,GrowthRate_print,DemogPeriods_print,DemogPeriodsNcte_print,DemogPeriodsNvar_print,Migrationmodel_print,MigrationRate_print,ConvDemes_print,HomogRec_print,SubsRate_print,outgroup_print,AminoacidFrequencies_print,SubsModel_print,AArateHetSites_print,AAPinv_print," -bsequences -c1 1 0 -y0 -:",ThisReplicate,sep="")

	write(paste("",ExecutionHeader,"",sep=""),"./ProteinEvolverABC_arguments.txt",append=T)
	Rho_sim <- 2 * 2 * PopSize * HomogRec * 604 
	HomogRec_print_Pop <- paste (" -r",Rho_sim,sep="")

	Theta_sim <- 2 * 2 * PopSize * SubsRate * 604 
	SubsRate_print_pop <- paste (" -u",Theta_sim,sep="")
	RHO[ThisReplicate]<-Rho_sim
	THETA[ThisReplicate]<-Theta_sim
	ExecutionHeader_Pop <- paste ("-n1",SampleSize_SequenceLength_print,PopSize_Hap_Dip_print,DatedTips_print,GenTime_print,GrowthRate_print,DemogPeriods_print,DemogPeriodsNcte_print,DemogPeriodsNvar_print,Migrationmodel_print,MigrationRate_print,ConvDemes_print,HomogRec_print_Pop,SubsRate_print_pop,outgroup_print,AminoacidFrequencies_print,SubsModel_print,AArateHetSites_print,AAPinv_print," -bsequences -c1 1 0 -y0 -:",ThisReplicate,sep="")

	write(paste("",ExecutionHeader_Pop,"",sep=""),"./ProteinEvolverABC_Pop_arguments.txt",append=T)

} # end replicates
	NumberOfProcessors <- 1
	NumberOfArguments <- 35
	bodyhereNew<-paste("more \"./ProteinEvolverABC_arguments.txt\" | xargs -n ",NumberOfArguments," \"./bin/ProteinEvolverProtABC1.2.0\"",sep="")
	write(paste("",bodyhereNew,"",sep=""),"./ProteinEvolverABC_Phase2.sh",append=T)


# Printing prior distributions for Recombination and Substitution rates 

figureName1<-paste("./Histogram_PriorRecombination.pdf",sep="")
pdf(figureName1)
par(mfrow = c(1,1))

if (NumberReplicates < 1000)	{
	hist (recombination_rate, breaks=NumberReplicates)
	} else {
	hist (recombination_rate, breaks=1000)
	}
dev.off()

figureName2<-paste("./Histogram_PriorSubstitution.pdf",sep="")
pdf(figureName2)
par(mfrow = c(1,1))

if (NumberReplicates < 1000)	{
	hist (substitution_rate, breaks=NumberReplicates)
	} else {
	hist (substitution_rate, breaks=1000)
	}
dev.off()

figureName4<-paste("./Histogram_PriorRho.pdf",sep="")
pdf(figureName4)
par(mfrow = c(1,1))

if (NumberReplicates < 1000)	{
	hist (RHO, breaks=NumberReplicates)
	} else {
	hist (RHO, breaks=1000)
	}
dev.off()

figureName5<-paste("./Histogram_PriorTheta.pdf",sep="")
pdf(figureName5)
par(mfrow = c(1,1))

if (NumberReplicates < 1000)	{
	hist (THETA, breaks=NumberReplicates)
	} else {
	hist (THETA, breaks=1000)
	}
dev.off()
