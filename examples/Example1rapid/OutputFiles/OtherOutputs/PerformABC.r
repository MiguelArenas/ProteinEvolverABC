# R script to perform the ABC analysis 
# 

# Load libraries
suppressPackageStartupMessages(library(abc))

##### GENERAL SETTINGS #####
nsstats <- 16
nparams <- 2
sstats.set <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
sstats.fname <- "./SSsimulations.csv"
params.fname <- "./PSimulations.csv"
target.fname <- "./SSRealData.ss"
abc.niter <- 250
abc.ntol <- 50
abc.method <- "loclinear"  #loclinear or rejection
abc.transf <- "logit"  #log, logit or none
abc.hcorr <- TRUE
multipage <- TRUE
summ1.fname <- "./Summary_ABC_rejection.csv"
summ2.fname <- "./Summary_ABC_loclinear.csv"
fig1.fname <- "./Posteriors_rejection"
fig2.fname <- "./Posteriors_loclinear"
fig3.fname <- "./Histograms_SStats"
fig4.fname <- "./PCA_SStats.pdf"
fig5.fname <- "./Scaterplots_SStatsVSParams"
############################

#read sstats data
sstats <- matrix(scan(sstats.fname,skip=1,nlines=abc.niter,sep=",",na.strings="--",what=double(),quiet=T),ncol=nsstats,byrow=T)
colnames(sstats) <- scan(sstats.fname,skip=0,nlines=1,sep=",",what="character",quiet=T)
cat(" Summary Statistics file read\n")

#read params data
params <- matrix(scan(params.fname,skip=1,nlines=abc.niter,sep=",",what=double(),quiet=T),ncol=nparams,byrow=T)
colnames(params) <- scan(params.fname,skip=0,nlines=1,sep=",",what="character",quiet=T)
cat(" Parameters file read\n")

#read target data
target <-as.numeric(scan(target.fname,what="raw",quiet=T)[-c(seq(1,2*nsstats-1,2))])

cat(" Sequences data file read\n")

#select sstats set
sstats<-sstats[,sstats.set]
target<-target[sstats.set]

#get rid of NA's and Inf's
sstats.na <- apply(sstats,1,function(x) all(!is.na(x)) && all(is.finite(x)))
if(abc.niter - sum(sstats.na) > 0){
    params<-params[sstats.na,]
    sstats<-sstats[sstats.na,]
    cat(paste("\n   Simulations whose summary statistics have NA or Inf values discarded: ",abc.niter - sum(sstats.na),"\n",sep=""))
}
if(dim(params)[1] < abc.ntol){
    abc.ntol <- dim(params)[1]
    cat("\n   Tolerance number is larger than number of simulations!\n")
    cat("\n   Tolerance set to number of simulations\n")
}

#get range of params
params.min <- apply(params,2,min)
params.max <- apply(params,2,max)
sstats.min <- apply(sstats,2,function(x) quantile(x,names=F,probs=0.001))
sstats.max <- apply(sstats,2,function(x) quantile(x,names=F,probs=0.999))

####### ABC ANALYSIS #######
abc.obj <- suppressWarnings(abc(target,params,sstats,abc.ntol/dim(params)[1],abc.method,hcorr=abc.hcorr,abc.transf,logit.bounds=cbind(params.min,params.max)))
if(abc.method == "loclinear" && any(is.na(abc.obj$adj.values))){
    cat("\n  Problems while performing the local linear regression!\n")
    cat("  Only presenting estimates from rejection step\n")
    cat("  Try increasing number of simulation and/or tolerance\n")
    abc.method <- "rejection"
}
cat("\n> ABC method performed, printing results ...\n")

#get estimates from rejection method
tab.rej <- summary(abc.obj,unadj=TRUE,intvl=0.95,print=FALSE)
write("Statistic,Rho,Theta",file=summ1.fname,app=F)
write.table(tab.rej,file=summ1.fname,sep=",",col.names=FALSE,app=T)
cat("\nEstimates from the rejection step:\n")
print(tab.rej)

#get estimates from loclinear method
if(abc.method == "loclinear"){
    tab.reg <- summary(abc.obj,unadj=FALSE,intvl=0.95,print=FALSE)
    write("Statistic,Rho,Theta",file=summ2.fname,app=F)
    write.table(tab.reg,file=summ2.fname,sep=",",col.names=FALSE,app=T)
    cat("\nEstimates from the local linear regression step:\n")
    print(tab.reg)
}

#plot posteriors from rejection method
niter <- ifelse(dim(params)[1]<1000,dim(params)[1],1000)
if(multipage)
    pdf(file=paste(fig1.fname,".pdf",sep=""))
for(i in 1:nparams){
    param.name <- colnames(params)[i]
    if(!multipage)
        pdf(file=paste(fig1.fname,"_",param.name,".pdf",sep=""))
    h.plot <- hist(abc.obj$unadj.values[,i],plot=F)
    maxy <- 1.2*max(h.plot$density)
    hist(abc.obj$unadj.values[,i],xlim=c(params.min[i],params.max[i]),ylim=c(0,maxy),freq=F,xlab=param.name,main="")
    posterior.d <- density(abc.obj$unadj.values[,i],from=params.min[i],to=params.max[i])
    lines(posterior.d,xlim=c(params.min[i],params.max[i]),col="blue")
    prior.d <- density(params[1:niter,i],from=params.min[i],to=params.max[i])
    lines(prior.d,xlim=c(params.min[i],params.max[i]),col="black")
    if(!multipage)
        invisible(dev.off())
}
if(multipage)
    invisible(dev.off())
cat("\n Plots of posterior distributions from the rejection step created\n")
#plot posteriors from loclinear method
if(abc.method == "loclinear"){
    if(multipage)
        pdf(file=paste(fig2.fname,".pdf",sep=""))
    for(i in 1:nparams){
        param.name <- colnames(params)[i]
        if(!multipage)
            pdf(file=paste(fig2.fname,"_",param.name,".pdf",sep=""))
        h.plot <- hist(abc.obj$adj.values[,i],plot=F)
        maxy <- 1.2*max(h.plot$density)
        hist(abc.obj$adj.values[,i],xlim=c(params.min[i],params.max[i]),ylim=c(0,maxy),freq=F,xlab=param.name,main="")
        posterior.d <- density(abc.obj$adj.values[,i],weight=abc.obj$weights/sum(abc.obj$weights),from=params.min[i],to=params.max[i])
        lines(posterior.d,xlim=c(params.min[i],params.max[i]),col="blue")
        prior.d <- density(params[1:niter,i],from=params.min[i],to=params.max[i])
        lines(prior.d,xlim=c(params.min[i],params.max[i]),col="black")
        if(!multipage)
            invisible(dev.off())
    }
    if(multipage)
        invisible(dev.off())
    cat(" Plots of posterior distributions from the local linear regression step created\n")
}
######### CHECK ABC ########
#check target on summstats
if(multipage)
    pdf(file=paste(fig3.fname,".pdf",sep=""))
for(i in 1:dim(sstats)[2]){
    if(!multipage)
        pdf(file=paste(fig3.fname,"_",colnames(sstats)[i],".pdf",sep=""))
    hist(abc.obj$ss[,i],xlim=c(sstats.min[i],sstats.max[i]),main="",xlab=colnames(sstats)[i])
    abline(v=target[i],lwd=2,lty=1,col="blue")
    if(!multipage)
        invisible(dev.off())
}
if(multipage)
    invisible(dev.off())

#check target on summstats 2
niter <- ifelse(dim(params)[1]<10000,dim(params)[1],10000)
pcaparams <- prcomp(sstats[1:niter,],scale=TRUE)
sstats.pca <- pcaparams$x
sstats.rej.pca <- t(apply(abc.obj$ss,1,function(x) ((t(x)-pcaparams$center)/pcaparams$scale)%*%pcaparams$rotation))
target.pca <- ((target-pcaparams$center)/pcaparams$scale)%*%pcaparams$rotation
pdf(file=fig4.fname)
plot(sstats.pca[1:niter,1],sstats.pca[1:niter,2],col="black",pch=21,xlab="PC1",ylab="PC2")
points(sstats.rej.pca[,1],sstats.rej.pca[,2],col="red",pch=21)
points(target.pca[1],target.pca[2],col = "black",pch=21)
points(target.pca[1],target.pca[2],col = "yellow",pch=20)
invisible(dev.off())

#check target on summstats vs params
niter <- ifelse(abc.ntol<1000,abc.ntol,1000)
if(multipage){
    pdf(file=paste(fig5.fname,".pdf",sep=""))
    oldpar<-par(mfrow=c(1,nparams),mar=c(4,2,4,1))
}
for(i in 1:dim(sstats)[2]){
    if(!multipage){
        pdf(file=paste(fig5.fname,"_",colnames(sstats)[i],".pdf",sep=""))
        oldpar<-par(mfrow=c(1,nparams),mar=c(4,2,4,1))
    }
    for(j in 1:nparams){
        param.name <- colnames(params)[j]
        plot(abc.obj$unadj.values[1:niter,j],abc.obj$ss[1:niter,i],pch=20,ylab="",ylim=c(sstats.min[i],sstats.max[i]),xlab=param.name)
        abline(h=target[i],lwd=1,lty=1,col="blue")
        mtext(colnames(sstats)[i], side=3, line=-2, outer=TRUE)
    }
    if(!multipage){
        par(oldpar)
        invisible(dev.off())
    }
}    
if(multipage){
    par(oldpar)
    invisible(dev.off())
}
cat(" Plots to assess ABC analysis integrity created\n")
############################

