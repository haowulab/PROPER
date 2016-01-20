########################################
## some graphical functions to visualize
## the power simulation results
########################################

add.axis1 <- function(y, strata) {
    tmp = axis(1, 1:length(strata), labels=FALSE)
    ypos = -diff(range(y, na.rm=TRUE))*.12
    ## the computation of ypos is tricky!!
    text(1:length(strata)-0.1, ypos, strata,
         srt=270+45, xpd=TRUE, adj=c(0,0))
}

############################################################
## plot power by strata. Stratified by mean expression
############################################################
plotPower <- function(powerOutput, cols=1:ncol(powerOutput$FD),
                      lty=1:ncol(powerOutput$power), main="Power",
                      ylab="Power", leg=TRUE, error.bar=TRUE) {

    ## average over simulations
    nsims = dim(powerOutput$power)[3]
    power = apply(powerOutput$power,c(1,2),mean, na.rm=TRUE)
    power.se = apply(powerOutput$power,c(1,2), sd, na.rm=TRUE) / sqrt(nsims)

    ## remove NAs, if any (this happens when stratifying by dispersion)
    ix.na = apply(power, 1, function(x) all(is.na(x)))
    power = power[!ix.na,]
    power.se = power.se[!ix.na,]
    strata=levels(cut(0,powerOutput$strata))
    strata = strata[!ix.na]

    ## plot
    matplot(power, type="l", lwd=2, col=cols, lty=lty,
            ylim=c(0, max(power, na.rm=TRUE)),xlim=c(0,nrow(power)+1),
            main=main, axes=FALSE, ylab="")
    ## draw error bars
    stratas = 1:dim(powerOutput$power)[1]
    power.lo = power - power.se*1.96
    power.hi = power + power.se*1.96

    if(error.bar) {
        for(j in 1:ncol(power))
            arrows(stratas, power.lo[,j], stratas, power.hi[,j], length=0.05, angle=90,
                   code=3, col=j, lty=1, lwd=1)
    }

    xlab <- switch(powerOutput$stratify.by,
                   expr = "Average count strata",
                   dispersion = "Dispersion strata"
                   )
    mtext(xlab, side=1, line=4)
    mtext(ylab, side=2, line=2.5)

    if(leg) {
        ll = paste0("SS = ", powerOutput$Nreps1, ", ", powerOutput$Nreps2)
        legend("bottomright",ll, col=cols,lty=lty, lwd=2)
    }
    grid()
    add.axis1(power, strata)
    axis(2);  box()
}


#####################################
## plot number of true discoveries
#####################################
plotPowerTD <- function(powerOutput,cols=1:ncol(powerOutput$TD),
                         lty=1:ncol(powerOutput$TD), main="",
                         ylab="# True Discoveries", leg=TRUE, error.bar=TRUE){
    nsims = dim(powerOutput$TD)[3]
    TD = apply(powerOutput$TD,c(1,2), mean, na.rm=TRUE)
    TD.se = apply(powerOutput$TD,c(1,2), sd, na.rm=TRUE) / sqrt(nsims)

    matplot(TD, type="l",col=cols, lty=lty, lwd=2, ylim=c(0,max(TD)+4),xlim=c(0,nrow(TD)+1),
            main=main, axes=FALSE, ylab="")
    TD.lo = TD - TD.se*1.96
    TD.hi = TD + TD.se*1.96
    stratas = 1:nrow(TD)
    if(error.bar) {
        for(j in 1:ncol(TD))
            arrows(stratas, TD.lo[,j], stratas, TD.hi[,j], length=0.05, angle=90,
                   code=3, col=j, lty=1, lwd=1)
    }

    xlab <- switch(powerOutput$stratify.by,
                   expr = "Average count strata",
                   dispersion = "Dispersion strata"
                   )
    mtext(xlab, side=1, line=4)
    mtext(ylab, side=2, line=2.5)

    strata = levels(cut(0, powerOutput$strata))
    add.axis1(TD, strata)
    axis(2);  box()
    grid()
    if(leg) {
        ll = paste0("SS = ", powerOutput$Nreps1, ", ", powerOutput$Nreps2)
        legend("topright",ll, col=c(cols),lty=lty, lwd=2)
    }
}

#####################################
## plot number of false discoveries
#####################################
plotPowerFD <- function(powerOutput,cols=1:ncol(powerOutput$FD),
                         lty=1:ncol(powerOutput$FD), main="",
                         ylab="# False Discoverie", leg=TRUE, error.bar=TRUE){

    nsims = dim(powerOutput$TD)[3]
    strata = levels(cut(0,powerOutput$strata))
    FD = apply(powerOutput$FD,c(1,2), mean, na.rm=TRUE)
    FD.se = apply(powerOutput$FD,c(1,2), sd, na.rm=TRUE) / sqrt(nsims)

    matplot(FD, type="l",col=cols, lty=lty, lwd=2, ylim=c(0,max(FD)*1.1),xlim=c(0,nrow(FD)+1),pch="1",
            main=main, axes=FALSE, ylab="")
    FD.lo = FD - FD.se*1.96
    FD.hi = FD + FD.se*1.96
    stratas = 1:nrow(FD)
    if(error.bar) {
        for(j in 1:ncol(FD))
            arrows(stratas, FD.lo[,j], stratas, FD.hi[,j], length=0.05, angle=90,
                   code=3, col=j, lty=1, lwd=1)
    }

    xlab <- switch(powerOutput$stratify.by,
                   expr = "Average count strata",
                   dispersion = "Dispersion strata"
                   )

    mtext(xlab, side=1, line=4)
    mtext(ylab, side=2, line=2.5)
    add.axis1(FD, strata);
    axis(2);  box(); grid()
    if(leg) {
        ll = paste0("SS = ", powerOutput$Nreps1, ", ", powerOutput$Nreps2)
        legend("topright",ll, col=cols, lty=lty, lwd=2)
    }
}


######################################
## plot FDR
######################################
plotFDR <- function(powerOutput,cols=1:ncol(powerOutput$FDR),
                     lty=1:ncol(powerOutput$FDR), main="",
                     ylab="FDR", leg=TRUE, error.bar=TRUE) {
    strata = levels(cut(0,powerOutput$strata))
    fdr=apply(powerOutput$FDR,c(1,2),mean, na.rm=TRUE)

    ## remove NAs, if any (this happens when stratifying by dispersion)
    ix.na = apply(fdr, 1, function(x) all(is.na(x)))
    fdr = fdr[!ix.na,]
    strata=levels(cut(0,powerOutput$strata))
    strata = strata[!ix.na]

    matplot(fdr, type="l",col=cols, lty=lty, lwd=2, ylim=c(0,max(fdr)*1.2),xlim=c(0,nrow(fdr)+1),
            main=main, axes=FALSE, ylab="")
    xlab <- switch(powerOutput$stratify.by,
                   expr = "Average count strata",
                   dispersion = "Dispersion strata"
                   )
    mtext(xlab, side=1, line=4)
    mtext(ylab, side=2, line=2.5)

    add.axis1(fdr, strata)
    axis(2);  box();  grid()
    if(leg) {
        ll = paste0("SS = ", powerOutput$Nreps1, ", ", powerOutput$Nreps2)
        legend("topright",ll, col=cols, lty=lty, lwd=2)
    }
}



#############################################
## plot FD cost
#############################################
plotFDcost <- function(powerOutput,cols=1:ncol(powerOutput$FD),
                        lty=1:ncol(powerOutput$FD), main="",
                        ylab="False discovery cost", leg=TRUE, error.bar=TRUE){
    strata=levels(cut(0,powerOutput$strata))

##     FD=apply(powerOutput$FD,c(1,2),mean, na.rm=TRUE)
##     TD=apply(powerOutput$TD,c(1,2),mean, na.rm=TRUE)
##     FDC = FD/TD
##     FDC[is.infinite(FDC)] = NA

    FDC.all = powerOutput$FD / powerOutput$TD
    FDC = apply(FDC.all, c(1,2), function(x) mean(x[is.finite(x)], na.rm=TRUE))
    FDC.se = apply(FDC.all, c(1,2), function(x) sd(x[is.finite(x)], na.rm=TRUE)) / sqrt(dim(FDC.all)[3])

    matplot(FDC,type="l", col=cols, lty=lty, lwd=2, ylim=c(0,max(FDC,na.rm=TRUE)*1.2),
            xlim=c(0,nrow(FDC)+1),
            main=main, axes=FALSE, ylab="")
    FDC.lo = FDC - FDC.se*1.96
    FDC.hi = FDC + FDC.se*1.96
    stratas = 1:nrow(FDC)
    if(error.bar) {
        for(j in 1:ncol(FDC))
            arrows(stratas, FDC.lo[,j], stratas, FDC.hi[,j], length=0.05, angle=90,
                   code=3, col=j, lty=1, lwd=1)
    }

    xlab <- switch(powerOutput$stratify.by,
                   expr = "Average count strata",
                   dispersion = "Dispersion strata"
                   )
    mtext(xlab, side=1, line=4)
    mtext(ylab, side=2, line=2.5)

    add.axis1(FDC, strata);
    axis(2);  box();  grid()

    if(leg) {
        ll = paste0("SS = ", powerOutput$Nreps1, ", ", powerOutput$Nreps2)
        legend("topright",ll, col=cols, lty=lty, lwd=2)
    }

}

##########################################
## plot type I error
##########################################
plotPowerAlpha <- function(powerOutput,cols=1:ncol(powerOutput$alpha),
                            lty=1:ncol(powerOutput$alpha),
                            main="", ylab="False positive rate", leg=TRUE, error.bar=TRUE) {
    strata = levels(cut(0,powerOutput$strata))
    alpha=apply(powerOutput$alpha,c(1,2),mean, na.rm=TRUE)
    alpha.se = apply(powerOutput$alpha,c(1,2), sd, na.rm=TRUE) / sqrt(dim(powerOutput$alpha)[3])

    matplot(alpha,type="l", col=cols, lty=lty, lwd=2, ylim=c(0,max(alpha,na.rm=TRUE))*1.2,
            xlim=c(0,nrow(alpha)+1),
            main=main, axes=FALSE, ylab="")
    alpha.lo = alpha - alpha.se*1.96
    alpha.hi = alpha + alpha.se*1.96
    stratas = 1:nrow(alpha)
    if(error.bar) {
        for(j in 1:ncol(alpha))
            arrows(stratas, alpha.lo[,j], stratas, alpha.hi[,j], length=0.05, angle=90,
                   code=3, col=j, lty=1, lwd=1)
    }

    xlab <- switch(powerOutput$stratify.by,
                   expr = "Average count strata",
                   dispersion = "Dispersion strata"
                   )
    mtext(xlab, side=1, line=4)
    mtext(ylab, side=2, line=2.5)

    add.axis1(alpha, strata);  axis(2);
    box(); grid()
    abline(h=powerOutput$powerlist$alpha.nominal)

    if(leg) {
        ll = paste0("SS = ", powerOutput$Nreps1, ", ", powerOutput$Nreps2)
        legend("topright",ll, col=cols, lty=lty, lwd=2)
    }
}

##########################################
## Plot histogram of power
##########################################
plotPowerHist <- function(powerOutput, simResult, main="Histogram of power", return=FALSE){
    Xcut = powerOutput$strata
    strata = levels(cut(0,Xcut))
    nsims = dim(powerOutput$TD)[3]
    N = length(powerOutput$Nreps1)
    table1 = table2 = matrix(NA,length(strata),nsims)
    for(i in 1:nsims){
        id = simResult$DEid[[i]]
        table1[,i] = as.vector(table(cut(simResult$xbar[,N,i],Xcut)))
        table2[,i] = as.vector(table(cut(simResult$xbar[id,N,i],Xcut)))
    }
    table1 = rowMeans(table1)
    table2 = rowMeans(table2)
    table2hist(table1,axes=FALSE, ylab="Counts", main=main)
    table2hist(table2,axes=FALSE, add=TRUE,col="blue")
    add.axis1(c(0,max(table1)), strata);
    axis(2)

    if(return){
        tmp = rbind(table1,table2)
        colnames(tmp) = strata
        rownames(tmp) = c("Number of genes", "Number of true DE genes")
        return(tmp)
    }
}

### make a histogram based on table. This is used by plot.powerHist function.
table2hist<-function(x,add=FALSE,axes=FALSE,...){
    hist1=hist(rnorm(1000), plot=FALSE)
    hist1$xname=""
    hist1$breaks=0:length(x)+.5
    hist1$counts=x;hist1$density=hist1$counts/sum(hist1$counts)
    plot(hist1, add=add, axes=FALSE, ...)
    if(axes){
        axis(1,at=1:length(x),labels=names(x),las=2)
        axis(2)
    }
    if (!add) box()
}

##########################################
## make all power plots in one figures
##########################################
plotAll <- function(powerOutput){
  par(mfrow=c(3,2),mai=c(.67,.5,.2,.1),mgp=c(0,.5,0))
    plotPower(powerOutput,main="");mtext("A",side=3,at=0)
    plotPowerTD(powerOutput);mtext("B",side=3,at=0)
    plotPowerFD(powerOutput);mtext("C",side=3,at=0)
    plotFDcost(powerOutput);mtext("D",side=3,at=0)
    plotFDR(powerOutput);mtext("E",side=3,at=0)
    plotPowerAlpha(powerOutput);mtext("F",side=3,at=0)
}


