######################################################
## A list of functions for power calculation
######################################################

## With simulation results, this function compute stratified power, etc.
comparePower <- function(simOutput, alpha.type=c("fdr","pval"), alpha.nominal=0.1,
                         stratify.by=c("expr", "dispersion"), strata,
                         filter.by=c("none", "expr"), strata.filtered=1,
                         target.by=c("lfc", "effectsize"), delta=0.5) {

    alpha.type = match.arg(alpha.type)
    stratify.by = match.arg(stratify.by)
    filter.by = match.arg(filter.by)
    target.by = match.arg(target.by)

    ## set strata if not given
    if(missing(strata)) {
        if(stratify.by == "expr")
            strata = c(0,10,2^(1:7)*10,Inf)
        else
            strata = c(0, c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1), Inf)
    }

    ## some general parameters
    Nreps1 = simOutput$Nreps1
    Nreps2 = simOutput$Nreps2
    ngenes = simOutput$sim.opts$ngenes
    sim.opts = simOutput$sim.opts
    DEids = simOutput$DEid
    lfcs = simOutput$lfcs
    nsims = dim(simOutput$pvalue)[3]

    ## initialize results
    ## determine dimension of results, for filtering
    nr = length(strata) - 1
    if(filter.by == "expr") {
        nr = nr - length(strata.filtered)
    }
    TD = FD = FDR = alpha = power = array(NA,dim=c(nr,length(Nreps1), nsims))
    power.marginal = alpha.marginal = FDR.marginal = matrix(NA,length(Nreps1), nsims)

    ## loop over simulation and replicates
    for(i in 1:nsims) {
        for(j in seq(along=Nreps1)) {
            nn1 = Nreps1[j]
            nn2 = Nreps2[j]
            ## get DE flags.
            DEid = simOutput$DEid[[i]]
            lfc = simOutput$lfcs[[i]]
            Zg = Zg2 = rep(0, ngenes)
            Zg[DEid] = 1
            ## find target (interesting) genes
            if(target.by == "lfc") {
                lfc = simOutput$lfcs[[i]]
                ix = abs(lfc) > delta
            } else if (target.by == "effectsize") {
                effectsize = simOutput$lfcs[[i]] / sqrt(exp(simOutput$sim.opts$lOD[DEid]))
                ix = abs(effectsize) > delta
            }
            Zg2[DEid[ix]] = 1

            ## stratification and filtering
            X.bar1 = simOutput$xbar[,j,i]
            ix.keep = which(X.bar1>0)
            if(stratify.by == "expr") { ## stratify by expression
                xgr = cut(X.bar1[ix.keep], strata)
                ## filtering. This only applies when stratifying by expression
                if(filter.by == "expr") { ## filter by expression levels (average counts)
                    lev = levels(xgr)
                    ix.keep = ix.keep[!(xgr %in% lev[strata.filtered])]
                    ## recut
                    xgr = cut(X.bar1[ix.keep], strata[-strata.filtered])
                }
            }
            else if(stratify.by == "dispersion") { ## stratify by dispersion
                xgr=cut(simOutput$sim.opts$lOD[ix.keep], log(strata))
                if(filter.by != "none") { ## filter
                    stop("Filtering only applies when stratifying by expression level.")
                }
            }

            ## get type I error
            if(alpha.type == "pval")
                x = simOutput$pvalue[ix.keep,j,i]
            else { ## use FDR.
                x = simOutput$fdr[ix.keep,j,i]
                if(filter.by != "none") { ## Need to recompute FDR after filtering. I'll just use BH
                    pval = simOutput$pvalue[ix.keep,j,i]
                    x = p.adjust(pval, "BH")
                }
            }

            ## update Zg flags after filtering
            Zg = Zg[ix.keep]
            Zg2 = Zg2[ix.keep]

            ## calculate stratified power-related quantities
            power00 = POWER1(x, p.crit=alpha.nominal, Zg, Zg2, xgr=xgr)
            TD[,j,i] = power00$TD
            FD[,j,i] = power00$FD
            alpha[,j,i] = power00$alpha
            alpha.marginal[j,i] = power00$alpha.marginal
            power[,j,i] = power00$power
            power.marginal[j,i] = power00$power.marginal
            FDR[,j,i] = power00$FDR
            FDR.marginal[j,i] = power00$FDR.marginal
        }
    }

    output <- list(TD=TD, FD=FD, FDR=FDR, alpha=alpha, power=power,
                   alpha.marginal=alpha.marginal, power.marginal=power.marginal, FDR.marginal=FDR.marginal,
                   ## below are input parameters
                   alpha.type=alpha.type, alpha.nominal=alpha.nominal,
                   stratify.by=stratify.by, strata=strata,
                   target.by=target.by, Nreps1=simOutput$Nreps1, Nreps2=simOutput$Nreps2,
                   delta=delta)

    output
}


########################################################
## Power function, compute the stratified power-realted quantities
##
## Return values:
## TD: number of True Discoveries in each stratum
## FD: number of False Discoveries in each stratum
## power: within strata, proportion of TD out of total DE
## alpha.nomial: cutoff of alpha on raw p-values
## alpha: empirical p-value in each stratum
## alpha.marginal: overall empirical p-value
########################################################

POWER1 <- function(p, p.crit, Zg, Zg2, xgr){
    ## p is input nominal p-value or FDR.
    ## alpha is cutoff of p
    ## !Zg is indicator for TN,  Zg2 is indicators for TP,
    ## xgr is grouping in covariate

    ix.D = p <= p.crit
    N = sum(ix.D) ## this is the total discovery
    N.stratified = tapply(ix.D, xgr, sum)

    ##  TD
    id.TP = Zg2==1
    TD = tapply(p[id.TP] <= p.crit, xgr[id.TP], sum)
    TD[is.na(TD)] = 0

    ##  FD
    id.TN = Zg==0
    FD = tapply(p[id.TN] <= p.crit, xgr[id.TN], sum)
    FD[is.na(FD)] = 0

    ## type I error
    alpha = as.vector(FD/table(xgr[id.TN]))
    alpha.marginal = sum(FD)/sum(id.TN)

    ## power
    power=as.vector(TD/table(xgr[id.TP]))
    power.marginal=sum(TD,na.rm=TRUE)/sum(id.TP)

    ## FDR
    FDR = FD / N.stratified
    FDR.marginal = sum(FD, na.rm=TRUE) / N

    list(TD=TD,FD=FD,alpha.nominal=p.crit,
         alpha=alpha, alpha.marginal=alpha.marginal,
         power=power, power.marginal=power.marginal,
         FDR=FDR, FDR.marginal=FDR.marginal)
}


###########################################################################
## summary the power calculation result.
## The result is an object from comparePower function.
###########################################################################
summaryPower <- function(powerOutput) {
    nn1 <- powerOutput$Nreps1
    nn2 <- powerOutput$Nreps2

    alpha.type <- powerOutput$alpha.type
    if(alpha.type == "pval") {
        alpha.nam <- "type I error"
        alpha.mar <- rowMeans(powerOutput$alpha.marginal)
    }  else {
        alpha.nam <- "FDR"
        alpha.mar <- rowMeans(powerOutput$FDR.marginal)
    }

    TD.avg <- colSums(apply(powerOutput$TD,c(1,2), mean, na.rm=TRUE))
    FD.avg <- colSums(apply(powerOutput$FD,c(1,2), mean, na.rm=TRUE))

    res <- cbind(nn1, nn2, powerOutput$alpha.nominal, alpha.mar, rowMeans(powerOutput$power.marginal),
                 TD.avg, FD.avg, FD.avg/TD.avg)
    colnames(res) <- c("SS1", "SS2",  paste(c("Nominal", "Actual"), alpha.nam),
                       "Marginal power", "Avg # of TD", "Avg # of FD", "FDC")

    print(signif(res,2))
    return(invisible(res))

}



#######################################################
## recompute power based on sequencing depth
#######################################################
power.seqDepth <- function(simResult, powerOutput,
                           depth.factor=c(0.2, 0.5, 1, 2, 5, 10)) {
    ## first get stratified power
    power.expr <- apply(powerOutput$power,c(1,2),mean, na.rm=TRUE)

    dd = dim(powerOutput$TD)
    nsims = dd[3]
    nSS = dd[2]

    res = matrix(0, nrow=length(depth.factor), ncol=nSS)
    rownames(res) = depth.factor
    target.by = powerOutput$target.by
    delta = powerOutput$delta
    for(idepth in 1:length(depth.factor)) {
        power.adj = matrix(0, nrow=nsims, ncol=nSS)
        for(isim in 1:nsims) {
            ## find biologically interesting DE genes.
            DEid = simResult$DEid[[isim]]
            lfc = simResult$lfcs[[isim]]
            if(target.by == "lfc") {
                ix = abs(lfc) > delta
            } else if (target.by == "effectsize") {
                effectsize = lfc / sqrt(exp(simResult$sim.opts$lOD[DEid]))
                ix = abs(effectsize) > delta
            }
            DEid2 = DEid[ix]

            for(iSS in 1:nSS) {
                X.bar1 = simResult$xbar[DEid2, iSS, isim] * depth.factor[idepth]
                K=which(X.bar1>0)
                xgr = cut(X.bar1[K], powerOutput$strata)
                tt = table(xgr)
                tt = tt/sum(tt)
                power.expr = powerOutput$power[,iSS,isim]
                power.adj[isim, iSS] = sum(tt*power.expr)
            }
        }
        res[idepth, ] = colMeans(power.adj)
    }

    colnames(res) = paste0("SS=", powerOutput$Nreps1, ",", powerOutput$Nreps2)
    print(signif(res,2))
    return(invisible(res))
}

