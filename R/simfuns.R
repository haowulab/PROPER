#####################################################
## functions for simulating RNA-seq count data.
## Take an option object and sample sizes.
#####################################################

simRNAseq <- function(simOptions, n1, n2) {
    ##  set.seed(simOptions$sim.seed)
    if(simOptions$design == "2grp")
        data <- simRNAseq.2grp(simOptions, n1, n2)

    data
}

#####################################################
## simulate 2-group RNA-seq data
#####################################################
simRNAseq.2grp <- function(simOptions, n1, n2) {

    ## make design vector
    design <- c(rep(0, n1), rep(1, n2))

    ## generate ID for DE genes
    ngenes = simOptions$ngenes
    p.DE = simOptions$p.DE
    DEid = sample(ngenes, round(ngenes*p.DE))

    ## generate lfc for all genes
    lfc = rep(0, ngenes)
    lfc[DEid] = simOptions$lfc

    ## generate mean expressions for all replicates
    lmeanExpr <- makeMeanExpr.2grp(simOptions$lBaselineExpr, lfc, n1, n2)

    ## generate counts
    allmu <- lmeanExpr
    phi <- exp(simOptions$lOD)
    x <- rnegbinom(length(allmu), allmu, phi)
    x <- matrix(x,ncol=n1+n2)

    ## return
    list(counts=x, designs=design, DEid=DEid, simOptions=simOptions)
}

###############################################################
## generate mean expression values in two group comparison.
###############################################################
makeMeanExpr.2grp <- function(lBaselineExpr, lfc, n1, n2) {
    ## generate mean expressions in two groups
    result <- matrix(exp(lBaselineExpr+c(rep(lfc/2,n1),rep(-lfc/2,n2))),
                     ncol=n1+n2)
    ##result <- sweep(result, 2, sizefactor, FUN="*")
    result
}


###############################################################
## run simulation and DE detection
###############################################################
runSims <- function(Nreps=c(3,5,7,10), Nreps2, nsims=100, sim.opts,
                    DEmethod=c("edgeR", "DSS", "DESeq", "DESeq2"),
                    verbose=TRUE) {

    DEmethod = match.arg(DEmethod)
    ## generate size factor if not given
    ##   if(missing(sizefactor)) {
    ##     sizefactor=rep(1, n1+n2)
    ##   } else if(!is.vector(sizefactor) | length(sizefactor) != (m.n1+m.n2) ) {
    ##     stop("sizefactor must be a vector of length n1+n2!\n")
    ##   }

    if(missing(Nreps2))
        Nreps2 = Nreps
    else {
        if(length(Nreps2) != length(Nreps))
            stop("Nreps and Nreps2 must be vectors of the same length")
    }
    n1 = max(Nreps)
    n2 = max(Nreps2)

    ## start simulation
    set.seed(sim.opts$sim.seed)
    pvalue = fdrs = xbar = array(NA,dim=c(sim.opts$ngenes,length(Nreps), nsims))
    DEids = lfcs = NULL
    for(i in 1:nsims) {
        if(verbose)
            cat("Simulation number", i, "\n")

        ## update the simulatino option. Reregenerate the DEid and lfc
        sim.opts$sim.seed = sim.opts$sim.seed + 1
        sim.opts = update.RNAseq.SimOptions.2grp(sim.opts)

        ## generate data
        ##sim.opts$sim.seed = sim.opts$sim.seed + 1
        dat.sim.big = simRNAseq(sim.opts, n1, n2)
        DEids[[i]] = dat.sim.big$DEid
        lfcs[[i]] = dat.sim.big$simOptions$lfc

        ##  for different sample sizes
        for(j in seq(along=Nreps)) {
            nn1 = Nreps[j]
            nn2 = Nreps2[j]

            ## take a subsample of the simulated counts
            idx = c(1:nn1, n1+(1:nn2))
            this.design = dat.sim.big$designs[idx]
            this.X = dat.sim.big$counts[,idx]
            this.simOpts = sim.opts
            ## filter out genes with all 0 counts
            ss = rowSums(this.X)
            ix.valid = ss>0
            this.X.valid = this.X[ix.valid,, drop=FALSE]

            ## create an object and pass into DE detection
            data0=list(counts=this.X.valid, designs=this.design)
            if (DEmethod == "edgeR")
                res1 = run.edgeR(data0)
            if (DEmethod == "DESeq")
                res1=run.DESeq(data0)
            if (DEmethod == "DSS")
                res1=run.DSS(data0)
            if (DEmethod == "DESeq2")
                res1=run.DESeq2(data0)

            ## store results. Be careful here about filtering
            pval = fdr = rep(1, nrow(this.X))
            X.bar1 = rep(0, nrow(this.X))
            pval[ix.valid] = res1[, "pval"]
            fdr[ix.valid] = res1[, "fdr"]
            sizeF = colSums(data0$count)
            sizeF = sizeF/median(sizeF) #size factor
            X.bar1[ix.valid] = rowMeans(sweep(data0$count,2,sizeF,FUN="/"))
            pvalue[,j,i] = pval
            fdrs[,j,i] = fdr
            xbar[,j,i]=X.bar1
        }
    }

    ## return
    list(pvalue=pvalue, fdrs=fdrs, xbar=xbar, DEid=DEids, lfcs=lfcs, Nreps1=Nreps, Nreps2=Nreps2, sim.opts=sim.opts)
}


###########################################################################
## update the simulation option object.
## This is to regenreate a new set of DE genes.
## They have different (regenreated) parameters.
###########################################################################
update.RNAseq.SimOptions.2grp <- function(sim.opts) {
    ## update DEid and lfc, but keep the others
    RNAseq.SimOptions.2grp(ngenes = sim.opts$ngenes,
                           lBaselineExpr=sim.opts$lBaselineExpr,
                           lOD=sim.opts$lOD,
                           p.DE=sim.opts$p.DE, lfc=sim.opts$lfc,
                           sim.seed=sim.opts$sim.seed)
}
