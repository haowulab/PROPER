################################################################
## interface for simulation options.
## Available options include:
## - ngenes: number of genes.
##
## - lBaselineExpr: Baseline expressions, in *log* scale,  for all genes. It could be:
##   (1) a constant, so that all genes share the same expression.
##   (2) a vector of length ngenes.
##   (3) a function taking a single parameter ngenes, so that we can sample from it.
##   (4) empirically sample from MAQC/Gilad/Cheung data.
##   By default it'll sample from Cheung data.
##
## - seqDepth: sequencing depth. Total number of reads for each experiment.
##   Assume all experiments have the same depth. This will be ignored if lBaselineExpr
##   is provided.
##
## - lOD: over-dispersion parameters, in *log* scale,  for all genes. It could be:
##   (1) a constant, so that all genes share the same OD.
##   (2) a vector of length ngenes.
##   (3) a function taking a single parameter ngenes, so that we can sample from it.
##   (4) empirically sample from MAQC/Gilad/Cheung data.
##   By default it'll sample from Cheung data.
##
## - p.DE: proportion of DE genes, 5% by default.
##
## - lfc: fold change, in *log* scale,  between two groups for DE genes. It could be:
##   (1) a constant, so that all DE genes share the same fold change.
##   (2) a vector of length ngenes*p.DE.
##   (3) a function taking a single parameter ngenes, so that we can sample from it.
##

### Need to do: add DEid as an option!!!
#######################################################

## return all current options
RNAseq.SimOptions.2grp <- function(ngenes=50000, seqDepth, lBaselineExpr, lOD,
                                   p.DE=0.05, lfc, sim.seed) {

    if(missing(sim.seed))
        sim.seed = 11111
    set.seed(sim.seed)
    param = NULL
    ## set Baseline expression
    if(missing(lBaselineExpr)) {
        if(missing(seqDepth)) {
            data(cheung, envir=environment())
            lBaselineExpr=sample(param$lmean, ngenes, replace=TRUE)
        } else { ## provide sequencing depth,
            lBaselineExpr=setBaselineExpr.seqDepth(seqDepth, ngenes)
        }
    } else {
        lBaselineExpr=setBaselineExpr(lBaselineExpr, ngenes)
    }

    ## set over dispersion parameter
    if(missing(lOD)) {
        data(cheung, envir=environment())
        lOD=sample(param$lOD, ngenes, replace=TRUE)
    } else {
        lOD=setOD(lOD, ngenes)
    }

    ## set up fold change for alternative genes
    nDE = round(ngenes*p.DE)
    if(missing(lfc)) {
        lfc = lfc.alt(nDE)
    }
    lfc = setFC(lfc, nDE)


    ##   ## size factor
    ##   m.n1=max(n1); m.n2=max(n2)
    ##   if(missing(sizefactor)) {
    ##     data(cheung)
    ##     sizefactor=rep(1,m.n1+m.n2)
    ##   } else if(!is.vector(sizefactor) | length(sizefactor) != (m.n1+m.n2) ) {
    ##     stop("sizefactor must be a vector of length n1+n2!\n")
    ##   }

    ##   ## generate mean expression
    ##   lmeanExpr <- makeMeanExpr.2grp(ngenes, DEid, lBaselineExpr, lfc, m.n1, m.n2, sizefactor)

    ## return
    list(ngenes=ngenes, p.DE=p.DE, lBaselineExpr=lBaselineExpr, lOD=lOD,
         lfc=lfc, sim.seed=sim.seed, design="2grp")
}


## function to set fold change for DE genes
## I'll let the lfc for null genes being 0 at this time,
## and only generate the lfc for alternative genes.
setFC <- function(input, nDEgenes) {

    if(is.vector(input)) { ## vector
        if(length(input)==1) { ## constant
            lfc = rep(input, nDEgenes)
        } else if (length(input)!=nDEgenes) { ## vector
            ## warning("The length of lfc doesn't equal to number of DE genes. I'll sample from it. \n")
            lfc = sample(input, nDEgenes, replace=TRUE)
        } else {
            lfc = input
        }
    } else if (is.function(input)) { # a function
        lfc=input(nDEgenes)
    }   else {
        stop("Unrecognized form of lfc!\n")
    }

    ##   lfc0 <- lfc.null(ngenes)
    ##   lfc0[DEid] <- lfc
    lfc
}


