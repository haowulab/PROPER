######### some utility functions for simulation

######################################################################
## function to set baseline expression
######################################################################
setBaselineExpr <- function(input, ngenes) {

    param = NULL
    if(is.numeric(input)) { ## numeric
        if(length(input)==1 & is.numeric(input)) { ## constant
            lmeanExpr=rep(input, ngenes)
        } else if (length(input)!=ngenes) { ## vector
            stop("The length of lmeanExpr doesn't equal to ngenes!\n")
        } else
        lmeanExpr = input
    } else if (is.function(input)) { # a function
        lmeanExpr=input(ngenes)
    } else if (is.character(input)) { ## a character. Sample from real data.
        ## Note I save the mean expr in log-scale so have to exp them.
        datatsets <- c("cheung", "gilad","maqc","bottomly")
        if(input %in% datatsets ) {
            eval(parse(text=paste0("data(",input,", envir=environment())")))
        } else {
            stop("Unrecognized string of lmeanExpr. It must be one of 'cheung', 'gilad', 'bottomly', or 'maqc'!\n")
        }
        lmeanExpr=sample(param$lmean, ngenes, replace=TRUE)
    }
    else {
        stop("Unrecognized form of lmeanExpr!\n")
    }
    lmeanExpr
}

######################################################################
## function to set baseline expression, given sequencing depth
######################################################################
setBaselineExpr.seqDepth <- function(seqDepth, ngenes) {
    GE.human = NULL
    data(GE.human, envir=environment())
    GE.sample = sample(GE.human, ngenes, replace=TRUE)
    ntotal = sum(GE.sample)
    p0 = seqDepth / ntotal
    lmeanExpr=log(GE.sample*p0)
    lmeanExpr
}

######################################################################
## function to set over dispersion parameter
######################################################################
setOD <- function(input, ngenes) {
    param <- NULL
    if(is.numeric(input)) {
        if(length(input)==1) { ## constant
            lOD=rep(input, ngenes)
        } else if (length(input)!=ngenes) { ## vector
            stop("The length of OD doesn't equal to ngenes!\n")
        } else
        lOD = input
    } else if (is.function(input)) { # a function
        lOD=input(ngenes)
    } else if (is.character(input)) { ## a character
        datatsets <- c("cheung", "gilad","maqc","bottomly")
        if(input %in% datatsets ) {
            eval(parse(text=paste0("data(",input,", envir=environment())")))
        } else {
            stop("Unrecognized string of lmeanExpr. It must be one of 'cheung', 'gilad', 'bottomly', or 'maqc'!\n")
        }
        lOD=sample(param$lOD, ngenes, replace=TRUE)
    }
    else {
        stop("Unrecognized form of lOD!\n")
    }
    lOD
}


##############################################################
## default function to generate null and
## alternative log fold change
##############################################################
lfc.null <- function(n0) {
    rnorm(n0, mean=0, sd=0.0)
}

lfc.alt <- function(nDE) {
    nDE1 = round(nDE/2); nDE2 = nDE - nDE1
    ##lfc = c(rnorm(nDE1, mean=1, sd=0.2), rnorm(nDE2, mean=-1, sd=0.2))
    ##lfc = c(runif(nDE1, 0.5, 3), runif(nDE2, -3, -0.5))
    lfc = rnorm(nDE, 0, 1.5)
    lfc
}


######################################################################
## generate negative binomial rv, given mu and phi (over-dispersion)
######################################################################
rnegbinom <- function (n, mu=1, phi=0.01){
    rpois(n, rgamma(n, shape=1/phi,scale=mu*phi))
}


