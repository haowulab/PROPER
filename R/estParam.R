#### added by Jean 2015-04
#### edited from Hao's 
require(edgeR)

####################################################
#### rowVars: variances per row for a matrix
####################################################
rowVars=function(x) {
  n0 <- ncol(x)
  EX <- rowMeans(x, na.rm=TRUE)
  vv <- rowSums(sweep(x,1,EX)^2)/(n0-1)
  vv
}


####################################################
#### est.param:
#### input: eset is a  or a matrix
####################################################

## function to estimate and return parameters, including:
## sequencing depth, sizefactors, log mean expressions, log OD
estParam <- function(X,type=1) {
  if (!class(X)%in%c("ExpressionSet","matrix"))
    stop("input X should be either an ExpressionSet or a matrix")
  if(!is.matrix(X))  X=exprs(X)

  if (!(type==1|type==2)) stop("type must be either 1 or 2")
  if(type==1){# use simple estimation with minimum assumption
    res=getDisp1(X)
  }
  
  if(type==2){ #use edgeR
    require(edgeR)
    res=getDisp2(X)
  }

  res
}
  

getDisp1=function(X,seed=2015){
  seqDepth=colSums(X)
  k=seqDepth/median(seqDepth)
  X2=sweep(X, 2, k, FUN="/")##pseudo Counts
    
  m=rowMeans(X2)
  v=rowVars(X2)
  phi.g0 = phi.g = (v-m)/m^2
  ## only keep those with good coverage
  Good=m>30 & rowMeans(X>0)>0.8
  phi.g0=phi.g0[Good]
  phi.g0=replace(phi.g0,is.na(phi.g0),.001)
  phi.g0=replace(phi.g0,phi.g0<.001,0.001)
  
  ## for those with unobserved dispersion, sample from phi.g0
  ii=(phi.g<=0.001) | is.na(phi.g)
  set.seed(seed)
  phi.g[ii]=sample(phi.g0, sum(ii), replace=TRUE)

  list(seqDepth=seqDepth,lmeans=log(m),lOD=log(phi.g))
}
##in edgeR, "The effective library size is then
##the original library size multiplied by the scaling factor."
## example: data(bottomlyCounts)
##
#param1=estParam(bottomly.eset[sub1,1:10])
getDisp2=function(X){
  X=DGEList(counts=X,lib.size=colSums(X))
  y=estimateCommonDisp(X)
  y= suppressWarnings(estimateTrendedDisp(y))
  y=estimateTagwiseDisp(y)
  lmean=y$AveLogCPM*log(2)+log(y$pseudo.lib.size/1e6)
  phi.g=y$tagwise.dispersion

  list(seqDepth=y$sample$lib.size, lmean=lmean, lOD=log(phi.g))
  
}
