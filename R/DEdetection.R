#############################################
## funcions to operate on simulated data,
## run DE detection and return results.
#############################################


#############################################
## Use DSS for DE detection.
#############################################
run.DSS <- function(dat)  {
    require(DSS)

    ## run DSS
    seqData=newSeqCountSet(dat$counts, dat$designs)
    seqData=estNormFactors(seqData)
    seqData=estDispersion(seqData)
    result=waldTest(seqData, 0, 1)

    ## construct results, ordered by the orignal order of genes.
    ix=result$geneIndex
    result=data.frame(geneIndex=ix, result[,c("pval", "fdr")])
    res=result[order(result$geneIndex),]
    res
}

######################################
## get pvalues and FDR from edgeR
######################################
run.edgeR <- function(dat) {
    require(edgeR)

    ## run edgeR
    d <- DGEList(counts=dat$counts, group=dat$designs)
    d <- calcNormFactors(d)
    d <- estimateCommonDisp(d)
    d <- estimateTagwiseDisp(d)
    fit.edgeR <- exactTest(d)
    pval.edgeR <- fit.edgeR$table$PValue
    a <- topTags(fit.edgeR, n=nrow(dat$counts))

    ## construct results, ordered by the orignal order of genes
    ix <- sort(fit.edgeR$table$PValue, index.return=TRUE)$ix
    result <- data.frame(geneIndex=ix, pval=a$table$PValue, fdr=a$table$FDR)
    res <- result[order(result$geneIndex),]
    res
}

## ######################################
## ## get pvalues from DEseq
## ######################################
## run.DESeq <- function(dat) {
##     require(DESeq)

##     ## run DESeq
##     cds <- newCountDataSet(dat$counts, dat$designs)
##     cds <- estimateSizeFactors( cds )
##     cds <- estimateDispersions(cds) ##, fitType="local")
##     fit.DEseq <- nbinomTest( cds, "0","1")
##     pval <- fit.DEseq$pval
##     pval[is.na(pval)] <- 1
##     fdr <- fit.DEseq$padj

##     ## construct results, ordered by the orignal order of genes
##     ix <- sort(pval, index.return=TRUE)$ix
##     result <- data.frame(geneIndex=ix, pval=pval[ix], fdr=fdr[ix])
##     res <- result[order(result$geneIndex),]
##     res
## }


######################################
## get pvalues from DEseq2
######################################
run.DESeq2 <- function(dat) {
    require(DESeq2)

    ## run DESeq2
    cond <- factor(dat$designs)
    dds <- DESeqDataSetFromMatrix(dat$counts, DataFrame(cond), ~ cond)
    dds <- DESeq(dds, quiet=TRUE)
    res <- results(dds)
    pval <- res$pvalue
    pval[is.na(pval)] <- 1
    fdr <- res$padj
    fdr[is.na(fdr)] <- 1

    ## construct results, ordered by the orignal order of genes
    ix <- sort(pval, index.return=TRUE)$ix
    result <- data.frame(geneIndex=ix, pval=pval[ix], fdr=fdr[ix])
    res <- result[order(result$geneIndex),]
    res
}


