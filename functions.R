library(Biobase); ## for featureData()
library(plyr);

getGeneArray <- function(rdsF, subjects) {
  ##----------------------------------------------------
  ## Function reads in the RDS file given. The microarray data is extracted
  ## and collapsed to gene level. Only the d0 columns for the requested 
  ## subjects are returned. 
  ##
  ## studyFileName = 'SDY123.rds' for example
  ## subjects = data.frame of subjectIDs you want d0 data for
  ##----------------------------------------------------
  
  
  if(!file.exists(rdsF)) stop(paste(rdsF,"not found in\n"));
  
  cat(paste("\nReading",rdsF,"\n"));
  Y3 <- readRDS(rdsF);
  
  ## Get geneSymbol to probeID mappings
  gs <- data.frame(
    probeIDs=as.character(featureData(Y3)$probeID),
    gs=as.character(featureData(Y3)$geneSymbol),
    stringsAsFactors=FALSE);
  
  ## Renaud edited the probeIDs to *NOT* include '_PM'. 
  ## We need to correct this in the 'gs' object in order to map the probes.
  if(length(grep("SDY269", rdsF)) > 0) gs$probeIDs=gsub("_PM", "", gs$probeIDs);
  
  ## Get probe-level normalized intensities.
  expr <- exprs(Y3);
  
  ## Get gene-level intensities
  m <- collapseArray(expr, gs);
  
  
  ## now select the columns for the subjects you want
  df <- data.frame(subjects, 
                   labs = paste(subjects$subjectID,"_d0",sep=""),
                   idx = 0,
                   newLabs = NA,
                   stringsAsFactors=FALSE);
  df$idx <- match(df$labs, colnames(m));
  df <- df[ !is.na(df$idx),]
  
  
  df$newLabs[ df$isTrueNR == 1 ] <- gsub("SUB..(.+)", "tNR_\\1", df$labs[ df$isTrueNR == 1 ])
  df$newLabs[ df$isTrueNR == 0 ] <- gsub("SUB..(.+)", "R_\\1", df$labs[ df$isTrueNR == 0 ])
  df <- df[ order(df$isTrueNR), ]
  
  if(sum(df$isTrueNR) == 0) {
    ## You don't have a single true non-responder in this data set, 
    ## so skip returning is expression matrix.
    ## You'll always get true responders but you need true non-responders
    ## for comparisons.
    m = NULL;
    return(m);
  }
  
  
  m <- m[, df$idx]; ## keep the subjects you are interested in
  colnames(m) <- df$newLabs; ## clearer labels for each subject
  
  return(m);
}



collapseArray <- function(mat, gs.df) {
  ##--------------------------------------------------------------------------
  ## Function collapses matrix 'mat' probe intensities to a single gene.
  ## When mulitple probes are found for same gene, the probe with the largest 
  ## average intensity is selected.
  ## 1 gene --> 1 probe 
  ##
  ## mat = matrix of expression values, rownames = probeIDs
  ## gs.df = geneSymbol-to-probeID data.frame
  ## geneCol = column index of gs.df with gene symbols
  ## probeCol = column index of gs.df with probeIDs
  ##--------------------------------------------------------------------------
  
#   mat = expr; gs.df = gs; ## debug
  
  ## get average intensity for each probe
  avgProbe = apply(mat, 1, mean);
  
  ## remove any probes that are 'NA' after taking their log
  avg.df <- data.frame(gs.df, avgProbe, stringsAsFactors=FALSE);
  i <- is.na(avgProbe);
  avg.df <- avg.df[!i,]
  
  ## remove any genes with nchar(gene_id) = 0 (ie: blanks)
  avg.df$nchar = sapply(avg.df$gs, nchar);
  avg.df <- avg.df[ avg.df$nchar > 1, ];
  i <- is.na(avg.df$gs);
  avg.df <- avg.df[!i,];
  rm(i);
  
  
  ## pick the probe for each gene with the *maximum* average intensity
  ## got this command from:
  ## http://stackoverflow.com/questions/9718711/selecting-rows-which-contain-daily-max-value-in-r
  x <- ddply(avg.df, ~gs,function(x){x[which.max(x$avgProbe),]})
  
  ## the probes listed in 'x' are the most intense probes for the genes.
  ret <- mat[x$probeIDs, ];
  rownames(ret) <- x$gs;
  
  return(ret);
}



runLimma <- function(curSDY, curGL) {
  ##--------------------------------------------------------------------------
  ## Function runs limma on the matrix provided for the selected gene list
  ##--------------------------------------------------------------------------
  require(limma);
  
  n.tNR <- length(grep("tNR_", colnames(curSDY)));
  n.R   <- length(grep("^R_", colnames(curSDY)));
  
  ## design matrix
  design <- model.matrix(~0 + as.factor(c(rep(0,n.R), rep(1,n.tNR))));
  colnames(design) <- c("R", "tNR");
  contrast <- makeContrasts(tNR - R, levels=design);
  
  ## get selected genes
  i <- match(curGL, rownames(curSDY));
  i <- i[!is.na(i)];
  m <- curSDY[ i, ]; 
  
  if(nrow(m) == 0) {
    ## none of the genes in the list were found in this data set, return NULL
    return(NULL);
  }
  
  ## calculate differential expression
  fit <- lmFit(m, design);
  fit2 <- contrasts.fit(fit, contrast);
  fit2 <- eBayes(fit2);
  res <- topTable(fit2, number=nrow(m), adjust.method="BH");
  res <- res[ order(res$adj.P.Val), ];
  
  return(res);
}


parseLimmaRes <- function(res, f, sdyID) {
  ##--------------------------------------------------------------------------
  ##  Function summarizes the results produced by limma for a given study
  ##--------------------------------------------------------------------------
  
  ## output data.frame
  ret <- data.frame(
    sdyID=sdyID,
    libCutOff=f,
    numProts=NA,
    nPval_0.01=0,
    nPval_0.05=0,
    nPval_0.10=0,
    nPval_0.15=0,
    nPval_0.20=0,
    stringsAsFactors=FALSE);
  
  ret$numProts=nrow(res);
  ret$nPval_0.01=sum(res$adj.P.Val <= 0.01);
  ret$nPval_0.05=sum(res$adj.P.Val <= 0.05);
  ret$nPval_0.10=sum(res$adj.P.Val <= 0.10);
  ret$nPval_0.15=sum(res$adj.P.Val <= 0.15);
  ret$nPval_0.20=sum(res$adj.P.Val <= 0.20);
  
  return(ret);
}


















