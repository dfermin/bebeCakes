library(utils);
library(Biobase); ## for featureData() 
library(ggplot2);

getSelectGenesArray <- function(studyFileName, GL) {
  ##----------------------------------------------------
  ## Function reads in the RDS file given. The microarray data is extracted
  ## and collapsed to gene level. From this, a new matrix of the selected 
  ## gene IDs is returned.
  ##
  ## studyFileName = 'SDY123.rds' for example
  ## GL = vector of gene IDs you want to extract
  ##----------------------------------------------------

  #studyFileName = "SDY269.rds"; GL = purvesh; ## debug;
  
  srcF = paste("../rds.files/",studyFileName,sep="");
  if(!file.exists(srcF)) stop(paste(studyFileName,"not found in ../rds.files\n"));
  
  cat(paste("\nReading",srcF,"\n"));
  Y3 <- readRDS(srcF);
  
  ## Get geneSymbol to probeID mappings
  gs <- data.frame(
    probeIDs=as.character(featureData(Y3)$probeID),
    gs=as.character(featureData(Y3)$geneSymbol),
    stringsAsFactors=FALSE);
  
  ## Renaud edited the probeIDs to *NOT* include '_PM'. 
  ## We need to correct this in the 'gs' object in order to map the probes.
  if(studyFileName == "SDY269.rds") gs$probeIDs=gsub("_PM", "", gs$probeIDs);
  
  ## reduce gs to be just the genes you want from GL
  ridx <- vector();
  for(g in GL) {
    pattern=paste("^",g,"$",sep=""); ## regex pattern: "/^geneID$/"
    j <- grep(pattern, gs$gs);
    if(length(j) > 0) ridx <- c(ridx, j);
  }
  ridx <- ridx[!is.na(ridx)];
  gs <- gs[ridx,];
  
  ## Get probe-level normalized intensities.
  expr <- exprs(Y3);
  
  ## Get gene-level intensities
  m <- collapseArray(expr, gs, geneCol=2, probeCol=1);
  
  return(m);
}



collapseArray <- function(mat, gs.df, geneCol=1, probeCol=2) {
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
  
  #geneCol=2; probeCol=1; mat = expr; gs.df = gs; ## debug
  
  
  ## remove any genes with nchar(gene_id) = 0 (ie: blanks)
  x <- sapply(gs.df[,geneCol], nchar);
  i <- which(x > 0)
  gs.df <- gs.df[i,];
  rm(i);
  
  Np <- nrow(mat); ## number of probes
  genes <- unique(gs.df[,geneCol]);
  Ng <- length(genes); ## number of unique genes
  
  ret <- matrix(data=NA, nrow=Ng, ncol=ncol(mat));
  rownames(ret) <- genes;
  colnames(ret) <- colnames(mat);
  
  gsMap <- data.frame(gs=genes, probeID=NA,stringsAsFactors=FALSE);
  
  cat("\nCollapsing probes...\n")
  PB <- txtProgressBar(min=1,max=Ng, width=50); ## initialize progress bar
  
  for(i in seq(Ng)) { ## loop over genes
    
    g <- genes[i];
    j <- base::which(gs.df[, geneCol] == g);
    k <- gs.df[j, probeCol];
    
    if( length(k) == 0 ) { next; } ## no probes matched
    
    ## none of the probes matched
    if( all( is.na( match(k,rownames(mat)) ) ) == TRUE ) { next; }
    
    ## if you got here, then at least 1 probe matched
    k <- k[!is.na(k)]
    tmp <- mat[ k, ];
    
    if(length(j) > 1) { ## multi-probe gene  
      avg <- apply(tmp,1,mean); ## get average intensity for each probe
      pid <- names(rev(sort(avg))[1]); ## probe with largest average intensity
      ret[g,] <- tmp[pid,];
      gsMap[i,] <- c(g,pid);
    } 
    else if(length(j) == 1) {
      ret[g,] <- tmp; ## 1 probe for this gene
      gsMap[i,] <- c(g,k);
    }
    setTxtProgressBar(PB,i);
  }
  close(PB);
  
  if( (nrow(ret) > 0) & ((all(is.na(ret))) == FALSE) ) {
    out <- data.frame(gsMap,ret, stringsAsFactors=FALSE);
    return(out);
  } else { return(NULL); }
}



makeGGscatter <- function(df, studyName, groupId, YLIMS) {
  ##--------------------------------------------------------------------------
  ## Function generate a ggplot2 scatterplot graphics object from the 
  ## data.frame (df) passed to it.
  ##--------------------------------------------------------------------------
  
  genes <- unique(df$gene);
  
  for(g in genes) {
    outF = paste("g_",g,"-",studyName,"_",grp,".pdf", sep="");
    pdf(file=outF, width=6, height=6);
    
    tmp <- df[ df$gene == g, ];
    
    med.bar <- data.frame(
      class=c("R","trueNR"),
      hline=c(median(tmp$expr[tmp$class == "R"]),
              median(tmp$expr[tmp$class == "trueNR"])));
    
    x <-  ggplot(tmp, aes(x=class, y=expr, fill=class)) + 
          geom_jitter(aes(colour=class), size=8, alpha=0.7) +
          xlab("") + ylab("") + ggtitle(g) + ylim(YLIMS) + 
          geom_errorbar(data=med.bar, aes(y=hline, ymax=hline, ymin=hline), size=1, colour="black") +
          theme(
            legend.position="none",
            axis.text.x=element_text(size=24,colour="black"),
            axis.text.y=element_text(size=24,colour="black"),
            plot.title=element_text(size=24)
          );
    print(x);  
    dev.off();
  }
}




makeGGscatterByGene <- function(df, G, groupId, YLIMS) {
  ##--------------------------------------------------------------------------
  ## Function generate a ggplot2 scatterplot graphics object from the 
  ## data.frame (df) passed to it.
  ##--------------------------------------------------------------------------
  
  
  
  outF = paste(G,"_",groupId,".pdf", sep="");
  pdf(file=outF, width=6, height=6);
  
  tmp <- df[, -4];
  
  med.bar <- data.frame(
    class=c("R","trueNR"),
    hline=c(median(tmp$expr[tmp$class == "R"]),
            median(tmp$expr[tmp$class == "trueNR"])));
  
  x <-  ggplot(tmp, aes(x=class, y=expr, fill=class)) + 
    geom_jitter(aes(colour=class), size=8, alpha=0.7) +
    xlab("") + ylab("") + ggtitle(G) + ylim(YLIMS) + 
    geom_errorbar(data=med.bar, aes(y=hline, ymax=hline, ymin=hline), size=1, colour="black") +
    theme(
      legend.position="none",
      axis.text.x=element_text(size=24,colour="black"),
      axis.text.y=element_text(size=24,colour="black"),
      plot.title=element_text(size=24)
    );
  print(x);  
  dev.off();
  
}
