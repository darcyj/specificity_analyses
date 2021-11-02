# Franzosa et al. 2019 analysis


# define ncats, the number of categories to analyze
  # basically looking at the top ncats categories with the most features within them
  ncats <- 60

# load packages
  library(vegan)
  library(specificity)

# read in otutable
  otutable <- read.delim("../../../data_sets/franzosa_ibd/microbiome.csv", sep=',', row.names=1, stringsAsFactors=F)

# read in metabolomics data and metabolomics metadata
  metabdata <- read.delim("../../../data_sets/franzosa_ibd/metabolites.csv", sep=",", row.names=1, stringsAsFactors=F)
  metabdata <- t(metabdata[8:nrow(metabdata),])
  metabdata <- apply(X=metabdata, MAR=2, FUN=as.numeric)
  metabmeta <- read.delim("../../../data_sets/franzosa_ibd/metaboliteMetadata.csv", sep=',', stringsAsFactors=F)

# check that otutable and metabdata match each other:
  all(rownames(metabdata) == rownames(otutable))

# need vegan for bray-curtis distance
  library(vegan)

# function to pull metabolite submatrix out of metabdata given a category in metabmeta
  get_metab_dm <- function(cat){
    metabs <- metabmeta$Metabolomic.Feature[which(metabmeta$Putative.Chemical.Class == cat)]
    vegdist(metabdata[, colnames(metabdata) %in% metabs])
  }

# table of how many features are within each category
  nfeatures <- table(metabmeta$Putative.Chemical.Class)
  # remove NA category
  nfeatures <- nfeatures[!grepl("N/A", names(nfeatures))]

# table of how many samples are within each category
  # (i.e. samples that aren't missing that cat altogether)
  nsamps <- sapply(X=names(nfeatures), FUN=function(cat){
    metabs <- metabmeta$Metabolomic.Feature[which(metabmeta$Putative.Chemical.Class == cat)]
    tab <- metabdata[, colnames(metabdata) %in% metabs, drop=FALSE]
    sum(rowSums(tab) > 0)
  })

# get top ncats, discarding those with missing samples
  cats2use <- sort( nfeatures[nsamps == nrow(otutable)], decreasing=TRUE )[1:ncats]



# run specificity for those cats, creating giant specs list
  specs_list <- lapply(X=names(cats2use), FUN=function(catname){
    phy_or_env_spec(otutable, get_metab_dm(catname), n_sim=500, n_cores=22, chunksize=44)
  })
  names(specs_list) <- names(cats2use)

# save
  save(list=c("cats2use", "specs_list"), file="specs_list_cats2use.rdata")

# table of number of significant species
  nsig <- sapply(X=specs_list, FUN=function(df){sum(df$Pval <= 0.05)})
  nsig <- sort(nsig)

# kmeans to visualize co-correlation
  # pre-calculate all beta div matrices
  library(parallel)
  dms <- mclapply(X=names(nsig), FUN=get_metab_dm, mc.cores=20)
  # make mantel correlation matrix
  nm <- length(names(nsig))
  cormat <- matrix(0, nrow=nm, ncol=nm)
  for(i in 1:nm){for(j in 1:nm){
    if(cormat[i,j] == 0){
      cormat[i,j] <- cormat[j,i] <- cor(
        x=dms[[i]],
        y=dms[[j]],
        use="complete"
      )
    }
  }}
  colnames(cormat) <- rownames(cormat) <- names(nsig)
  pdf("correlation_histogram.pdf")
  hist(cormat[lower.tri(cormat)])
  dev.off()
  # kmeans for multiple ks
  dm2 <- 1 - abs(cormat)
  ks <- 2:(nrow(dm2) -1)
  k_rsquared <- sapply(X=ks, FUN=function(k){
    res <- kmeans(dm2, centers=k)
    return(res$betweenss / res$totss)
  })

  # distances to corner
  dist2corner <- mapply(
    FUN=function(x,y){
      yscaled <- y / (1-min(k_rsquared))
      xscaled <- x / length(ks)
      sqrt(((1-yscaled)^2) + xscaled^2)
    },
    x=ks,
    y=k_rsquared
  )
  bestk <- ks[which.min(dist2corner)]

  pdf("kmeans_rsquareds.pdf")
  plot(k_rsquared ~ ks, xlab="Number of clusters (k)", ylab="R-squared", pch=20, cex=3)
  points(x=bestk, y=k_rsquared[ks==bestk], pch=20, cex=3, col="red")
  text(x=bestk, y=k_rsquared[ks==bestk] + 0.05, labels=paste0("k=", bestk), 
    col="red", adj=c(0.5,0.5))
  dev.off()

  kclust <- kmeans(dm2, centers=bestk)$cluster
  capture.output(kclust[order(kclust, decreasing = T)], file="k_clusters.txt")

# generate colors for kclust
  set.seed(5432)
  cols2use <- sample(rainbow(bestk, start=0.3))
  kcols <- rep("black", length(kclust))
  for(i in 1:bestk){
    kcols[kclust==i] <- cols2use[i]
  }

# barplot sorted low to hi
  par(mar=c(12,4,1,1))
  centers <- barplot(nsig, ylab="Count of significantly specific features", 
    las=2, cex.names=0.515, col=kcols)

# barplot sorted by cluster/clustermax, low2hi within cluster
  dec <- TRUE
  maxcount <- sapply(X=1:bestk, FUN=function(k){max(nsig[kclust==k])})
  kclust_order <- order(order(maxcount, decreasing=dec))
  nsig2 <- c()
  for(i in 1:bestk){
    k <- which(kclust_order == i)
    namesk <- names(kclust[kclust==k])
    countsk <- sapply(X=namesk, FUN=function(nm){nsig[names(nsig)==nm]})
    names(countsk) <- namesk
    nsig2 <- c(nsig2, countsk <- sort(countsk, decreasing=dec))
  }
  kcols2 <- sapply(X=names(nsig2), FUN=function(nm){cols2use[kclust[names(kclust)==nm]]})

  # make spacing
  spacesize <- 0.4
  binsizes <- table(factor(kcols2, levels=unique(kcols2)))
  spacing <- c()
  for(b in binsizes){ spacing <- c(spacing, spacesize, rep(0, b-1)) }


  # draw plot
  pdf("barplot_sorted_kmeans.pdf")
  par(mar=c(12,4,1,1))
  par(mfrow=c(2,1))
  centers2 <- barplot(nsig2, ylab="Significant features", 
    las=2, cex.names=0.515, col=kcols2, space=spacing)
  depth <- -4
  wid <- 2
  nextbar <- FALSE
  par(lend=1) # forflat but lines
  for(bar in 1:length(centers2)){
    xpos <- centers2[bar]
    barcol <- kcols2[bar]
    if(bar==1){
      segments(x0=xpos, y0=0, x1=xpos, y1=depth, col=barcol, xpd = TRUE, lwd=wid)
      lastxpos <- xpos
    }else if(nextbar){
      segments(x0=xpos, y0=0, x1=xpos, y1=depth, col=barcol, xpd = TRUE, lwd=wid)
      lastxpos <- xpos
      nextbar <- FALSE
    }else if((bar==length(centers2)) || (kcols2[bar] != kcols2[bar+1])){
      segments(x0=xpos, y0=0, x1=xpos, y1=depth, col=barcol, xpd = TRUE, lwd=wid)
      segments(x0=lastxpos, y0=depth, x1=xpos, y1=depth, col=barcol, xpd = TRUE, lwd=wid)
      nextbar <- TRUE
    }
  }
  dev.off()

