# Franzosa et al. 2019 analysis


# load packages
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

# function to pull metabolite submatrix out of metabdata given a category in metabmeta
  get_metab_dm <- function(cat){
    metabs <- metabmeta$Metabolomic.Feature[which(metabmeta$Putative.Chemical.Class == cat)]
    dist(metabdata[, colnames(metabdata) %in% metabs])
    # change vegdist to dist to do euclidean distance instead
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
  cats2use <- sort( nfeatures[nsamps == nrow(otutable)], decreasing=TRUE )
  # discard cats with fewer than 5 member features
  cats2use <- cats2use[cats2use > 10]
  
# run specificity for those cats, creating giant specs list; or not
  specs_list <- lapply(X=names(cats2use), FUN=function(catname){
    phy_or_env_spec(otutable, get_metab_dm(catname), n_sim=500, n_cores=22, 
      denom_type="sim_center", chunksize=44)
  })
  names(specs_list) <- names(cats2use)
  save(list=c("cats2use", "specs_list"), file="specs_list_cats2use_euc.rdata")

# table of number of significant species
  nsig <- sapply(X=specs_list, FUN=function(df){sum(df$Pval <= 0.05)})
  nsig <- sort(nsig)

# organize nfeatures and nsig, compare
  nfeatures2 <- sapply(X=names(nsig), FUN=function(nm){nfeatures[names(nfeatures) == nm]})
  pdf("diagnostic_nmetabolites_vs_nsig.pdf")
  plot(nsig ~ log(nfeatures2), xlab="log Number of metabolomic features in matrix",
    ylab="number of significant bacterial species")
  dev.off()


# draw plots
  pdf("barplot_sorted_all.pdf")
  par(mar=c(12,4,1,1))
  par(mfrow=c(2,1))
  centers2 <- barplot(nsig, ylab="# Sig. Species", 
    las=2, cex.names=0.4)
  dev.off()

  pdf("barplot_sorted_nonzero.pdf")
  nsig_nonzero <- nsig[nsig>0]
  par(mar=c(12,4,1,1))
  par(mfrow=c(2,1))
  centers2 <- barplot(nsig_nonzero, ylab="# Sig. Species", 
    las=2, cex.names=1)
  dev.off()

  pdf("barplot_sorted_nonzero_angle.pdf")
  barplot_rotatex <- function(x, labels_vec, labels_cex, labels_angle, ...) {
      plt <- barplot(x, xaxt="n", ...)
      text(plt, par("usr")[3], labels = labels_vec, srt = labels_angle, adj = c(1,1), xpd = TRUE, cex=labels_cex) 
  }
  par(mar=c(12,4,1,1))
  barplot_rotatex(nsig_nonzero, names(nsig_nonzero), 0.9, 45, ylab="Number of significant species")
  dev.off()


# correlations
  # pre-calculate all beta div matrices
  library(parallel)
  dms <- mclapply(X=names(nsig), FUN=get_metab_dm, mc.cores=20)
  # make mantel correlation matrix
  nm <- length(names(nsig))
  metab_cormat <- matrix(0, nrow=nm, ncol=nm)
  rownames(metab_cormat) <- colnames(metab_cormat) <- names(nsig)
  for(i in 1:nm){for(j in 1:nm){
    if(metab_cormat[i,j] == 0){
      metab_cormat[i,j] <- metab_cormat[j,i] <- cor(
        x=dms[[i]],
        y=dms[[j]],
        use="complete"
      )
    }
  }}

  range( metab_cormat[lower.tri(metab_cormat)])

  # also do correlations for Spec
  nm <- length(names(nsig))
  spec_cormat <- matrix(0, nrow=nm, ncol=nm)
  for(i in 1:nm){for(j in 1:nm){
    if(i==j){
      spec_cormat[i,j] <- 1
    }else if(spec_cormat[i,j] == 0){
      iname <- names(nsig)[[i]]
      jname <- names(nsig)[[j]]
      spec_cormat[i,j] <- spec_cormat[j,i] <- cor(
        x=specs_list[[iname]]$Spec,
        y=specs_list[[jname]]$Spec,
        use="complete"
      )
    }
  }}

  range( spec_cormat[lower.tri(spec_cormat)])


  pdf("diagnostic_metabolite_corr_vs_spec_corr.pdf")
  plot(
    x=as.dist(metab_cormat),
    y=as.dist(spec_cormat),
    xlab="Metabolomic Mantel correlation",
    ylab="Spec correlation"
  )
  dev.off()

# which metabolite classes were strongly correlated?
  # transform metab_cormat into edge list
  # m - a square distance matrix (with row+colnames) or a dist
  dm2edgelist <- function(m){
    if("dist" %in% class(m)){ m <- as.matrix(m) }
    colidx <- sapply(X=1:nrow(m), FUN=function(x){rep(x, ncol(m))})
    rowidx <- t(colidx)[lower.tri(t(colidx))]
    colidx <- colidx[lower.tri(colidx)]
    data.frame(
      value=mapply(FUN=function(i,j){m[i,j]}, i=rowidx, j=colidx),
      namei=sapply(X=rowidx, FUN=function(i){rownames(m)[i]}),
      namej=sapply(X=colidx, FUN=function(j){colnames(m)[j]}),
      stringsAsFactors=FALSE
    )
  }

  metab_cor_df <- dm2edgelist(metab_cormat)
  # sort it
  metab_cor_df <- metab_cor_df[order(metab_cor_df$value, decreasing=TRUE), ] 
  write.table(file="metab_cor_df.tsv", x=metab_cor_df, sep="\t", quote=F, row.names=F)

  # only in figure:
  metab_cor_df_fig <- metab_cor_df[(metab_cor_df$namei %in% names(nsig_nonzero)) & (metab_cor_df$namej %in% names(nsig_nonzero)), ]
  write.table(file="metab_cor_df_fig.tsv", x=metab_cor_df_fig, sep="\t", quote=F, row.names=F)


# is anything strongly correlated with bile acids?
  biles <- "Bile acids, alcohols and derivatives"
  biles_cor_df <- metab_cor_df[metab_cor_df$namei==biles | metab_cor_df$namej==biles, ]
  max(biles_cor_df$value)
  biles_cor_df[ which.max(biles_cor_df$value),]

# what bacteria are significantly specific to biles
  specs_biles <- specs_list[[biles]]
  specs_biles <- specs_biles[specs_biles$Pval <= 0.05,]
  specs_biles[order(specs_biles$Spec),]

# aggregate all results into matrix and look for cdiff specificity
  sigmat <- sapply(X=specs_list, FUN=function(x){x$Pval <= 0.05})
  rownames(sigmat) <- rownames(specs_list[[1]])
  which( sigmat[grep("difficile", rownames(sigmat)), ] )
