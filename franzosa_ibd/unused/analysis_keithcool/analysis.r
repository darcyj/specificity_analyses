# Franzosa et al. 2019 analysis
# load packages
  library(vegan)
  library(specificity)

# read in otutable
  otutable <- read.delim("../../data_sets/franzosa_ibd/microbiome.csv", sep=',', row.names=1, stringsAsFactors=F)

# read in metabolomics data and metabolomics metadata
  metabdata <- read.delim("../../data_sets/franzosa_ibd/metabolites.csv", sep=",", row.names=1, stringsAsFactors=F)
  metabdata <- t(metabdata[8:nrow(metabdata),])
  metabdata <- apply(X=metabdata, MAR=2, FUN=as.numeric)
  metabmeta <- read.delim("../../data_sets/franzosa_ibd/metaboliteMetadata.csv", sep=',', stringsAsFactors=F)

# check that otutable and metabdata match each other:
  all(rownames(metabdata) == rownames(otutable))

# what metabolite categories are worth looking at? Need to have lots of members
  mtab <- table(metabmeta$Putative.Chemical.Class)
  mtab <- mtab[order(mtab, decreasing=TRUE)]
  write.table(x=data.frame(category=names(mtab), n=as.numeric(mtab)), 
    file="metab_categories_table.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# function to pull metabolite submatrix out of metabdata given a category in metabmeta
  get_metab_mat <- function(cat){
    metabs <- metabmeta$Metabolomic.Feature[which(metabmeta$Putative.Chemical.Class == cat)]
    output <- metabdata[, colnames(metabdata) %in% metabs]
  }

# subset otutable to only include otus in 10+ samples
  otutable <- prop_abund(otutable)
  otutable <- occ_threshold(otutable, 10)

# do specificity analysis
  nc <- 20 # number of cpu cores
  specs_list <- list()
  specs_list$"Long-chain fatty acids" <- phy_or_env_spec(
    abunds_mat=otutable,
    env=vegdist(get_metab_mat("Long-chain fatty acids")),
    n_cores=nc
  )
  specs_list$"Sphingolipids" <- phy_or_env_spec(
    abunds_mat=otutable,
    env=vegdist(get_metab_mat("Sphingolipids")),
    n_cores=nc
  )
  specs_list$"Bile acids/alcohols/derivs" <- phy_or_env_spec(
    abunds_mat=otutable,
    env=vegdist(get_metab_mat("Bile acids, alcohols and derivatives")),
    n_cores=nc
  )
  specs_list$"Flavonoids" <- phy_or_env_spec(
    abunds_mat=otutable,
    env=vegdist(get_metab_mat("Flavonoids")),
    n_cores=nc
  )
  specs_list$"Cholesteryl esters" <- phy_or_env_spec(
    abunds_mat=otutable,
    env=vegdist(get_metab_mat("Cholesteryl esters")),
    n_cores=nc
  )

# colors
  speccols <- c("red", "blue", "forestgreen", "orange", "purple")

# make plot
  pdf("specs_plot.pdf", useDingbats=F)
  plot_specs_violin(specs_list, label_cex=0.7, cols=speccols)
  dev.off()

  pdf("specs_plot_ranged.pdf", useDingbats=F)
  plot_specs_violin(specs_list, label_cex=0.7, cols=speccols, minval=-0.6, maxval=0.6)
  dev.off()


  pdf("pairwise_specs_plot.pdf", useDingbats=F)
  plot_pairwise_spec(specs_list)
  dev.off()

# save
  save(list=ls(), file="specs_save.rdata")

# record version info
  capture.output(file="specificity_version.txt", capture.output(packageVersion("specificity")))

# make interactive app
  library(specificity.shiny)
  

