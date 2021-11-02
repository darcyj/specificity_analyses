#!/usr/bin/env Rscript


# emp analysis

# read in data
	library(data.table)
	library(specificity)
	library(ape)
	metadata <- read.delim("../../../data_sets/emp/clean_map_global16S.tsv", 
		sep="\t", header=T,	stringsAsFactors=FALSE)

# format metadata
	message("formatting metadata")
	metadata <- metadata[metadata$empo_1 == "Free-living", ]
	metadata <- metadata[!is.na(metadata$empo_3), ]
	metadata <- metadata[!is.na(metadata$empo_2), ]
	metadata <- metadata[!is.na(metadata$empo_1), ]
	metadata$empo_3[metadata$empo_3=="Sediment(saline)"] <- "Sediment (saline)"
	metadata$empo_3[metadata$empo_3=="Water(saline)"] <- "Water (saline)"

	prep4phy <- function(x){
		x <- gsub(")", "", x)
		x <- gsub("\\(", "", x)
		x <- gsub(" ", "_", x)
		return(x)
	}
	metadata$empo_3 <- prep4phy(metadata$empo_3)
	metadata$empo_2 <- prep4phy(metadata$empo_2)
	metadata$empo_1 <- prep4phy(metadata$empo_1)

# make a smaller otutable, in bash. commented out because run once.
	# message("subsetting otutable")
	# message("making smaller otu gzip")
	# read in OTU table as a lit of columns via CUT, first reading in sample names only
	# otutable_fp <- "../../../data_sets/emp/global16S_subset_fixed.tsv.gz"
	# frcmd <- paste0("gunzip -c ", otutable_fp, " | head -n 10")
	# cols <- colnames(fread(cmd=frcmd, sep="\t"))
	# cols2get <- c(1, which(cols %in% metadata$X.SampleID))
	# frcmd2 <- paste0("gunzip -c ", otutable_fp, " | cut -f", paste(cols2get, collapse=","),
	# 	" | sed '2,3d' | gzip >2read.tsv.gz")
	# system("rm 2read.tsv.gz")
	# system(frcmd2)

# read in smaller otutable per above
	message("reading subset otutable")
	otutable <- fread("2read.tsv.gz", sep="\t")
	otuids <- as.vector(otutable[[1]])
	otutable <- as.matrix(otutable[,2:ncol(otutable)])
	rownames(otutable) <- otuids

# remove samples with fewer than 5000 observations
	message("removing low count samples and calculating colsums")
	# can't do colSums(otutable), uses too much ram! Have to chunk it!
	colchunks <- split(x=1:ncol(otutable), f=rep(1:100, 
		each=ceiling(ncol(otutable)/100))[1:ncol(otutable)] )
	otutable_colsums <- list()
	pb <- txtProgressBar(min=0, max=length(colchunks), style=3)
	for(i in 1:length(colchunks)){
		otutable_colsums[[i]] <- sapply(X=colchunks[[i]], FUN=function(x){sum(otutable[,x])})
		garbage <- gc()
		setTxtProgressBar(pb, i)
	}
	close(pb)
	otutable_colsums <- unlist(otutable_colsums)
	otutable <- otutable[ , otutable_colsums >= 5000]
	otutable_colsums <- otutable_colsums[otutable_colsums >= 5000]
	rm(colchunks, garbage, pb)

# remove species with occupancy lower than 2% of 12700=250
	message("removing low occupancy OTUs")
	# can't do rowSums(otutable>0), uses too much ram! Have to chunk it!
	rowchunks <- split(x=1:nrow(otutable), f=rep(1:100, 
		each=ceiling(nrow(otutable)/100))[1:nrow(otutable)] )
	otutable_occs <- list()
	pb <- txtProgressBar(min=0, max=length(rowchunks), style=3)
	for(i in 1:length(rowchunks)){
		otutable_occs[[i]] <- sapply(X=rowchunks[[i]], FUN=function(x){sum(otutable[x,] > 0)})
		garbage <- gc()
		setTxtProgressBar(pb, i)
	}
	close(pb)
	otutable_occs <- unlist(otutable_occs)
	otutable <- otutable[otutable_occs>=250, ]
	rm(rowchunks, garbage, pb, otutable_occs)

# otutable is hopefully small enough to transpose now, and prop abund using otutable_colsums
	message("transposing and prop-abunding otutable")
	otutable <- t(otutable)
	# some code to show this works (remember, otutable_colsums are actually rows now)
	# testmat <- matrix(12, nrow=3, ncol=3)
	# testmat / 1:3
	otutable <- otutable / otutable_colsums

# organize metadata
	message("subsetting and organizing metadata")
	metadata <- metadata[metadata$X.SampleID %in% rownames(otutable),]
	metadata <- metadata[ sapply(X=1:nrow(otutable), FUN=function(i){
		which(metadata$X.SampleID == rownames(otutable)[i])}), ]
	message("following must be true (check sort):")
	message(all(rownames(otutable) == metadata$X.SampleID))

# make some plots
	message("making plots")
	# chart of empo3 categories
	pdf("sampling_barplot.pdf")
	par(mar=c(7,12,4,2))
	barplot(table(metadata$empo_3), horiz=T, las=2, xlab="number of samples")
	dev.off()
	ontophy <- read.tree(text=onto2nwk( data.frame(
			metadata$empo_1, 
			metadata$empo_2, 
			metadata$empo_3)
	))
	pdf("ontology_as_phylogeny.pdf")
	plot(ontophy, show.node.label=T)
	dev.off()

# saving data
	saveRDS(object = otutable, file="otutable.RDS")
	saveRDS(object = metadata, file="metadata.RDS")
	saveRDS(object = ontophy, file="ontophy.RDS")

