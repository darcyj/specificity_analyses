#!/usr/bin/env Rscript


# emp analysis

# set up env
	library(data.table)
	library(specificity)
	library(ape)

# read in metadata and do some formatting
	message("formatting metadata")
	metadata <- read.delim("../../../data_sets/emp/clean_map_global16S.tsv", 
		sep="\t", header=T,	stringsAsFactors=FALSE)

	prep4phy <- function(x){
		x <- gsub(")", "", x)
		x <- gsub("\\(", "", x)
		x <- gsub(" ", "_", x)
		return(x)
	}
	metadata$empo_3 <- prep4phy(metadata$empo_3)
	metadata$empo_2 <- prep4phy(metadata$empo_2)
	metadata$empo_1 <- prep4phy(metadata$empo_1)

# raw otutable is WAY too large to transpose or to do rowsums
	# (species are rows). strategy is to read it in in chunks. 
	otutable_fp <- "../../../data_sets/emp/global16S_subset_fixed.tsv.gz"
	# first, gotta get the total number of lines.
	otutab_nlines <- 309469
	# otutab_nlines <- system(paste0("gunzip -c ", otutable_fp, " | wc -l"), intern=TRUE)
	# now in chunks of 10000 lines, read it in and calculate occupancy and col sums
	# need to start at 4th line, because first three are colnames and some
	# comment lines
	# make gz connection, then read in "trash" lines to skip them.
	otutab_conn <- gzfile(otutable_fp, open="r") # open is NOT OPTIONAL, will fail if not used
	header <- unlist(strsplit(readLines(otutab_conn, n=1), split="\t"))
	all_sampleids <- header[-1]
	# skip first 3 lines since they are junk
	trash <- readLines(otutab_conn, n=3)
	rm(trash)
	# read in chunk-by-chunk
	chunksize <- 10000
	col_sums <- list()
	n_chunks <- ceiling((otutab_nlines-3)/chunksize)
	pb <- txtProgressBar(min=0, max=n_chunks, style=3)
	for(i in 1:n_chunks){
		txt <- readLines(otutab_conn, chunksize)
		chunkdata <- fread( text=txt, sep="\t", drop=1 )
		col_sums[[i]] <- colSums(chunkdata)
		setTxtProgressBar(pb, i)
	}
	close(pb)
	close(otutab_conn)
	# total col sums is just sum of col_sums
	col_sums <- Reduce("+", col_sums)

# create column mask - boolean vector describing which columns to keep
	# discard columns that are not well sampled
	col_mask <- col_sums >= 5000
	occ_thresh <- 500

# Read in only desired data:
	# samples with depth >= 5000
	# otus occupancy > 500 (just in good samples, above)
	otutab_conn <- gzfile(otutable_fp, open="r") # open is NOT OPTIONAL, will fail if not used
	# skip first 3 lines since they are junk
	trash <- readLines(otutab_conn, n=3)
	rm(trash)
	chunksize <- 10000
	table_pieces <- list()
	n_chunks <- ceiling((otutab_nlines-3)/chunksize)
	pb <- txtProgressBar(min=0, max=n_chunks, style=3)
	for(i in 1:n_chunks){
		# raw text read in
		txt <- readLines(otutab_conn, chunksize)
		# transform to matrix, and get otuids
		chunkdata <- as.matrix(fread( text=txt, sep="\t", drop=1 ))
		otunames <- as.vector(gsub(pattern="\t.*", replacement="", x=txt))
		# occupancy of each OTU, but only within samples of approprate depth
		occs <- rowSums( chunkdata[, col_mask] > 0 )
		# store approprate subset of chunkdata and name the OTUs
		table_pieces[[i]] <- chunkdata[occs>=occ_thresh, col_mask, drop=FALSE]
		rownames(table_pieces[[i]]) <- otunames[occs>=occ_thresh]
		setTxtProgressBar(pb, i)
	}
	close(pb)
	close(otutab_conn)

	# prop abund each chunk using col_sums
	# matrix m divided by vector v m/v proceeds columnwise
	# so all values in m[1,] are divided by v[1]
	# since each matrix in table_pieces has samples as columns,
	# they must be transposed first and THEN divided.
	col_sums_masked <- col_sums[col_mask]
	for(i in 1:length(table_pieces)){
		table_pieces[[i]] <- t(table_pieces[[i]]) / col_sums_masked
	}

	# combine into matrix
	otutable <- do.call("cbind", table_pieces)

# add sampleids back to otutable
	rownames(otutable) <- all_sampleids[col_mask]

# save otutable
	saveRDS(otutable, file="otutable.rds")

# format metadata to match up with otutable
	metadata <- metadata[metadata$X.SampleID %in% rownames(otutable),]
	metadata <- metadata[ sapply(X=1:nrow(otutable), FUN=function(i){
		which(metadata$X.SampleID == rownames(otutable)[i])}), ]
	# check
	all( rownames(otutable) == metadata$X.SampleID )


# fix EMPO3 stuff because EMP people did a BAD JOB
	sort( unique(metadata$empo_3) )
	metadata$empo_3[metadata$empo_3=="Watersaline"] <- "Water_saline"
	metadata$empo_3[metadata$empo_3=="Sedimentsaline"] <- "Sediment_saline"

	table(metadata$empo_3)
	table(metadata$empo_2)
	table(metadata$empo_1)

	a <- ape::read.tree(text=onto2nwk(data.frame(
		empo1=metadata$empo_1, 
		empo2=metadata$empo_2, 
		empo3=metadata$empo_3)
	))
	plot(a)


# save metadata
	saveRDS(metadata, file="metadata.rds")

