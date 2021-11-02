#!/usr/bin/env Rscript


# emp analysis

# read in data
library(data.table)
library(specificity)
library(ape)

message("formatting metadata")
metadata <- read.delim("../../data_sets/emp/clean_map_global16S.tsv", 
	sep="\t", header=T,	stringsAsFactors=FALSE)

# get 5000 samps, with sample intensity inversely proportional to empo3
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

weights <- sapply(X=metadata$empo_3, FUN=function(x){1/sum(metadata$empo_3==x)})
set.seed(12345)
metadata_less <- metadata[sample(1:nrow(metadata), size=5000, prob=weights), ]

# message("making smaller otu gzip")
# read in OTU table as a lit of columns via CUT, first reading in sample names only
# otutable_fp <- "../../data_sets/emp/global16S_subset_fixed.tsv.gz"
# frcmd <- paste0("gunzip -c ", otutable_fp, " | head -n 10")
# cols <- colnames(fread(cmd=frcmd, sep="\t"))
# cols2get <- c(1, which(cols %in% metadata_less$X.SampleID))
# frcmd2 <- paste0("gunzip -c ", otutable_fp, " | cut -f", paste(cols2get, collapse=","),
# 	" | sed '1,2d' | gzip >2read.tsv.gz")
# system("rm 2read.tsv.gz")
# system(frcmd2)

message("reading in smaller otu gzip")
otutable <- fread("2read.tsv.gz", sep="\t")
otuids <- as.vector(otutable[[1]])
otutable <- as.matrix(otutable[,2:ncol(otutable)])
rownames(otutable) <- otuids

# remove zero count OTUs and otus with less than 3 observations, since they
# make table big and won't throw off prop abund TOO MUCH
message("otutable dims before otus >=3 and samps >= 5000 seqs operation:")
dim(otutable)
otutable <- otutable[rowSums(otutable) >= 3,]
otutable <- otutable[,colSums(otutable) >= 5000]
message("otutable dims after:")
dim(otutable)

# prop abund
message("prop abunding otutable")
otutable <- t(otutable)
otutable <- prop_abund(otutable)

# drop low occupancy otus
otutable <- occ_threshold(otutable, threshold=30)

# organize metadata
message("subsetting and organizing metadata")
metadata_less <- metadata[metadata$X.SampleID %in% rownames(otutable),]
metadata_less <- metadata_less[order(metadata_less$X.SampleID), ]
otutable <- otutable[order(rownames(otutable)), ]
message("following must be true (check sort):")
message(all(rownames(otutable) == metadata_less$X.SampleID))

message("making plots")
# chart of empo3 categories
pdf("sampling_barplot.pdf")
par(mar=c(7,12,4,2))
barplot(table(metadata_less$empo_3), horiz=T, las=2, xlab="number of samples")
dev.off()
ontophy <- read.tree(text=onto2nwk( data.frame(
		metadata_less$empo_1, 
		metadata_less$empo_2, 
		metadata_less$empo_3)
))
pdf("ontology_as_phylogeny.pdf")
plot(ontophy, show.node.label=T)
dev.off()


# otutable2 <- otutable[,1:200]


message("running specificity")
specs <- phy_or_env_spec( otutable, hosts_phylo=ontophy, 
	hosts=metadata_less$empo_3, 
	n_sim=100, p_method="gamma_fit", n_cores=10)
write.table(specs, file="emp_specs_results.txt", sep='\t', quote=FALSE)
save(list="specs", file="specs.rdata")
 
