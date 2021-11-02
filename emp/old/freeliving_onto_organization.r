# emp analysis

# read in data
library(data.table)
library(specificity)
library(ape)
source("onto2nwk.r")
metadata <- read.delim("../../data_sets/emp/clean_map_global16S.tsv", 
	sep="\t", header=T,	stringsAsFactors=FALSE)


# get 3000 samps, with sample intensity inversely proportional to empo3
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
metadata3k <- metadata[sample(1:nrow(metadata), size=3000, prob=weights), ]

# chart of empo3 categories
pdf("sampling_barplot.pdf")
par(mar=c(7,12,4,2))
barplot(table(metadata3k$empo_3), horiz=T, las=2, xlab="number of samples")
dev.off()

ontophy <- read.tree(text=onto2nwk(data.frame(metadata$empo_1, metadata$empo_2, metadata$empo_3)))
pdf("ontology_as_phylogeny.pdf")
plot(ontophy, show.node.label=T)
dev.off()

# read in OTU table as a lit of columns via CUT, first reading in sample names only
# otutable_fp <- "../../data_sets/emp/global16S_subset_fixed.tsv.gz"
# frcmd <- paste0("gunzip -c ", otutable_fp, " | head -n 10")
# cols <- colnames(fread(cmd=frcmd, sep="\t"))
# cols2get <- c(1, which(cols %in% metadata3k$X.SampleID))
# frcmd2 <- paste0("gunzip -c ", otutable_fp, " | cut -f", paste(cols2get, collapse=","),
# 	" | sed '1,2d' | gzip >2read.tsv.gz")
# system(frcmd2)
otutable <- fread("2read.tsv.gz", sep="\t")

# remove low occupancy OTUs
occ <- rowSums(otutable >0)
otutable <- otutable[ occ>=30, ]

metadata3k <- metadata3k[metadata3k$X.SampleID %in% colnames(otutable), ]
metadata3k <- metadata3k[order(metadata3k$X.SampleID), ]

# otutable into matrix and sort
rn <- otutable[[1]]
otutable <- as.matrix(otutable[,2:ncol(otutable)])
rownames(otutable) <- rn
otutable <- otutable[,order(colnames(otutable))]

all(metadata3k$X.SampleID == colnames(otutable))


# get ready for specificity
otutable <- t(otutable)

save(file="freeliving_inputs.rdata", list=c("otutable", "metadata3k", "ontophy"))



