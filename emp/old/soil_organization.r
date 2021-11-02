# emp analysis

# read in data
library(data.table)
# otutable <- fread("../../data_sets/emp/global16S_subset_fixed.tsv.gz", sep="\t")
metadata <- read.delim("../../data_sets/emp/clean_map_global16S.tsv", sep="\t", header=T)

# format data and reduce
metadata <- metadata[metadata$sample_type == "soil",]

# rn <- otutable[[1]]
# cols2keep <- colnames(otutable)[colnames(otutable) %in% metadata$X.SampleID]
# otutable <- otutable[, cols2keep, with=FALSE]
# rownames(otutable) <- rn
# fwrite(x=otutable, file="otutable_soil.tsv", sep="\t", row.names=T)
otutable <- fread("otutable_soil.tsv.gz", sep="\t", header=T)
rn <- otutable[[1]]
otutable <- otutable[, 1:=NULL]
otutable <- as.matrix(otutable)
rownames(otutable) <- rn
rm(rn)

occ <- rowSums(otutable >0)
otutable <- otutable[occ>=30, ]

metadata <- metadata[metadata$X.SampleID %in% colnames(otutable), ]
metadata <- metadata[order(metadata$X.SampleID), ]
otutable <- otutable[,order(colnames(otutable))]

all(metadata$X.SampleID == colnames(otutable))


# get ready for specificity
otutable <- t(otutable)

save(file="soil_inputs.rdata", list=c("otutable", "metadata"))


