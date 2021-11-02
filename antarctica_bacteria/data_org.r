# antarctica data organization
	library(specificity)
	library(ape)
	source("funs.r")

# read in otutables
	bac_otutab <- read.table("../../data_sets/antarctica/16S/final_filtering/otutable.txt",
		row.names=1, header=T, sep='\t', stringsAsFactors=FALSE)
	euk_otutab <- read.table("../../data_sets/antarctica/18S/final_filtering/otutable.txt",
		row.names=1, header=T, sep='\t', stringsAsFactors=FALSE)

# read in tax
	euk_tax <- read.delim("../../data_sets/antarctica/18S/final_filtering/tax.txt", 
		header=T, sep='\t', stringsAsFactors=FALSE)
	bac_tax <- read.delim("../../data_sets/antarctica/16S/gg_tax/repset_tax_assignments.txt",
		header=F, sep='\t', stringsAsFactors=FALSE)
	colnames(bac_tax) <- c("ASVID", "tax", "e-val", "gg_id")

# read in gg bac tree
	gg_tree <- ape::read.tree("../../data_sets/antarctica/16S/gg_tax/99_otus_unannotated.tree")

# get rid of blanks (they don't agree, not contamination)
	bac_otutab <- bac_otutab[!grepl(x=rownames(bac_otutab), pattern="BLANK"),]
	euk_otutab <- euk_otutab[!grepl(x=rownames(euk_otutab), pattern="BLANK"),]

# check sorting
	all(rownames(bac_otutab) == rownames(euk_otutab))
	# TRUE - good work, dada2

# make euk subtables
	subtable_tax <- function(tax2get, otutable=euk_otutab, tax=euk_tax){
		otus2keep <- NULL
		for(taxon in tax2get){
			otus2keep <- c(otus2keep, tax[,1][grepl(pattern=taxon, x=tax$tax)])
		}
		out <- otutable[, colnames(otutable) %in% otus2keep]
		return(out)
	}

	euk_otutab <- prop_abund(euk_otutab)

	euk_otutab_fungi <- subtable_tax("Fungi")
	euk_otutab_algae <- subtable_tax("phyta")
	euk_otutab_meta <- subtable_tax("Metazoa")
	euk_otutab_cerc <- subtable_tax("Cercozoa")
	euk_otutab_cilio <- subtable_tax("Ciliophora")

# check for samples with 0s
	min(rowSums(euk_otutab_fungi > 0))
	min(rowSums(euk_otutab_algae > 0))
	min(rowSums(euk_otutab_meta > 0))
	min(rowSums(euk_otutab_cerc > 0))
	min(rowSums(euk_otutab_cilio > 0))

# make dists
	library(vegan)
	# aitch <- function(x){dist(ilr(x))}
	dist_fungi <- vegdist(prop_abund(euk_otutab_fungi))
	dist_algae <- vegdist(prop_abund(euk_otutab_algae))
	dist_meta <- vegdist(prop_abund(euk_otutab_meta))
	dist_cerc <- vegdist(prop_abund(euk_otutab_cerc))
	dist_cilio <- vegdist(prop_abund(euk_otutab_cilio))
	dist_alleuks <- vegdist(prop_abund(euk_otutab))

# read in metadata
	metadata <- read.table("../../data_sets/antarctica/16S/metadata.txt", 
		sep='\t', header=T, stringsAsFactors=FALSE)
	all(metadata[,1] == rownames(bac_otutab))
	# TRUE!

# fix metadata names, esp since specificity calculation is scale-invariant
	colnames(metadata)[colnames(metadata) == "N_percent"] <- "N"
	colnames(metadata)[colnames(metadata) == "C_percent"] <- "C"
	colnames(metadata)[colnames(metadata) == "P_mg.kg"] <- "P"

# calculate hole area
	metadata$area <- pi * metadata$ns_diam * metadata$ew_diam

# make geographic distance matrix
	dist_geo <- as.dist(distcalc(lat=metadata$LAT, lng=metadata$LNG))

# make otutable that only includes otus in10 or more samps
	bac_otutab <- prop_abund(bac_otutab)
	bac_otutab_ovr10 <- bac_otutab[, colSums(bac_otutab > 0)>=10]

# figure out tree garbage
	# for each tip name in tree, figure out if any of my asvs were assigned to it. if so,
	# replace gg tip with polytomy including all asvids.
	# first, get rid of tips not in bac_tax$gg_id
	gg_tree <- keep.tip(phy=gg_tree, tip=bac_tax$gg_id[bac_tax$gg_id != "None"])
	# now do the rest
	for(tn in gg_tree$tip.label){
		# which asvs were assigned to it?
		asv_hits <- bac_tax$ASVID[bac_tax$gg_id == tn]
		gg_tree <- replace_tip_with_polytomy(gg_tree, tn, asv_hits)
	}

# save organized data
	save(list=c("bac_otutab_ovr10", "metadata", "dist_geo", "dist_fungi",
		"dist_algae", "dist_meta", "dist_cerc", "dist_cilio", 
		"dist_alleuks", "gg_tree"), file="organized_data.rdata")
