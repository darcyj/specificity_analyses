# set up
	library(specificity)
	library(ape)
	otutable <- readRDS("otutable.rds")
	metadata <- readRDS("metadata.rds")

# remove "unknown" samples from metadata
	metadata <- metadata[ ! metadata$empo_1 %in% c("Control", "unknown") , ]

# remove "NA" samples from metadata
	metadata <- metadata[!is.na(metadata$empo_3), ]
	metadata <- metadata[!is.na(metadata$empo_2), ]
	metadata <- metadata[!is.na(metadata$empo_1), ]


# downsample for more even sampling across tips
	barplot(table(metadata$empo_3), horiz=T, las=2, xlab="number of samples")
	abline(v=600, col="red")
	downlvl <- 600

	set.seed(12345)
	for(e3 in unique(metadata$empo_3)){
		if(sum(metadata$empo_3==e3) >= downlvl){
			e3tokeep <- sample(which(metadata$empo_3 == e3), size=downlvl)
			tokeep <- metadata$empo_3 != e3
			tokeep[e3tokeep] <- TRUE
			metadata <- metadata[tokeep, ]
			rm(tokeep, e3tokeep)
		}
	}
	barplot(table(metadata$empo_3), horiz=T, las=2, xlab="number of samples")

# subset otutable to match and check
	otutable <- otutable[rownames(otutable) %in% metadata$X.SampleID, ]
	all(rownames(otutable) == metadata$X.SampleID)

# make ontology/phylogeny figure
	metadata$empo_3_n <- metadata$empo_3
	for(e3 in unique(metadata$empo_3)){
		metadata$empo_3_n[metadata$empo_3==e3] <- paste0(e3, "_<", sum(metadata$empo_3==e3), ">")
	}

	ontophy <- ape::read.tree(text=onto2nwk( data.frame(
			metadata$empo_1, 
			metadata$empo_2, 
			metadata$empo_3_n)
	))
	pdf("ontology_as_phylogeny.pdf")
	plot(ontophy, show.node.label=T)
	dev.off()

# do specificity
	message("running specificity")
	ga_par <- get_ga_defaults()
	ga_par$maxiters <- 1000
	specs <- phy_or_env_spec( otutable, hosts_phylo=ontophy, hosts=metadata$empo_3_n, 
		n_sim=500, n_cores=20, chunksize=100)

# make plot
	pdf("spec_plot.pdf")
	plot_specs_violin(list("Earth Microbiome Project ontology"=specs), minval=-1, maxval=0.2)
	dev.off()

# save stuff for later
	write.table(specs, file="emp_specs_results.txt", sep='\t', quote=FALSE)
	save(list=c("metadata", "ontophy", "otutable", "specs"), file="workspace.rdata")

# get taxonomy
	# turn them into fasta file; seqIDs are the ROW of specs the sequence matches
	fasta <- unlist(lapply(X=1:nrow(specs), FUN=function(i){
		c(paste0(">", i), rownames(specs)[i])
	}))
	writeLines(text=fasta, con="specs_seqs.fasta")
	# use qiime to get taxonomy, from another terminal:
	# # conda activate qiime1
	# # parallel_assign_taxonomy_blast.py -i specs_seqs.fasta -o specs_tax -O 20
	tax <- read.table("specs_tax/specs_seqs_tax_assignments.txt", 
		stringsAsFactors=FALSE, sep='\t', header=FALSE)
	# tax[[1]] is which row of specs it corresponds to, so re-order that way
	tax <- tax[order(tax[[1]]),]
	# check sorting:
	all(tax[[1]] == 1:nrow(specs))

# save again
	save(list=c("metadata", "ontophy", "otutable", "specs", "tax"), file="workspace.rdata")

# function to plot an OTU
	plotOTU <- function(otuidx, logotu=F, text3=""){
		categories <- metadata$empo_3_n
		otuvec <- as.vector(otutable[,otuidx])
		otuvec_nz <- otuvec > 0
		if(logotu){
			xrange <- log10(range(otuvec[otuvec_nz]))
			otuvec <- log10(otuvec)
		}else{
			xrange <- range(otuvec)
		}
		otutax <- tax[tax[[1]]==otuidx, 2]
		ntip=length(ontophy$tip.label)
		par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(0,0,1.5,0))
		# plot tree (first panel)
		plot(ontophy, cex=1)
		# plot points (second panel)
		plot(0, type="n", xlim=xrange, ylim=c(1, ntip),
			xaxt="n", yaxt="n", bty="n")
		for(i in 1:ntip){
			cat_i <- ontophy$tip.label[i]
			vals_i <- otuvec[categories == cat_i & otuvec_nz]
			points(x=vals_i, y=rep(i, length(vals_i)))
		}
		mtext(paste(otuidx, "-", rownames(specs)[otuidx]), outer = TRUE, cex = 0.6)
		mtext(paste0("Spec=", specs$Spec[otuidx], "; P=", round(specs$Pval[otuidx], 3), "; " ,otutax), outer = TRUE, cex = 0.6, line=-0.6)
		mtext(text3, outer=TRUE, tex=0.6, line=-1.2)
	}

# what are some widely dispersed OTUs?
	spec_0_idxs <- ((1:nrow(specs))[order(abs(0-specs$Spec))])[1:10]
	pdf("near_zero_specs_species.pdf")
	trash <- sapply(X=spec_0_idxs, FUN=plotOTU)
	dev.off()


# what are some specific otus around Spec=-0.6
	spec_0.6_idxs <- ((1:nrow(specs))[order(abs(-0.6 - specs$Spec))])[1:10]
	pdf("near_-0.6_specs_species.pdf")
	trash <- sapply(X=spec_0.6_idxs, FUN=plotOTU)
	dev.off()


# what are some specific otus around Spec=-0.9
	spec_0.9_idxs <- ((1:nrow(specs))[order(abs(-0.9 - specs$Spec))])[1:10]
	pdf("near_-0.9_specs_species.pdf")
	trash <- sapply(X=spec_0.9_idxs, FUN=plotOTU)
	dev.off()


# what are some specific otus around Spec=-1
	spec_1_idxs <- ((1:nrow(specs))[order(abs(-1 - specs$Spec))])[1:10]
	pdf("near_-1_specs_species.pdf")
	trash <- sapply(X=spec_1_idxs, FUN=plotOTU)
	dev.off()

# look at pseudomonas sp
	pseudomonas_idxs <- which( grepl(tax[[2]], pattern="pseudomonas", ignore.case=T) )
	pdf("pseudomonas_species.pdf")
	trash <- sapply(X=pseudomonas_idxs, FUN=plotOTU)
	dev.off()




