# specificity analysis of antarctic cryoconite bacteria

library(specificity)
library(dunn.test)
load("organized_data.rdata")

specs_list <- list()
specs_list$N <- phy_or_env_spec(bac_otutab_ovr10, env=metadata$N, n_sim=1000, n_cores=10)
specs_list$P <- phy_or_env_spec(bac_otutab_ovr10, env=metadata$P, n_sim=1000, n_cores=10)
# removed pH because it didn't show much
# specs_list$pH <- phy_or_env_spec(bac_otutab_ovr10, env=metadata$pH, n_sim=1000, n_cores=10)
specs_list$Geo <- phy_or_env_spec(bac_otutab_ovr10, env=dist_geo, n_sim=1000, n_cores=10)
# removed Fungi since it's highly correlated with algae and also methodologically not independent
# mainly to simplify figure for discussion in paper
# specs_list$Fungi <- phy_or_env_spec(bac_otutab_ovr10, env=dist_fungi, n_sim=1000, n_cores=10)
specs_list$Algae <- phy_or_env_spec(bac_otutab_ovr10, env=dist_algae, n_sim=1000, n_cores=10)
save(list="specs_list", file="specs_list.rdata")

pdf("specs_plot.pdf")
plot_specs_violin(specs_list, cols=c("red", "blue", "black", "forestgreen"))
dev.off()

pdf("specs_cor_plot.pdf")
plot_pairwise_spec(specs_list)
dev.off()

capture.output(file="specificity_version.txt", capture.output(packageVersion("specificity")))

# deal with taxonomy
options(stringsAsFactors=FALSE)
taxonomy <- read.delim("../../data_sets/antarctica/16S/gg_tax/repset_tax_assignments.txt", header=F, sep='\t', stringsAsFactors=FALSE)
taxonomy <- taxonomy[taxonomy[[1]] %in% rownames(specs_list[[1]]),]
tax_split <- lapply(X=taxonomy[[2]], FUN=function(x){
	x <- unlist(strsplit(x, split="; "))
	if(length(x) == 7){
		x <- gsub("^.__", "", x)
	}else if(tolower(x) == "no blast hit"){
		x <- c(x, rep("", 6))
	}
	names(x) <- c("k", "p", "c", "o", "f", "g", "s")
	x
})
tax_split <- do.call("rbind", tax_split)
taxonomy <- data.frame(id=taxonomy[[1]], eval=taxonomy[[2]], hit=taxonomy[[3]], tax_split)
rownames(taxonomy) <- NULL

# save stuff
saveRDS(specs_list, "sl.RDS")
saveRDS(taxonomy, "fd.RDS")





# mantel correlations
library(vegan)
mantel(dist(metadata$N), dist(metadata$p))
mantel(dist(metadata$N), dist_geo)
mantel(dist(metadata$N), dist_algae)

mantel(dist(metadata$p), dist_geo)
mantel(dist(metadata$p), dist_algae)

mantel(dist_algae, dist_geo)





# load back in (done later)
load("organized_data.rdata")
taxonomy <- readRDS("fd.RDS")
specs_list <- readRDS("sl.RDS")
library("dunn.test")

# classify each OTU as overabundant in glaciers.
# not a generic function - just works on 3-level cats
classify_glacier <- function(featurevec, catsvec=metadata$glacier, drawplot=FALSE){
	means_table <- aggregate(featurevec, by=list(catsvec), FUN=mean)
	means_table <- means_table[order(means_table[[2]], decreasing=TRUE),]
	dunn <- dunn.test(featurevec, catsvec)
	# A is most abundant lvl, B middle, C low.
	# bools to simplify logic
	AoverB_sig <- dunn$P[grepl(means_table[1,1], dunn$comparisons) & grepl(means_table[2,1], dunn$comparisons)] < 0.05
	AoverC_sig <- dunn$P[grepl(means_table[1,1], dunn$comparisons) & grepl(means_table[3,1], dunn$comparisons)] < 0.05
	BoverC_sig <- dunn$P[grepl(means_table[2,1], dunn$comparisons) & grepl(means_table[3,1], dunn$comparisons)] < 0.05
	A <- means_table[1,1]; B <- means_table[2,1]; C <- means_table[3,1]
	# build output
	if(AoverB_sig+AoverC_sig+BoverC_sig==0){
		output <- ""
	}else if(!AoverB_sig){
		# A=B
		output <- paste0(c(A, rep("*", sum(AoverB_sig, AoverC_sig)), "=", B, rep("*", sum(BoverC_sig))), collapse="")
	}else{
		# A>B
		if(BoverC_sig){
			output <- paste0(c(A, rep("*", sum(AoverB_sig, AoverC_sig)), ">", B, "*"), collapse="")
		}else{
			output <- paste0(c(A, rep("*", sum(AoverB_sig, AoverC_sig))), collapse="")
		}
	}
	return(output)
}
taxonomy$glacierclass <- sapply(X=taxonomy$id, FUN=function(id){
	classify_glacier(bac_otutab_ovr10[,colnames(bac_otutab_ovr10)==id])
})



# make shiny
# specs_list <- readRDS("sl.RDS")
# taxonomy <- readRDS("fd.RDS")
library(specificity.shiny)
make_specs_app(sl=specs_list, fd=taxonomy, 1, "antarctica_specs_app")
shiny::runApp("antarctica_specs_app")
rsconnect::deployApp("antarctica_specs_app")

