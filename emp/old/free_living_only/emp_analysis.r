library(specificity)

otutable <- readRDS("otutable.RDS")
metadata <- readRDS("metadata.RDS")
ontophy <- readRDS("ontophy.RDS")

specs <- phy_or_env_spec( otutable, hosts_phylo=ontophy, 
	hosts=metadata$empo_3, n_sim=100, p_method="gamma_fit", 
	n_cores=20, chunksize=100)


write.table(specs, file="emp_specs_results.txt", sep='\t', quote=FALSE)
save(list="specs", file="specs.rdata")
 
# get taxonomic information for generalist species
