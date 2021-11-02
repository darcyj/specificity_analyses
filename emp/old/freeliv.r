# specificity analysis for free-living critters

load("freeliving_inputs.rdata")
library(specificity)

otutable <- prop_abund(otutable)

envmat <- tree2mat(ontophy, x=metadata3k$empo_3, n_cores=5)

spec <- phy_or_env_spec(otutable, 
	env=envmat,
	n_sim=100, p_method="gamma_fit",
	center="mode", n_cores=2
)



 
