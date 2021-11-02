#!/usr/bin/env Rscript


# read in data
library(specificity)
library(ape)
source("../onto2nwk.r")
library("Rcpp")
sourceCpp("../rao_1spv2.cpp")
source("../phy_or_env_spec_new.r")

message("getting data")
data(endophyte)
attach(endophyte)

otutable <- prop_abund(otutable)
otutable <- otutable[, colSums(otutable>0)>10]


message("running specificity")
specs <- phy_or_env_spec2( otutable, env=metadata$Elevation, 
	n_sim=100, p_method="gamma_fit", n_cores=10)
write.table(specs, file="emp_specs_results.txt", sep='\t', quote=FALSE)
save(list="specs", file="specs.rdata")
