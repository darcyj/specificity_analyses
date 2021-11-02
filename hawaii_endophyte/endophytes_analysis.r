# Analysis of Hawaii Endophyte data
# load specificity package and set working dir
  setwd("~/Documents/first_author/specificity/analyses/hawaii_endophyte/")
  library("specificity")

# load endophytes data
  attach(endophyte)
  names(endophyte)

# get rid of species that are observed too infrequently to analyze.
  dim(otutable)
  otutable <- prop_abund(otutable)
  otutable_over10 <- occ_threshold(otutable, threshold = 10)

# create geographic distance matrix
  geo_distmat <- distcalc(metadata$Lat, metadata$Lon, rownames(metadata))

# run various specificity analyses
  specs_list <- list()
  nc <- 20

  # The following variables are just column vectors from `metadata`. 
  specs_list$NDVI <- phy_or_env_spec(otutable_over10, env=metadata$NDVI, n_cores=nc)
  specs_list$Elevation <- phy_or_env_spec(otutable_over10, env=metadata$Elevation, n_cores=nc)
  specs_list$Evapotranspiration <- phy_or_env_spec(otutable_over10, 
    env=metadata$Evapotranspiration, n_cores=nc)
  specs_list$Rainfall <- phy_or_env_spec(otutable_over10, env=metadata$Rainfall, n_cores=nc)
  specs_list$Host <- phy_or_env_spec(otutable_over10, hosts=metadata$PlantGenus, 
    hosts_phylo=supertree, n_cores=nc)
  specs_list$Geography <- phy_or_env_spec(otutable_over10, env=geo_distmat, 
    n_cores=nc)

# save specificity results
  save(list="specs_list", file="specs_results.rdata")
