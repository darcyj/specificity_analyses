# load endophytes data
  library(specificity)
  attach(endophyte)
  names(endophyte)

# get rid of species that are observed too infrequently to analyze.
  dim(otutable)
  otutable <- prop_abund(otutable)
  otutable_over10 <- occ_threshold(otutable, threshold = 10)




  zscore <- function(emp, sim){(emp - mean(sim))/sd(sim)}
  test_plant_genera <- function(otuid){
    y <- otutable_over10[, colnames(otutable_over10) == otuid]
    y01 <- y>0
    plantfactor <- factor(metadata$PlantGenus)
    plantemp <- table(plantfactor[y01])
    plantperms <- replicate(1000, table(plantfactor[sample(y01)]))
    plantpvals <- sapply(X=1:nrow(plantperms), FUN=function(i){
      expected_mean <- mean(plantperms[i,])
      expected_sd <- sd(plantperms[i,])
      if(plantemp[i]==0){
        list(P=1, Z=0, obs=0, sim_mean=expected_mean, sim_sd=expected_sd)
      }else{
        list(
          P=specificity::p_from_perms_or_gfit( emp=plantemp[i], perm=plantperms[i,] ),
          Z=zscore(plantemp[i], plantperms[i,]),  sim_mean=expected_mean, sim_sd=expected_sd
        )
      }
    })
    planttable <- suppressWarnings( do.call("rbind", plantpvals))
    rownames(planttable) <- levels(plantfactor)
    return(planttable)
  }
  plot_otu <- function(otuid){
    y <- otutable_over10[, colnames(otutable_over10) == otuid]
    par(mfrow=c(2,2))
    plot(y ~ metadata$Elevation, pch=20, ylab=otuid, xlab="Elevation (masl)")
    taxstr<- taxonomy$tax[taxonomy$OTUID == otuid]
    mtext(side=1, line=5, taxstr, cex=0.5, adj=c(0,0), at=0)
    island_nums <- as.numeric( factor(metadata$Island))
    island_labs <- levels(factor(metadata$Island))
    plot(y ~ island_nums, xaxt="n", xlab="Island", pch=20)
    axis(1, labels=island_labs, at=1:length(island_labs))
    plot(y ~ metadata$Evapotranspiration, pch=20, ylab=otuid, xlab="Evapotranspiration (mm/year)")
    plot(y ~ metadata$NDVI, pch=20, ylab=otuid, xlab="NDVI (index)")


  }


  pdf("plotted_otus.pdf")
    plot_otu("ASV14995")
    plot_otu("ASV1089")
    plot_otu("ASV1122")
    plot_otu("ASV1222")
    plot_otu("ASV126")
    plot_otu("ASV12876")
    plot_otu("ASV1328")
    plot_otu("ASV14664")
    plot_otu("ASV16269")
    plot_otu("ASV363")
    plot_otu("ASV610")
    plot_otu("ASV14995")
    plot_otu("ASV13159")
  dev.off()

  plot_otu("ASV8414")



# What plant hosts?
  asvname <- "ASV14995"
  table(metadata$PlantGenus[ otutable[ , colnames(otutable) == asvname] > 0 ])


# plot geo, testing
  library(raster)
  library(rgdal)
  elev <- raster("Elevation_all/Hawaii_WGS84_USGS13arcsec.txt")
  elev2 <- elev
  elev2[elev2 >0] <- NA
  plot(elev, legend=FALSE)
  plot(elev2, col="blue", add=TRUE, legend=FALSE)

# plot geo, function
  plot_otu_geo <- function(otuid, ...){
    otuvec <- otutable_over10[ , colnames(otutable_over10) == otuid, drop=TRUE]
    x <- metadata$Lon[otuvec > 0]
    y <- metadata$Lat[otuvec > 0]
    plot(elev, legend=FALSE, main=otuid, ...)
    plot(elev2, col="blue", add=TRUE, legend=FALSE)
    points(x=x, y=y, pch=20, col="red")
  }

  pdf("plotted_otus_geo.pdf")
  plot_otu_geo("ASV14995", xlim=c(-156,-154.5), ylim=c(19,20))
  plot_otu_geo("ASV1089")
  plot_otu_geo("ASV1122")
  plot_otu_geo("ASV1222")
  plot_otu_geo("ASV126")
  plot_otu_geo("ASV12876")
  plot_otu_geo("ASV1328")
  plot_otu_geo("ASV14664")
  plot_otu_geo("ASV16269")
  plot_otu_geo("ASV363")
  plot_otu_geo("ASV610")
  plot_otu_geo("ASV14995")
  plot_otu_geo("ASV13159")
  dev.off()
  