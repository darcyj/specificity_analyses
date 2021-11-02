# Visualization of Hawaii Endophyte specificity

# load analysis results
  library(specificity.shiny)
  library(specificity)
  load("specs_results.rdata")

# visualizing aggregate specificity trends
  speccols <- c("forestgreen", "red", "orange", "blue", "purple", "black")
  # plot_specs_violin(specs_list, label_cex=0.7, cols=speccols)
  pdf("specs.pdf")
  plot_specs_violin(specs_list, label_cex=0.7, cols=speccols)
  dev.off()

# Correlation between specificities
  pdf("pairwise_specs.pdf")
  plot_pairwise_spec(specs_list)
  dev.off()

# Record specificity version
  # packageVersion("specificity")
  capture.output(file="specificity_version.txt", capture.output(packageVersion("specificity")))

# Make shiny app
  appname <- "endophytes_specs_app"
  if(file.exists(appname)){ system(paste0("rm -r ", appname)) }
  make_specs_app(sl=specs_list, fd=endophyte$taxonomy, 1, appname)
  # shiny::runApp(appname)
  rsconnect::deployApp(appname)
