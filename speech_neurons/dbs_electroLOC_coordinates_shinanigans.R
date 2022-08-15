# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list packages to use
pkgs <- c( "tidyverse", "dplyr", # data wrangling
           "ggplot2", "patchwork", # plotting
           "R.matlab" # reda MATLAB files
)

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# create folders for figures, tables and models if they do not already exist
# prints NULL if the folder already exists, prints TRUE if it was newly created
#sapply( c("models","figures","tables"), function(i) if( !dir.exists(i)) dir.create(i) )


# ----------- read and wrangle the data -----------

# read all electrode reconstruction files into a nested list structure
d0 <- list()
for ( i in dir( path="data" ) ) {
  for ( j in dir( path = paste0("data/",i) ) ) {
    # because 'j' is a file name with '.mat' at the end I use gsub() to extract patients' name only
    # (i.e., characters before the dot)
    d0[[i]][[gsub("\\..*","",j)]] <- readMat( paste0("data/",i,"/",j) )$reco
  }
}
