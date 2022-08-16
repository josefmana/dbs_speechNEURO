# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list packages to use
pkgs <- c( "tidyverse", "dplyr", # data wrangling
           "ggplot2", "patchwork", # plotting
           "R.matlab" # read MATLAB files
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

# since readMat does not copy field names, let's add them manually
for ( i in names(d0) ) {
  for ( j in names(d0[[i]]) ) {
    
    # level 1
    names(d0[[i]][[j]]) <- c("props","native","scrf","mni")
    
    # level 2
    names(d0[[i]][[j]]$props) <- c(NA,NA,"elmodel","manually_corrected")
    for ( k in c("native","scrf","mni") ) {
      # need to go via ifelse() cause mni is ordered differently
      names(d0[[i]][[j]][[k]]) <- case_when( k == "mni" ~ c("coords_mm","markers","trajectory"),
                                             k != "mni" ~ c("coords_mm","trajectory","markers")
                                             )
      
      # level 3
      for ( l in c("coords_mm","trajectory") ) {
        # need to work around missing electrodes which are coded as NULL in this structure
        # I assume that at least one electrode is always present
        names(d0[[i]][[j]][[k]][[l]]) <- case_when( is.null(d0[[i]][[j]][[k]][[l]][[1]]) ~ c(NA,"left"), # label only left if right is missing
                                                    is.null(d0[[i]][[j]][[k]][[l]][[2]]) ~ c("right",NA), # label only right if left is missing
                                                    !list(NULL) %in% d0[[i]][[j]][[k]][[l]] ~ c("right","left") # label both if neither is missing
                                                    )
        
      }
      # now name markers that have somewhat different structure than coordinates and trajectories
      # similarly to above I need to work around missing electrodes
      names(d0[[i]][[j]][[k]]$markers) <- case_when( is.null(d0[[i]][[j]][[k]]$markers[[1]]) ~ c( rep(NA,4), paste0("left",c("head","tail","x","y")) ),
                                                     is.null(d0[[i]][[j]][[k]]$markers[[5]]) ~ c( paste0("right_",c("head","tail","x","y")), rep(NA,4) ),
                                                     !list(NULL) %in% d0[[i]][[j]][[k]]$markers ~ c( paste0("right_",c("head","tail","x","y")), paste0("left",c("head","tail","x","y")) )
                                                     )

    }
  }
}
