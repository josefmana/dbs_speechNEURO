# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# where the reconstructions are?
reco.dir <- "reconstructions"

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


# ----------- read and wrangle the data -----------

# read all electrode reconstruction files into a single list
d0 <- lapply( dir(path=reco.dir),
              function(i) readMat( paste(reco.dir,i,sep="/") )$reco
              ) %>%
  # because 'i' is a file name with '.mat' at the end I use gsub() to extract patients' name only
  # (i.e., characters before the dot)
  `names<-`( gsub( "\\..*", "", dir(path=reco.dir) ) )


# since readMat does not copy field names, let's add them manually
for ( i in names(d0) ) {
  
  # level 1
  names(d0[[i]]) <- c("props","native","scrf","mni")
  
  # level 2
  names(d0[[i]]$props) <- c( ifelse( is.null(d0[[i]]$props[[1]]), NA, "elmodel_right" ),
                             ifelse( is.null(d0[[i]]$props[[2]]), NA, "manually_corrected_right" ),
                             ifelse( is.null(d0[[i]]$props[[3]]), NA, "elmodel_left" ),
                             ifelse( is.null(d0[[i]]$props[[4]]), NA, "manually_corrected_left" )
                             )
  
  for ( j in c("native","scrf","mni") ) {
    
    # need to go via ifelse() cause mni is ordered differently
    names(d0[[i]][[j]]) <- case_when( j == "mni" ~ c("coords_mm","markers","trajectory"),
                                      j != "mni" ~ c("coords_mm","trajectory","markers")
    )
    
    # level 3
    for ( k in c("coords_mm","trajectory") ) {
      # need to work around missing electrodes which are coded as NULL in this structure
      # I assume that at least one electrode is always present
      names(d0[[i]][[j]][[k]]) <- case_when( is.null(d0[[i]][[j]][[k]][[1]]) ~ c(NA,"left"), # label only left if right is missing
                                             is.null(d0[[i]][[j]][[k]][[2]]) ~ c("right",NA), # label only right if left is missing
                                             !list(NULL) %in% d0[[i]][[j]][[k]] ~ c("right","left") # label both if neither is missing
                                             )
    }
    
    # name markers that have somewhat different structure than coordinates and trajectories
    # similarly to above I need to work around missing electrodes
    names(d0[[i]][[j]]$markers) <- case_when( is.null(d0[[i]][[j]]$markers[[1]]) ~ c( rep(NA,4), paste0("left_",c("head","tail","x","y")) ),
                                              is.null(d0[[i]][[j]]$markers[[5]]) ~ c( paste0("right_",c("head","tail","x","y")), rep(NA,4) ),
                                              !list(NULL) %in% d0[[i]][[j]]$markers ~ c( paste0("right_",c("head","tail","x","y")), paste0("left_",c("head","tail","x","y")) )
                                              )
    
    # level 4
    # label axes in all of contacts nad trajectories
    for ( k in c("right","left") ) {
      for ( l in c("coords_mm","trajectory") ) {
        
        if ( is.null(d0[[i]][[j]][[l]][[k]][[1]]) ) next # skip empty data frames (i.e., no electrode locations)
        else d0[[i]][[j]][[l]][[k]][[1]] <- d0[[i]][[j]][[l]][[k]][[1]] %>% `colnames<-`( c("x","y","z") ) # otherwise rename columns
        
      }
      
      # for contacts, rename rows as well to label each contact position by it's own index
      if ( is.null(d0[[i]][[j]]$coords_mm[[k]][[1]]) ) next
      else d0[[i]][[j]]$coords_mm[[k]][[1]] <- d0[[i]][[j]]$coords_mm[[k]][[1]] %>% `rownames<-`(
        # using the following coding convention for the contacts:
        # (i) ventral = smaller, dorsal = bigger,
        # (ii) the same on both sides (the side specific nomenclature can be recovered from the final table by combination of electrode model, hemisphere and contact codes)
        # (iii) Medtronic 3389 electrodes begin with "m_", St. Jude 6172 electrodes begin with "s_"
        if( d0[[i]]$props[[paste0("elmodel_",k)]][1,1] == "Medtronic 3389" ) paste0( "m_", 0:3 )
        # since PaCER does not differentiate directional electrodes in St. Jude (unless in MNI space where it does!)
        # need to add special treatment for coordinates that were not manually corrected
        else if( d0[[i]]$props[[paste0("elmodel_",k)]][1,1] == "St. Jude Directed 6172 (short)" & grepl("pacer",i) & j != "mni" ) paste0( "s_", 1:4 )
        else if( d0[[i]]$props[[paste0("elmodel_",k)]][1,1] == "St. Jude Directed 6172 (short)" & grepl("pacer",i) & j == "mni" ) paste0( "s_", c(1,paste0(2,c("a","b","c")), paste0(3,c("a","b","c")), 4 ) )
        else if( d0[[i]]$props[[paste0("elmodel_",k)]][1,1] == "St. Jude Directed 6172 (short)" & !grepl("pacer",i) ) paste0( "s_", c(1,paste0(2,c("a","b","c")), paste0(3,c("a","b","c")), 4 ) )
      )
    }
    
    # due to their different format, rename markers' columns separately
    for ( k in names(d0[[i]][[j]]$markers) ) {
      if ( is.null(d0[[i]][[j]]$markers[[k]]) ) next
      else d0[[i]][[j]]$markers[[k]] <- d0[[i]][[j]]$markers[[k]] %>% `colnames<-`( c("x","y","z") )
      
    }
  }
}


# ----------- put the data into a single file -----------

# prepare a data frame containing all 'marker' variables (i.e., those that say what any number in outcome
# variables is supposed to mean)
# for now, using the strategy of creating a maximum table (i.e., creating rows for each patient/electrode combination)
# and will trim the table to only those electrodes that do exist in the data set later
# also for now, this code is specific to the problem at hand (i.e., 'speech neurons')
df <- expand.grid( c("x","y","z"), # space axes
                   c( paste0("m_",0:3), paste0("s_",1:4), paste0("s_",c("2a","2b","2c","3a","3b","3c")) ), # contact ids
                   c("left","right"), # hemisphere
                   c("native","scrf","mni"), # space
                   names(d0) # patients
                   ) %>%
  `colnames<-`( c("axis","contact","hemisphere","space", "pid") ) %>%
  select( ncol(.):1 ) # reverse columns' order

# fill-in df by the data from d0
for ( i in names(d0) ) {
  for ( j in c("right","left") ) {
  
    # fill-in electrode model and manually corrected indicator
    # will also print which electrodes are missing
    tryCatch( df[ with( df, pid == i & hemisphere == j ), "elmodel" ] <- d0[[i]]$props[[paste0("elmodel_",j)]][1,1], error = function(e) { print( paste(i,j,"missing",sep=", ") ) } )
    tryCatch( df[ with( df, pid == i & hemisphere == j ), "manually_corrected"] <- d0[[i]]$props[[paste0("manually_corrected_",j)]][1,1], error = function(e) { print( paste(i,j,"missing",sep=", ") ) } )
    
  }
  
  # next loop through all spaces to fill-in contact coordinates, trajectories and markers
  for ( j in c("native","scrf","mni") ) {
    for ( k in c("x","y","z") ) { # loop through axes
      for ( l in names(d0[[i]][[j]]$coords_mm) ) { # loop through hemispheres
        
        if ( is.na(l) ) next
        else {
          
          # fill-in contact coordinates in a long format
          for ( m in rownames(d0[[i]][[j]]$coords_mm[[l]][[1]]) ) df[ with( df, pid == i & space == j & hemisphere == l & axis == k & contact == m ) , "coords_mm" ] <- d0[[i]][[j]]$coords_mm[[l]][[1]][m,k]
          
          # fill-in trajectories in a wide format
          for ( m in 1:nrow(d0[[i]][[j]]$trajectory[[l]][[1]]) ) df[ with( df, pid == i & space == j ) , paste("traj",l,k,m,sep="_") ] <- d0[[i]][[j]]$trajectory[[l]][[1]][m,k]
          
        }
      }
      
      # fill-in markers in a wide format
      for ( l in names(d0[[i]][[j]]$markers) ) {
        if ( is.na(l) ) next
        else df[ with( df, pid == i & space == j ) , paste(m,k,sep="_") ] <- d0[[i]][[j]]$markers[[l]][1,k]
        
      }
    }
  }
}

# in df keep only electrodes that exist (get rid of all dummy rows)
df <- df[ complete.cases(df$coords_mm), ]

# create a folder for tables if it doesn't already exist
if( !dir.exists("tables") ) dir.create( "tables" )

# save df as csv
write.table( df , "tables/dbs_speech_neurons_lead_coords.csv", sep = ",", row.names = F )
