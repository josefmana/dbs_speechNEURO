# This little one will extract electrode coordinates as localized via Lead-DBS and print the result into a neat csv file.

# list packages to be used
pkgs <- c("rstudioapi", # setting working directory via RStudio API
          "dplyr", "tidyverse", # data wrangling
          "R.matlab" ) # reading MatLab files

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# set working directory (works in RStudio only)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# write down some important parameters
d.dir <- "_no_github/imaging" # Where the MRI data is?

# extract IDs of subjects to compare
id <- list.files( d.dir, recursive = F )
#id <- read.csv( "data/tabs/subjects_to_compare.csv" ) %>% t() %>% as.character()


# ---- extract electrode trajectories ----

# prepare a folder for trajectories
if( !dir.exists("data/traj") ) dir.create( "data/traj" )

# loop through patients to copy their raw ea_reconstruction files
for ( i in id ) {
  # prepare a folder for patient "i"
  if( !dir.exists( paste0("data/traj/", i) ) ) dir.create( paste0("data/traj/", i) )
  # copy the ea_reconstruction file
  file.copy( from = paste(d.dir, i, "ea_reconstruction.mat", sep = "/" ),
             to = paste0("data/traj/", i, "/ea_reconstruction.mat") )
}

# next read the MatLab files
d.traj <- lapply( setNames(id,id) , function(i)
  tryCatch( readMat( paste( d.dir, i, "ea_reconstruction.mat", sep = "/") )$reco,
            error = function(e) print( paste0(i, " reconstruction missing") )
            )
  )

# deal with patients tha have only one electrode localised
for ( i in "IPN243" ) {
  
  # fill-in electrode model
  d.traj[[i]][[1]][ , , 1]$elmodel <- NA
  d.traj[[i]][[1]][ , , 1]$manually.corrected <- NA
  
  # fill-in the others
  for ( j in 2:4 ) {
    
    d.traj[[i]][[j]][[1]][[1]][[1]] <- matrix( rep(NA,12), nrow = 4 )
    
    if ( j == 4 ) {
      
      d.traj[[i]][[j]][[3]][[1]][[1]] <- matrix( rep(NA,150), nrow = 50 )
      for ( k in 1:4 ) d.traj[[i]][[j]][[2]][ , , 1][[k]] <- matrix( rep(NA,3), nrow = 1 )
      
    } else {
        
        d.traj[[i]][[j]][[2]][[1]][[1]] <- matrix( rep(NA,150), nrow = 50 )
        for ( k in 1:4 ) d.traj[[i]][[j]][[3]][ , , 1][[k]] <- matrix( rep(NA,3), nrow = 1 )
      
    }
  }
}


# loop through the files to get them into a nice format
for ( i in names(d.traj) ) {
  
  if ( is.list(d.traj[[i]]) ) {
    
    # level 1
    d.traj[[i]] <- d.traj[[i]][ , , 1] # explicitly name each facet of level 1
    
    # level 2
    # electrode model - re-format to a data.frame and name rows by hemisphere
    d.traj[[i]]$props <- with( d.traj[[i]], rbind.data.frame( props[ , , 1], props[ , , 2] ) ) %>% `rownames<-`( c("right","left") )
    
    # get rid of the "electrode" sub-list, dunno what it is but I believe it contains info for visualisation
    d.traj[[i]]$electrode <- NULL
    
    # loop through lead coordinates in different space
    for ( j in c("native","scrf","mni") ) {
      
      # if any of these is missing, print a message and march on
      if ( !exists( j, where = d.traj[[i]] ) ) {
        print( paste0( i, ", ", j, " missing!" ) )
        next
      }
      
      # name the sub-lists in each space accordingly
      d.traj[[i]][[j]] <- d.traj[[i]][[j]][ , , 1]
      
      # level 3
      # label contact coordinates and lead trajectories
      for ( k in c("coords.mm","trajectory") ) d.traj[[i]][[j]][[k]] <- list( right = d.traj[[i]][[j]][[k]][[1]][[1]] %>% `colnames<-`( c("x","y","z") ),
                                                                              left = d.traj[[i]][[j]][[k]][[2]][[1]]  %>% `colnames<-`( c("x","y","z") ) )
      
      # next deal with markers
      d.traj[[i]][[j]]$markers <- list( right = do.call( rbind.data.frame , d.traj[[i]][[j]]$markers[ , , 1] ) %>% `colnames<-`( c("x","y","z") ),
                                        left = do.call( rbind.data.frame , d.traj[[i]][[j]]$markers[ , , 2] ) %>% `colnames<-`( c("x","y","z") ) )
        
    }
  }
}


# ---- prepare a data frame with contacts' coordinates ----

# prepare a list to be filled-in with the coordinates
df <- list()

# drop patients with no data
for ( i in names(d.traj) ) if( !is.list(d.traj[[i]]) ) d.traj[[i]] <- NULL

# loop through patients
for ( i in names(d.traj) ) {
  # prepare a temporary sub-list for each patient
  df[[i]] <- list()
  
  # loop through spaces
  for ( j in c("native","scrf","mni") ) {
    
    # if any of these is missing, print a message and march on
    if ( !exists( j, where = d.traj[[i]] ) ) {
      print( paste0( i, ", ", j, " missing!" ) )
      next
    }
    
    # in each space create a data.frame with parameters
    df[[i]][[j]] <- do.call( rbind.data.frame, d.traj[[i]][[j]]$coords.mm ) %>%
      add_column( contact = sub( ".*\\.", "", rownames(.) ), .before = 1 ) %>% # add contact number (smaller = more ventral)
      add_column( hemisphere = sub( "\\..*", "", rownames(.) ), .before = 1 ) %>% # add hemisphere
      add_column( elmodel = NA, .after = "hemisphere" ) %>% # prepare a column for electrode model
      add_column( space = j , .before = 1 ) # add space
  }
  
  # merge all data frames for each patient
  df[[i]] <- do.call( rbind.data.frame, df[[i]] ) %>% add_column( id = i, .before = 1 )
  
}

# collapse to a single long table and tidy the rownames
df <- do.call( rbind.data.frame, df ) %>% `rownames<-`( 1:nrow(.) )

# finishing touches - add the electrode model
for ( i in 1:nrow(df) ) df$elmodel[i] <- with( df, d.traj[[id[i]]]$props[ hemisphere[i], "elmodel" ] )

# save df as csv
write.table( df , "data/dbs_speech_neurons_lead_coords.csv", sep = ",", row.names = F )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/dbs_coordEXTRACT_session_info.txt" )
