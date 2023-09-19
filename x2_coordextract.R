# This little one will extract electrode coordinates as localized via Lead-DBS and print the result into a neat csv file.


################################################################################################
# ---- THIS SCRIPT IS SUPPOSED TO BE RUN ONLY AFTER ELECTRODES WERE LOCALISED IN LEAD-DBS ---- #
################################################################################################


# list packages to be used
pkgs <- c("rstudioapi", "dplyr", "tidyverse", "R.matlab" )

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# set working directory (works in RStudio only)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# write down some important parameters
d.mri <- "_nogithub/mri/lead" # Where the MRI data is?
d.rec <- "_nogithub/mri/recon" # Where are we going to put single patients' reconstructions in?
d.out <- "_nogithub/coords" # Where are we going to put extracted coordinates table in?

# extract IDs of subjects to compare
id <- list.files( d.mri, recursive = F ) %>% subset( . ,  grepl( "IPN", . ) )
#id <- read.csv( "data/tabs/subjects_to_compare.csv" ) %>% t() %>% as.character()


# ---- extract electrode reconstructions ----

# prepare a folder for trajectories
if( !dir.exists( d.rec ) ) dir.create( d.rec )

# loop through patients to copy their raw ea_reconstruction files
for ( i in id ) file.copy( from = paste( d.mri, i, "ea_reconstruction.mat", sep = "/" ),
                           to = paste0( d.rec, "/", i, ".mat" ) )

# next read the MatLab files
recon <- lapply( setNames(id,id),
                 function(i)
                   tryCatch( readMat( paste0( d.rec, "/", i, ".mat" ) )$reco,
                             error = function(e) print( paste0(i, " reconstruction missing") )
                             )
                 )

# deal with patients tha have only one electrode localised
for ( i in "IPN243" ) {
  
  # fill-in electrode model
  recon[[i]][[1]][ , , 1]$elmodel <- NA
  recon[[i]][[1]][ , , 1]$manually.corrected <- NA
  
  # fill-in the others
  for ( j in 2:4 ) {
    
    recon[[i]][[j]][[1]][[1]][[1]] <- matrix( rep(NA,12), nrow = 4 )
    
    if ( j == 4 ) {
      
      recon[[i]][[j]][[3]][[1]][[1]] <- matrix( rep(NA,150), nrow = 50 )
      for ( k in 1:4 ) recon[[i]][[j]][[2]][ , , 1][[k]] <- matrix( rep(NA,3), nrow = 1 )
      
    } else {
        
        recon[[i]][[j]][[2]][[1]][[1]] <- matrix( rep(NA,150), nrow = 50 )
        for ( k in 1:4 ) recon[[i]][[j]][[3]][ , , 1][[k]] <- matrix( rep(NA,3), nrow = 1 )
      
    }
  }
}

# loop through the files to get them into a nice format
for ( i in names(recon) ) {
  
  if ( is.list(recon[[i]]) ) {
    
    # level 1
    recon[[i]] <- recon[[i]][ , , 1] # explicitly name each facet of level 1
    
    # level 2
    # electrode model - re-format to a data.frame and name rows by hemisphere
    recon[[i]]$props <- with( recon[[i]], rbind.data.frame( props[ , , 1], props[ , , 2] ) ) %>% `rownames<-`( c("right","left") )
    
    # get rid of the "electrode" sub-list, dunno what it is but I believe it contains info for visualisation
    recon[[i]]$electrode <- NULL
    
    # loop through lead coordinates in different space
    for ( j in c("native","scrf","mni") ) {
      
      # if any of these is missing, print a message and march on
      if ( !exists( j, where = recon[[i]] ) ) {
        print( paste0( i, ", ", j, " missing!" ) )
        next
      }
      
      # name the sub-lists in each space accordingly
      recon[[i]][[j]] <- recon[[i]][[j]][ , , 1]
      
      # level 3
      # label contact coordinates and lead trajectories
      for ( k in c("coords.mm","trajectory") ) recon[[i]][[j]][[k]] <- list( right = recon[[i]][[j]][[k]][[1]][[1]] %>% `colnames<-`( c("x","y","z") ),
                                                                             left = recon[[i]][[j]][[k]][[2]][[1]]  %>% `colnames<-`( c("x","y","z") ) )
      
      # next deal with markers
      recon[[i]][[j]]$markers <- list( right = do.call( rbind.data.frame , recon[[i]][[j]]$markers[ , , 1] ) %>% `colnames<-`( c("x","y","z") ),
                                       left = do.call( rbind.data.frame , recon[[i]][[j]]$markers[ , , 2] ) %>% `colnames<-`( c("x","y","z") ) )
        
    }
  }
}


# ---- prepare a data frame with contacts' coordinates ----

# prepare a list to be filled-in with the coordinates
df <- list()

# drop patients with no data
for ( i in names(recon) ) if( !is.list(recon[[i]]) ) recon[[i]] <- NULL

# loop through patients
for ( i in names(recon) ) {
  # prepare a temporary sub-list for each patient
  df[[i]] <- list()
  
  # loop through spaces
  for ( j in c("native","scrf","mni") ) {
    
    # if any of these is missing, print a message and march on
    if ( !exists( j, where = recon[[i]] ) ) {
      print( paste0( i, ", ", j, " missing!" ) )
      next
    }
    
    # in each space create a data.frame with parameters
    df[[i]][[j]] <- do.call( rbind.data.frame, recon[[i]][[j]]$coords.mm ) %>%
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
for ( i in 1:nrow(df) ) df$elmodel[i] <- with( df, recon[[id[i]]]$props[ hemisphere[i], "elmodel" ] )

# save df as csv
write.table( df , paste0( d.out, "/coords_perm.csv" ) , sep = ",", row.names = F, quote = F )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/coordextract.txt" )
