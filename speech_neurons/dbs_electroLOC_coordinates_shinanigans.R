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

# get rid of the "electroCOORD_" label in each folder's name
d0 <- d0 %>% `names<-`( sub(".*_", "", names(.) ) )

# since readMat does not copy field names, let's add them manually
for ( i in names(d0) ) {
  for ( j in names(d0[[i]]) ) {
    
    # level 1
    names(d0[[i]][[j]]) <- c("props","native","scrf","mni")
    
    # level 2
    names(d0[[i]][[j]]$props) <- c( ifelse( is.null(d0[[i]][[j]]$props[[1]]), NA, "elmodel_right" ),
                                    ifelse( is.null(d0[[i]][[j]]$props[[2]]), NA, "manually_corrected_right" ),
                                    ifelse( is.null(d0[[i]][[j]]$props[[3]]), NA, "elmodel_left" ),
                                    ifelse( is.null(d0[[i]][[j]]$props[[4]]), NA, "manually_corrected_left" )
                                    )
    
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
      
      # name markers that have somewhat different structure than coordinates and trajectories
      # similarly to above I need to work around missing electrodes
      names(d0[[i]][[j]][[k]]$markers) <- case_when( is.null(d0[[i]][[j]][[k]]$markers[[1]]) ~ c( rep(NA,4), paste0("left_",c("head","tail","x","y")) ),
                                                     is.null(d0[[i]][[j]][[k]]$markers[[5]]) ~ c( paste0("right_",c("head","tail","x","y")), rep(NA,4) ),
                                                     !list(NULL) %in% d0[[i]][[j]][[k]]$markers ~ c( paste0("right_",c("head","tail","x","y")), paste0("left_",c("head","tail","x","y")) )
                                                     )
      
      # level 4
      # label axes in all of contacts nad trajectories
      for ( l in c("right","left") ) {
        for ( m in c("coords_mm","trajectory") ) {
          
          if ( is.null(d0[[i]][[j]][[k]][[m]][[l]][[1]]) ) next # skip empty data frames (i.e., no electrode locations)
          else d0[[i]][[j]][[k]][[m]][[l]][[1]] <- d0[[i]][[j]][[k]][[m]][[l]][[1]] %>% `colnames<-`( c("x","y","z") ) # otherwise rename columns
        
        }
        
        # for contacts, rename rows as well to label each contact position by it's own index
        if ( is.null(d0[[i]][[j]][[k]]$coords_mm[[l]][[1]]) ) next
        else d0[[i]][[j]][[k]]$coords_mm[[l]][[1]] <- d0[[i]][[j]][[k]]$coords_mm[[l]][[1]] %>% `rownames<-`(
          # using the following coding convention for the contacts:
          # (i) ventral = smaller, dorsal = bigger,
          # (ii) the same on both sides (the side specific nomenclature can be recovered from the final table by combination of electrode model, hemisphere and contact codes)
          # (iii) Medtronic 3389 electrodes begin with "m_", St. Jude 6172 electrodes begin with "s_"
          if( d0[[i]][[j]]$props[[paste0("elmodel_",l)]][1,1] == "Medtronic 3389" ) paste0( "m_", 0:3 )
          # since PaCER does not differentiate directional electrodes in St. Jude (unless in MNI space where it does!)
          # need to add special treatment for coordinates that were not manually corrected
          else if( d0[[i]][[j]]$props[[paste0("elmodel_",l)]][1,1] == "St. Jude Directed 6172 (short)" & grepl("pacer",i) & k != "mni" ) paste0( "s_", 1:4 )
          else if( d0[[i]][[j]]$props[[paste0("elmodel_",l)]][1,1] == "St. Jude Directed 6172 (short)" & grepl("pacer",i) & k == "mni" ) paste0( "s_", c(1,paste0(2,c("a","b","c")), paste0(3,c("a","b","c")), 4 ) )
          else if( d0[[i]][[j]]$props[[paste0("elmodel_",l)]][1,1] == "St. Jude Directed 6172 (short)" & !grepl("pacer",i) ) paste0( "s_", c(1,paste0(2,c("a","b","c")), paste0(3,c("a","b","c")), 4 ) )
        )
      }
      
      # due to their different format, rename markers' columns separately
      for ( l in names(d0[[i]][[j]][[k]]$markers) ) {
        if ( is.null(d0[[i]][[j]][[k]]$markers[[l]]) ) next
        else d0[[i]][[j]][[k]]$markers[[l]] <- d0[[i]][[j]][[k]]$markers[[l]] %>% `colnames<-`( c("x","y","z") )
        
      }
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
                   c("mac","pacer","win"), # type of pre-processing
                   union( names(d0$mac), union( names(d0$pacer), names(d0$win) ) ) ) %>% # selecting all patients including across types of pre-processing
  `colnames<-`( c("axis","contact","hemisphere","space", "processing","pid") ) %>%
  select( ncol(.):1 ) # reverse columns' order

# fill-in df by the data from d0
for ( i in names(d0) ) {
  for ( j in names(d0[[i]]) ) {
    for ( k in c("right","left") ) {
      
      # fill-in electrode model and manually corrected indicator
      # will also print which electrodes are missing
      tryCatch( df[ with( df, processing == i & pid == j & hemisphere == k ), "elmodel" ] <- d0[[i]][[j]]$props[[paste0("elmodel_",k)]][1,1], error = function(e) { print( paste(j,i,k,"missing",sep=",") ) } )
      tryCatch( df[ with( df, processing == i & pid == j & hemisphere == k ), "manually_corrected"] <- d0[[i]][[j]]$props[[paste0("manually_corrected_",k)]][1,1], error = function(e) { print( paste(j,i,k,"missing",sep=",") ) } )
      
    }
    
    # next loop through all spaces to fill-in contact coordinates, trajectories and markers
    for ( k in c("native","scrf","mni") ) {
      for ( l in c("x","y","z") ) { # loop through axes
        for ( m in names(d0[[i]][[j]][[k]]$coords_mm) ) { # loop through hemispheres
          
          if ( is.na(m) ) next
          else {
            
            # fill-in contact coordinates in a long format
            for ( n in rownames(d0[[i]][[j]][[k]]$coords_mm[[m]][[1]]) ) df[ with( df, processing == i & pid == j & space == k & hemisphere == m & axis == l & contact == n ) , "coords_mm" ] <- d0[[i]][[j]][[k]]$coords_mm[[m]][[1]][n,l]
            
            # fill-in trajectories in a wide format
            for ( n in 1:nrow(d0[[i]][[j]][[k]]$trajectory[[m]][[1]]) ) df[ with( df, processing == i & pid == j & space == k ) , paste("traj",m,l,n,sep="_") ] <- d0[[i]][[j]][[k]]$trajectory[[m]][[1]][n,l]
            
          }
        }
        
        # fill-in markers in a wide format
        for ( m in names(d0[[i]][[j]][[k]]$markers) ) {
          if ( is.na(m) ) next
          else df[ with( df, processing == i & pid == j & space == k ) , paste(n,l,sep="_") ] <- d0[[i]][[j]][[k]]$markers[[m]][1,l]
          
        }
      }
    }
  }
}

# in df keep only electrodes that exist (get rid of all dummy rows)
df <- df[ complete.cases(df$coords_mm), ]

# save df as csv
write.table( df , "data/dbs_electroLOC_speech_neurons_coords.csv", sep = ",", row.names = F )
