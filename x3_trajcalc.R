# This is a script for calculating exploratory microelectrodes trajectories from surgery protocols in DBS operation (in the Leksell ram
# space) and then transforming to patients' raw native T1 MRI space ("raw_anat_t1.nii") and Lead-DBS native T1 space ("anat_t1.nii")


######################################################################################################
# ---- THIS SCRIPT IS SUPPOSED TO BE RUN ONLY AFTER MRI WAS COREGISTERED/NORMALISED IN LEAD-DBS ---- #
######################################################################################################


# list packages to be used
pkgs <- c("rstudioapi", "dplyr", "tidyverse", "pracma", "rmatio", "RNifti" ) # pracma is for cross product calculation

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# set working directory (works in RStudio only)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# Where do the data live?
indir <- "_nogithub/prots/sums" # summaries (a.k.a. 'inputs directory')
mrdir <- "_nogithub/mri/lead" # MRIs processed in Lead-DBS (a.k.a. 'MR directory')

# create folders for outcomes and session info
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c( "sess", "_nogithub", "_nogithub/coords" ), function(i) if( !dir.exists(i) ) dir.create(i) )

# read coordinates of target, entry, AC, PC a MS
for ( i in c("dic","leks") ) assign(
  
  # read DICOM and Leksell coordinates separately
  paste0( "d.", i ), read.csv( paste0( indir, "/", i, ".csv" ) ) %>%
    
    # calculate vectors showing midcommisural point (MC) as well as AC-to-PC and AC-to-MS lines
    mutate( mc = (ac + pc) / 2, acpc = pc-ac, acms = ms-ac, pcms = ms-pc )

)

# read values of distances to be computed and compared and angles to be used for contact location estimation
# as well as algorithms used for normalisation vi Lead-DBS
d.dists <- read.csv( paste0( indir, "/dists.csv" ), sep = "," )
d.angs <- read.csv( paste0( indir, "/angs.csv" ), sep = "," )
d.algs <- read.csv( paste0( mrdir, "/norms.csv"), sep = "," )

# some spaces (in millimeters)
sp = 2 # space between electrodes in the gun
distd = 0.75 # space between center of distal permanent contact away and fringe of this contact
proxd = 6.75 # space between center of proximal permanent contact away and fringe of this contact


# ---- rotate in Leksell space ----

# first prepare function for rotating along x axis (by beta, i.e., 90°-minus-ring_angle) and y axis (alpha, i.e., 90°-minus-arc_angle
# resulting in positive angles for left- whereas negative for right-side electrodes)
# note that this procedure depends on mounting and is only valid for lateral_right
rot <- function( alpha = NA, beta = NA, v = c(0,0,1) ) {
  
  # re-scale from degrees to radians
  for ( i in c("alpha","beta") ) assign( i, pi * ( get(i)/180 ) )
  
  # prepare rotation matrices
  Rx <- matrix( c( 1, 0, 0, 0, cos(beta), -sin(beta), 0, sin(beta), cos(beta) ), nrow = 3, byrow = T ) # rotation along x-axis
  Ry <- matrix( c( cos(-alpha), 0, sin(-alpha), 0, 1, 0, -sin(-alpha), 0, cos(-alpha) ), nrow = 3, byrow = T ) # rotation along y-axis
                                                                                                       # sign reversing the angle because
                                                                                                       # DICOM is LPS, Leksell is LAI
  
  # rotate 
  t <- c(Rx %*% c(Ry %*% v) )
  return(t)
  
}

# prepare vectors that are to be rotated
v0 <- lapply( seq(3,-6,-.5), # loop through -3 to 6 mm vertically with 0.5 mm steps
              function(i)
                # first prepare all combinations of x and y coordinates shifted by between-contact-space sp
                # then keep only those that are orthogonal in the xy plane
                # finally, for each combination of x and y add z spanning -3 to +6 mm spaced by 0.5 mm
                # Z is negative because Leksell space is coded in LAI
                expand.grid( Ux = c(0,sp,-sp), Uy = c(0,sp,-sp) ) %>% filter( ( abs(Ux) + abs(Uy) ) < 2*sp ) %>% add_column( Uz = i )
              ) %>%
  
  # pull all vectors to a single file
  do.call( rbind.data.frame, . ) %>%
  
  # label the contacts (c = central, a = anterior, p = posterior, l = left, r = right)
  mutate( cont = paste0( case_when( (Ux == 0 & Uy == 0) ~ "c", Ux == sp ~ "l", Ux == -sp ~ "r", Uy == -sp ~ "p", Uy == sp ~ "a" ), Uz ), .before = 1 )

# prepare a data frame for shifting coordinates of each patient
shift <- with( d.angs, lapply(
  
  # loop through all the patients
  setNames( id, id ), function(i)
    
    # bind the template created before with patient-specific rotations
    cbind.data.frame(
      v0, # the template
      sapply(
        
        # calculate rotated coordinates for each vector separately
        1:nrow(v0), function(j) rot( alpha = ( arc_angle[id == i] - 90 ), beta = ( 90 - ring_angle[id == i] ), v = c( t( v0[j,-1] ) ) )
      
        # re-format and tidy-up for each patient separately
      ) %>% t() %>% `colnames<-`( paste0( "S", c("x","y","z") ) )
    
      # shift patient's target coordinates in Leksell space
    ) %>% mutate( Ax = ( Lx[id == i] + Sx ), Ay = ( Ly[id == i] + Sy ), Az = ( Lz[id == i] + Sz ) )
  
  # merge data sets across patients
  ) %>% do.call( rbind.data.frame, . )

# add patient identificator
) %>% mutate( id = sub( "\\..*", "", rownames(.) ), .before = 1 )

# re-label left/right electrodes to lateral/medial
for ( i in unique(shift$id) ) {
  
  # for patients with left-sided electrode left -> lateral and right -> medial
  if ( !( i %in% with( d.angs, id[ side == "right" ] ) ) ) shift[ shift$id == i, "cont" ] <- shift[ shift$id == i, "cont" ] %>% sub( "r", "m", . )
  
  # for patients with right-sided electrode left -> medial and right -> lateral
  else shift[ shift$id == i, "cont" ] <- shift[ shift$id == i, "cont" ] %>% sub( "l", "m", . ) %>% sub( "r", "l", . )
  
}


# --- mapping from Leksell to DICOM space ----

# prepare a function for computing transformation matrix from Leksell to DICOM space
l2d_trans <- function( i = NA, L = leks, D = dic ) {
  
  # Leksell space coordinates of AC, PC, MS and target (Source)
  s <- with( L, cbind( rbind( ac[ id == i ], pc[ id == i ], ms[ id == i ], target[ id == i ] ) ) )
  
  # extract DICOM coordinates for AC, PC and MS points
  d <- t( D[ D$id == i, c("ac","pc","ms",'target') ] )
  
  # compute the transformation matrix
  # use homogeneous coordinates and solve
  # s * tx = d
  A <- t( solve( cbind(s,1) , cbind(d,1) ) )
  return(A)
  
}

# compute transformation matrices for each patient
trans <- lapply( setNames(d.angs$id, d.angs$id), function(i) l2d_trans( i = i, L = d.leks, D = d.dic ) )

# check the correspondence with DICOM space by comparing entry DICOM coordinates with entry transformed Leksell coordinates
sapply( d.angs$id, function(i) {
  
  t <- c( trans[[i]] %*% c( d.leks[ d.leks$id == i, "entry" ], 1 ) )[-4] # entry in DICOM coordinates according to the transformation
  d <- t - d.dic[ d.dic$id == i, "entry" ] # target in DICOM coordinates according to the protocol
  return( sqrt( d %*% d ) ) # return a distance between estimated entry DICOM coordinates and protocol target DICOM coordinates

} ) %>% round(2)

# transform shifted coordinates from Leksell to DICOM space
shift[ , c("Dx","Dy","Dz") ] <- sapply( 1:nrow(shift), function(i) c( trans[[ shift[i,"id"] ]] %*% c( t( shift[ i, c("Ax","Ay","Az") ] ), 1 ) )[-4] ) %>% t()

# because the real DICOM (as read by Slicer) has flipped x and y axes with respect to protocols, flip them now
shift[ , c("Dx","Dy") ] <- -shift[ , c("Dx","Dy") ]


# ---- mapping from raw_anat_t1 to anat_t1 space ----

# add "anat_t1.nii" coordinates via transformation matrix in Lead MRI folders
shift[ , c("Nx","Ny","Nz") ] <- sapply( 1:nrow(shift), function(i) {
  
  # start by extracting ids of patients with MRIs processed and raw dicom coordinates from the row
  pats <- list.files( mrdir, recursive = F )
  x <- c( t( shift[ i , c("Dx","Dy","Dz") ] ) , 1 )
  
  # if the row belongs to a patient without MRI processed, fill with NAs
  if ( !( shift$id[i] %in% pats ) ) return( rep(NA,3) )
  
  # otherwise transform raw_anat_t1 to anat_t1 coordinates
  else {
    
    # extract transformation matrix and compute the transformation
    t <- solve( read.mat( paste( mrdir, shift$id[i], "ea_precoregtransformation.mat", sep = "/" ) )$tmat )
    x <- t %*% x
    
    # return the result
    return( x[-4] )
    
  }
  
} ) %>% t()

# loop through patients and add a column with algorithm used for normalisation (to be used in further steps in Matlab)
for ( i in unique(shift$id) ) shift[ shift$id == i, "norm_algo" ] <- with( d.algs, algo[ id == i ] )

# save the outcomes
write.table(
  
  # keep only raw_anat_t1 and anat_t1 coordinates and rename columns
  shift[ , grepl("id|cont|D|N|algo", names(shift) ) ] %>% `colnames<-`( c( "id", "cont", paste0("raw_",c("x","y","z")), paste0( "nat_",c("x","y","z")), "norm_algo" ) ),
  "_nogithub/coords/coords_expl.csv", sep = ",", row.names = F, quote = F
  
)


# ---- prepare the distal and proximal contacts estimation ----

# this part will use rot() and l2d_trans() on contacts derived from surgery protocols
# the result will be a table including coordinates of expected centers of the permanent electrode
# these will be later compared to Lead-DBS reconstructions of the permanent electrode

# read the final target for permanent electrode according to protocols
d.targ <- read.csv( paste0( indir, "/explorprots.csv" ) , sep = "," )

# begin by preparing contact IDs that correspond with the most distal and most proximal contact of permanent electrodes
d.targ <- d.targ %>% mutate(
  
  # extract ditinct contact information (microelectrodes and depth chosen)
  cont = substr( chosen_one, 1, 1 ), # the position chosen
  targ = as.numeric( gsub( "\\D", "", chosen_one) ) %>% ifelse( . > 10, ./10, . ) %>% ifelse( grepl( "-", chosen_one), -. , . ),
  
  # add centers of distal and proximal contacts based on
  # https://www.manualslib.com/manual/1723343/Medtronic-Dbs-3389.html?page=7#manual
  # these correspond to Uz coordinate of v0 above
  dist = targ - distd, # the center of distal contact of permanent electrode is 0.75 mm "above" the distal fringe
  prox = targ - proxd, # the center of proximal contact of permanent electrode is 6.75 mm "above" the distal fringe
  
  # based on target microelectrode (targ) add Ux and Uy coordinates
  Ux = case_when( cont == "c" ~ 0, cont == "a" ~ 0, cont == "p" ~ 0, cont == "l" ~ ifelse( side == "sin", sp, -sp), cont == "m" ~ ifelse( side == "sin", -sp, sp) ),
  Uy = case_when( cont == "c" ~ 0, cont == "a" ~ sp, cont == "p" ~ -sp, cont == "l" ~ 0, cont == "m" ~ 0 )

) %>%
  
  # pivot such that there are rows per contact (not patient)
  pivot_longer( cols = c("dist","prox"), names_to = "contact", values_to = "Uz" ) %>% as.data.frame()

# rotate the vector
d.targ <- d.targ %>%
  
  # add columns with rotated vectors that are:
  # beginning in point c(0,0,0) (labelled c("Sx","Sy",Sz")) and
  # shifted to Leksell space (labelled c("Ax","Ay","Az"))
  cbind( . , sapply( 1:nrow(d.targ), function(i) {
    
    # extract patient id for current row
    pat = d.targ[i,"id"]
    
    # calculate rotated vectors
    with( d.angs,
          rot( alpha = ( arc_angle[ id == pat ] - 90 ),
               beta = ( 90 - ring_angle[ id == pat ] ),
               v = with( d.targ, c( Ux[i], Uy[i], Uz[i] ) )
          ) %>%
            
            # shift vectors (defined by c("Sx","Sy",Sz")) to Leksell space
            t() %>% as.data.frame() %>% `colnames<-`( paste0( "S", c("x","y","z") ) ) %>%
            mutate( Ax = ( Lx[id == pat] + Sx ), Ay = ( Ly[id == pat] + Sy ), Az = ( Lz[id == pat] + Sz ) )
          
    )
    
  } ) %>% t()
  
) %>% mutate_if( is.list, as.numeric ) # formatting shinanigans

# transform shifted coordinates from Leksell to DICOM space
d.targ[ , c("Dx","Dy","Dz") ] <- sapply( 1:nrow(d.targ), function(i) c( trans[[ d.targ[i,"id"] ]] %*% c( t( d.targ[ i, c("Ax","Ay","Az") ] ), 1 ) )[-4] ) %>% t()

# because the real DICOM (as read by Slicer) has flipped x and y axes with respect to protocols, flip them now
d.targ[ , c("Dx","Dy") ] <- -d.targ[ , c("Dx","Dy") ] 

# add "anat_t1.nii" coordinates via transformation matrix in Lead MRI folders
d.targ[ , c("Nx","Ny","Nz") ] <- sapply( 1:nrow(d.targ), function(i) {
  
  # start by extracting ids of patients with MRIs processed and raw dicom coordinates from the row
  pats <- list.files( mrdir, recursive = F )
  x <- c( t( d.targ[ i , c("Dx","Dy","Dz") ] ) , 1 )
  
  # if the row belongs to a patient without MRI processed, fill with NAs
  if ( !( shift$id[i] %in% pats ) ) return( rep(NA,3) )
  
  # otherwise transform raw_anat_t1 to anat_t1 coordinates
  else {
    
    # extract transformation matrix and compute the transformation
    t <- solve( read.mat( paste( mrdir, d.targ$id[i], "ea_precoregtransformation.mat", sep = "/" ) )$tmat )
    x <- t %*% x
    
    # return the result
    return( x[-4] )
    
  }
  
} ) %>% t()

# loop through patients and add a column with algorithm used for normalisation (to be used in further steps in Matlab)
for ( i in unique(d.targ$id) ) d.targ[ d.targ$id == i, "norm_algo" ] <- with( d.algs, algo[ id == i ] )

# save the outcomes
write.table(
  
  # keep only raw_anat_t1 and anat_t1 coordinates and rename columns
  d.targ[ , grepl("id|contact|D|N|algo", names(d.targ) ) ] %>%
    `colnames<-`( c( "id", "side","cont", paste0("raw_",c("x","y","z")), paste0( "nat_",c("x","y","z")), "norm_algo" ) ),
  
  # select output file
  "_nogithub/coords/coords_check.csv", sep = ",", row.names = F, quote = F

)


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/trajcalc.txt" )
