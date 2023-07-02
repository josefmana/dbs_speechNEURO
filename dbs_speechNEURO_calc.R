# This is a script for extracting microelectrodes trajectories from surgery protocols in DBS operation (in the Leksell ram space) and
# then transforming to patients T1 MRI space

# list packages to be used
pkgs <- c("rstudioapi", # setting working directory via RStudio API
          "dplyr", "tidyverse", # data wrangling
          "pracma" # cross product calculation
          )

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# set working directory (works in RStudio only)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# create folders for figures, tables and session info
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("figs", "tabs", "sess"), function(i) if( !dir.exists(i) ) dir.create(i) )

# read coordinates of the target, the entry AC, PC a MS
for ( i in c("dic","leks") ) assign(
  # read DICOM and Leksell coordinates separately
  i, read.csv( paste0( "data/dbs_speechNEURO_", i, ifelse( i == "dic", "om", "ell" ), ".csv" ) ) %>%
    # calculate vectors showing midcommisural point (MC) as well as AC-to-PC and AC-to-MS lines
    mutate( mc = (ac + pc) / 2, acpc = pc-ac, acms = ms-ac, pcms = ms-pc )
)

# read values of distances to be computed and compared and angles to be used for contact location estimation
dists <- read.csv( "data/dbs_speechNEURO_distances.csv", sep = "," )
angs <- read.csv( "data/dbs_speechNEURO_angles.csv", sep = "," )

# space between contacts in the gun (in millimeters)
sp = 1


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
  t <- c(Ry %*% c(Rx %*% v) )
  return(t)
  
}

# prepare vectors that are to be rotated
v0 <- lapply( seq(-3,6,.5), # loop through -3 to 6 mm vertically with 0.5 mm steps
              function(i)
                # first prepare all combinations of x and y coordinates shifted by between-contact-space sp
                # then keep only those that are orthogonal in the xy plane
                # finally, for each combination of x and y add z spanning -3 to +6 mm spaced by 0.5 mm
                # Z is negative because Leksell space is coded in LAI
                expand.grid( Ux = c(0,sp,-sp), Uy = c(0,sp,-sp) ) %>% filter( ( abs(Ux) + abs(Uy) ) < 2 ) %>% add_column( Uz = -i )
              ) %>%
  # pull all vectors to a single file
  do.call( rbind.data.frame, . ) %>%
  # label the contacts (c = central, a = anterior, p = posterior, l = left, r = right)
  mutate( cont = paste0( case_when( (Ux == 0 & Uy == 0) ~ "c", Ux == 1 ~ "l", Ux == -1 ~ "r", Uy == -1 ~ "p", Uy == 1 ~ "a" ), -Uz ), .before = 1 )

# prepare a data frame for shifting coordinates of each patient
shift <- with( angs, lapply(
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

# compute transformation matrixes for each patient
trans <- lapply( setNames(angs$id, angs$id), function(i) l2d_trans( i = i, L = leks, D = dic ) )

# check the correspondence with DICOM space by comparing entry DICOM coordinates with entry transformed Leksell coordinates
sapply( angs$id, function(i) {
  t <- c( trans[[i]] %*% c( leks[ leks$id == i, "entry" ], 1 ) )[-4] # entry in DICOM coordinates according to the transformation
  d <- t - dic[ dic$id == i, "entry" ] # target in DICOM coordinates according to the protocol
  return( sqrt( d %*% d ) ) # return a distance between estimated entry DICOM coordinates and protocol target DICOM coordinates
} ) %>% round(2)

# transform shifted coordinates from Leksell to DICOM space
shift[ , c("Dx","Dy","Dz") ] <- sapply( 1:nrow(shift), function(i) c( trans[[ shift[i,"id"] ]] %*% c( t( shift[ i, c("Ax","Ay","Az") ] ), 1 ) )[-4] ) %>% t()

# save the outcomes
write.table( shift, "tabs/dbs_speechNEURO_dicom_traj.csv", sep = ",", row.names = F, quote = F )


# ---- validate that changing bases by shifting sagital plane to ACPCMS plane reproduces values from the protocol ----

# prepare a function that will return distance of the target from a plane defined by its normal and one point
dist_val <- function( N = NA, p = NA, t = NA ) {
  
  # start by computing plane's shift (the "d" parameter in the equation of plane in 3D: "ax + by + cz + d = 0")
  d <- -( N %*% p ) # it is a scalar product of the normal and the point by the way
  
  # now we can easily compute the displacement
  # inspired by this one https://math.stackexchange.com/questions/2874812/orthogonal-projection-of-a-point-to-plane
  return( ( ( (N %*% t) + d ) / sqrt( N %*% N ) ) %>% abs() )
    
}

# compute all the distances
dists <- dists %>%
  # fill-in the comp column
  mutate(
    # check the dist column to know what to compute
    comp = case_when(
      dist == "acpc" ~ sapply( id, function(i) with( dic[dic$id == i, ], sqrt(acpc %*% acpc) ) ), # ACPC distance
      dist %in% c("left","right") ~ sapply( id, function(i) with( dic[dic$id == i, ], dist_val( N = cross( acpc, acms ), p = mc, t = target ) ) ), # lateral displacement (x axis)
      dist %in% c("ant","post") ~ sapply( id, function(i) with( dic[dic$id == i, ], dist_val( N = acpc, p = mc, t = target ) ) ), # anterior/posterior displacement (y axis)
      dist %in% c("inf","sup") ~ sapply( id, function(i) with( dic[dic$id == i, ], dist_val( N = cross( acpc, cross( acpc, acms ) ), p = mc, t = target ) ) ) # inferior/superior displacement (z axis)
    ),
    # add difference of protocol-minus-computed for a good measure
    diff = prot-comp
  )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/dbs_speechNEURO_calc.txt" )
