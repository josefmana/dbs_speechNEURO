# This is a script for calculating basis changes in pre-surgery MRI native space with respect to peri-operatively used ACPC space.

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

# read coordinates of the target, the entry AC, PC a MS
coords <- read.csv( "data/dbs_speechNEURO_protocols.csv", sep = "," ) %>% # coordinates in MRI space
  # calculate vectors showing midcommisural point (MC) as well as AC-to-PC and AC-to-MS lines
  mutate( mc = (ac + pc) / 2, acpc = pc-ac, acms = ms-ac, pcms = ms-pc )

# read values of sitances to be computed and compared
dists <- read.csv( "data/dbs_speechNEURO_distances.csv", sep = "," )


v

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
      dist == "acpc" ~ sapply( id, function(i) with( coords[coords$id == i, ], sqrt(acpc %*% acpc) ) ), # ACPC distance
      dist %in% c("left","right") ~ sapply( id, function(i) with( coords[coords$id == i, ], dist_val( N = cross( acpc, acms ), p = mc, t = target ) ) ), # lateral displacement (x axis)
      dist %in% c("ant","post") ~ sapply( id, function(i) with( coords[coords$id == i, ], dist_val( N = acpc, p = mc, t = target ) ) ), # anterior/posterior displacement (y axis)
      dist %in% c("inf","sup") ~ sapply( id, function(i) with( coords[coords$id == i, ], dist_val( N = cross( acpc, cross( acpc, acms ) ), p = mc, t = target ) ) ) # inferior/superior displacement (z axis)
    ),
    # add difference of protocol-minus-computed for a good measure
    diff = prot-comp
  )


# ---- get preparations of the central electrodes in ACPCMS space ----

# prepare a function for computing standardized normal vectors
U_comp <- function( N = NA ) return( N / c( sqrt(N %*% N) ) )

# prepare a data set for the standardized normal vectors for the acpcms coordinates
acpc_space <- with( coords, data.frame( id = id, axis = axis ) ) %>%
  # add standardized normals
  mutate(
    Ux = sapply( unique(id), function(i) with( coords[coords$id == i, ], U_comp( cross( acpc, acms ) ) ) ) %>% c(),
    Uy = sapply( unique(id), function(i) with( coords[coords$id == i, ], U_comp( acpc ) ) ) %>% c(),
    Uz = -sapply( unique(id), function(i) with( coords[coords$id == i, ], U_comp( cross( acpc, cross( acpc, acms ) ) ) ) ) %>% c(),
    # add transoformed coordinates of the target
    t0 = sapply( unique(id), function(i) solve( cbind( Ux[id==i], Uy[id==i], Uz[id==i] ), coords[ coords$id == i, "target" ] ) ) %>% c()
  )

# add columns for the other contacts
# each one is spaced 2 milimeters above (+2 on Uz)
for( i in 1:4 ) {
  for ( j in 1:nrow(acpc_space) ) acpc_space[ j, paste0("t",i) ] <- with( acpc_space, ifelse( axis[j] == "z", t0[j] + 2*i, t0[j] ) )
}

# checking the shift along Uz makes sense, the shifted coordinates of t1-t4 should be shifted by the right amount in DICOM space as well
check_shift <- function( pat = NA, t = NA, acpc = acpc_space, coords = coords ) {
  
  # prepare vectors to be compared (both in the DICOM space)
  v <- ( as.matrix( acpc[ acpc$id == pat, c("Ux","Uy","Uz") ] ) %*% acpc[ acpc$id == pat, t ] ) %>% c() # transform tX coordinates from ACPC to DICOM space
  c <- coords[ coords$id == pat, "target" ] # extract coordinates of the target in DICOM space
  diff <-  v - c # difference between v and c represents a vector starting from the target (t0) and ending at selected contact (tX)
  
  # compute length of the diff vector (i.e., square-root of it's self-dot product)
  # should yield 2 for t1, 4 for t2, 6 for t3, and 8 for t4
  return( sqrt( diff %*% diff ) )

}

# check distances in DICOM space for each patient/contact pair
sapply( unique(coords$id), function(i) sapply( paste0("t",1:4), function(j) check_shift( pat = i, t = j, acpc = acpc_space, coords = coords ) ) )


# ---- rotate the central electrodes in ACPCMS space ----

# use the Rodrigues' rotation formula to get coordinates of rotated phantom electrodes calculated in previous step
# see https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula


