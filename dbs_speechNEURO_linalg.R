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


# ---- validate that changing bases by shifting saggital plane to ACPCMS plane reproduces values from the protocol ----

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
    # check the dist colmn to know what to compute
    comp = case_when(
      dist == "acpc" ~ sapply( id, function(i) with( coords[coords$id == i, ], sqrt(acpc %*% acpc) ) ), # ACPC distance
      dist %in% c("left","right") ~ sapply( id, function(i) with( coords[coords$id == i, ], dist_val( N = cross( acpc, acms ), p = mc, t = target ) ) ), # lateral displacement
      dist %in% c("ant","post") ~ sapply( id, function(i) with( coords[coords$id == i, ], dist_val( N = acpc, p = mc, t = target ) ) ), # anterior/posterior displacement
      dist %in% c("inf","sup") ~ sapply( id, function(i) with( coords[coords$id == i, ], dist_val( N = cross( acpc, cross( acpc, acms ) ), p = mc, t = target ) ) ) # inferior/superior displacement
    ),
    # add difference of protocol-minus-computed for a good measure
    diff = prot-comp
  )

