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

# write down coordinates of the target, the entry AC, PC a MS
coords <- data.frame( target = c(10.91, -39.91, -18.38), entry = c(39.66, -77.13, 57.35),
                      ac = c(-0.97, -53.21, -12.13), pc = c(-1.04, -29.51, -13.87), ms = c(-0.39, -24.40, 55.94),
                      row.names = c("x","y","z")
                      )

# calculate vectors showing midcommisural point (MC) as well as AC-to-PC and AC-to-MS lines
coords <- coords %>% mutate( mc = (ac + pc) / 2, acpc = pc-ac, acms = ms-ac, pcms = ms-pc )

# list values to be computed
dists <- data.frame( prot = c(23.76, 11.96, 10.07, 5.15), comp = NA, row.names = c("acpc","lat_left","ant","inf") )


# ---- validate that changing bases by shifting saggital plane to ACPCMS plane reproduces values from the protocol ----

# prepare a function that will return distance of the target from a plane defined by its normal and one point
dist_val <- function( N = NA, p = NA, t = NA ) {
  
  # start by computing plane's shift (the "d" parameter in the equation of plane in 3D: "ax + by + cz + d = 0")
  d <- -( N %*% p ) # it is a scalar product of the normal and the point by the way
  
  # now we can easily compute the displacement
  # inspired by this one https://math.stackexchange.com/questions/2874812/orthogonal-projection-of-a-point-to-plane
  return( ( ( (N %*% t) + d ) / sqrt( N %*% N ) ) %>% abs() )
    
}

# compute the distances
dists <- dists %>%
  # fill-in the comp column
  mutate(
    # calculate them one-by-one like a good boy
    comp = with( coords, c( sqrt(acpc %*% acpc), # ACPC line
                            dist_val( N = cross( acpc, acms ), p = ac, t = target ), # lateral displacement
                            dist_val( N = acpc, p = pc, t = target), # anterior displacement
                            dist_val( N = cross( cross( acpc, acms ), acpc), p = pc, t = target ) # inferior displacement
                            ) ),
    # add difference of protocol-minus-computed for a good measure
    diff = prot-comp
  )


