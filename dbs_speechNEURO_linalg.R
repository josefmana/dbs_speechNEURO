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
coords <- data.frame( x = c(10.91,39.66,-0.97,-1.04,-0.39), # right/left
                      y = c(-39.91,-77.13,-53.21,-29.51,-24.4), # anterior/posterior
                      z = c(-18.38,57.35,-12.13,-13.87,55.94), # superior/inferior
                      row.names = c("target","entry","ac","pc","ms")
                      )

# calculate vectors showing AC-to-PC and AC-to-MS lines
coords <- coords %>% t() %>% as.data.frame() %>% mutate( acpc = pc-ac, acms = ms-ac, pcms = ms-pc ) %>% t() %>% as.data.frame()


# ---- find out rotation angles to transform xy plane to AC-PC-MS plane ----

# inspired by this https://math.stackexchange.com/questions/2249307/orientation-of-a-3d-plane-using-three-points
# get the normal (via cross product)
N <- cross( as.numeric( coords[ "acpc", ] ), as.numeric( coords[ "acms", ] ) )

# get unit normal vector of the plane (i.e., normalize N)
U <- ( N / sqrt(N %*% N) ) %>% `names<-`( c("x","y","z") )

# calculate the angles from x, y and z axes
for( i in c("x","y","z") ) assign( paste0("theta_",i), asin( U[i] ) )

# prepare rotation matrices
R_x <- matrix( c( 1, 0, 0, 0, cos(theta_x), -sin(theta_x), 0, sin(theta_x), cos(theta_x) ), nrow = 3, ncol = 3, byrow = T )
R_y <- matrix( c( cos(theta_y), 0, sin(theta_y), 0, 1, 0, -sin(theta_y), 0, cos(theta_y) ), nrow = 3, ncol = 3, byrow = T )
R_z <- matrix( c( cos(theta_z), -sin(theta_z), 0, sin(theta_z), cos(theta_z), 0, 0, 0, 1 ), nrow = 3, ncol = 3, byrow = T )


# ---- compute equation of the plane that goes through AC, PC and MS ----

# prepare an object for plane parameters
sag_plane <- c( a = NA, b = NA, c = NA, d = NA )

# get the normal via cross product of acpc and acms vectors
sag_plane[1:3] <- N

# get "d" by putting AC on the plane
sag_plane[4] <- -( N %*% as.numeric( t(coords)[,"ac"] ) )


# ---- calculate some distances as a validation ----

# compute distance from the target to the ACPC line
sqrt( as.numeric( t(coords)[,"acpc"] ) %*% as.numeric( t(coords)[,"acpc"] ) )

# left displacement of the target w.r.t. a plane defined by AC, PC and MS
# inspired by this one https://math.stackexchange.com/questions/2874812/orthogonal-projection-of-a-point-to-plane
( N %*% as.numeric( t(coords)[,"target"] ) + sag_plane["d"] ) /
  sqrt( N %*% N )

# anterior displacement of the target w.r.t. a plane that has ACPC as the norm and includes PC
# compute the plane
cor_plane <- c( with( as.data.frame( t(coords) ), acpc ), -( as.numeric( t(coords)[,"acpc"] ) %*% as.numeric( t(coords)[,"pc"] ) ) )

# check the distance anterior displacement
( as.numeric( t(coords)[,"acpc"] ) %*% as.numeric( t(coords)[,"target"] ) + cor_plane[4] ) /
  sqrt( as.numeric( t(coords)[,"acpc"] ) %*% as.numeric( t(coords)[,"acpc"] ) )

# inferior displacement of the target w.r.t. a plane orthogonal to both previous planes
ax_plane <- c ( cross( N, as.numeric( t(coords)[,"acpc"] ) ),
                -( cross( N, as.numeric( t(coords)[,"acpc"] ) ) %*% as.numeric( t(coords)[,"pc"] ) )
                )

# check the distance
( ax_plane[1:3] %*% as.numeric( t(coords)[,"target"] ) + ax_plane[4] ) /
  sqrt( ax_plane[1:3] %*% ax_plane[1:3] )
