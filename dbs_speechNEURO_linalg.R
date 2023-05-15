# This is a script for calculating basis changes in pre-surgery MRI native space with respect to ACPCMS space.

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

# read values of ditances to be computed and compared and angles to be used for contact location estimation
dists <- read.csv( "data/dbs_speechNEURO_distances.csv", sep = "," )
angs <- read.csv( "data/dbs_speechNEURO_angles.csv", sep = "," )

# space between contacts in millimeters
sp = 1


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


# ---- rotate the central electrodes in ACPCMS space ----

# first prepare function for rotating along x axis (by beta, i.e., anterior angle) and y axis (alpha, i.e., lateral angle)
rot <- function( alpha = NA, beta = NA ) {
  
  # rescale from degrees to radians
  for ( i in c("alpha","beta") ) assign( i, pi * ( get(i)/180 ) )
  
  # prepare rotation matrixes
  Rx <- matrix( c( 1, 0, 0, 0, cos(beta), -sin(beta), 0, sin(beta), cos(beta) ), nrow = 3, byrow = T ) # rotation along x-axis
  Ry <- matrix( c( cos(alpha), 0, sin(alpha), 0, 1, 0, -sin(alpha), 0, cos(alpha) ), nrow = 3, byrow = T ) # rotation along y-axis
  
  # prepare a unit vector in the z-axis direction
  v <- c(0,0,1)
  
  # add rotated vectors and normalize the resulting one
  t <- c( c(Rx %*% v) + c(Ry %*% v) )
  t <- t / sqrt( c( t %*% t ) )
  return(t)

}

# prepare a data frame that includes the shifting vector for each patient (based on their alpha and beta angles)
# to be applied in the ACPCMS space
shift <- with( angs, sapply(
  
  # loop through IDs in the angles data set
  id, function(i)
    rot( alpha = ifelse( side[id == i] == "left", lat_angle[id == i], -lat_angle[id == i]), # rotating to the right incur negative sign
         beta = ant_angle[id == i]
    )
  )
) %>%
  
  # tidy-up the results
  as.data.frame() %>% add_column( axis = c("x","y","z"), .before = 1 ) %>%
  pivot_longer( cols = -1, names_to = "id", values_to = "shift" )


# ---- get preparations for extraction of central electrodes in ACPCMS space ----

# prepare a function for computing standardized normal vectors
U_comp <- function( N = NA ) return( N / c( sqrt(N %*% N) ) )

# prepare a data set of the standardized normal vectors of the acpcms space
acpc <- with( dic, data.frame( id = id, axis = axis ) ) %>%
  # add standardized normals
  mutate(
    Ux = sapply( unique(id), function(i) with( dic[dic$id == i, ], U_comp( cross( acpc, acms ) ) ) ) %>% c(),
    Uy = sapply( unique(id), function(i) with( dic[dic$id == i, ], U_comp( acpc ) ) ) %>% c(),
    Uz = -sapply( unique(id), function(i) with( dic[dic$id == i, ], U_comp( cross( acpc, cross( acpc, acms ) ) ) ) ) %>% c(),
    # add transformed coordinates of the target
    c0 = sapply( unique(id), function(i) solve( cbind( Ux[id==i], Uy[id==i], Uz[id==i] ), dic[dic$id == i, "target" ] ) ) %>% c()
  )

# add columns for the other contacts
# shift by a multiplicative of the shifting vector from shift data frame
for( i in 1:4 ) for ( j in 1:nrow(acpc) ) acpc[ j, paste0("c",i) ] <- with( acpc, c0[j] + sp * i * shift[ shift$id == id[j]  & shift$axis == axis[j], "shift" ] )

# checking the shift along Uz makes sense, the shifted coordinates of t1-t4 should be shifted by the right amount in DICOM space as well
check_shift <- function( pat = NA, t = NA, acpc = acpc, coords = dic ) {
  
  # prepare vectors to be compared (both in the DICOM space)
  v <- ( as.matrix( acpc[ acpc$id == pat, c("Ux","Uy","Uz") ] ) %*% acpc[ acpc$id == pat, t ] ) %>% c() # transform tX coordinates from ACPC to DICOM space
  c <- coords[ coords$id == pat, "target" ] # extract coordinates of the target in DICOM space
  diff <-  v - c # difference between v and c represents a vector starting from the target (t0) and ending at selected contact (tX)
  
  # compute length of the diff vector (i.e., square-root of it's self-dot product)
  # should yield x * sp for each tx
  return( sqrt( diff %*% diff ) )
  
}

# check distances in DICOM space for each patient/contact pair
sapply( unique(dic$id), function(i) sapply( paste0("c",1:4), function(j) check_shift( pat = i, t = j, acpc = acpc, coords = dic ) ) )

# add DICOM coordinates of all 
for ( i in 0:4 ) for ( j in unique(dic$id) ) dic[ dic$id == j, paste0("c",i) ] <- ( as.matrix( acpc[ acpc$id == j, c("Ux","Uy","Uz") ] ) %*% acpc[ acpc$id == j, paste0("c",i) ] ) %>% c()

# save the outcome
write.table( dic, "tabs/dbs_speechNEURO_contacts.csv", sep = ",", quote = F, row.names = F )


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/dbs_speechNEURO_calc.txt" )
