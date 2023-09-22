# This is a script that run the following sanity checks:
# 1) compare contact coordinates inferred from surgery protocols vs coordinates reconstructed via Lead-DBS


###################################################################################################
# ---- THIS SCRIPT IS SUPPOSED TO BE RUN ONLY AFTER NATIVE-TO-MNI TRANSFORMATION WAS APPLIED ---- #
###################################################################################################


# The following scripts ought to be run before this one:
# x2_coordextract.R
# x3_trajcalc.R
# x4_nattomni.m

# list packages to be used
pkgs <- c("rstudioapi", "dplyr", "tidyverse" )

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# set working directory (works in RStudio only)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# where the data is?
indir <- "_nogithub/coords" # directory containing inputs

# read the data
for ( i in c("perm","check") ) assign( paste0( "d.", i ), read.csv( paste0( indir, "/coords_", i, ".csv" ), sep = "," ) )


# ---- in-house functions ---

# angle between two vectors
ang <- function( x, y ) {
  
  dp <- x %*% y # dot (inner) product
  nx <- sqrt( x %*% x ) # norm of vector x
  ny <- sqrt( y %*% y ) # norm of vector y
  
  theta <- acos( dp / (nx * ny ) ) # final angle in radians
  return( as.numeric( ( theta / pi ) * 180 ) ) # return the results in degrees
}


# ---- pre-processing ----

# prepare a common data set for angles as well as distances calculatio
d0 <- d.perm %>%
  
  # keep only native and MNI coordinates
  filter( space != "scrf" ) %>%
  
  # reformat some variables so the files mix well
  mutate(
    source = "perm", # label the source of estimation ('perm' for permanent electrode)
    side = case_when( hemisphere == "left" ~ "sin", hemisphere == "right" ~ "dex" ), # re-code according to d.targ
    space = ifelse( space == "native", "nat", space ), # re-code according to d.targ
    # keep only the most proximal and distal contacts (the others won't work for St. Jude)
    cont = case_when(
      contact == 1 ~ "dist",
      contact == 4 * ifelse( elmodel == "Medtronic 3389", 1, 2 ) ~ "prox"
    )
  ) %>%
  
  # keep only variables of interest and rows with data
  select( id, side, cont, source, x, y, z, space ) %>%
  filter( complete.cases(cont) ) %>%
  
  # re-format to a format ready to joining d.targ
  pivot_wider(
    id_cols = c("id","cont","side","source"),
    values_from = c("x","y","z"),
    names_from = space,
    names_glue = "{space}_{.value}"
  ) %>%
  
  # add d.targ and get rid of raw coordinates (which are not part of Lead-DBS output)
  full_join( d.check %>% mutate( source = "check", .after = side ) ) %>%
  mutate( ctqual = factor( ifelse( id %in% paste0( "IPN", c(243,244,248,256) ), "low", "ok" ), levels = c("ok","low"), ordered = T ) ) %>% # label low-qualty CT subjects
  select( -norm_algo, -contains("raw") ) %>%
  
  # put 'perm' and 'check' data next to each other to locate irrelevant rows (side not evaluated)
  pivot_wider(
    names_from = source,
    values_from = contains("_"),
    id_cols = c("id","cont","side","ctqual")
  ) %>%
  
  # keep only cases (sides) used for the experiment
  filter( complete.cases(nat_x_check) ) %>%
  
  # prepare unique rows for patient/space/source combinations
  pivot_longer(
    cols = contains("_"),
    names_to = c("space","axis","source"),
    names_pattern = "(.*)_(.*)_(.*)",
    values_to = "estimate"
  )

# calculate Eucledian distance between 'perm' and 'check' estimates
d1 <- d0  %>%
  
  # return coordinates to their own distinct columns
  pivot_wider(
    names_from = c("source","axis"),
    values_from = "estimate",
    id_cols = c("id","cont","side","space","ctqual")
  ) %>%
  
  # calculate the distances
  mutate(
    x_diff = perm_x - check_x,
    y_diff = perm_y - check_y,
    z_diff = perm_z - check_z,
    dist_mm = sqrt( x_diff^2 + y_diff^2 + z_diff^2 )
  )

# calculate angles w.r.t. YZ plane (alpha) and XZ plane (beta)
d2 <- d0 %>%
  
  # return coordinates to their own distinct columns
  pivot_wider(
    names_from = c("cont","axis"),
    values_from = "estimate",
    id_cols = c("id","side","space","source","ctqual")
  ) %>%
  
  # calculate the angles
  mutate(
    x = - prox_x + dist_x,
    y = - prox_y + dist_y,
    z = - prox_z + dist_z,
    alpha = sapply( 1:nrow(.), function(i) 90-ang( c(x[i],y[i],z[i]), c(1,0,0) ) ),
    beta = sapply( 1:nrow(.), function(i) 90-ang( c(x[i],y[i],z[i]), c(0,-1,0) ) )
  ) %>%
  
  # put the angles side-to-side
  pivot_wider(
    names_from = "source",
    values_from = c("alpha","beta"),
    id_cols = c("id","side","space","ctqual")
  ) %>%
  
  # compute angle differences
  mutate(
    alpha_diff = alpha_perm - alpha_check,
    beta_diff = beta_perm - beta_check
  )
  