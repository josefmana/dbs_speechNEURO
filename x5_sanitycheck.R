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


# ---- pre-processing ----

# prepare the permanent electrodes data set
df <- d.perm %>%
  filter( space != "scrf" ) %>%
  mutate(
    source = "perm",
    side = case_when( hemisphere == "left" ~ "sin", hemisphere == "right" ~ "dex" ),
    space = ifelse( space == "native", "nat", space ),
    cont = case_when(
      contact == 1 ~ "dist",
      contact == 4 * ifelse( elmodel == "Medtronic 3389", 1, 2 ) ~ "prox"
    )
  ) %>%
  select( id, side, cont, source,x, y, z, space ) %>%
  filter( complete.cases(cont) ) %>%
  pivot_wider( id_cols = c("id","cont","side","source"),
               values_from = c("x","y","z"),
               names_from = space,
               names_glue = "{space}_{.value}"
  ) %>%
  
  full_join( d.check %>% mutate( source = "check", .after = side ) )
