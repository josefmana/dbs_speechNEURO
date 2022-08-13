setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# list required packages into a character object
pkgs <- c(
  "dplyr", # for objects manipulation
  "tidyverse", # for objects manipulation when dplyr ain't enough
  "brms", # fitting models
  "ggplot" # for plotting
)

# load required packages
for ( i in pkgs ) {
  if ( !require( i , character.only = T ) ){
    # it's important to have the 'character.only = T' command here, not remember why though
    install.packages(i)
    library(i)
  }
}

# create folders "models", "figures" and "tables" to store results in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply(
  c("models", "figures", "tables"), # folders to be created
  function(i)
    if( !dir.exists(i) ) dir.create(i)
)

# read the data
d <- read.csv( "dbs_coord.csv", sep = ",")

# calculate raw differences to get an idea
d1 <- d %>%
  pivot_wider( names_from = "rater", values_from = c("x","y","z") ) %>%
  mutate(x_dif = x_a - x_j,
         y_dif = y_a - y_j,
         z_dif = z_a - z_j)

# continue with MNI coordinates only as native space is very weird
d.mni <- d %>% slice( which(space == "mni") ) %>%
  mutate( rater = case_when( rater == "a" ~ "Andrej", rater == "j" ~ "Josef"),
          contact = ifelse( manuf == "medtronic", paste0("k", contact), contact) ) %>%
  mutate_if( is.character, as.factor )


# ----------- multivatiate multilevel models for MNI coordinates  -----------

# first set some contrast
contrasts(d.mni$rater) <- contr.sum(2)/2 # Andrej = 0.5, Josef = -0.5
contrasts(d.mni$hemisphere) <- contr.sum(2)/2 # left = 0.5, right = -0.5

# list for models gathering
m <- list()

# fist model with rater as a random-effect (i.e., generalizable to new raters)
m$rand <- brm(
  bf( mvbind(x,y,z) ~ 1 + hemisphere + (1 | rater ) + (1 | contact/manuf ) + (1 | subject ) ) + set_rescor(T),
  family = gaussian(), data = d.mni, chains = 4, iter = 2000, warmup = 500,
  cores = 4, seed = 87542, control = list(adapt_delta = .99), file = "models/m1_random_raters"
)

# second model with rater as a fixed-effect for a direct though non-genralizable comparison
m$fix <- brm(
  bf( mvbind(x,y,z) ~ 1 + hemisphere + rater + (1 | contact/manuf ) + (1 | subject ) ) + set_rescor(T),
  family = gaussian(), data = d.mni, chains = 4, iter = 2000, warmup = 500,
  cores = 4, seed = 87542, control = list(adapt_delta = .99), file = "models/m2_fixed_raters"
)

# make it simpler by considering fixed position of each contact (ventral to dorsal)
# instead of manufacturer/contact (random) pairs
d.mni <- d.mni %>%
  mutate(
    loc = case_when(
      ( manuf == "medtronic" & contact %in% c("k0", "k8") ) | (manuf == "st_jude" & contact %in% c("1","9") ) ~ "c1",
      ( manuf == "medtronic" & contact %in% c("k1", "k9") ) | (manuf == "st_jude" & contact %in% c("2a","2b","2c","10a","10b","10c") ) ~ "c2",
      ( manuf == "medtronic" & contact %in% c("k2", "k10") ) | (manuf == "st_jude" & contact %in% c("3a","3b","3c","11a","11b","11c") ) ~ "c3",
      ( manuf == "medtronic" & contact %in% c("k3", "k11") ) | (manuf == "st_jude" & contact %in% c("4","12") ) ~ "c4"
    ) %>% as.ordered
  )

# fit a model with rater as a fixed-effect and electrode location as a monotonic predictor
m$mono <- brm(
  bf( mvbind(x,y,z) ~ 1 + hemisphere + rater + mo(loc) + (1 | subject ) ) + set_rescor(T),
  family = gaussian(), data = d.mni, chains = 4, iter = 2000, warmup = 500,
  cores = 4, seed = 87542, control = list(adapt_delta = .99), file = "models/m3_monotonic_location"
)

# fit a model with monotic location * rater interaction
m$inter <- brm(
  bf( mvbind(x,y,z) ~ 1 + hemisphere + rater * mo(loc) + (1 | subject ) ) + set_rescor(T),
  family = gaussian(), data = d.mni, chains = 4, iter = 2000, warmup = 500,
  cores = 4, seed = 87542, control = list(adapt_delta = .99), file = "models/m4_monotonic_interaction"
)

# full fixed-effects (i.e., 'descriptive') with (almost) all interactions included 
m$desc <- brm(
  bf( mvbind(x,y,z) ~ 1 + hemisphere * rater * subject + mo(loc) ) + set_rescor(T),
  family = gaussian(), data = d.mni, chains = 4, iter = 2000, warmup = 500,
  cores = 4, seed = 87542, control = list(adapt_delta = .99), file = "models/m5_full_fixed"
)

# add PSIS-LOO to the models
for ( i in names(m) ) {
  for ( j in c("x","y","z") ) m[[i]] <- add_criterion( m[[i]], criterion = "loo", resp = j)
}

# compare the models
loo( m$rand, m$fix, m$mono, m$inter, m$desc, resp = "x" )
loo( m$rand, m$fix, m$mono, m$inter, m$desc, resp = "y" )
loo( m$rand, m$fix, m$mono, m$inter, m$desc, resp = "z" )


# ----------- visualization  -----------

# prepare hemisphere-specific means for centering
means <- sapply( c("r","l") , function(i)
  sapply( c("x","y","z") , function(j)
    mean(d.mni[d.mni$hemisphere == i, j ] )
    )
  )

# prepare data frame for plotting
d.mni.img <- d.mni %>%
  # center all coordinates
  mutate( x = case_when( hemisphere == "l" ~ x - means["x","l"], hemisphere == "r" ~ x - means["x","r"] ),
          y = case_when( hemisphere == "l" ~ y - means["y","l"], hemisphere == "r" ~ y - means["y","r"] ),
          z = case_when( hemisphere == "l" ~ z - means["z","l"], hemisphere == "r" ~ z - means["z","r"] )
          ) %>%
  # change to a suitable format for plotting
  pivot_longer( cols = c("x","y","z"), names_to = "axis", values_to = "coord") %>%
  pivot_wider( values_from = "coord", names_from = "rater" ) %>%
  mutate( Subject = subject,
          `Manuf.` = case_when(manuf == "medtronic" ~ "Medtronic", manuf == "st_jude" ~ "St.Jude"),
          Hemisphere = case_when(hemisphere == "l" ~ "Left", hemisphere == "r" ~ "Right")
          )

# plot it
d.mni.img %>%
  ggplot( aes(x = Andrej, y = Josef, color = Subject ) ) +
  geom_abline( intercept = 0, slope = 1 , size = 1.25 ) +
  geom_abline( intercept = -1, slope = 1, linetype = "dotted", size = 1.1) +
  geom_abline( intercept = +1, slope = 1, linetype = "dotted", size = 1.1) +
  geom_jitter( size = 3 ) +
  facet_grid( axis ~ Hemisphere ) +
  theme_classic( base_size = 18)

# save the figure
ggsave( "figures/comparisons in MNI space.png", width = 1.5 * 10.1, height = 1.5 * 5.39, dpi = "retina" )

# plot histograms
d1.mni %>%
  pivot_longer( cols = paste0( c("x","y","z"), "_dif"), values_to = "difference (mm)", names_to = "axis" ) %>%
  mutate( Hemisphere = case_when(hemisphere == "l" ~ "Left", hemisphere == "r" ~ "Right"),
          axis = substr( axis, 1, 1) ) %>%
  ggplot( aes(x = `difference (mm)`, fill = stat( abs(x) > 1 ) ) ) +
  geom_histogram() +
  geom_vline( xintercept = -1, linetype = "dotted", size = 1.1 ) +
  geom_vline( xintercept = +1, linetype = "dotted", size = 1.1 ) +
  labs( y = "") +
  facet_grid( axis ~ Hemisphere ) +
  scale_fill_manual( values = c("grey69", "skyblue") ) +
  theme_bw( base_size = 18) +
  theme( legend.position = "none" )

# save the histograms
ggsave( "figures/histogram of differences in MNI space.png", width = 1.5 * 10.1, height = 1.5 * 5.39, dpi = "retina" )


# ----------- some redundant stats  -----------

# prepare and array for correlations and t-test
t <- array( data = NA, dim = c(2, 3, 2), dimnames = list( c("l", "r"), c("x","y","z"), c("cor","t.test") ) )

# loop through hemispheres and axes
for ( i in c("l","r") ) {
  for ( j in c("x","y","z") ) {
    # (Pearson`s) correlation
    # pre-compute
    cor <- cor.test( d.mni.img[d.mni.img$hemisphere == i & d.mni.img$axis == j,]$Andrej,
                     d.mni.img[d.mni.img$hemisphere == i & d.mni.img$axis == j,]$Josef )
    # fill-in
    t[i,j,"cor"] <- paste0(
      "r = ", sprintf("%.2f", round( cor$estimate, 2 )), ", (p ",
      ifelse( cor$p.value < .001, "< 0.001", paste0("= ", sprintf( "%.3f", round(cor$p.value, 3) ) ) ), ")"
    )
    # (paired) t-test
    # pre-compute
    tt <- t.test( d.mni.img[d.mni.img$hemisphere == i & d.mni.img$axis == j,]$Andrej,
                  d.mni.img[d.mni.img$hemisphere == i & d.mni.img$axis == j,]$Josef,
                  paired = T )
    # fill-in
    t[i,j,"t.test"] <- paste0(
      "dif. = ", sprintf("%.2f", round( tt$estimate, 2 )), ", (95% CI [",
      sprintf( "%.2f", round(tt$conf.int[[1]], 2) ), ", ", sprintf( "%.2f", round(tt$conf.int[[2]], 2) ), "])"
    )
  }
}

# save as tables
write.table( t(t[,,"cor"]), "tables/pearson_correlations.csv", sep = ";", col.names = c("left","right"), row.names = T )
write.table( t(t[,,"t.test"]), "tables/paired_t-tests.csv", sep = ";", col.names = c("left","right"), row.names = T )
