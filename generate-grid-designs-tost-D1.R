## Generating regular grid and optimal (spacing) grid designs for Tost, Mongolia, under spatially-varying detection and density
## To be compared with designs in Figure 3 of the paper.

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(secr)
library(secrdesign)
library(raster)
library(gdistance)
library(oSCR)

# user parameters
lambda0 <- 1  
dens_per_100km2 <- 2 # mean animal density per 100km2, SLs are ~1
D <- dens_per_100km2 / 10000
sigma <- 3000

buffer <- 0
dt <- "count"

# load secr masks -- note that these have been constructed to have buffer = 0 (hard boundary to study area)
load("data/Tost.RData")
fullmask <- TostMask 

# reduce resolution of mesh so have fewer possible camera locations
cellsize <- 2/3 * sigma
red_factor <- cellsize[1] / attr(fullmask, "spacing")
if ((trunc(red_factor) - red_factor) != 0) stop("Check spacing, causing non-integer reduction factor for mesh")
mask <- secr::raster(fullmask, "stdGC") %>% raster::aggregate(fact = red_factor, fun = mean) 
mask_df <- data.frame(coordinates(mask), stdGC = matrix(mask)) %>% filter(!is.na(stdGC))
mask <- read.mask(data = mask_df)
plot(mask)

# going to create grids that are centred around either highest D or highest lambda0 cell.
# D ~ stdGC, lambda0 ~ -x, but either relationship (coeff) can be pos or neg
# so, get the cells with the highest or lowest values of stdGC and x

# mask point with highest covariate value (highest D when b_ac > 0)
highD_posbac <- mask_df[which.max(mask_df$stdGC), ]
# mask point with lowest covariate value (highest D when b_ac < 0)
highD_negbac <- mask_df[which.min(mask_df$stdGC), ]
# mask point with highest covariate value (highest lambda0 when b_d > 0, remembering lam ~ -x)
highlam0_posbd <- mask_df[which.min(mask_df$x), ]
# mask point with lowest covariate value (highest lambda when b_d < 0)
highlam0_negbd <- mask_df[which.max(mask_df$x), ]

# set up parameter values in same order used to generate optimal designs, so can associate a regular
# grid to each optimized design

b_ac <- c(-1, -1, 1, 1, 3, 3, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
b_d <- c(0, 0, 0, 0, 0, 0, -0.75, -0.75, 0.75, 0.75, 1.5, 1.5, -0.75, -0.75, 0.75, 0.75, 1.5, 1.5)
nT <- rep(c(40,60), times = 9)

################################################
### Regular grids centred on 'best' (highest D, highest lambda0) spot
################################################

grid_traps_all <- data.frame(nT = as.integer(), sigma = as.numeric(), lambda0 = as.numeric(), b_ac = as.numeric(), b_det = as.numeric(), buffer = as.numeric(), dt = as.character(),
                             x = as.numeric(), y = as.numeric(), trap_id = as.integer())

opt_grid_traps_all <- data.frame(nT = as.integer(), sigma = as.numeric(), lambda0 = as.numeric(), b_ac = as.numeric(), b_det = as.numeric(), buffer = as.numeric(), dt = as.character(),
                                 x = as.numeric(), y = as.numeric(), trap_id = as.integer())

for(i in 1:18){
  
  bac <- b_ac[i]
  bd <- b_d[i]
  n <- nT[i]
  
  # choose the appropriate mask cell to center the grid at
  if(bac > 0){
    grid_center <- highD_posbac
  } else if (bac < 0) {
      grid_center <- highD_negbac 
  } else if (bd > 0) {
    grid_center <- highlam0_posbd
  } else if (bd < 0) {
    grid_center <- highlam0_negbd 
  } else stop("Check b_ac or b_det conditions")
  
  ##########################
  ### Regular 2 sigma grid
  ##########################
  
  # hacky way to make a polygon just bigger than boundary of mask points
  # used later to decide which randomly generated detectors are on the mask and so allowed
  sap_1 <- mask_df %>% st_as_sf(coords = c("x", "y")) %>% st_buffer(dist = cellsize, endCapStyle = "SQUARE") %>% st_union() 
  sap_1 <- sap_1 %>% st_buffer(dist = -cellsize * 0.4, endCapStyle = "SQUARE")
  
  # place a grid over the area, with cells X sigma apart
  my_grid <- st_make_grid(sap_1, cellsize = c(2 * sigma, 2 * sigma), what = "centers") %>%
    st_intersection(sap_1)
  grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
  
  # choose a subset of nT of these, starting from a random trap and then choosing nT nearest neighbours
  all_grid_traps <- list()
  
  grid_traps <- grid_traps_full %>% mutate(d2c = e2dist(grid_traps_full, grid_center[c(1,2)]))
  xy_rand <- grid_traps[which.min(grid_traps$d2c), ]
  grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
  grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 1)
  grid_traps <- grid_traps %>% dplyr::select(-d2c, -dist2pt)
  
  grid_traps_all <- rbind(grid_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, b_ac = bac, b_det = bd, buffer = 0, dt = dt, grid_traps))
  
  #plot(mask)
  #plot(read.traps(data = grid_traps, detector = "count"), add = TRUE)
  
  ###############
  ### regular grid with optimal spacing
  ###############
  
  grid_traps <- read.traps(data = grid_traps, detector = dt)
  
  # compute optimal spacing for previous grid 
  dd <- optimalSpacing(D = D, traps = grid_traps, detectpar = list(lambda0 = lambda0, sigma = sigma), noccasions = 1)$rotRSE$optimum.spacing
  
  # sometimes optimal spacing is too big to fit desired number of traps in; if so, reduce dd until it is
  n_opt_grid <- 0
  red_dd <- 0.9
  while(n_opt_grid < n){
    
    dd <- dd * red_dd
    
    # place a grid over the area, with cells optimally spaced
    my_grid <- st_make_grid(sap_1, cellsize = c(dd, dd), what = "centers") %>% st_intersection(sap_1)
    opt_grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
    
    # choose a subset of nT of these points, starting from a random trap and then choosing nT nearest neighbours
    opt_grid_traps <- opt_grid_traps_full %>% mutate(d2c = e2dist(opt_grid_traps_full, grid_center[c(1,2)]))
    xy_rand <- opt_grid_traps[which.min(opt_grid_traps$d2c), ]
    opt_grid_traps$dist2pt <- (opt_grid_traps$x - xy_rand$x)^2 + (opt_grid_traps$y - xy_rand$y)^2
    opt_grid_traps <- opt_grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 1)
    opt_grid_traps <- opt_grid_traps %>% dplyr::select(-d2c, -dist2pt)
    
    n_opt_grid <- nrow(opt_grid_traps)
    
  }
  
  #plot(mask)
  #plot(read.traps(data = opt_grid_traps, detector = "count"), add = TRUE)
  
  opt_grid_traps_all <- rbind(opt_grid_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, b_ac = bac, b_det = bd, buffer = 0, dt = dt, opt_grid_traps))

  }

save(mask_df, grid_traps_all, opt_grid_traps_all, file = "output/new/Tost_examples_nonopt_D1_new.Rdata")
