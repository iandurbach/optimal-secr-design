## Generating regular grid and optimal (spacing) grid designs for Tost, Mongolia, under uniform density
## To be compared with designs in Figure 2 of the paper.

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(stringr)
library(purrr)
library(secr)
library(secrdesign)
library(oSCR)
library(raster)
library(gdistance)
library(kofnGA)

# user parameters
lambda0 <- 1  
dens_per_100km2 <- 2 # mean animal density per 100km2, SLs are ~1
D <- dens_per_100km2 / 10000
sigma <- 3000
b_ac <- 1

nT <- 60
buffer <- 4 * sigma
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

# create a buffered mask (the mask used for the designs) from the mask above
buffer_mult <- ceiling(buffer / cellsize)
newmask_df <- expand.grid(x = seq(from = min(mask_df$x) - buffer_mult * cellsize, 
                                  to = max(mask_df$x) + buffer_mult * cellsize, 
                                  by = cellsize),
                          y = seq(from = min(mask_df$y) - buffer_mult * cellsize, 
                                  to = max(mask_df$y) + buffer_mult * cellsize, 
                                  by = cellsize)) 
newmask_df <- newmask_df %>% 
  mutate(keep = (apply(e2dist(as.matrix(newmask_df), as.matrix(mask_df)), 1, min) <= buffer)) %>%
  filter(keep == TRUE) %>% dplyr::select(-keep)
newmask <- read.mask(data = newmask_df)
plot(newmask,axes=T)

# create grid of possible trap locations

cellsize <- 2/3 * sigma
red_factor <- cellsize[1] / attr(fullmask, "spacing")
if ((trunc(red_factor) - red_factor) != 0) stop("Check spacing, causing non-integer reduction factor for mesh")
alltraps <- secr::raster(fullmask, "stdGC") %>% raster::aggregate(fact = red_factor, fun = mean) 
alltraps_df <- data.frame(coordinates(alltraps), stdGC = matrix(alltraps)) %>% filter(!is.na(stdGC))
alltraps <- as.matrix(alltraps_df)[,c(1,2)]
plot(read.mask(data = alltraps_df), col = 'blue', add = TRUE, axes =T)              

### going to generate En, Er, CV under different designs in this survey area, varying nT, lambda0, sigma, buffer

################################################
### To compare with min(n,r) designs in Figure 2: uniform D, uniform habitat use
################################################

grid_traps_all <- data.frame(nT = as.integer(), sigma = as.numeric(), beta0 = as.numeric(), buffer = as.numeric(), dt = as.character(),
                             x = as.numeric(), y = as.numeric(), trap_id = as.integer())
opt_grid_traps_all <- data.frame(nT = as.integer(), sigma = as.numeric(), beta0 = as.numeric(), buffer = as.numeric(), dt = as.character(),
                                 x = as.numeric(), y = as.numeric(), trap_id = as.integer())

# Fig 1 a-c, d-f: varying number of cameras

for(n in c(20,40,60)){
  
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
  
  grid_traps <- grid_traps_full
  set.seed(700)
  xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
  grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
  grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 1)
  
  grid_traps_all <- rbind(grid_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, grid_traps))
  
  
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
    opt_grid_traps <- opt_grid_traps_full
    set.seed(700)
    xy_rand <- opt_grid_traps[sample(1:nrow(opt_grid_traps), 1), ]
    opt_grid_traps$dist2pt <- (opt_grid_traps$x - xy_rand$x)^2 + (opt_grid_traps$y - xy_rand$y)^2
    opt_grid_traps <- opt_grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 1)
    
    n_opt_grid <- nrow(opt_grid_traps)
    
  }
  
  opt_grid_traps_all <- rbind(opt_grid_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, opt_grid_traps))
  
}

for(n in c(20,40,60)){
  
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
  
  grid_traps <- grid_traps_full
  set.seed(222)
  xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
  grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
  grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 2)
  
  grid_traps_all <- rbind(grid_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, grid_traps))
  
  
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
    opt_grid_traps <- opt_grid_traps_full
    set.seed(222)
    xy_rand <- opt_grid_traps[sample(1:nrow(opt_grid_traps), 1), ]
    opt_grid_traps$dist2pt <- (opt_grid_traps$x - xy_rand$x)^2 + (opt_grid_traps$y - xy_rand$y)^2
    opt_grid_traps <- opt_grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 2)
    
    n_opt_grid <- nrow(opt_grid_traps)
    
  }
  
  opt_grid_traps_all <- rbind(opt_grid_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, opt_grid_traps))
  
}

# Fig 1 g-i: with zero buffer

buffer <- 0 # doesn't change anything, just says what pars were used to get these results

for(n in c(20, 40, 60)){
  
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
  
  grid_traps <- grid_traps_full
  set.seed(700)
  xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
  grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
  grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 1)
  
  grid_traps_all <- rbind(grid_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, grid_traps))
  
  
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
    opt_grid_traps <- opt_grid_traps_full
    set.seed(700)
    xy_rand <- opt_grid_traps[sample(1:nrow(opt_grid_traps), 1), ]
    opt_grid_traps$dist2pt <- (opt_grid_traps$x - xy_rand$x)^2 + (opt_grid_traps$y - xy_rand$y)^2
    opt_grid_traps <- opt_grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 1)
    
    n_opt_grid <- nrow(opt_grid_traps)
    
  }
  
  opt_grid_traps_all <- rbind(opt_grid_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, opt_grid_traps))
  
}

# Fig 1 j-l: vary beta0/lambda0

nT <- 60
buffer <- 4 * sigma  # doesn't change anything, just says what pars were used to get these results

for(lam in c(0.5*lambda0, 1.5*lambda0, 2*lambda0)){
  
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
  
  grid_traps <- grid_traps_full
  set.seed(700)
  xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
  grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
  grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= nT) %>% mutate(trap_id = 1)
  
  grid_traps_all <- rbind(grid_traps_all, cbind(nT = nT, lambda0 = lam, sigma = sigma, buffer = buffer, dt = dt, grid_traps))
  
  
  ###############
  ### regular grid with optimal spacing
  ###############
  
  grid_traps <- read.traps(data = grid_traps, detector = dt)
  
  # compute optimal spacing for previous grid 
  dd <- optimalSpacing(D = D, traps = grid_traps, detectpar = list(lambda0 = lam, sigma = sigma), noccasions = 1)$rotRSE$optimum.spacing
  
  # sometimes optimal spacing is too big to fit desired number of traps in; if so, reduce dd until it is
  n_opt_grid <- 0
  red_dd <- 0.9
  while(n_opt_grid < nT){
    
    dd <- dd * red_dd
    
    # place a grid over the area, with cells optimally spaced
    my_grid <- st_make_grid(sap_1, cellsize = c(dd, dd), what = "centers") %>% st_intersection(sap_1)
    opt_grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
    
    # choose a subset of nT of these points, starting from a random trap and then choosing nT nearest neighbours
    opt_grid_traps <- opt_grid_traps_full
    set.seed(700)
    xy_rand <- opt_grid_traps[sample(1:nrow(opt_grid_traps), 1), ]
    opt_grid_traps$dist2pt <- (opt_grid_traps$x - xy_rand$x)^2 + (opt_grid_traps$y - xy_rand$y)^2
    opt_grid_traps <- opt_grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= nT) %>% mutate(trap_id = 1)
    
    n_opt_grid <- nrow(opt_grid_traps)

  }
  
  opt_grid_traps_all <- rbind(opt_grid_traps_all, cbind(nT = nT, lambda0 = lam, sigma = sigma, buffer = buffer, dt = dt, opt_grid_traps))
  
}

# Fig 1 m-o: vary sigma (buffer is always 4 * sigma)

for(s in c(0.5*sigma, 1.5*sigma, 2*sigma)){
  
  buffer = 4 * s 
  
  ##########################
  ### Regular 2 sigma grid
  ##########################
  
  # hacky way to make a polygon just bigger than boundary of mask points
  # used later to decide which randomly generated detectors are on the mask and so allowed
  sap_1 <- mask_df %>% st_as_sf(coords = c("x", "y")) %>% st_buffer(dist = cellsize, endCapStyle = "SQUARE") %>% st_union() 
  sap_1 <- sap_1 %>% st_buffer(dist = -cellsize * 0.4, endCapStyle = "SQUARE")
  
  # 60 detectors only fit on if sigma < 3000, so reduce sigma if necessary
  sm <- min(s, 3000)
  
  # place a grid over the area, with cells X sigma apart
  my_grid <- st_make_grid(sap_1, cellsize = c(2 * sm, 2 * sm), what = "centers") %>%
    st_intersection(sap_1)
  grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
  
  # choose a subset of nT of these, starting from a random trap and then choosing nT nearest neighbours
  all_grid_traps <- list()
  
  grid_traps <- grid_traps_full
  set.seed(700)
  xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
  grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
  grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= nT) %>% mutate(trap_id = 1)
  
  grid_traps_all <- rbind(grid_traps_all, cbind(nT = nT, lambda0 = lambda0, sigma = s, buffer = buffer, dt = dt, grid_traps))
  
  ###############
  ### regular grid with optimal spacing
  ###############
  
  grid_traps <- read.traps(data = grid_traps, detector = dt)
  
  # compute optimal spacing for previous grid 
  dd <- optimalSpacing(D = D, traps = grid_traps, detectpar = list(lambda0 = lambda0, sigma = s), noccasions = 1)$rotRSE$optimum.spacing
  
  # sometimes optimal spacing is too big to fit desired number of traps in; if so, reduce dd until it is
  n_opt_grid <- 0
  red_dd <- 0.9
  while(n_opt_grid < nT){
    
    dd <- dd * red_dd
    
    # place a grid over the area, with cells optimally spaced
    my_grid <- st_make_grid(sap_1, cellsize = c(dd, dd), what = "centers") %>% st_intersection(sap_1)
    opt_grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
    
    # choose a subset of nT of these points, starting from a random trap and then choosing nT nearest neighbours
    opt_grid_traps <- opt_grid_traps_full
    set.seed(700)
    xy_rand <- opt_grid_traps[sample(1:nrow(opt_grid_traps), 1), ]
    opt_grid_traps$dist2pt <- (opt_grid_traps$x - xy_rand$x)^2 + (opt_grid_traps$y - xy_rand$y)^2
    opt_grid_traps <- opt_grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= nT) %>% mutate(trap_id = 1)
    
    n_opt_grid <- nrow(opt_grid_traps)
    
  }
  
  opt_grid_traps_all <- rbind(opt_grid_traps_all, cbind(nT = nT, lambda0 = lambda0, sigma = s, buffer = buffer, dt = dt, opt_grid_traps))
  
}

# Fig 1 p-r: other kinds of detectors 
# note: these use different values for sigma, lambda0 to show differences between designs more clearly

for(d in c("multi", "proximity", "count")){
  
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
  
  grid_traps <- grid_traps_full
  set.seed(700)
  xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
  grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
  grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 1)
  
  grid_traps_all <- rbind(grid_traps_all, cbind(nT = nT, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = d, grid_traps))
  
  
  ###############
  ### regular grid with optimal spacing
  ###############
  
  grid_traps <- read.traps(data = grid_traps, detector = d)
  
  # compute optimal spacing for previous grid 
  dd <- optimalSpacing(D = D, traps = grid_traps, detectpar = list(lambda0 = lambda0/5, sigma = sigma), noccasions = 5)$rotRSE$optimum.spacing
  
  # sometimes optimal spacing is too big to fit desired number of traps in; if so, reduce dd until it is
  n_opt_grid <- 0
  red_dd <- 0.9
  while(n_opt_grid < nT){
    
    dd <- dd * red_dd
    
    # place a grid over the area, with cells optimally spaced
    my_grid <- st_make_grid(sap_1, cellsize = c(dd, dd), what = "centers") %>% st_intersection(sap_1)
    opt_grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
    
    # choose a subset of nT of these points, starting from a random trap and then choosing nT nearest neighbours
    opt_grid_traps <- opt_grid_traps_full
    set.seed(700)
    xy_rand <- opt_grid_traps[sample(1:nrow(opt_grid_traps), 1), ]
    opt_grid_traps$dist2pt <- (opt_grid_traps$x - xy_rand$x)^2 + (opt_grid_traps$y - xy_rand$y)^2
    opt_grid_traps <- opt_grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= nT) %>% mutate(trap_id = 1)
    
    n_opt_grid <- nrow(opt_grid_traps)
    
  }
  
  opt_grid_traps_all <- rbind(opt_grid_traps_all, cbind(nT = nT, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = d, opt_grid_traps))
  
}

save(mask_df, grid_traps_all, opt_grid_traps_all, file = "output/Tost_examples_nonopt_D0_new.Rdata")
