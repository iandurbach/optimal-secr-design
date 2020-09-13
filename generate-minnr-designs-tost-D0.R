## Generating min(n,r) designs for Tost, Mongolia, under uniform density
## These are the designs in Figure 2 of the paper.

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(stringr)
library(secr)
library(secrdesign)
library(oSCR)
library(raster)
library(kofnGA)

source("oSCR/scrdesignGAenrm.R")
source("oSCR/scrdesignOFenrm.R")

# https://www.r-bloggers.com/three-ways-to-call-cc-from-r/
source("oSCR/LambdaL.R")
source("oSCR/utils_for_enrmL.R")
dyn.load("oSCR/mysecrdesign.so")

# user parameters
lambda0 <- 1  # beta0 = log(lambda0) = log(K * 'p0')
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

### going to generate En, Er, CV under different designs in this survey area, varying nT, beta0, sigma, buffer

################################################
### Results for Figure 2: uniform D, uniform habitat use
################################################

opt_traps_all <- data.frame(nT = as.integer(), sigma = as.numeric(), lambda0 = as.numeric(), buffer = as.numeric(), dt = as.character(),
                            x = as.numeric(), y = as.numeric(), trap_id = as.integer())

all_masks <- list()
mask_cnt <- 1

# Fig 1 a-c, d-f: varying number of cameras

for(n in c(20,40,60)){
  
  # store the mask used
  all_masks[[mask_cnt]] <- newmask
  mask_cnt <- mask_cnt + 1
  
  mnr <- scrdesignGAenrm(statespace = newmask,
                         alltraps = alltraps,
                         ntraps = n,
                         beta0 = log(lambda0),
                         sigma = sigma,
                         D = D,
                         occasions = 1,
                         detector = "count",
                         ngen = 50,
                         popsize = 1000,
                         crit = 3,
                         seed = 700)
  
  # extract camera locations
  mnr_traps <- as.data.frame(mnr$optimaltraps) %>% mutate(trap_id = 1)
  opt_traps_all <- rbind(opt_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, mnr_traps))
  
}

for(n in c(20,40,60)){
  
  # store the mask used
  all_masks[[mask_cnt]] <- newmask
  mask_cnt <- mask_cnt + 1
  
  mnr <- scrdesignGAenrm(statespace = newmask,
                         alltraps = alltraps,
                         ntraps = n,
                         beta0 = log(lambda0),
                         sigma = sigma,
                         D = D,
                         occasions = 1,
                         detector = "count",
                         ngen = 50,
                         popsize = 1000,
                         crit = 3,
                         seed = 222)
  
  # extract camera locations
  mnr_traps <- as.data.frame(mnr$optimaltraps) %>% mutate(trap_id = 2)
  opt_traps_all <- rbind(opt_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, mnr_traps))
  
}

# Fig 1 g-i: with zero buffer

buffer <- 0 # doesn't change anything, just says what pars were used to get these results

for(n in c(20, 40, 60)){
  
  # store the mask used
  all_masks[[mask_cnt]] <- mask
  mask_cnt <- mask_cnt + 1
  
  mnr <- scrdesignGAenrm(statespace = mask,
                         alltraps = alltraps,
                         ntraps = n,
                         beta0 = log(lambda0),
                         sigma = sigma,
                         D = D,
                         occasions = 1,
                         detector = "count",
                         ngen = 50,
                         popsize = 1000,
                         crit = 3,
                         seed = 222)
  
  # extract camera locations
  mnr_traps <- as.data.frame(mnr$optimaltraps) %>% mutate(trap_id = 1)
  opt_traps_all <- rbind(opt_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, mnr_traps))
  
}

# Fig 1 j-l: vary beta0/lambda0

nT <- 60
buffer <- 4 * sigma  # doesn't change anything, just says what pars were used to get these results

for(lam in c(0.5*lambda0, 1.5*lambda0, 2*lambda0)){
  
  # store the mask used
  all_masks[[mask_cnt]] <- newmask
  mask_cnt <- mask_cnt + 1
  
  mnr <- scrdesignGAenrm(statespace = newmask,
                         alltraps = alltraps,
                         ntraps = nT,
                         beta0 = log(lam),
                         sigma = sigma,
                         D = D,
                         occasions = 1,
                         detector = "count",
                         ngen = 50,
                         popsize = 1000,
                         crit = 3,
                         seed = 222)
  
  # extract camera locations
  mnr_traps <- as.data.frame(mnr$optimaltraps) %>% mutate(trap_id = 1)
  opt_traps_all <- rbind(opt_traps_all, cbind(nT = nT, lambda0 = lam, sigma = sigma, buffer = buffer, dt = dt, mnr_traps))
  
}

# Fig 1 m-o: vary sigma (buffer is always 4 * sigma)

for(s in c(0.5*sigma, 1.5*sigma, 2*sigma)){
  
  buffer <- 4 * s
  
  # create statespace and traps objects
  buffer_mult <- ceiling(buffer / cellsize)
  newmask_df_s <- expand.grid(x = seq(from = min(mask_df$x) - buffer_mult * cellsize, 
                                      to = max(mask_df$x) + buffer_mult * cellsize, 
                                      by = cellsize),
                              y = seq(from = min(mask_df$y) - buffer_mult * cellsize, 
                                      to = max(mask_df$y) + buffer_mult * cellsize, 
                                      by = cellsize)) 
  newmask_df_s <- newmask_df_s %>% 
    mutate(keep = (apply(e2dist(as.matrix(newmask_df_s), as.matrix(mask_df)), 1, min) <= buffer)) %>%
    filter(keep == TRUE) %>% dplyr::select(-keep)
  
  mask_s <- read.mask(data = newmask_df_s)
  
  # store the mask used
  all_masks[[mask_cnt]] <- mask_s
  mask_cnt <- mask_cnt + 1
  
  mnr <- scrdesignGAenrm(statespace = mask_s,
                         alltraps = alltraps,
                         ntraps = nT,
                         beta0 = log(lambda0),
                         sigma = s,
                         D = D,
                         occasions = 1,
                         detector = "count",
                         ngen = 50,
                         popsize = 1000,
                         crit = 3,
                         seed = 222)
  
  # extract camera locations
  mnr_traps <- as.data.frame(mnr$optimaltraps) %>% mutate(trap_id = 1)
  opt_traps_all <- rbind(opt_traps_all, cbind(nT = nT, lambda0 = lambda0, sigma = s, buffer = buffer, dt = dt, mnr_traps))
  
}

# Fig 1 p-r: other kinds of detectors 
# note: these use different values for sigma, beta0 to show differences between designs more clearly

for(d in c("multi", "proximity", "count")){
  
  # store the mask used
  all_masks[[mask_cnt]] <- newmask
  mask_cnt <- mask_cnt + 1
  
  mnr <- scrdesignGAenrm(statespace = newmask,
                         alltraps = alltraps,
                         ntraps = nT,
                         beta0 = log(lambda0/5),
                         sigma = sigma,
                         D = D,
                         occasions = 5,
                         detector = d,
                         ngen = 50,
                         popsize = 1000,
                         crit = 3,
                         seed = 222)
  
  # extract camera locations
  mnr_traps <- as.data.frame(mnr$optimaltraps) %>% mutate(trap_id = 1)
  opt_traps_all <- rbind(opt_traps_all, cbind(nT = nT, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = d, mnr_traps))
  
}

save(all_masks, alltraps_df, mask_df, maskbuffer_df, newmask_df, opt_traps_all, file = "output/new/Tost_examples_all_D0_new.Rdata")