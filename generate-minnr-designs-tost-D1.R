## Generating min(n,r) designs for Tost, Mongolia, under spatially-varying detection and density
## These are the designs in Figure 3 of the paper.

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
### Figure 2: non-uniform D, uniform habitat use
################################################

opt_traps_D1_all <- data.frame(nT = as.integer(), sigma = as.numeric(), lambda0 = as.numeric(), b_ac = as.numeric(), b_det = as.numeric(), buffer = as.numeric(), dt = as.character(),
                               x = as.numeric(), y = as.numeric(), trap_id = as.integer())

all_masks_D1 <- list()
mask_cnt <- 1

buffer <- 0

for(ba in c(-1,1,3)){
  
  # D = population density animals / hectare; may be scalar or vector of length nrow(mask)
  # convert user specified dens_per_100km2 into dens per cell on the mask
  D_per_mask_cell <- dens_per_100km2 / 100 * (attr(mask, "area") / 100)
  # make density on the mask a function of some covariate (b_ac = 0 for constant density)
  covariates(mask)$Dac <- exp(ba * as.numeric(scale(covariates(mask)$stdGC)))
  # now standardize the Dac so that its mean is = D_per_mask_cell, or equivalently that its sum is = Dpmc * N mask cells
  Dcov_for_sim <- covariates(mask)$Dac / sum(covariates(mask)$Dac) * (D_per_mask_cell * length(mask$x))
  # turn this into equivalent density per hectate, for use with Enrm (secr functions want per ha values)
  Dcov_for_sim_ha <- Dcov_for_sim / attr(mask, "area")
  
  for(nT in c(40,60)){
    
    # store the mask used
    all_masks_D1[[mask_cnt]] <- mask
    mask_cnt <- mask_cnt + 1
    
    mnr <- scrdesignGAenrm(statespace = mask,
                           alltraps = alltraps,
                           ntraps = nT,
                           lambda0 = lambda0,
                           sigma = sigma,
                           D = Dcov_for_sim_ha,
                           occasions = 1,
                           detector = "count",
                           ngen = 50,
                           popsize = 1000,
                           crit = 3,
                           pen_wt = 100,
                           seed = 700)
    
    # extract camera locations
    mnr_traps <- as.data.frame(mnr$optimaltraps) %>% mutate(trap_id = 1)
    opt_traps_D1_all <- rbind(opt_traps_D1_all, cbind(nT = nT, lambda0 = lambda0, sigma = sigma, b_ac = ba, b_det = 0, buffer = buffer, dt = dt, mnr_traps))
    
  }
  
}

### detection function covariates

for(ba in c(0,1)){
  
  # D = population density animals / hectare; may be scalar or vector of length nrow(mask)
  # convert user specified dens_per_100km2 into dens per cell on the mask
  D_per_mask_cell <- dens_per_100km2 / 100 * (attr(mask, "area") / 100)
  # make density on the mask a function of some covariate (b_ac = 0 for constant density)
  covariates(mask)$Dac <- exp(ba * as.numeric(scale(covariates(mask)$stdGC)))
  # now standardize the Dac so that its mean is = D_per_mask_cell, or equivalently that its sum is = Dpmc * N mask cells
  Dcov_for_sim <- covariates(mask)$Dac / sum(covariates(mask)$Dac) * (D_per_mask_cell * length(mask$x))
  # turn this into equivalent density per hectate, for use with Enrm (secr functions want per ha values)
  Dcov_for_sim_ha <- Dcov_for_sim / attr(mask, "area")
  
  for(bd in c(-0.75,0.75,1.5)){
    
    # detection function parameters
    db1 <- bd
    dLam <- exp(as.numeric(db1 * scale(-alltraps[,1])))
    dLam <- dLam / sum(dLam) * (lambda0 * nrow(alltraps)) 
    
    for(nT in c(40,60)){
      
      # store the mask used
      all_masks_D1[[mask_cnt]] <- mask
      mask_cnt <- mask_cnt + 1
      
      mnr <- scrdesignGAenrm(statespace = mask,
                             alltraps = alltraps,
                             ntraps = nT,
                             lambda0 = dLam,
                             sigma = sigma,
                             D = Dcov_for_sim_ha,
                             occasions = 1,
                             detector = "count",
                             ngen = 50,
                             popsize = 1000,
                             crit = 3,
                             pen_wt = 100,
                             seed = 700)
      
      # extract camera locations
      mnr_traps <- as.data.frame(mnr$optimaltraps) %>% mutate(trap_id = 1)
      opt_traps_D1_all <- rbind(opt_traps_D1_all, cbind(nT = nT, lambda0 = lambda0, sigma = sigma, b_ac = ba, b_det = bd, buffer = buffer, dt = dt, mnr_traps))
      
    }
  }
}

save(all_masks_D1, alltraps_df, opt_traps_D1_all, file = "output/new/Tost_examples_nonuniD_new.Rdata")