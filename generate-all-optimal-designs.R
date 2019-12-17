## 

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

source("oSCR/SCRdesign_mcvd.R")
source("oSCR/Qfn_mcvd.R")

# user parameters
cellsize <- 2000
sigma <- 6000
buffer <- 3 * sigma
ndesigns <- 2
dens_per_100km2 <- 1 # mean animal density per 100km2, SLs are ~1
dt <- "count"

# load secr masks -- note that these have been constructed to have buffer = 0 (hard boundary to study area)
load("data/Tost.RData")
mask <- TostMask 

# reduce resolution of mesh so have fewer possible camera locations
red_factor <- cellsize[1] / attr(mask, "spacing")
if ((trunc(red_factor) - red_factor) != 0) stop("Check spacing, causing non-integer reduction factor for mesh")
mask <- secr::raster(mask, "stdGC") %>% raster::aggregate(fact = red_factor, fun = mean) 
mask_df <- data.frame(coordinates(mask), stdGC = matrix(mask)) %>% filter(!is.na(stdGC))
mask <- read.mask(data = mask_df)

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

save(mask_df, file = "output/TostMask-for-plotting.RData")

### going to generate En, Er, CV under different designs in this survey area, varying nT, beta0, sigma, buffer

################################################
### Results for Figure 1: uniform D, uniform habitat use
################################################

optEnrm <- data.frame(nT = as.integer(), En = as.numeric(), Er = as.numeric(), Em = as.numeric())
opt_traps_all <- data.frame(nT = as.integer(), sigma = as.numeric(), beta0 = as.numeric(), buffer = as.numeric(), dt = as.character(),
                            x = as.numeric(), y = as.numeric(), trap_id = as.integer())

# Fig 1 a-c, d-f: varying number of cameras

beta0 <- -0.8  # beta0 = log(lambda0) = log(K * 'p0')

# statespace = all mask points, including any buffer region
statespace <- newmask_df[,1:2] %>% as.matrix()
# traps = subset of mask points that are potential camera locations
my_traps <- mask_df[,1:2] %>% as.matrix()

for(nT in c(20,40,60)){
  
  optSCR <- SCRdesign(statespace = statespace,
                      all.traps = my_traps,
                      ntraps = nT, # number of cameras available
                      ndesigns = ndesigns, # number of random starting points
                      beta0 = beta0,
                      sigma = sigma, # SCR pars (log(lambda0) and usual sigma)
                      D_per_mask_cell = dens_per_100km2 / 100 * (attr(mask, "area") / 100))
  
  opt_traps_nT <- optSCR$Xlst %>%
    purrr::reduce(rbind) %>% data.frame() %>%
    mutate(trap_id = rep(1:ndesigns, each = nT))
  
  opt_traps_all <- rbind(opt_traps_all, cbind(nT = nT, beta0 = beta0, sigma = sigma, buffer = buffer, dt = dt, opt_traps_nT))
  
}

# Fig 1 g-i: with zero buffer

beta0 <- -0.8
sigma <- 6000
buffer <- 0 # doesn't change anything, just says what pars were used to get these results

statespace <- mask_df[,1:2] %>% as.matrix()
my_traps <- mask_df[,1:2] %>% as.matrix()

for(nT in c(20, 40, 60)){
  
  optSCR <- SCRdesign(statespace = statespace,
                      all.traps = my_traps,
                      ntraps = nT, # number of cameras available
                      ndesigns = ndesigns, # number of random starting points
                      beta0 = beta0,
                      sigma = sigma, # SCR pars (log(lambda0) and usual sigma)
                      D_per_mask_cell = dens_per_100km2 / 100 * (attr(mask, "area") / 100))
  
  opt_traps_nT <- optSCR$Xlst %>%
    purrr::reduce(rbind) %>% data.frame() %>%
    mutate(trap_id = rep(1:ndesigns, each = nT))
  
  opt_traps_all <- rbind(opt_traps_all, cbind(nT = nT, beta0 = beta0, sigma = sigma, buffer = buffer, dt = dt, opt_traps_nT))
  
  
}

# Fig 1 j-l: vary beta0/lambda0

nT <- 40
buffer <- 3 * sigma  # doesn't change anything, just says what pars were used to get these results

statespace <- newmask_df[,1:2] %>% as.matrix()
my_traps <- mask_df[,1:2] %>% as.matrix()

for(beta0 in c(-0.113,-0.4,-1.5)){
  
  optSCR <- SCRdesign(statespace = statespace,
                      all.traps = my_traps,
                      ntraps = nT, # number of cameras available
                      ndesigns = ndesigns, # number of random starting points
                      beta0 = beta0,
                      sigma = sigma, # SCR pars (log(lambda0) and usual sigma)
                      D_per_mask_cell = dens_per_100km2 / 100 * (attr(mask, "area") / 100))
  
  opt_traps_nT <- optSCR$Xlst %>%
    purrr::reduce(rbind) %>% data.frame() %>%
    mutate(trap_id = rep(1:ndesigns, each = nT))
  
  opt_traps_all <- rbind(opt_traps_all, cbind(nT = nT, beta0 = beta0, sigma = sigma, buffer = buffer, dt = dt, opt_traps_nT))
  
  
}

# Fig 1 m-o: vary sigma (buffer is always 3 * sigma)

nT <- 40
beta0 <- -0.8

for(sigma in c(3000, 9000, 12000)){
  
  buffer <- 3 * sigma
  
  # create statespace and traps objects
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
  
  statespace <- newmask_df[,1:2] %>% as.matrix()
  my_traps <- mask_df[,1:2] %>% as.matrix()
  
  optSCR <- SCRdesign(statespace = statespace,
                      all.traps = my_traps,
                      ntraps = nT, # number of cameras available
                      ndesigns = ndesigns, # number of random starting points
                      beta0 = beta0,
                      sigma = sigma, # SCR pars (log(lambda0) and usual sigma)
                      D_per_mask_cell = dens_per_100km2 / 100 * (attr(mask, "area") / 100))
  
  opt_traps_nT <- optSCR$Xlst %>%
    purrr::reduce(rbind) %>% data.frame() %>%
    mutate(trap_id = rep(1:ndesigns, each = nT))
  
  opt_traps_all <- rbind(opt_traps_all, cbind(nT = nT, beta0 = beta0, sigma = sigma, buffer = buffer, dt = dt, opt_traps_nT))
  
}

# Fig 1 p-r: other kinds of detectors 
# note: these use different values for sigma, beta0 to show differences between designs more clearly

sigma <- 2000
buffer <- 3 * sigma
beta0 <- -0.1  # beta0 = log(lambda0) = log(K * 'p0')
nT <- 30

# create statespace and traps objects
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

statespace <- newmask_df[,1:2] %>% as.matrix()
my_traps <- mask_df[,1:2] %>% as.matrix()

for(dt in c("multi", "proximity", "count")){
  
  set.seed(123)
  optSCR <- SCRdesign(statespace = statespace,
                      all.traps = my_traps,
                      ntraps = nT, # number of cameras available
                      ndesigns = ndesigns, # number of random starting points
                      beta0 = beta0,
                      sigma = sigma, # SCR pars (log(lambda0) and usual sigma)
                      D_per_mask_cell = dens_per_100km2 / 100 * (attr(mask, "area") / 100),
                      detector = dt,
                      occasions = 5)
  
  opt_traps_nT <- optSCR$Xlst %>%
    purrr::reduce(rbind) %>% data.frame() %>%
    mutate(trap_id = rep(1:ndesigns, each = nT))
  
  opt_traps_all <- rbind(opt_traps_all, cbind(nT = nT, beta0 = beta0, sigma = sigma, buffer = buffer, dt = dt, opt_traps_nT))
  
}

save(mask_df, opt_traps_all, file = "output/Tost_examples_all_D0.Rdata")

################################################
### Figure 2: non-uniform D, uniform habitat use
################################################

optEnrm <- data.frame(nT = as.integer(), En = as.numeric(), Er = as.numeric(), Em = as.numeric())
opt_traps_all <- data.frame(nT = as.integer(), b_ac = as.numeric(), x = as.numeric(), y = as.numeric(), trap_id = as.integer())

statespace <- mask_df[,1:2] %>% as.matrix()
my_traps <- mask_df[,1:2] %>% as.matrix()

for(b_ac in c(-1,1,3)){
  
  # D = population density animals / hectare; may be scalar or vector of length nrow(mask)
  # convert user specified dens_per_100km2 into dens per cell on the mask
  D_per_mask_cell <- dens_per_100km2 / 100 * (attr(mask, "area") / 100)
  # make density on the mask a function of some covariate (b_ac = 0 for constant density)
  if (cov_name == "stdGC") {
    covariates(mask)$Dac <- exp(b_ac * as.numeric(scale(covariates(mask)$stdGC)))
  } else if (cov_name == "x"){
    covariates(mask)$Dac <- exp(b_ac * as.numeric(scale(mask$x)))
  } else {stop("Covariate choice not supported")}
  # now standardize the Dac so that its mean is = D_per_mask_cell, or equivalently that its sum is = Dpmc * N mask cells
  Dcov_for_sim <- covariates(mask)$Dac / sum(covariates(mask)$Dac) * (D_per_mask_cell * length(mask$x))
  # turn this into equivalent density per hectate, for use with Enrm (secr functions want per ha values)
  Dcov_for_sim_ha <- Dcov_for_sim / attr(mask, "area")
  
  for(nT in c(20,40,60)){
    
    # 4) with approximately optimal design
    optSCR <- SCRdesign(statespace = statespace,
                        all.traps = my_traps,
                        ntraps = nT, # number of cameras available
                        ndesigns = ndesigns, # number of random starting points
                        beta0 = beta0,
                        sigma = sigma, # SCR pars (log(lambda0) and usual sigma)
                        D_per_mask_cell = Dcov_for_sim)
    
    opt_traps_nT <- optSCR$Xlst %>%
      purrr::reduce(rbind) %>% data.frame() %>%
      rename(x = X, y = Y) %>%
      mutate(trap_id = rep(1:ndesigns, each = nT))
    
    for(i in 1:ndesigns){
      opt_traps_i <- opt_traps_nT %>% dplyr::filter(trap_id == i)
      opt_traps <- read.traps(data = opt_traps_i, detector = "count")
      optEnrm <- rbind(optEnrm,
                       c(nT, b_ac, Enrm(D = Dcov_for_sim_ha, traps = opt_traps, mask = mask, detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1)) %>% t() %>% data.frame())
    }
    
    opt_traps_all <- rbind(opt_traps_all, cbind(nT = nT, b_ac = b_ac, opt_traps_nT))
    
  }
  
}

save(mask_df, optEnrm, opt_traps_all, file = "output/Tost_examples_nonuniD.Rdata")

# mask_df %>% ggplot(aes(x = x, y = y)) + geom_tile(aes(fill = stdGC), colour = "black", size = 0.5) +
#   geom_point(data = opt_traps_all, inherit.aes = FALSE, aes(x = x, y = y), size = 2, colour = "red") +
#   facet_grid(nT ~ trap_id) + scale_fill_viridis_c() + coord_equal() + theme_bw()

################################################
### Figure 4: uniform D, non-uniform habitat use
################################################

optEnrm <- data.frame(nT = as.integer(), En = as.numeric(), Er = as.numeric(), Em = as.numeric())
opt_traps_all <- data.frame(nT = as.integer(), alpha2 = as.numeric(), x = as.numeric(), y = as.numeric(), trap_id = as.integer())

for(nT in c(20,40,60)){
  for(alpha2 in c(-1, 1, 3)){
    
    set.seed(123)
    optSCR <- SCRdesign_biObj(statespace = statespace,
                              all.traps = my_traps,
                              ntraps = nT, # number of cameras available
                              ndesigns = ndesigns, # number of random starting points
                              beta0 = beta0,
                              sigma = sigma, # SCR pars (log(lambda0) and usual sigma)
                              D_per_mask_cell = dens_per_100km2 / 100 * (attr(mask, "area") / 100),
                              noneuc_costs = covariates(mask)$stdGC,
                              alpha2 = alpha2)
    
    opt_traps_nT <- optSCR$Xlst %>%
      purrr::reduce(rbind) %>% data.frame() %>%
      rename(x = X, y = Y) %>%
      mutate(trap_id = rep(1:ndesigns, each = nT))
    
    for(i in 1:ndesigns){
      opt_traps_i <- opt_traps_nT %>% dplyr::filter(trap_id == i)
      opt_traps <- read.traps(data = opt_traps_i, detector = "count")
      optEnrm <- rbind(optEnrm,
                       c(nT, alpha2, Enrm(D = Dcov_for_sim_ha, traps = opt_traps, mask = mask, detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1)) %>% t() %>% data.frame())
    }
    
    opt_traps_all <- rbind(opt_traps_all, cbind(nT = nT, alpha2 = alpha2, opt_traps_nT))
    
  }
}


save(mask_df, optEnrm, opt_traps_all, file = "output/Tost_examples_noneuc.Rdata")

# mask_df %>% ggplot(aes(x = x, y = y)) + geom_tile(aes(fill = stdGC), colour = "black", size = 0.5) +
#   geom_point(data = opt_traps_nT, inherit.aes = FALSE, aes(x = x, y = y), size = 2, colour = "red") +
#   scale_fill_viridis_c() + coord_equal() + theme_bw()

########################
### Table 1: same as above, just running different combinations of b_ac, nT, beta0
### and recording En, Er for OPTIMAL and OTHER (GRID, ACTUAL SURVEY, GRID + OPT.SP) designs
########################

## first for optimal designs

optEnrm <- data.frame(nT = as.integer(), b_ac = as.numeric(), beta0 = as.numeric(), En = as.numeric(), Er = as.numeric(), Em = as.numeric())
opt_traps_all <- data.frame(nT = as.integer(), b_ac = as.numeric(), beta0 = as.numeric(), x = as.numeric(), y = as.numeric(), trap_id = as.integer())

statespace <- mask_df[,1:2] %>% as.matrix()
my_traps <- mask_df[,1:2] %>% as.matrix()

for(b_ac in c(0, 1)){
  
  # D = population density animals / hectare; may be scalar or vector of length nrow(mask)
  # convert user specified dens_per_100km2 into dens per cell on the mask
  D_per_mask_cell <- dens_per_100km2 / 100 * (attr(mask, "area") / 100)
  # make density on the mask a function of some covariate (b_ac = 0 for constant density)
  if (cov_name == "stdGC") {
    covariates(mask)$Dac <- exp(b_ac * as.numeric(scale(covariates(mask)$stdGC)))
  } else if (cov_name == "x"){
    covariates(mask)$Dac <- exp(b_ac * as.numeric(scale(mask$x)))
  } else {stop("Covariate choice not supported")}
  # now standardize the Dac so that its mean is = D_per_mask_cell, or equivalently that its sum is = Dpmc * N mask cells
  Dcov_for_sim <- covariates(mask)$Dac / sum(covariates(mask)$Dac) * (D_per_mask_cell * length(mask$x))
  # turn this into equivalent density per hectate, for use with Enrm (secr functions want per ha values)
  Dcov_for_sim_ha <- Dcov_for_sim / attr(mask, "area")
  
  for(beta0 in c(-1.5, -0.8, -0.4)){
    
    for(nT in c(20,30,40,50,60)){
      
      # 4) with approximately optimal design
      optSCR <- SCRdesign(statespace = statespace,
                          all.traps = my_traps,
                          ntraps = nT, # number of cameras available
                          ndesigns = ndesigns, # number of random starting points
                          beta0 = beta0,
                          sigma = sigma, # SCR pars (log(lambda0) and usual sigma)
                          D_per_mask_cell = Dcov_for_sim)
      
      opt_traps_nT <- optSCR$Xlst %>%
        purrr::reduce(rbind) %>% data.frame() %>%
        rename(x = X, y = Y) %>%
        mutate(trap_id = rep(1:ndesigns, each = nT))
      
      for(i in 1:ndesigns){
        opt_traps_i <- opt_traps_nT %>% dplyr::filter(trap_id == i)
        opt_traps <- read.traps(data = opt_traps_i, detector = "count")
        optEnrm <- rbind(optEnrm,
                         c(nT, b_ac, beta0, Enrm(D = Dcov_for_sim_ha, traps = opt_traps, mask = mask, detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1)) %>% t() %>% data.frame())
      }
      
      opt_traps_all <- rbind(opt_traps_all, cbind(nT = nT, b_ac = b_ac, beta0 = beta0, opt_traps_nT))
      
    }
    
  }
  
}

save(mask_df, optEnrm, opt_traps_all, file = "output/Tost_Enrm_opt.Rdata")

## now for non-optimal/ other designs

# disgracefully hacky way to make a polygon just bigger than boundary of mask points
# used later to decide which randomly generated detectors are on the mask are so allowed
sap_1 <- mask_df %>% st_as_sf(coords = c("x", "y")) %>% st_buffer(dist = cellsize, endCapStyle = "SQUARE") %>% st_union() 
sap_1 <- sap_1 %>% st_buffer(dist = -cellsize * 0.4, endCapStyle = "SQUARE")

# get past trap locations
all_traps <- TostCams

### going to generate En, Er, CV under different designs in this survey area, varying different SCR parameters

survey_traps_all <- data.frame(nT = as.integer(), b_ac = as.numeric(), beta0 = as.numeric(), x = as.numeric(), y = as.numeric(), trap_id = as.integer())
grid_traps_all <- data.frame(nT = as.integer(), b_ac = as.numeric(), beta0 = as.numeric(), x = as.numeric(), y = as.numeric(), trap_id = as.integer())
opt_grid_traps_all <- data.frame(nT = as.integer(), b_ac = as.numeric(), beta0 = as.numeric(), x = as.numeric(), y = as.numeric(), trap_id = as.integer())

surveyEnrm <- data.frame(nT = as.integer(), b_ac = as.numeric(), beta0 = as.numeric(), En = as.numeric(), Er = as.numeric(), Em = as.numeric())
gridEnrm <- data.frame(nT = as.integer(), b_ac = as.numeric(), beta0 = as.numeric(), En = as.numeric(), Er = as.numeric(), Em = as.numeric())
optgridEnrm <- data.frame(nT = as.integer(), b_ac = as.numeric(), beta0 = as.numeric(), dd = as.numeric(), En = as.numeric(), Er = as.numeric(), Em = as.numeric())

for(b_ac in c(0, 1)){
  
  # D = population density animals / hectare; may be scalar or vector of length nrow(mask)
  # convert user specified dens_per_100km2 into dens per cell on the mask
  D_per_mask_cell <- dens_per_100km2 / 100 * (attr(mask, "area") / 100)
  # make density on the mask a function of some covariate (b_ac = 0 for constant density)
  if (cov_name == "stdGC") {
    covariates(mask)$Dac <- exp(b_ac * as.numeric(scale(covariates(mask)$stdGC)))
  } else if (cov_name == "x"){
    covariates(mask)$Dac <- exp(b_ac * as.numeric(scale(mask$x)))
  } else {stop("Covariate choice not supported")}
  # now standardize the Dac so that its mean is = D_per_mask_cell, or equivalently that its sum is = Dpmc * N mask cells
  Dcov_for_sim <- covariates(mask)$Dac / sum(covariates(mask)$Dac) * (D_per_mask_cell * length(mask$x))
  # turn this into equivalent density per hectate, for use with Enrm (secr functions want per ha values)
  Dcov_for_sim_ha <- Dcov_for_sim / attr(mask, "area")
  
  for(beta0 in c(-1.5, -0.8, -0.4)){
    
    for(nT in c(20,30,40,50,60)){
      
      # 1) using a subset of actual traps used for survey
      
      # get all the traps and put into df
      my_t <- all_traps %>% filter(str_detect(Mountain, region_name)) %>% filter(!is.na(Longitude)) %>% filter(!is.na(Latitude))
      my_t <- st_as_sf(my_t, coords = c("Longitude", "Latitude"), crs = 4326) %>%
        st_transform(crs = "+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs")
      traps_all <- st_coordinates(my_t) %>% as.data.frame() %>% rename(x = X, y = Y)
      
      # choose a subset of nT of these, starting from a random trap and then choosing nT nearest neighbours
      for(i in 1:nsims){
        traps <- traps_all
        idd <- sample(1:nrow(traps), 1)
        xy_rand <- traps_all[idd, ]
        traps$dist2pt <- (traps$x - xy_rand$x)^2 + (traps$y - xy_rand$y)^2
        traps <- traps %>% filter(rank(dist2pt) <= nT)
        traps <- read.traps(data = traps, detector = "count")
        
        survey_traps_all <- rbind(survey_traps_all, cbind(nT = nT, b_ac = b_ac, beta0 = beta0, x = traps$x, y = traps$y, trap_id = i))
        
        surveyEnrm <- rbind(surveyEnrm,
                            c(nT, b_ac, beta0, Enrm(D = Dcov_for_sim_ha, traps = traps, mask = mask, detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1)) %>% t() %>% data.frame())
        
      }
      
      # 2) with a grid design X sigma apart
      
      # place a grid over the area, with cells X sigma apart
      my_grid <- st_make_grid(sap_1, cellsize = c(1 * sigma, 1 * sigma), what = "centers") %>%
        st_intersection(sap_1)
      grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
      
      # choose a subset of nT of these, starting from a random trap and then choosing nT nearest neighbours
      all_grid_traps <- list()
      for(i in 1:nsims){
        grid_traps <- grid_traps_full
        xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
        grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
        grid_traps <- grid_traps %>% filter(rank(dist2pt) <= nT)
        grid_traps <- read.traps(data = grid_traps, detector = "count")
        
        all_grid_traps[[i]] <- grid_traps # for use in optimalSpacing, needs to be secr traps obj
        
        grid_traps_all <- rbind(grid_traps_all, cbind(nT = nT, b_ac = b_ac, beta0 = beta0, x = grid_traps$x, y = grid_traps$y, trap_id = i))
        
        gridEnrm <- rbind(gridEnrm,
                          c(nT, b_ac, beta0, Enrm(D = Dcov_for_sim_ha, traps = grid_traps, mask = mask, detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1)) %>% t() %>% data.frame())
        
      }
      
      # 3) with a optimal spaced grid design
      
      # for each of the partial grids above, try work out optimalSpacing and (if it exists) then
      # run another survey with this spacing
      for(i in 1:nsims){
        
        # try compute optimal spacing for that grid i
        dd <- NULL
        attempt <- 0
        while( is.null(dd) && attempt <= 1 ) {
          attempt <- attempt + 1
          try(
            dd <- optimalSpacing(D = mean(Dcov_for_sim), traps = all_grid_traps[[i]], detectpar = list(lambda0 = exp(beta0), sigma = sigma),
                                 noccasions = 1)$rotRSE$optimum.spacing
          )
        }
        
        # if optimal spacing is available
        if (!is.null(dd)) {
          
          # place a grid over the area, with cells optimally spaced
          my_grid <- st_make_grid(sap_1, cellsize = c(dd, dd), what = "centers") %>%
            st_intersection(sap_1)
          opt_grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
          
          # choose a subset of nT of these points, starting from a random trap and then choosing nT nearest neighbours
          opt_grid_traps <- opt_grid_traps_full
          xy_rand <- opt_grid_traps[sample(1:nrow(opt_grid_traps), 1), ]
          opt_grid_traps$dist2pt <- (opt_grid_traps$x - xy_rand$x)^2 + (opt_grid_traps$y - xy_rand$y)^2
          opt_grid_traps <- opt_grid_traps %>% filter(rank(dist2pt) <= nT)
          opt_grid_traps <- read.traps(data = opt_grid_traps, detector = "count")
          
          opt_grid_traps_all <- rbind(opt_grid_traps_all, cbind(nT = nT, b_ac = b_ac, beta0 = beta0, x = opt_grid_traps$x, y = opt_grid_traps$y, trap_id = i))
          
          optgridEnrm <- rbind(optgridEnrm,
                               c(nT, b_ac, beta0, dd, Enrm(D = Dcov_for_sim_ha, traps = opt_grid_traps, mask = mask, detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1)) %>% t() %>% data.frame())
          
        }
      }
      
      print(c(b_ac, beta0, nT))
      
    }
    
  }
  
}

save(surveyEnrm, gridEnrm, optgridEnrm, survey_traps_all, grid_traps_all, opt_grid_traps_all,
     file = "output/Tost_Enrm_nonopt.Rdata")
