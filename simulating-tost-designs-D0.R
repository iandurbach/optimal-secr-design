######## 
### simulation check of all designs
########

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(stringr)
library(purrr)
library(secr)
library(secrdesign)

source("fn-simulating-cvD.R")

# load results
load("output/Tost_examples_all_D0.Rdata")  
load("output/Tost_examples_nonopt_D0.Rdata")  

# put ids
ids <- c(rep(1,20), rep(2,40), rep(3, 60), rep(4,20), rep(5,40), rep(6, 60), 
         rep(7,20), rep(8,40), rep(9, 60), rep(10,60), rep(11,60), rep(12,60),
         rep(13,60), rep(14,60), rep(15,60), rep(16,60), rep(17,60), rep(18,60))
opt_traps_all <- opt_traps_all %>% mutate(id = ids)
grid_traps_all <- grid_traps_all %>% mutate(id = ids)
opt_grid_traps_all <- opt_grid_traps_all %>% mutate(id = ids)

nsim <- 1000

mnr_sim <- list()
for(i in 1:18){
  
  D <- 2/10000
  mask <- all_masks[[i]]
  mnr_traps <- opt_traps_all %>% filter(id == i)
  lambda0 <- ifelse(i < 16, mnr_traps$lambda0[1], mnr_traps$lambda0[1]/5)
  sigma <- mnr_traps$sigma[1]
  nocc <- ifelse(i < 16, 1, 5)
  dt <- mnr_traps$dt[1]
  traps <- read.traps(data = data.frame(x = mnr_traps$x, y = mnr_traps$y), detector = dt)
  mnr_sim[[i]] <- simulating_cvD(grids = traps, masks = mask, D = D, lambda0 = lambda0, sigma = sigma, noccasions = nocc, my_pop = NULL, new_pop_each_rep = TRUE, sim.popn.Ndist = "poisson", nrepl = nsim, seed = 123)
  enrm_mnr <- Enrm(D = D, traps, mask, list(lambda0 = lambda0, sigma = sigma), noccasions = nocc)
  if(i == 1){ 
    mnr_sim_sum <- mnr_sim[[i]]$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "mnr_mod", trap_id = i, En = enrm_mnr[1], Er = enrm_mnr[2], Em = enrm_mnr[3])
  } else {
    mnr_sim_sum <- rbind(mnr_sim_sum, mnr_sim[[i]]$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "mnr_mod", trap_id = i, En = enrm_mnr[1], Er = enrm_mnr[2], Em = enrm_mnr[3]))
  }
  
}

#save(mnr_sim_sum, mnr_sim, file="output/mnr-res.RData")

grid_sim <- list()
for(i in 1:18){
  
  D <- 2/10000
  mask <- all_masks[[i]]
  grid_traps <- grid_traps_all %>% filter(id == i)
  lambda0 <- ifelse(i < 16, grid_traps$lambda0[1], grid_traps$lambda0[1]/5)
  sigma <- grid_traps$sigma[1]
  nocc <- ifelse(i < 16, 1, 5)
  dt <- as.character(grid_traps$dt[1])
  traps <- read.traps(data = data.frame(x = grid_traps$x, y = grid_traps$y), detector = dt)
  grid_sim[[i]] <- simulating_cvD(grids = traps, masks = mask, D = D, lambda0 = lambda0, sigma = sigma, noccasions = nocc, my_pop = NULL, new_pop_each_rep = TRUE, sim.popn.Ndist = "poisson", nrepl = nsim, seed = 123)
  enrm_grid <- Enrm(D = D, traps, mask, list(lambda0 = lambda0, sigma = sigma), noccasions = nocc)
  if(i == 1){ 
    grid_sim_sum <- grid_sim[[i]]$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "grid_mod", trap_id = i, En = enrm_grid[1], Er = enrm_grid[2], Em = enrm_grid[3])
  } else {
    grid_sim_sum <- rbind(grid_sim_sum, grid_sim[[i]]$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "grid_mod", trap_id = i, En = enrm_grid[1], Er = enrm_grid[2], Em = enrm_grid[3]))
  }
  
}

#save(grid_sim_sum, grid_sim, file="output/grid-res.RData")

opt_grid_sim <- list()
for(i in 10:12){
  
  D <- 2/10000
  mask <- all_masks[[i]]
  opt_grid_traps <- opt_grid_traps_all %>% filter(id == i)
  lambda0 <- ifelse(i < 16, opt_grid_traps$lambda0[1], opt_grid_traps$lambda0[1]/5)
  sigma <- opt_grid_traps$sigma[1]
  nocc <- ifelse(i < 16, 1, 5)
  dt <- as.character(opt_grid_traps$dt[1])
  traps <- read.traps(data = data.frame(x = opt_grid_traps$x, y = opt_grid_traps$y), detector = dt)
  opt_grid_sim[[i]] <- simulating_cvD(grids = traps, masks = mask, D = D, lambda0 = lambda0, sigma = sigma, noccasions = nocc, my_pop = NULL, new_pop_each_rep = TRUE, sim.popn.Ndist = "poisson", nrepl = nsim, seed = 123)
  enrm_opt_grid <- Enrm(D = D, traps, mask, list(lambda0 = lambda0, sigma = sigma), noccasions = nocc)
  if(i == 1){ 
    opt_grid_sim_sum <- opt_grid_sim[[i]]$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "opt_grid_mod", trap_id = i, En = enrm_opt_grid[1], Er = enrm_opt_grid[2], Em = enrm_opt_grid[3])
  } else {
    opt_grid_sim_sum <- rbind(opt_grid_sim_sum, opt_grid_sim[[i]]$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "opt_grid_mod", trap_id = i, En = enrm_opt_grid[1], Er = enrm_opt_grid[2], Em = enrm_opt_grid[3]))
  }
  
}

save(opt_grid_sim_sum, opt_grid_sim, file="output/no_commit/opt_grid-res.RData")

save(grid_sim_sum, opt_grid_sim_sum, mnr_sim_sum, file = "output/all-res-D0-simsum-only.RData")
