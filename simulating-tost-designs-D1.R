######## 
### simulation check of all D1 designs
########

library(dplyr)
library(tidyr)
library(oSCR)
library(secr)
library(secrdesign)

source("fn-simulating-cvD.R")

# updates Enrm to accommdate detector covariates
# https://www.r-bloggers.com/three-ways-to-call-cc-from-r/
source("oSCR/LambdaL.R")
source("oSCR/utils_for_enrmL.R")
dyn.load("oSCR/mysecrdesign.so")

# gives error if detector type is a factor, check!

# load results
load("output/new/Tost_examples_nonuniD_new.Rdata")  
load("output/new/Tost_examples_nonopt_D1_new.Rdata")  

# put ids
ids <- c(rep(1,40), rep(2,60), rep(3, 40), rep(4,60), rep(5,40), rep(6, 60), 
         rep(7,40), rep(8,60), rep(9, 40), rep(10,60), rep(11,40), rep(12,60),
         rep(13,40), rep(14,60), rep(15,40), rep(16,60), rep(17,40), rep(18,60))
opt_traps_all <- opt_traps_D1_all %>% mutate(id = ids)
grid_traps_all <- grid_traps_all %>% mutate(id = ids)
opt_grid_traps_all <- opt_grid_traps_all %>% mutate(id = ids)

all_masks <- all_masks_D1
rm(all_masks_D1, opt_traps_D1_all)

nsim <- 500

# set up parameter values in same order used to generate optimal designs, so can associate possibly spatially
# varying D and lambda0 with each combination

b_ac <- c(-1, -1, 1, 1, 3, 3, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
b_d <- c(0, 0, 0, 0, 0, 0, -0.75, -0.75, 0.75, 0.75, 1.5, 1.5, -0.75, -0.75, 0.75, 0.75, 1.5, 1.5)
nT <- rep(c(40,60), times = 9)
dens_per_100km2 <- 2

# arguments for secr fit 
model.args_list <- list()
for(i in 1:6){model.args_list[[i]] <- formula(D ~ stdGC)}
for(i in 7:12){model.args_list[[i]] <- formula(D ~ stdGC)}
for(i in 13:18){model.args_list[[i]] <- list(formula(D ~ stdGC), formula(lambda0 ~ mnsx))}

Dcov_for_sim_ha <- list()
for(i in 1:18){
  
  # densities for sim.popn
  if(b_ac[i] == 0){ 
    Dcov_for_sim_ha[[i]] <- dens_per_100km2 / 10000 
  } else {
    # D = population density animals / hectare; may be scalar or vector of length nrow(mask)
    # convert user specified dens_per_100km2 into dens per cell on the mask
    D_per_mask_cell <- dens_per_100km2 / 100 * (attr(all_masks[[i]], "area") / 100)
    # make density on the mask a function of some covariate (b_ac = 0 for constant density)
    Dac <- exp(b_ac[i] * as.numeric(scale(covariates(all_masks[[i]])$stdGC)))
    # now standardize the Dac so that its mean is = D_per_mask_cell, or equivalently that its sum is = Dpmc * N mask cells
    Dcov_for_sim <- Dac / sum(Dac) * (D_per_mask_cell * length(all_masks[[i]]$x))
    # turn this into equivalent density per hectate, for use with secr
    Dcov_for_sim_ha[[i]] <- Dcov_for_sim / attr(all_masks[[i]], "area") 
  }
  
}

mnr_sim <- list()
lambda0_list <- list()
for(i in 1:18){
  
  D <- 2/10000
  mask <- all_masks[[i]]
  mnr_traps <- opt_traps_all %>% filter(id == i)
  # detector covariate part: possible vector lambda0 for sim.capthist
  lambda0 <- mnr_traps$lambda0[1]
  if(b_d[i] == 0) { 
    lambda0_list[[i]] <- lambda0
  } else {
    # recalculate detection lambda0 for each possible trap location
    db1 <- b_d[i]
    dLam <- exp(as.numeric(db1 * scale(-alltraps_df$x)))
    dLam <- dLam / sum(dLam) * (lambda0 * nrow(alltraps_df)) 
    # which possible trap location (in alltraps_df) is closest to each of the actual traps
    closest_locs <- apply(e2dist(mnr_traps[,c("x","y")], alltraps_df[,c("x","y")]), 1, which.min)
    lambda0_list[[i]] <- dLam[closest_locs]
  }
  sigma <- mnr_traps$sigma[1]
  nocc <- 1
  dt <- as.character(mnr_traps$dt[1])
  traps <- read.traps(data = data.frame(x = mnr_traps$x, y = mnr_traps$y), detector = dt)
  covariates(traps)$mnsx <- scale(-mnr_traps$x) # add detector covariate
  mnr_sim[[i]] <- simulating_cvD(grids = traps, masks = mask, D = Dcov_for_sim_ha[[i]], lambda0 = lambda0_list[[i]], sigma = sigma, noccasions = nocc, model.args = model.args_list[[i]], my_pop = NULL, new_pop_each_rep = TRUE, sim.popn.Ndist = "poisson", nrepl = nsim, seed = 123)
  # if detector covariates, use EnrmL
  if(length(lambda0_list[[i]]) == 1){
    enrm_mnr <- Enrm(D = Dcov_for_sim_ha[[i]], traps, mask, list(lambda0 = lambda0_list[[i]], sigma = sigma), noccasions = nocc)
  } else {
    enrm_mnr <- EnrmL(D = Dcov_for_sim_ha[[i]], traps, mask, list(lambda0 = lambda0_list[[i]], sigma = sigma), noccasions = nocc)
  }
  if(i == 1){ 
    mnr_sim_sum <- mnr_sim[[i]]$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "mnr_mod", trap_id = i, En = enrm_mnr[1], Er = enrm_mnr[2], Em = enrm_mnr[3])
  } else {
    mnr_sim_sum <- rbind(mnr_sim_sum, mnr_sim[[i]]$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "mnr_mod", trap_id = i, En = enrm_mnr[1], Er = enrm_mnr[2], Em = enrm_mnr[3]))
  }
  
}

#save(mnr_sim_sum, mnr_sim, Dcov_for_sim_ha, lambda0_list, file="output/new/mnr-res-D1.RData")

grid_sim <- list()
for(i in 1:18){
  
  D <- 2/10000
  mask <- all_masks[[i]]
  grid_traps <- grid_traps_all %>% filter(id == i)
  # detector covariate part: possible vector lambda0 for sim.capthist
  lambda0 <- grid_traps$lambda0[1]
  if(b_d[i] == 0) { 
    lambda0_list[[i]] <- lambda0
  } else {
    # recalculate detection lambda0 for each possible trap location
    db1 <- b_d[i]
    dLam <- exp(as.numeric(db1 * scale(-alltraps_df$x)))
    dLam <- dLam / sum(dLam) * (lambda0 * nrow(alltraps_df)) 
    # which possible trap location (in alltraps_df) is closest to each of the actual traps
    closest_locs <- apply(e2dist(grid_traps[,c("x","y")], alltraps_df[,c("x","y")]), 1, which.min)
    lambda0_list[[i]] <- dLam[closest_locs]
  }
  sigma <- grid_traps$sigma[1]
  nocc <- 1
  dt <- as.character(grid_traps$dt[1])
  traps <- read.traps(data = data.frame(x = grid_traps$x, y = grid_traps$y), detector = dt)
  covariates(traps)$mnsx <- scale(-grid_traps$x) # add detector covariate
  grid_sim[[i]] <- simulating_cvD(grids = traps, masks = mask, D = Dcov_for_sim_ha[[i]], lambda0 = lambda0_list[[i]], sigma = sigma, noccasions = nocc,  model.args = model.args_list[[i]], my_pop = NULL, new_pop_each_rep = TRUE, sim.popn.Ndist = "poisson", nrepl = nsim, seed = 123)
  # if detector covariates, use EnrmL
  if(length(lambda0_list[[i]]) == 1){
    enrm_grid <- Enrm(D = Dcov_for_sim_ha[[i]], traps, mask, list(lambda0 = lambda0_list[[i]], sigma = sigma), noccasions = nocc)
  } else {
    enrm_grid <- EnrmL(D = Dcov_for_sim_ha[[i]], traps, mask, list(lambda0 = lambda0_list[[i]], sigma = sigma), noccasions = nocc)
  }
  if(i == 1){ 
    grid_sim_sum <- grid_sim[[i]]$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "grid_mod", trap_id = i, En = enrm_grid[1], Er = enrm_grid[2], Em = enrm_grid[3])
  } else {
    grid_sim_sum <- rbind(grid_sim_sum, grid_sim[[i]]$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "grid_mod", trap_id = i, En = enrm_grid[1], Er = enrm_grid[2], Em = enrm_grid[3]))
  }
  
}

#save(grid_sim_sum, grid_sim,  Dcov_for_sim_ha, lambda0_list, file="output/new/grid-res-D1.RData")

opt_grid_sim <- list()
for(i in 1:18){
  
  D <- 2/10000
  mask <- all_masks[[i]]
  opt_grid_traps <- opt_grid_traps_all %>% filter(id == i)
  # detector covariate part: possible vector lambda0 for sim.capthist
  lambda0 <- opt_grid_traps$lambda0[1]
  if(b_d[i] == 0) { 
    lambda0_list[[i]] <- lambda0
  } else {
    # recalculate detection lambda0 for each possible trap location
    db1 <- b_d[i]
    dLam <- exp(as.numeric(db1 * scale(-alltraps_df$x)))
    dLam <- dLam / sum(dLam) * (lambda0 * nrow(alltraps_df)) 
    # which possible trap location (in alltraps_df) is closest to each of the actual traps
    closest_locs <- apply(e2dist(opt_grid_traps[,c("x","y")], alltraps_df[,c("x","y")]), 1, which.min)
    lambda0_list[[i]] <- dLam[closest_locs]
  }
  sigma <- opt_grid_traps$sigma[1]
  nocc <- 1
  dt <- as.character(opt_grid_traps$dt[1])
  traps <- read.traps(data = data.frame(x = opt_grid_traps$x, y = opt_grid_traps$y), detector = dt)
  covariates(traps)$mnsx <- scale(-opt_grid_traps$x) # add detector covariate
  opt_grid_sim[[i]] <- simulating_cvD(grids = traps, masks = mask, D = Dcov_for_sim_ha[[i]], lambda0 = lambda0_list[[i]], sigma = sigma, noccasions = nocc,  model.args = model.args_list[[i]], my_pop = NULL, new_pop_each_rep = TRUE, sim.popn.Ndist = "poisson", nrepl = nsim, seed = 123)
  # if detector covariates, use EnrmL
  if(length(lambda0_list[[i]]) == 1){
    enrm_opt_grid <- Enrm(D = Dcov_for_sim_ha[[i]], traps, mask, list(lambda0 = lambda0_list[[i]], sigma = sigma), noccasions = nocc)
  } else {
    enrm_opt_grid <- EnrmL(D = Dcov_for_sim_ha[[i]], traps, mask, list(lambda0 = lambda0_list[[i]], sigma = sigma), noccasions = nocc)
  }
  if(i == 1){ 
    opt_grid_sim_sum <- opt_grid_sim[[i]]$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "opt_grid_mod", trap_id = i, En = enrm_opt_grid[1], Er = enrm_opt_grid[2], Em = enrm_opt_grid[3])
  } else {
    opt_grid_sim_sum <- rbind(opt_grid_sim_sum, opt_grid_sim[[i]]$mod_summary %>% filter(E.N < 4 * true.N) %>% mutate(design = "opt_grid_mod", trap_id = i, En = enrm_opt_grid[1], Er = enrm_opt_grid[2], Em = enrm_opt_grid[3]))
  }
  
}

#save(opt_grid_sim_sum, opt_grid_sim,  Dcov_for_sim_ha, lambda0_list, file="output/new/opt_grid-res-D1.RData")

#save(mnr_sim_sum, grid_sim_sum, opt_grid_sim_sum, file="output/new/all-res-D1-simsum-only.RData")

###

library(TSP)
library(kableExtra)

load("output/new/all-res-D1-simsum-only.RData")

# combine
sim_sum <- rbind(grid_sim_sum, opt_grid_sim_sum, mnr_sim_sum)

# add shortest paths

opt_traps_all <- opt_traps_all %>% group_by(id) %>% mutate(shortest_cycle = tour_length(solve_TSP(ETSP(data.frame(x,y)))),
                                                           shortest_path = tour_length(solve_TSP(insert_dummy(TSP(dist(data.frame(x,y))))))) %>% ungroup()

grid_traps_all <- grid_traps_all %>% group_by(id) %>% mutate(shortest_cycle = tour_length(solve_TSP(ETSP(data.frame(x,y)))),
                                                             shortest_path = tour_length(solve_TSP(insert_dummy(TSP(dist(data.frame(x,y))))))) %>% ungroup()
opt_grid_traps_all <- opt_grid_traps_all %>% group_by(id) %>% mutate(shortest_cycle = tour_length(solve_TSP(ETSP(data.frame(x,y)))),
                                                                     shortest_path = tour_length(solve_TSP(insert_dummy(TSP(dist(data.frame(x,y))))))) %>% ungroup()

all_shortest_paths <- rbind(opt_traps_all %>% mutate(design = "mnr_mod"), 
                            grid_traps_all %>% mutate(design = "grid_mod"), 
                            opt_grid_traps_all %>% mutate(design = "opt_grid_mod")) %>% 
  group_by(design, id) %>% summarize(shortest_path = first(shortest_path)) %>% rename(trap_id = id) %>% ungroup()


# compare

all_sum <- sim_sum %>% group_by(trap_id, design) %>% summarize(approx_cv = 1 / sqrt(min(mean(En), mean(Er))),
                                                               sim_cv = sd(E.N) / mean(E.N),
                                                               n = mean(n), En = mean(En),
                                                               r = mean(r), Er = mean(Er),
                                                               m = mean(moves), Em = mean(Em),
                                                               ecva = sd(esa) / mean(esa),
                                                               esa = mean(esa),
                                                               cvN = mean(cvN), cva = mean(cva), cvD = mean(cvD),
                                                               nsim = n()) %>% ungroup() %>%
  left_join(all_shortest_paths, by = c("trap_id", "design"))

all_sum %>% 
  mutate(b_ac = rep(b_ac, each = 3), b_d = rep(b_d, each = 3), nT = rep(nT, each = 3), cv_a = round(100*ecva,0),
         sim_cv = round(100*sim_cv,0), approx_cv = round(100*approx_cv,0), shortest_path = round(shortest_path/1000,0)) %>%
  mutate(design = ifelse(design == "mnr_mod", "amnr_mod", design)) %>% arrange(design) %>% 
  filter(nT == 40) %>%
  dplyr::select(trap_id, b_ac, b_d, design, sim_cv, approx_cv, cv_a) %>%
  pivot_wider(names_from = "design", values_from = c("sim_cv", "approx_cv", "cv_a")) %>%
  #dplyr::select(trap_id, b_ac, b_d, mnr_mod, grid_mod, opt_grid_mod) %>%
  kable("latex", caption = "Group Rows", booktabs = T) %>% 
  kable_styling() %>%
  pack_rows("", 1, 3) %>%
  pack_rows("", 4, 6) %>%
  pack_rows("", 7, 9) 

all_sum %>% 
  mutate(b_ac = rep(b_ac, each = 3), b_d = rep(b_d, each = 3), nT = rep(nT, each = 3), cv_a = round(100*ecva,0),
         sim_cv = round(100*sim_cv,0), approx_cv = round(100*approx_cv,0), shortest_path = round(shortest_path/1000,0)) %>%
  mutate(design = ifelse(design == "mnr_mod", "amnr_mod", design)) %>% arrange(design) %>% 
  filter(nT == 40) %>%
  dplyr::select(trap_id, b_ac, b_d, design, shortest_path) %>%
  pivot_wider(names_from = "design", values_from = c("shortest_path")) %>%
  #dplyr::select(trap_id, b_ac, b_d, mnr_mod, grid_mod, opt_grid_mod) %>%
  kable("latex", caption = "Group Rows", booktabs = T) %>% 
  kable_styling() %>%
  pack_rows("", 1, 3) %>%
  pack_rows("", 4, 6) %>%
  pack_rows("", 7, 9) 
