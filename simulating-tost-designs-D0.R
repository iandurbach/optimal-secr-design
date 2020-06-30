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
load("output/new/Tost_examples_all_D0_new.Rdata")  
load("output/new/Tost_examples_nonopt_D0_new.Rdata")  

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

save(mnr_sim_sum, mnr_sim, file="output/new/mnr-res.RData")

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

save(grid_sim_sum, grid_sim, file="output/new/grid-res.RData")

opt_grid_sim <- list()
for(i in 1:18){
  
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

save(opt_grid_sim_sum, opt_grid_sim, file="output/new/opt_grid-res.RData")

save(grid_sim_sum, opt_grid_sim_sum, mnr_sim_sum, file = "output/new/all-res-D0-simsum-only.RData")

###

load("output/new/all-res-D0-simsum-only.RData")

# combine
sim_sum <- rbind(grid_sim_sum, opt_grid_sim_sum, mnr_sim_sum)

# compare

all_sum <- sim_sum %>% group_by(trap_id, design) %>% summarize(approx_cv = 1 / sqrt(min(mean(En), mean(Er))),
                                                               sim_cv = sd(E.N) / mean(E.N),
                                                               n = mean(n), En = mean(En),
                                                               r = mean(r), Er = mean(Er),
                                                               m = mean(moves), Em = mean(Em),
                                                               ecva = sd(esa) / mean(esa),
                                                               esa = mean(esa),
                                                               cvN = mean(cvN), cva = mean(cva), cvD = mean(cvD),
                                                               nsim = n()) 

library(kableExtra)
all_sum %>% mutate(design = ifelse(design == "mnr_mod", "amnr_mod", design)) %>% arrange(design) %>%
  mutate(En = round(En,1), Er = round(Er, 1), 
         nn = paste0(round(n,1)," (",round(En,1),")"), rr = paste0(round(r,1)," (",round(Er,1),")"), 
         cv_a = round(100*ecva,0),
         sim_cv = round(100*sim_cv,0), approx_cv = round(100*approx_cv,0), shortest_path = round(shortest_path/1000,0)) %>% 
  dplyr::select(trap_id, design, sim_cv, approx_cv, cv_a) %>%
  pivot_wider(names_from = "design", values_from = c("sim_cv", "approx_cv", "cv_a")) %>%
  #dplyr::select(trap_id, mnr_mod, grid_mod, opt_grid_mod) %>%
  kable("latex", caption = "Group Rows", booktabs = T) %>% 
  kable_styling() %>%
  pack_rows("", 1, 3) %>%
  pack_rows("", 4, 6) %>%
  pack_rows("", 7, 9) %>%
  pack_rows("", 10, 12) %>%
  pack_rows("", 13, 15) %>%
  pack_rows("", 16, 18)

# add shortest paths

library(TSP)
opt_traps_all <- opt_traps_all %>% group_by(id) %>% mutate(shortest_cycle = tour_length(solve_TSP(ETSP(data.frame(x,y)))),
                                                           shortest_path = tour_length(solve_TSP(insert_dummy(TSP(dist(data.frame(x,y))))))) %>% ungroup()
grid_traps_all <- grid_traps_all %>% dplyr::select(-dist2pt) %>% group_by(id) %>% mutate(shortest_cycle = tour_length(solve_TSP(ETSP(data.frame(x,y)))),
                                                                                         shortest_path = tour_length(solve_TSP(insert_dummy(TSP(dist(data.frame(x,y))))))) %>% ungroup()
opt_grid_traps_all <- opt_grid_traps_all %>% dplyr::select(-dist2pt) %>% group_by(id) %>% mutate(shortest_cycle = tour_length(solve_TSP(ETSP(data.frame(x,y)))),
                                                                                                 shortest_path = tour_length(solve_TSP(insert_dummy(TSP(dist(data.frame(x,y))))))) %>% ungroup()

all_shortest_paths <- rbind(opt_traps_all %>% mutate(design = "mnr_mod"), 
                            grid_traps_all %>% mutate(design = "grid_mod"), 
                            opt_grid_traps_all %>% mutate(design = "opt_grid_mod")) %>% 
  group_by(design, id) %>% summarize(shortest_path = first(shortest_path)) %>% rename(trap_id = id) %>% ungroup()

all_sum <-  all_sum %>%
  left_join(all_shortest_paths, by = c("trap_id", "design"))

all_sum %>% mutate(design = ifelse(design == "mnr_mod", "amnr_mod", design)) %>% arrange(design) %>%
  mutate(En = round(En,1), Er = round(Er, 1), 
         nn = paste0(round(n,1)," (",round(En,1),")"), rr = paste0(round(r,1)," (",round(Er,1),")"), 
         sim_cv = round(100*sim_cv,0), approx_cv = round(100*approx_cv,0), shortest_path = round(shortest_path/1000,0)) %>% 
  dplyr::select(trap_id, design, shortest_path) %>%
  pivot_wider(names_from = "design", values_from = c("shortest_path")) %>%
  #dplyr::select(trap_id, mnr_mod, grid_mod, opt_grid_mod) %>%
  kable("latex", caption = "Group Rows", booktabs = T) %>% 
  kable_styling() %>%
  pack_rows("", 1, 3) %>%
  pack_rows("", 4, 6) %>%
  pack_rows("", 7, 9) %>%
  pack_rows("", 10, 12) %>%
  pack_rows("", 13, 15) %>%
  pack_rows("", 16, 18)


