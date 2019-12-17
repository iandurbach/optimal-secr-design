## 

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(secr)
library(oSCR)
library(secrdesign)

source("oSCR/SCRdesign_mcvd.R")
source("oSCR/Qfn_mcvd.R")

# input file contains data frames with mask locations, and all potential camera locations
load("data/TTostExample.RData")

# design parameters
ndesigns <- 1
nT <- 30

# assumptions about animal density and movement
sigma <- 4000
beta0 <- -0.1
dens_per_100km2 <- 2 

# statespace = all mask points
statespace <- meshlocs[,1:2] %>% as.matrix()
# traps = subset of mask points that are potential camera locations
traps <- traplocs[,1:2] %>% as.matrix()

# calculate area of mask cell 
# just specify if known, need to convert dens_per_100km2 into D_per_mask_cell
secr_mask <- read.mask(data = meshlocs)
area_of_mask_cell <- attr(secr_mask, "area") # in hectares 

# generate optimal design
optSCR <- SCRdesign(statespace = statespace,
                    all.traps = traps,
                    ntraps = nT, # number of cameras available
                    ndesigns = ndesigns, # number of random starting points
                    beta0 = beta0,
                    sigma = sigma, # SCR pars (log(lambda0) and usual sigma)
                    D_per_mask_cell = dens_per_100km2 / 100 * (area_of_mask_cell / 100),
                    occasions = 1,
                    detector = "count")

# extract camera locations
opt_traps_nT <- optSCR$Xlst %>%
  purrr::reduce(rbind) %>% data.frame() %>%
  mutate(trap_id = rep(1:ndesigns, each = nT))

# choose one of the runs to plot, might need to rename variables
opt_traps <- opt_traps_nT %>% filter(trap_id == 1)

# plot
opt_traps %>%
  ggplot(aes(x = x, y = y)) + 
  geom_tile(data = traplocs, fill = "gray80", colour = "black") +
  geom_point(size = 1, colour = "red") + 
  scale_fill_viridis_c() + coord_equal()

# do some checks on the design

my_mask <- read.mask(data = meshlocs)
my_traps <- read.traps(data = opt_traps, detector = "count") 
Dcov_for_sim_ha <- dens_per_100km2 / 10000 # assumes uniform D
lambda0 <- exp(beta0)

# simulating directly from secr
all_mod <- list()
all_ch <- list()
mod_summary <- data.frame(j = as.integer(), 
                          E.N = as.numeric(), R.N = as.numeric(), E.Nlcl = as.numeric(), R.Nlcl = as.numeric(), 
                          E.Nucl = as.numeric(), R.Nucl = as.numeric(), true.N = as.numeric(),
                          n = as.numeric(), detections = as.numeric(), dets_visited = as.numeric())
cnt <- 1
for(j in 1:10){ # run more when doing for real!
  
  simulated_points_Dcov <- sim.popn(D = Dcov_for_sim_ha, 
                                    core = my_mask, 
                                    model2D = "IHP",
                                    Ndist = "fixed")
  
  ch <- sim.capthist(traps = my_traps, pop = simulated_points_Dcov, 
                     noccasions = 1,
                     detectpar = list(lambda0 = lambda0, sigma = sigma), 
                     detectfn = "HHN")
  
  # make starting values for secr.fit
  startvals <- list(D = mean(Dcov_for_sim_ha), lambda0 = lambda0, sigma = sigma)
  
  mod <- try(secr.fit(capthist = ch, mask = my_mask, model = list(D ~ 1), 
                      detectfn = "HHN", start = startvals))
  
  if(class(mod) != "try-error"){
    
    rN <- try(region.N(mod))
    
    if(class(rN) != "try-error"){
      
      all_mod[[cnt]] <- mod
      all_ch[[cnt]] <- ch
      mod_summary <- rbind(mod_summary, 
                           data.frame(j = j, 
                                      E.N = region.N(mod)$estimate[1], R.N = region.N(mod)$estimate[2],
                                      E.Nlcl = region.N(mod)$lcl[1], R.Nlcl = region.N(mod)$lcl[2], 
                                      E.Nucl = region.N(mod)$ucl[1], R.Nucl = region.N(mod)$ucl[2], 
                                      true.N = nrow(simulated_points_Dcov),
                                      n = summary(ch)[[4]][1,1], detections = summary(ch)[[4]][6,1], dets_visited = summary(ch)[[4]][7,1]))
      
      cnt <- cnt + 1
      
    }
    
    
  }
  
}

# some checks
all_mod[[1]]
region.N(all_mod[[1]])
mod_summary 

# using secrdesign
scen <- make.scenarios(D = Dcov_for_sim_ha, detectfn = "HHN", lambda0 = lambda0, sigma = sigma, noccasions = 1, popindex = 1)
sims <- run.scenarios(100, scenarios = scen, trapset = my_traps, maskset = my_mask, pop.args = list(Ndist = "fixed"), fit = FALSE)
summary(sims)
sims_full <- run.scenarios(10, scenarios = scen, trapset = my_traps, maskset = my_mask, pop.args = list(Ndist = "fixed"), fit = TRUE)
stats <- select.stats(sims_full, parameter = "D", statistics = c("estimate", "lcl", "ucl", "RB", "RSE", "COV"))
stats$output
