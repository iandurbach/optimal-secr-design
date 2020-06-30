library(dplyr)
library(sf)
library(secr)
library(secrdesign)
library(oSCR)
library(stringr)
library(ggplot2)

# load example designs for non-Uniform D, euclidean
load("output/new/Tost_examples_nonuniD_new.Rdata")  
# load the non-optimized arrays
load("output/new/Tost_examples_nonopt_D1_new.Rdata")  

opt_traps_all <- opt_traps_D1_all %>% mutate(method = "mnr")
grid_traps_all <- grid_traps_all %>% mutate(method = "grid")
opt_grid_traps_all <- opt_grid_traps_all %>% mutate(method = "opt_grid")

# only keep those with non-uniform density
opt_traps_all <- opt_traps_all %>% filter(b_det == 0)
grid_traps_all <- grid_traps_all %>% filter(b_det == 0)
opt_grid_traps_all <- opt_grid_traps_all %>% filter(b_det == 0)

mask <- read.mask(data = mask_df)

traps <- opt_traps_all %>% filter(trap_id == 1) %>%
  mutate(ID = paste0(nT, "_", b_ac)) %>% dplyr::select(ID,x,y)
ID <- traps$ID
traps <- split(traps[,-1], traps$ID)

# average D
dens_per_100km2 <- 2 # mean animal density per 100km2, SLs are ~1
D_per_mask_cell <- dens_per_100km2 / 100 * (attr(mask, "area") / 100)
D_per_ha <- dens_per_100km2 / 10000

# detection function parameters
lambda0 = 1
beta0 = log(lambda0)
g0 = 1-exp(-lambda0)
sigma <- 3000

Enrm_results <- c()

for(i in 1:length(traps)){
  
  # get the number of detectors from the name of the list element
  nT <- as.numeric(str_split(names(traps)[i], pattern = "_", simplify = TRUE)[1])  
  # get the value of the covariate from the name of the list element
  b_ac <- as.numeric(str_split(names(traps)[i], pattern = "_", simplify = TRUE)[2])
  # make density on the mask a function of some covariate (b_ac = 0 for constant density)
  covariates(mask)$Dac <- exp(b_ac * as.numeric(scale(covariates(mask)$stdGC)))
  # now standardize the Dac so that its mean is = D_per_mask_cell, or equivalently that its sum is = Dpmc * N mask cells
  Dcov_for_sim <- covariates(mask)$Dac / sum(covariates(mask)$Dac) * (D_per_mask_cell * length(mask$x))
  # turn this into equivalent density per hectate, for use with Enrm (secr functions want per ha values)
  Dcov_for_sim_ha <- Dcov_for_sim / attr(mask, "area")
  
  # optimal trap placements
  my_opt_traps <- read.traps(data = traps[[i]], detector = "count")
  # grid trap placement 
  this_nT <- nT; this_b_ac <- 1; this_beta0 <- beta0
  my_grid_traps <- grid_traps_all %>% dplyr::filter(nT == this_nT, b_ac == this_b_ac, beta0 == this_beta0)
  ID <- my_grid_traps$trap_id
  my_grid_traps <- split(data.frame(x = my_grid_traps$x, y = my_grid_traps$y), ID)
  # opt_grid trap placement 
  my_opt_grid_traps <- opt_grid_traps_all %>% dplyr::filter(nT == this_nT, b_ac == this_b_ac, beta0 == this_beta0)
  ID <- my_opt_grid_traps$trap_id
  my_opt_grid_traps <- split(data.frame(x = my_opt_grid_traps$x, y = my_opt_grid_traps$y), ID)

  seqprop <- seq(from = 0.5, to = 1.5, length.out = 11)
  # now evaluate Enrm under a potentially INCORRECT guess for lambda0 
  
  guess_beta <- log(seq(from = 0.5, to = 1.5, length.out = 11) * lambda0)
  
  for(j2 in 1:length(guess_beta)){
    
    j <- guess_beta[j2]
    
    # opt traps
    opt_res <- c(error = seqprop[j2], nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 1, method = 1, k = 1,
                 Enrm(D = Dcov_for_sim_ha, traps = my_opt_traps, mask = mask, detectpar = list(lambda0 = exp(j), sigma = sigma), noccasions = 1))
    # grid traps
    grid_res <- c()
    for(k in 1:length(my_grid_traps)){
      grid_res <- rbind(grid_res, c(error = seqprop[j2], nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 1, method = 2, k = k,
                                    Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_grid_traps[[k]], detector = "count"), mask = mask, 
                                         detectpar = list(lambda0 = exp(j), sigma = sigma), noccasions = 1)))
    }
    # opt grid traps
    opt_grid_res <- c()
    for(k in 1:length(my_opt_grid_traps)){
      opt_grid_res <- rbind(opt_grid_res, c(error = seqprop[j2], nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 1, method = 3, k = k,
                                            Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_opt_grid_traps[[k]], detector = "count"), mask = mask, 
                                                 detectpar = list(lambda0 = exp(j), sigma = sigma), noccasions = 1)))
    }

    Enrm_results <- rbind(Enrm_results, opt_res, grid_res, opt_grid_res)
    
  }
  
  # and now evaluate Enrm under a potentially INCORRECT guess for sigma 
  
  guess_sigma <- seq(from = 0.5, to = 1.5, length.out = 11) * sigma
  
  for(j2 in 1:length(guess_sigma)){
    
    j <- guess_sigma[j2]
    
    # opt traps
    opt_res <- c(error = seqprop[j2], nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 2, method = 1, k = 1,
                 Enrm(D = Dcov_for_sim_ha, traps = my_opt_traps, mask = mask, detectpar = list(lambda0 = exp(beta0), sigma = j), noccasions = 1))
    # grid traps
    grid_res <- c()
    for(k in 1:length(my_grid_traps)){
      grid_res <- rbind(grid_res, c(error = seqprop[j2], nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 2, method = 2, k = k,
                                    Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_grid_traps[[k]], detector = "count"), mask = mask, 
                                         detectpar = list(lambda0 = exp(beta0), sigma = j), noccasions = 1)))
    }
    # opt grid traps
    opt_grid_res <- c()
    for(k in 1:length(my_opt_grid_traps)){
      opt_grid_res <- rbind(opt_grid_res, c(error = seqprop[j2], nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 2, method = 3, k = k,
                                            Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_opt_grid_traps[[k]], detector = "count"), mask = mask, 
                                                 detectpar = list(lambda0 = exp(beta0), sigma = j), noccasions = 1)))
    }

    Enrm_results <- rbind(Enrm_results, opt_res, grid_res, opt_grid_res)
    
  }
  
  # and now evaluate Enrm under a potentially INCORRECT guess for b_ac 
  
  guess_b_ac <- seq(from = 0.5, to = 1.5, length.out = 11) * b_ac
  
  for(j2 in 1:length(guess_b_ac)){
    
    j <- guess_b_ac[j2]
    
    # make density on the mask a function of some covariate (b_ac = 0 for constant density)
    covariates(mask)$Dac2 <- exp(j * as.numeric(scale(covariates(mask)$stdGC)))
    # now standardize the Dac so that its mean is = D_per_mask_cell, or equivalently that its sum is = Dpmc * N mask cells
    Dcov_for_sim <- covariates(mask)$Dac2 / sum(covariates(mask)$Dac2) * (D_per_mask_cell * length(mask$x))
    # turn this into equivalent density per hectate, for use with Enrm (secr functions want per ha values)
    Dcov_for_sim_ha <- Dcov_for_sim / attr(mask, "area")
    
    # opt traps
    opt_res <- c(error = seqprop[j2], nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 3, method = 1, k = 1,
                 Enrm(D = Dcov_for_sim_ha, traps = my_opt_traps, mask = mask, detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1))
    # grid traps
    grid_res <- c()
    for(k in 1:length(my_grid_traps)){
      grid_res <- rbind(grid_res, c(error = seqprop[j2], nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 3, method = 2, k = k,
                                    Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_grid_traps[[k]], detector = "count"), mask = mask, 
                                         detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1)))
    }
    # opt grid traps
    opt_grid_res <- c()
    for(k in 1:length(my_opt_grid_traps)){
      opt_grid_res <- rbind(opt_grid_res, c(error = seqprop[j2], nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 3, method = 3, k = k,
                                            Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_opt_grid_traps[[k]], detector = "count"), mask = mask, 
                                                 detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1)))
    }

    Enrm_results <- rbind(Enrm_results, opt_res, grid_res, opt_grid_res)
    
  }
  
}

save(Enrm_results, file = "output/new/posthoc-check-robustness-to-misspec.RData")
