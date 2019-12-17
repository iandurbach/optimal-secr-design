library(dplyr)
library(sf)
library(secr)
library(secrdesign)
library(oSCR)
library(stringr)
library(ggplot2)

# load example designs for non-Uniform D, euclidean
load("output/Tost_examples_nonuniD.Rdata")  
# load all the non-optimized arrays
load("output/Tost_many_D0_for_dbns.Rdata")  

grid_traps_all <- all_traps %>% filter(method == "grid")
opt_grid_traps_all <- all_traps %>% filter(method == "opt_grid")
survey_traps_all <- all_traps %>% filter(method == "survey")

mask <- read.mask(data = mask_df)

traps <- opt_traps_all %>% filter(trap_id == 1) %>%
  mutate(ID = paste0(nT, "_", b_ac)) %>% dplyr::select(ID,x,y)
ID <- traps$ID
traps <- split(traps[,-1], traps$ID)

# average D
dens_per_100km2 <- 1 # mean animal density per 100km2, SLs are ~1
D_per_mask_cell <- dens_per_100km2 / 100 * (attr(mask, "area") / 100)
D_per_ha <- dens_per_100km2 / 10000

# detection function parameters
beta0 <- -0.8
lambda0 = exp(beta0)
g0 = 1-exp(-lambda0)
sigma <- 6000

Enrm_results <- c()

for(i in 1:9){
  
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
  # survey trap placement (not run if nT > 40 cos survey only used 40)
  if(nT <= 40){
    my_survey_traps <- survey_traps_all %>% dplyr::filter(nT == this_nT, b_ac == this_b_ac, beta0 == this_beta0)
    ID <- my_survey_traps$trap_id
    my_survey_traps <- split(data.frame(x = my_survey_traps$x, y = my_survey_traps$y), ID)
  }
  
  # now evaluate Enrm under a potentially INCORRECT guess for lambda0 
  
  guess_beta <- seq(from = 0.5, to = 1.5, length.out = 11) * beta0
  
  for(j in guess_beta){
    
    # opt traps
    opt_res <- c(error = j / beta0, nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 1, method = 1, k = 1,
                 Enrm(D = Dcov_for_sim_ha, traps = my_opt_traps, mask = mask, detectpar = list(lambda0 = exp(j), sigma = sigma), noccasions = 1))
    # grid traps
    grid_res <- c()
    for(k in 1:length(my_grid_traps)){
      grid_res <- rbind(grid_res, c(error = j / beta0, nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 1, method = 2, k = k,
                                    Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_grid_traps[[k]], detector = "count"), mask = mask, 
                                         detectpar = list(lambda0 = exp(j), sigma = sigma), noccasions = 1)))
    }
    # opt grid traps
    opt_grid_res <- c()
    for(k in 1:length(my_opt_grid_traps)){
      opt_grid_res <- rbind(opt_grid_res, c(error = j / beta0, nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 1, method = 3, k = k,
                                            Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_opt_grid_traps[[k]], detector = "count"), mask = mask, 
                                                 detectpar = list(lambda0 = exp(j), sigma = sigma), noccasions = 1)))
    }
    # survey traps
    survey_res <- c()
    if(nT <= 40){
      for(k in 1:length(my_survey_traps)){
        survey_res <- rbind(survey_res, c(error = j / beta0, nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 1, method = 4, k = k,
                                          Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_survey_traps[[k]], detector = "count"), mask = mask, 
                                               detectpar = list(lambda0 = exp(j), sigma = sigma), noccasions = 1)))
      }
    }
    
    Enrm_results <- rbind(Enrm_results, opt_res, grid_res, opt_grid_res, survey_res)
    
  }
  
  # and now evaluate Enrm under a potentially INCORRECT guess for sigma 
  
  guess_sigma <- seq(from = 0.5, to = 1.5, length.out = 11) * sigma
  
  for(j in guess_sigma){
    
    # opt traps
    opt_res <- c(error = j / sigma, nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 2, method = 1, k = 1,
                 Enrm(D = Dcov_for_sim_ha, traps = my_opt_traps, mask = mask, detectpar = list(lambda0 = exp(beta0), sigma = j), noccasions = 1))
    # grid traps
    grid_res <- c()
    for(k in 1:length(my_grid_traps)){
      grid_res <- rbind(grid_res, c(error = j / sigma, nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 2, method = 2, k = k,
                                    Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_grid_traps[[k]], detector = "count"), mask = mask, 
                                         detectpar = list(lambda0 = exp(beta0), sigma = j), noccasions = 1)))
    }
    # opt grid traps
    opt_grid_res <- c()
    for(k in 1:length(my_opt_grid_traps)){
      opt_grid_res <- rbind(opt_grid_res, c(error = j / sigma, nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 2, method = 3, k = k,
                                            Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_opt_grid_traps[[k]], detector = "count"), mask = mask, 
                                                 detectpar = list(lambda0 = exp(beta0), sigma = j), noccasions = 1)))
    }
    # survey traps
    survey_res <- c()
    if(nT <= 40){
      for(k in 1:length(my_survey_traps)){
        survey_res <- rbind(survey_res, c(error = j / sigma, nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 2, method = 4, k = k,
                                          Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_survey_traps[[k]], detector = "count"), mask = mask, 
                                               detectpar = list(lambda0 = exp(beta0), sigma = j), noccasions = 1)))
      }
    }
    
    Enrm_results <- rbind(Enrm_results, opt_res, grid_res, opt_grid_res, survey_res)
    
  }
  
  # and now evaluate Enrm under a potentially INCORRECT guess for b_ac 
  
  guess_b_ac <- seq(from = 0.5, to = 1.5, length.out = 11) * b_ac
  
  for(j in guess_b_ac){
    
    # make density on the mask a function of some covariate (b_ac = 0 for constant density)
    covariates(mask)$Dac2 <- exp(j * as.numeric(scale(covariates(mask)$stdGC)))
    # now standardize the Dac so that its mean is = D_per_mask_cell, or equivalently that its sum is = Dpmc * N mask cells
    Dcov_for_sim <- covariates(mask)$Dac2 / sum(covariates(mask)$Dac2) * (D_per_mask_cell * length(mask$x))
    # turn this into equivalent density per hectate, for use with Enrm (secr functions want per ha values)
    Dcov_for_sim_ha <- Dcov_for_sim / attr(mask, "area")
    
    # opt traps
    opt_res <- c(error = j / b_ac, nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 3, method = 1, k = 1,
                 Enrm(D = Dcov_for_sim_ha, traps = my_opt_traps, mask = mask, detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1))
    # grid traps
    grid_res <- c()
    for(k in 1:length(my_grid_traps)){
      grid_res <- rbind(grid_res, c(error = j / b_ac, nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 3, method = 2, k = k,
                                    Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_grid_traps[[k]], detector = "count"), mask = mask, 
                                         detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1)))
    }
    # opt grid traps
    opt_grid_res <- c()
    for(k in 1:length(my_opt_grid_traps)){
      opt_grid_res <- rbind(opt_grid_res, c(error = j / b_ac, nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 3, method = 3, k = k,
                                            Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_opt_grid_traps[[k]], detector = "count"), mask = mask, 
                                                 detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1)))
    }
    # survey traps
    survey_res <- c()
    if(nT <= 40){
      for(k in 1:length(my_survey_traps)){
        survey_res <- rbind(survey_res, c(error = j / b_ac, nT = nT, b_ac = b_ac, beta0 = beta0, sigma = sigma, varying = 3, method = 4, k = k,
                                          Enrm(D = Dcov_for_sim_ha, traps = read.traps(data = my_survey_traps[[k]], detector = "count"), mask = mask, 
                                               detectpar = list(lambda0 = exp(beta0), sigma = sigma), noccasions = 1)))
      }
    }
    
    Enrm_results <- rbind(Enrm_results, opt_res, grid_res, opt_grid_res, survey_res)
    
  }
  
}

save(Enrm_results, file = "output/posthoc-check-robustness-to-misspec.RData")
