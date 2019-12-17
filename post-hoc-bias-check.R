library(dplyr)
library(sf)
library(secr)
library(secrdesign)
library(oSCR)
library(stringr)
library(raster)
library(gdistance)

source("oSCR/Qfn_mcvd.R")

##### for designs assuming non-uniform density (uniform habitat use)

# load example designs 
load("output/Tost_examples_nonuniD.Rdata")  

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
lambda0 <- exp(beta0)
g0 <- 1-exp(-lambda0)
sigma <- 6000

results <- data.frame(i = as.integer(), j = as.integer(), 
                      E.N = as.numeric(), R.N = as.numeric(), E.Nlcl = as.numeric(), R.Nlcl = as.numeric(), 
                      E.Nucl = as.numeric(), R.Nucl = as.numeric(), true.N = as.numeric())

for(i in 1:length(traps)){
  
  # get the value of the covariate from the name of the list element
  b_ac <- as.numeric(str_split(names(traps)[i], pattern = "_", simplify = TRUE)[2])
  # make density on the mask a function of some covariate (b_ac = 0 for constant density)
  covariates(mask)$Dac <- exp(b_ac * as.numeric(scale(covariates(mask)$stdGC)))
  # now standardize the Dac so that its mean is = D_per_mask_cell, or equivalently that its sum is = Dpmc * N mask cells
  Dcov_for_sim <- covariates(mask)$Dac / sum(covariates(mask)$Dac) * (D_per_mask_cell * length(mask$x))
  # turn this into equivalent density per hectate, for use with Enrm (secr functions want per ha values)
  Dcov_for_sim_ha <- Dcov_for_sim / attr(mask, "area")
  
  my_traps <- read.traps(data = traps[[i]], detector = "count")
  
  for(j in 1:100){
    
    simulated_points_Dcov <- sim.popn(D = Dcov_for_sim_ha, 
                                      core = mask, 
                                      model2D = "IHP",
                                      Ndist = "fixed")
    
    ch <- sim.capthist(traps = my_traps, pop = simulated_points_Dcov, 
                       noccasions = 1,
                       detectpar = list(lambda0 = lambda0, sigma = sigma), 
                       detectfn = "HHN")
    summary(ch)
    
    # make starting values for secr.fit
    startvals <- list(D = mean(Dcov_for_sim_ha), lambda0 = lambda0, sigma = sigma)
    
    mod <- try(secr.fit(capthist = ch, mask = mask, model = list(D ~ stdGC), 
                        detectfn = "HHN", start = startvals))
    
    rmod <- try(region.N(mod))
    if((class(mod) != "try-error") & (class(rmod) != "try-error")){
      
      results <- rbind(results, 
                       data.frame(i = i, j = j, 
                                  E.N = region.N(mod)$estimate[1], R.N = region.N(mod)$estimate[2],
                                  E.Nlcl = region.N(mod)$lcl[1], R.Nlcl = region.N(mod)$lcl[2], 
                                  E.Nucl = region.N(mod)$ucl[1], R.Nucl = region.N(mod)$ucl[2], 
                                  true.N = nrow(simulated_points_Dcov),
                                  n = summary(ch)[[4]][1,1], detections = summary(ch)[[4]][6,1], dets_visited = summary(ch)[[4]][7,1]))
      
    }
    
  }
  
}

##### for designs assuming non-uniform habitat use (uniform density)

# load example designs 
load("output/Tost_examples_noneuc.Rdata")  

# when generating noneuc designs I standardised the LC distances by dividing by 
# a constant "divby" that made the distribution of the LC distances similar to the 
# distribution of the Euclidean distances, so that the same sigma parameter could be
# used in both. Now I want to use package secr, where you can only adjust sigma, so now
# I need to recalculate those divby factors (because I didn't store them) 

# user parameters
cellsize <- 2000
sigma <- 6000
beta0 <- -0.8  # beta0 = log(lambda0) = log(K * 'p0')
dens_per_100km2 <- 1 # mean animal density per 100km2, SLs are ~1

# load secr mask
load("data/Tost.RData")
mask <- TostMask 

# reduce resolution of mesh so have fewer possible camera locations
red_factor <- cellsize[1] / attr(mask, "spacing")
if ((trunc(red_factor) - red_factor) != 0) stop("Check spacing, causing non-integer reduction factor for mesh")
mask <- secr::raster(mask, "stdGC") %>% raster::aggregate(fact = red_factor, fun = mean) 
mask_df <- data.frame(coordinates(mask), stdGC = matrix(mask)) %>% filter(!is.na(stdGC))
mask <- read.mask(data = mask_df)

# mask and traps
statespace <- mask_df %>% st_as_sf(coords = c("x", "y"))  %>% st_coordinates() %>% as.matrix()
all.traps <- statespace

traps <- opt_traps_all %>% filter(trap_id == 1) %>%
  mutate(ID = paste0(nT, "_", alpha2)) %>% dplyr::select(ID,x,y)
ID <- traps$ID
traps <- split(traps[,-1], traps$ID)

opt_noneuc_EnEr <- c()
for(i in 1:100){
  
  # get the value of the covariate from the name of the list element
  alpha2 <- as.numeric(str_split(names(traps)[i], pattern = "_", simplify = TRUE)[2])
  
  D_per_mask_cell = dens_per_100km2 / 100 * (attr(mask, "area") / 100)
  noneuc_costs = covariates(mask)$stdGC
  
  cov <- rasterFromXYZ(cbind(statespace, noneuc_costs))
  cost <- exp(alpha2 * cov)
  trans <- transition(cost, transitionFunction=function(x) (mean(x)), direction = 16)
  trans <- geoCorrection(trans)
  # noneuc changes sigma and I want a way of automatically scaling lcd so that don't need to specify 
  # another sigma-noneuc par.
  # do this by calculate a factor that scales the lcd between all posstraps and all mask cells so that 
  # the dbn of lcds has the same minimum and pr[x>sigma] as the euclidean distances 
  pp_base <- e2dist(all.traps, statespace) # euc distance between all cells and traps
  lcd_base <- costDistance(trans, all.traps, statespace) 
  qe <- sum(pp_base < sigma) / length(pp_base)  # percentile of euc dist dbn = sigma
  pp0 <- pp_base - min(pp_base)
  lcd0 <- lcd_base - min(lcd_base)
  divby <- as.numeric(quantile(lcd0, probs = qe)) / as.numeric(quantile(pp0, probs = qe))
  addback <- min(pp_base) - (min(lcd_base) / divby) # scaled lcd <- (lcd - a) / divby + b where a = min(lcd) and b = min(pp)
  
  Qres <- Qfn(X = as.matrix(traps[[i]]), S = statespace, N = 100, sigma = sigma, beta0 = beta0, D_per_mask_cell = D_per_mask_cell, 
              transitions = trans, divby = divby, addback = addback, occasions = 1, detector = "count") 
  
  # 
  # # doesn't matter if you divide distance by divby or multiply sigma by divby, get same answer
  # Qfn(X = as.matrix(traps[[i]]), S = statespace, N = 100, sigma = sigma * divby, beta0 = beta0, D_per_mask_cell = D_per_mask_cell, 
  #      transitions = trans, divby = 1, addback = addback, occasions = 1, detector = "count) 
  # 
  opt_noneuc_EnEr <- rbind(opt_noneuc_EnEr, data.frame(i = i, divby = divby, addback = addback, En = -Qres[4], Er = -Qres[5]))
  
}

# mean function for non-euc distance calcs
mymean <- function(x) exp(mean(log(x))) 

# standard least cost distance function
# least cost distance function
myLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = mymean, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

results_noneuc <- data.frame(i = as.integer(), j = as.integer(), 
                             E.N = as.numeric(), R.N = as.numeric(), E.Nlcl = as.numeric(), R.Nlcl = as.numeric(), 
                             E.Nucl = as.numeric(), R.Nucl = as.numeric(), true.N = as.numeric(),
                             n = as.numeric(), detections = as.numeric(), dets_visited = as.numeric())

for(i in 1:length(traps)){
  
  # get the value of the covariate from the name of the list element
  b_con <- as.numeric(str_split(names(traps)[i], pattern = "_", simplify = TRUE)[2])
  # make conductance a function of some covariate 
  covariates(mask)$noneuc <- (b_con * as.numeric((covariates(mask)$stdGC)))
  
  my_traps <- read.traps(data = traps[[i]], detector = "count")
  
  for(j in 1:100){
    
    covariates(mask)$noneuc <- exp(b_con * as.numeric((covariates(mask)$stdGC)))
    
    simulated_points_Dcov <- sim.popn(D = D_per_ha, 
                                      core = mask, 
                                      model2D = "IHP",
                                      Ndist = "poisson")
    
    
    ch <- sim.capthist(traps = my_traps, 
                       pop = simulated_points_Dcov, 
                       userdist = myLCdist, 
                       noccasions = 1,
                       detectpar = list(lambda0 = lambda0, sigma = sigma * opt_noneuc_EnEr$divby[i]), 
                       detectfn = "HHN")
    
    # make starting values for secr.fit
    startvals <- list(D = D_per_ha, lambda0 = lambda0, sigma = sigma)
    
    mod <- secr.fit(capthist = ch, mask = mask, model = list(D ~ 1, noneuc ~ stdGC -1), 
                    detectfn = "HHN", start = startvals, details = list(userdist = myLCdist))
    
    summary(mod)
    region.N(mod)
    results_noneuc <- rbind(results_noneuc, 
                            data.frame(i = i, j = j, 
                                       E.N = region.N(mod)$estimate[1], R.N = region.N(mod)$estimate[2],
                                       E.Nlcl = region.N(mod)$lcl[1], R.Nlcl = region.N(mod)$lcl[2], 
                                       E.Nucl = region.N(mod)$ucl[1], R.Nucl = region.N(mod)$ucl[2], 
                                       true.N = nrow(simulated_points_Dcov),
                                       n = summary(ch)[[4]][1,1], detections = summary(ch)[[4]][6,1], dets_visited = summary(ch)[[4]][7,1]))
    
  }
  
}

# save
save(results, results_noneuc, file = "output/posthoc-check-bias-nonuniD-noneuc.RData")

