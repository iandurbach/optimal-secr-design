########################################################################################
## Efford, M. G. Efford, M. G. Fast evaluation of study designs for spatially explicit
##     capture--recapture. Ecological Applications.

## R code to compare simulated and approximate RSE(D-hat)

########################################################################################

runsim <- function (sp.ratio = seq(0.25,1.25,0.5), nx = 6, ny = 6, nrepeats = 1,
                    D = 10, lambda0 = 0.2, sigma = 20, nocc = 5, detector = 'proximity',
                    nrepl = 100, model.args = list(D ~ 1), seed = 123, ...) {
  
  Dorig <- D
  
  # mask is defined by the area around a 2sigma regular grid (since we won't look at spacings > 2s)
  grids_2s <- read.traps(data = expand.grid(x = seq(from = 500, by = 2 * sigma, length.out = 1.5*nx),
                                            y = seq(from = 500, by = 2 * sigma, length.out = 1.5*ny)), detector = detector)
  masks <- make.mask(grids_2s, spacing = sigma / 2, type = 'trapbuffer', buffer = 4*sigma)
  # add covariates
  covariates(masks) <- data.frame(D1 = D, 
                                  D2 = exp(0.5 * as.numeric(scale(sqrt(masks$x * masks$y)))),
                                  D3 = exp(-0.5 * as.numeric(scale(sqrt((masks$x-mean(masks$x))^2 + (masks$y-mean(masks$y))^2))))) 
  
  # define the "all points" grid -- demarcates possible trap locations
  all_allgrids <- list()
  for(i in 1:length(sp.ratio)){ 
    all_allgrids[[i]] <- read.traps(data = expand.grid(x = seq(from = min(grids_2s$x), to = max(grids_2s$x), by = sp.ratio[i] * sigma),
                                                       y = seq(from = min(grids_2s$y), to = max(grids_2s$y), by = sp.ratio[i] * sigma)), detector = detector)
  }
  
  # for each possible all points grid, find the optimal detector locations (a trap array)
  
  if(grepl("D ~ 1", model.args)) { D <- Dorig }
  if(grepl("D ~ D2", model.args)) { D <- covariates(masks)$D2 / sum(covariates(masks)$D2) * (Dorig * length(masks$x)) }
  if(grepl("D ~ D3", model.args)) { D <- covariates(masks)$D3 / sum(covariates(masks)$D3) * (Dorig * length(masks$x)) }
  
  # statespace = all mask points
  statespace <- as.matrix(masks[,1:2])
  # traps = potential camera locations
  all_allgrids_mat <- lapply(all_allgrids, as.matrix)
  
  # generate optimal designs
  all_grids <- lapply(all_allgrids_mat, 
                      SCRdesignEAenrm, 
                      statespace = as.matrix(masks[,1:2]),
                      #all.traps = grids[,1:2] %>% as.matrix(),
                      ntraps = nx * ny, # number of cameras available
                      ndesigns = 1, # number of random starting points
                      beta0 = log(lambda0),
                      sigma = sigma, 
                      D = D, # per mask cell not per ha!
                      occasions = nocc,
                      detector = detector,
                      crit = 6, 
                      use.secr = TRUE)
  
  
  sims <- list()
  for(i in 1:length(all_grids)){
    
    grids <- read.traps(data = as.data.frame(all_grids[[i]]$Xlst[[1]]), detector = detector)
    #masks <- all_masks[[i]]
    
    print(nrow(grids))
    set.seed(seed)
    
    # simulating directly from secr
    
    if(grepl("D ~ 1", model.args)) { D <- Dorig }
    if(grepl("D ~ D2", model.args)) { D <- covariates(masks)$D2 / sum(covariates(masks)$D2) * (Dorig * length(masks$x)) }
    if(grepl("D ~ D3", model.args)) { D <- covariates(masks)$D3 / sum(covariates(masks)$D3) * (Dorig * length(masks$x)) }
    
    mod_summary <- data.frame(j = as.integer(), 
                              D = as.numeric(), se.D = as.numeric(),
                              E.N = as.numeric(), R.N = as.numeric(), 
                              E.Nlcl = as.numeric(), R.Nlcl = as.numeric(), 
                              E.Nucl = as.numeric(), R.Nucl = as.numeric(), 
                              true.N = as.numeric(), n = as.numeric(),
                              detections = as.numeric(), 
                              dets_visited = as.numeric(),
                              esa = as.numeric(),
                              cvN = as.numeric(),
                              cva = as.numeric(),
                              cvD = as.numeric(),
                              maskarea = as.numeric(),
                              ncells_lt1s = as.numeric(), ncells_lt1.5s = as.numeric(),
                              ncells_lt2s = as.numeric(), ncells_lt4s = as.numeric())
    cnt <- 1
    for(j in 1:nrepl){ # run more when doing for real!
      
      if(length(D) > 1) { 
        my_pop <- sim.popn(D = D, core = masks, buffer = 0, model2D = "IHP") 
      } else { my_pop <- sim.popn(D = D, core = masks, buffer = 0) }
      
      ch <- sim.capthist(traps = grids, pop = my_pop, 
                         noccasions = nocc,
                         detectpar = list(lambda0 = lambda0, sigma = sigma), 
                         detectfn = "HHN")
      
      # make starting values for secr.fit
      startvals <- list(D = mean(D), lambda0 = lambda0, sigma = sigma)
      
      masks_grid <- make.mask(traps = grids, buffer = 4 * sigma, spacing = sigma / 2, type = "trapbuffer")

      mod <- try(secr.fit(capthist = ch, mask = masks_grid, model = model.args, 
                          detectfn = "HHN", start = startvals, trace = FALSE))
      
      if(class(mod) != "try-error"){
        
        rN <- try(region.N(mod))
        
        if(class(rN) != "try-error"){
          
          mod_summary <- rbind(mod_summary, 
                               data.frame(j = j, 
                                          D = coef(mod)[,"beta"][1],
                                          se.D = coef(mod)[,"SE.beta"][1],
                                          E.N = region.N(mod)$estimate[1], 
                                          R.N = region.N(mod)$estimate[2],
                                          E.Nlcl = region.N(mod)$lcl[1], 
                                          R.Nlcl = region.N(mod)$lcl[2], 
                                          E.Nucl = region.N(mod)$ucl[1], 
                                          R.Nucl = region.N(mod)$ucl[2], 
                                          true.N = nrow(my_pop),
                                          n = summary(ch)[[4]][2,"Total"], 
                                          detections = summary(ch)[[4]][6,"Total"], 
                                          dets_visited = summary(ch)[[4]][7,"Total"],
                                          esa = esa(mod)[1],
                                          cvN = derived(mod)[2,5],
                                          cva = derived(mod)[2,6],
                                          cvD = derived(mod)[2,7],
                                          maskarea = maskarea(masks),
                                          ncells_lt1s = nrow(make.mask(traps = grids, buffer = 1 * sigma, spacing = sigma / 2, type = "trapbuffer")),
                                          ncells_lt1.5s = nrow(make.mask(traps = grids, buffer = 1.5 * sigma, spacing = sigma / 2, type = "trapbuffer")),
                                          ncells_lt2s = nrow(make.mask(traps = grids, buffer = 2 * sigma, spacing = sigma / 2, type = "trapbuffer")),
                                          ncells_lt4s = nrow(masks_grid)))
          
          cnt <- cnt + 1
          
        }
        
        
      }
      
    }
    
    # can be fitting errors, remove any dodgy obs (this is subjective)
    # mod_summary <- mod_summary %>% filter(E.N < 4 * true.N)
    
    mod_summary$r <- mod_summary$detections - mod_summary$n
    
    # add theoretical values
    th <- Enrm(D = D, grids, masks, list(lambda0 = lambda0, sigma = sigma), noccasions = nocc)
    mod_summary$th_En <- th[1]
    mod_summary$th_Er <- th[2]
    mod_summary$th_Em <- th[3]
    mod_summary$th_cvD <- 1 / sqrt(min(th[1], th[2]))
    
    # add input variables
    mod_summary$nT <- nrow(grids)
    mod_summary$sigma <- sigma
    mod_summary$lambda0 <- lambda0
    mod_summary$dt <- detector
    mod_summary$sp.ratio <- summary(grids)$spacing / sigma
    
    sims[[i]] <- mod_summary
    
  }
  
  do.call(rbind, sims)
  
}

########################################################################################

## setup

library(secrdesign)
library(oSCR)
source("oSCR/SCRdesignEAenrm.R")
spr <- seq(0.25, 3.25, 0.5)
nrepl <- 100

#############################
## Simulations

## 2-D array
sim66_D1 <- runsim(sp.ratio = spr, nx = 6, ny = 6, nrepl = nrepl, model.args = list(D ~ 1), seed = 123)
sim66_D2 <- runsim(sp.ratio = spr, nx = 6, ny = 6, nrepl = nrepl, model.args = list(D ~ D2), seed = 123)
sim66_D3 <- runsim(sp.ratio = spr, nx = 6, ny = 6, nrepl = nrepl, model.args = list(D ~ D3), seed = 123)
sim88_D1 <- runsim(sp.ratio = spr, nx = 8, ny = 8, nrepl = nrepl, model.args = list(D ~ 1), seed = 123)
sim88_D2 <- runsim(sp.ratio = spr, nx = 8, ny = 8, nrepl = nrepl, model.args = list(D ~ D2), seed = 123)
sim88_D3 <- runsim(sp.ratio = spr, nx = 8, ny = 8, nrepl = nrepl, model.args = list(D ~ D3), seed = 123)
sim1010_D1 <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ 1), seed = 123)
sim1010_D2 <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ D2), seed = 123)
sim1010_D3 <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ D3), seed = 123)

sim1010.05_D2   <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ D2), lambda0 = 0.05, seed = 123)
sim1010.05_D3   <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ D3), lambda0 = 0.05, seed = 123)
sim1010.1_D2   <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ D2), lambda0 = 0.1, seed = 123)
sim1010.1_D3   <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ D3), lambda0 = 0.1, seed = 123)

simlist <- list(g66 = rbind(sim66_D1,sim66_D2,sim66_D3), 
                g88 = rbind(sim66_D1,sim88_D2,sim88_D3), 
                g1010 = rbind(sim66_D1,sim1010_D2,sim1010_D3),
                g1010.2 = rbind(sim1010.05_D2,sim1010.05_D3,sim1010.1_D3,sim1010.1_D3))

# save(simlist, file = "output/approxchecks/simlist-allgridspacing-fullsim.Rdata")