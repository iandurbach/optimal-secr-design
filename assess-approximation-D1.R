### Comparing simulated and approximate CV(D) for two spatially-varying density surfaces
### Code is adapted from Efford and Boulanger (2019) MEE, for spatially-varying D.

runsim <- function (sp.ratio = seq(0.25,3.25,0.25), nx = 6, ny = 6, nrepeats = 1,
                    D = 10, lambda0 = 0.2, sigma = 20, nocc = 5, detector = 'proximity',
                    nrepl = 100, model.args = list(D ~ 1), seed = 123, ...) {
  
  all_grids <- lapply(sp.ratio * sigma, make.grid, nx = nx, ny = ny, detector = detector, originxy = c(500,500))
  all_masks <- lapply(all_grids, make.mask, spacing = sigma/2, type = 'trapbuffer',  buffer = 4*sigma)
  
  # add covariates
  for(i in 1:length(all_masks)){
    covariates(all_masks[[i]]) <- data.frame(D1 = D, 
                                         D2 = exp(0.5 * as.numeric(scale(sqrt(all_masks[[i]]$x * all_masks[[i]]$y)))),
                                         D3 = exp(-0.5 * as.numeric(scale(sqrt((all_masks[[i]]$x-mean(all_masks[[i]]$x))^2 + (all_masks[[i]]$y-mean(all_masks[[i]]$y))^2))))) 
  }
  
  sims <- list()
  Dorig <- D
  for(i in 1:length(all_masks)){
    
    grids <- all_grids[[i]]
    masks <- all_masks[[i]]
    
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
                              dets_visited = as.numeric())
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
      
      mod <- try(secr.fit(capthist = ch, mask = masks, model = model.args, 
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
                                          dets_visited = summary(ch)[[4]][7,"Total"]))
          
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
spr <- seq(0.25, 3.25, 0.5)
nrepl <- 100

#############################
## Simulations

## Uniform D
sim66_D1 <- runsim(sp.ratio = spr, nx = 6, ny = 6, nrepl = nrepl, model.args = list(D ~ 1), seed = 123)
sim88_D1 <- runsim(sp.ratio = spr, nx = 8, ny = 8, nrepl = nrepl, model.args = list(D ~ 1), seed = 123)
sim1010_D1 <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ 1), seed = 123)
sim1010.05_D1 <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ 1), lambda0 = 0.05, seed = 123)
sim1010.1_D1 <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ 1), lambda0 = 0.1, seed = 123)

## Spatially-varying D
sim66_D2 <- runsim(sp.ratio = spr, nx = 6, ny = 6, nrepl = nrepl, model.args = list(D ~ D2), seed = 123)
sim66_D3 <- runsim(sp.ratio = spr, nx = 6, ny = 6, nrepl = nrepl, model.args = list(D ~ D3), seed = 123)
sim88_D2 <- runsim(sp.ratio = spr, nx = 8, ny = 8, nrepl = nrepl, model.args = list(D ~ D2), seed = 123)
sim88_D3 <- runsim(sp.ratio = spr, nx = 8, ny = 8, nrepl = nrepl, model.args = list(D ~ D3), seed = 123)
sim1010_D2 <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ D2), seed = 123)
sim1010_D3 <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ D3), seed = 123)
sim1010.05_D2   <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ D2), lambda0 = 0.05, seed = 123)
sim1010.05_D3   <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ D3), lambda0 = 0.05, seed = 123)
sim1010.1_D2   <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ D2), lambda0 = 0.1, seed = 123)
sim1010.1_D3   <- runsim(sp.ratio = spr, nx = 10, ny = 10, nrepl = nrepl, model.args = list(D ~ D3), lambda0 = 0.1, seed = 123)

simlist <- list(g66 = rbind(sim66_D1,sim66_D2,sim66_D3), 
                g88 = rbind(sim88_D1,sim88_D2,sim88_D3), 
                g1010 = rbind(sim1010_D1,sim1010_D2,sim1010_D3), 
                g1010.1 = rbind(sim1010.1_D1,sim1010.1_D2,sim1010.1_D3), 
                g1010.05 = rbind(sim1010.05_D1,sim1010.05_D2,sim1010.05_D3))

#save(simlist, file = "output/approxchecks/simlist-fullsim.Rdata")

