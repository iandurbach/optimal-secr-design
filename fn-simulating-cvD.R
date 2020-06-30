simulating_cvD = function(grids, masks, D, lambda0, sigma, noccasions, model.args = list(D ~ 1), my_pop = NULL, new_pop_each_rep = TRUE, sim.popn.Ndist = "poisson", secr.binomN = 0, secr.dbn = "poisson", nrepl, seed = NULL){
  
  if(!is.null(seed)){ set.seed(seed = seed) }
  
  if(length(D) > 1) { covariates(masks)$Dcov <- D }
  
  if(is.null(my_pop) & !new_pop_each_rep){
    if(length(D) > 1) {
      my_pop <- sim.popn(D = D, core = masks, buffer = 0, model2D = "IHP", Ndist = sim.popn.Ndist)
    } else { my_pop <- sim.popn(D = D, core = masks, buffer = 0, Ndist = sim.popn.Ndist) }
  }
  
  # simulating directly from secr
  all_mod <- list()
  all_ch <- list()
  all_pop <- list()
  mod_summary <- data.frame(j = as.integer(), 
                            D = as.numeric(), se.D = as.numeric(),
                            E.N = as.numeric(), R.N = as.numeric(), 
                            E.Nlcl = as.numeric(), R.Nlcl = as.numeric(), 
                            E.Nucl = as.numeric(), R.Nucl = as.numeric(), 
                            true.N = as.numeric(), n = as.numeric(),
                            m = as.numeric(),
                            detections = as.numeric(), 
                            dets_visited = as.numeric(),
                            esa = as.numeric(),
                            cvN = as.numeric(),
                            cva = as.numeric(),
                            cvD = as.numeric(),
                            maskarea = as.numeric())
  
  cnt <- 1
  for(j in 1:nrepl){ # run more when doing for real!
    
    if(j %% 100 == 0) { print(j) }
    
    if(new_pop_each_rep){
      if(length(D) > 1) {
        my_pop <- sim.popn(D = D, core = masks, buffer = 0, model2D = "IHP", Ndist = sim.popn.Ndist)
      } else { my_pop <- sim.popn(D = D, core = masks, buffer = 0, Ndist = sim.popn.Ndist) }
    }
    
    #plot(masks_red)
    #plot(grids, add=T)
    #plot(my_pop, add=T)
    
    ch <- sim.capthist(traps = grids, pop = my_pop, 
                       noccasions = noccasions,
                       detectpar = list(lambda0 = lambda0, sigma = sigma), 
                       detectfn = "HHN")
    
    # make starting values for secr.fit
    startvals <- list(D = mean(D), lambda0 = mean(lambda0), sigma = sigma)
    
    # NOTE REDUCE MASK TO SAVE TIME
    if((length(D) == 1) & (length(lambda0) == 1)){
      masks_red <- make.mask(traps = grids, buffer = 4 * sigma, spacing = sigma / 2, type = 'trapbuffer')
    } else { masks_red <- masks }
    
    mod <- try(secr.fit(capthist = ch, mask = masks_red, model = model.args, binomN = secr.binomN,
                        detectfn = "HHN", start = startvals, trace = FALSE, details = list(distribution = secr.dbn)))
    
    if(class(mod) != "try-error"){
      
      rN <- try(region.N(mod))
      
      if(class(rN) != "try-error"){
        
        all_pop[[cnt]] <- my_pop
        all_mod[[cnt]] <- mod
        all_ch[[cnt]] <- ch
        mod_summary <- rbind(mod_summary, 
                             data.frame(j = j, 
                                        D = coef(mod)[,"beta"][1],
                                        se.D = coef(mod)[,"SE.beta"][1],
                                        E.N = region.N(mod, region = masks)$estimate[1], 
                                        R.N = region.N(mod, region = masks)$estimate[2],
                                        E.Nlcl = region.N(mod, region = masks)$lcl[1], 
                                        R.Nlcl = region.N(mod, region = masks)$lcl[2], 
                                        E.Nucl = region.N(mod, region = masks)$ucl[1], 
                                        R.Nucl = region.N(mod, region = masks)$ucl[2], 
                                        true.N = nrow(my_pop),
                                        n = summary(ch)[[4]][2,"Total"], 
                                        moves = summary(ch, moves = TRUE, terse = TRUE)[5],
                                        detections = summary(ch)[[4]][6,"Total"], 
                                        dets_visited = summary(ch)[[4]][7,"Total"],
                                        esa = esa(mod)[1],
                                        cvN = derived(mod)[2,5],
                                        cva = derived(mod)[2,6],
                                        cvD = derived(mod)[2,7],
                                        maskarea = maskarea(masks)))
        
        cnt <- cnt + 1
        
      }
      
      
    }
    
  }
  
  mod_summary$r <- mod_summary$detections - mod_summary$n
  
  return(list(all_mod = all_mod, all_ch = all_ch, all_pop = all_pop, mod_summary = mod_summary))
}