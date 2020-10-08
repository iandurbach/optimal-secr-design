SCRdesignEAenrm = function (statespace = NULL, all.traps = NULL, clusters = NULL, 
                      fix = NULL, clust.mids = clust.mids, ntraps = 9, ndesigns = 10, 
                      nn = 19, lambda0 = exp(-0.6), sigma = 2, N = 100, 
                      occasions = 1, detector = NULL, D_per_mask_cell = 1, 
                      noneuc_cov = NULL, alpha2 = 1, crit = 8, use.secr = FALSE,
                      pen_wt = 0, pen_gridsigma = 2)  
{
  g_n234 <- c(0,0,0)
  if(use.secr){
    statespace <- read.mask(data = data.frame(statespace))
    
    #pop <- sim.popn(D = D_per_mask_cell, core = statespace)
    
    # find distribution of trap spacings on a close to regular grid, to ensure later optimized grid has spaced
    # enough detectors sufficiently far apart to get low var(sigma)
    
    # hacky way to make a polygon just bigger than boundary of possible detector points
    sp <- diff(sort(unique(all.traps[,1])))[1] 
    pg <- data.frame(all.traps) %>% st_as_sf(coords = c("x", "y")) %>% st_buffer(dist = sp, endCapStyle = "SQUARE") %>% st_union() 
    pg <- pg %>% st_buffer(dist = -sp * 0.4, endCapStyle = "SQUARE")
    
    # place a grid over the area, with cells X sigma apart
    my_grid <- st_make_grid(pg, cellsize = c(pen_gridsigma * sigma, pen_gridsigma * sigma), what = "centers") %>% st_intersection(pg)
    grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
    
    # choose a subset of nT of these, starting from a random trap and then choosing nT nearest neighbours
    all_grid_traps <- list()
    
    grid_traps <- grid_traps_full
    xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
    grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
    grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= ntraps)
    grid_traps <- read.traps(data = grid_traps, detector = "count")
    
    # find out how many detector pairs are between 2-3 and 3-4 sigma apart; these are important for sigma estimation         
    td <- e2dist(grid_traps, grid_traps) 
    n2ish <- sum((td <= 2.499 * sigma) & (td > 1.499 * sigma))
    n3ish <- sum((td <= 3.499 * sigma) & (td > 2.499 * sigma))
    n4ish <- sum((td <= 4.449 * sigma) & (td > 3.499 * sigma))
    
    # set reasonable targets for optimized grids
    g_n2ish <- 0 # n2ish * 0.8
    g_n3ish <- n3ish * 1
    g_n4ish <- n4ish * 1
    
    g_n234 <- c(g_n2ish, g_n3ish, g_n4ish)
    
    Qfn <- function (X, S, N = 100, sigma, lambda0, 
                     occasions, detector, D_per_mask_cell, ##
                     transitions = NULL, divby = 1, addback = 0, g_n234, pen_wt)  ##
    {
      # penalty for too clustered
       td <- e2dist(X, X) 
       n2ish <- sum((td <= 2.499 * sigma) & (td > 1.499 * sigma))
       n3ish <- sum((td <= 3.499 * sigma) & (td > 2.499 * sigma))
       n4ish <- sum((td <= 4.449 * sigma) & (td > 3.499 * sigma))
       pen <- pen_wt * (max(0, g_n234[1] - n2ish) + max(0, g_n234[2] - n3ish) + max(0, g_n234[3] - n4ish))
      traps <- read.traps(data = data.frame(X), detector = detector)
      enrm <- Enrm(D = D_per_mask_cell, traps = traps, mask = S, noccasions = occasions, detectpar = list(lambda0 = lambda0[1], sigma = sigma))
      # ch <- sim.capthist(traps = traps, pop = pop, noccasions = 5, detectpar = list(lambda0 = lambda0, sigma = sigma), detectfn = 14)
      # mod <- secr.fit(capthist = ch, mask = masks, model = list(D ~ 1), method = "none",
      #                 detectfn = "HHN", start = c(log(D), log(lambda0), log(sigma)), trace = FALSE)
      # 
      # cva <- derived(mod)[2,6]
      # if(is.na(cva)|is.nan(cva)){ cva <- 2}
      # c(0, 0,  0, -enrm[1], -enrm[3], -enrm[2], sqrt(1/enrm[1] + cva^2))
      c(0, 0,  0, -enrm[1], -enrm[2], pen-min(enrm[1],enrm[2]))
    }
  }
  if (is.null(statespace))
    stop("Must supply a 'statespace' object (coordinates of the study area)")
  if (is.null(all.traps)) 
    stop("Must supply a 'statespace' object (coordinates of the study area)")
  if (is.null(detector))
    stop("Must supply a 'detector' object")
  if (!(detector %in% c("count", "multi", "proximity")))
    stop("Detector must be one of 'count', 'multi', 'proximity'")
  if((occasions>1) & (length(lambda0)==1)){
    lambda0 <- rep(lambda0, occasions)
    warning("Single lambda0 value given for >1 occasion, assuming lambda0 is per occasion")
  }
  # if noneuc then work out transition matrix with conductances
  if(!is.null(noneuc_cov)){
    cov <- rasterFromXYZ(cbind(statespace, noneuc_cov))
    cost <- exp(alpha2 * cov)
    trans <- transition(cost, transitionFunction=function(x) (mean(x)), direction = 16)
    trans <- geoCorrection(trans)
    # noneuc changes sigma so automatically scaling lcd so that don't need to specify 
    # another sigma-noneuc par.
    # do this by calculating a factor that scales the lcd between all possibe detectors and all 
    # mask cells so that the dbn of lcds has the same minimum and pr[x>sigma] as the euclidean distances 
    pp_base <- e2dist(all.traps, statespace) # euc distance between all cells and traps
    lcd_base <- costDistance(trans, all.traps, statespace) 
    qe <- sum(pp_base < sigma) / length(pp_base)  # percentile of euc dist dbn = sigma
    pp0 <- pp_base - min(pp_base)
    lcd0 <- lcd_base - min(lcd_base)
    divby <- as.numeric(quantile(lcd0, probs = qe)) / as.numeric(quantile(pp0, probs = qe))
    addback <- min(pp_base) - (min(lcd_base) / divby) # scaled lcd <- (lcd - a) / divby + b where a = min(lcd) and b = min(pp)
  } else { 
    trans <- NULL
    divby <- 1
    addback <- 0
  }
  ngrid <- nrow(statespace)
  Dlist <- list()
  Qhistory <- NULL
  Qdesign <- NULL
  if (is.null(clusters)) {
    Cd <- round(e2dist(all.traps, all.traps), 8)
    NN2 <- NN <- matrix(0, nrow = nrow(Cd), ncol = ncol(Cd))
    for (i in 1:nrow(Cd)) {
      xx <- Cd[i, ]
      NN[i, ] <- (xx > 0 & xx <= sort(xx)[nn])
      NN2[i, ] <- (xx > 0 & xx <= sort(xx)[3])
    }
  }
  if (!is.null(clusters)) {
    Cd <- round(e2dist(clust.mids[, -1], clust.mids[, -1]), 
                8)
    NN2 <- NN <- matrix(0, nrow = nrow(Cd), ncol = ncol(Cd))
    for (i in 1:nrow(Cd)) {
      xx <- Cd[i, ]
      NN[i, ] <- (xx > 0 & xx <= sort(xx)[nn])
      NN2[i, ] <- (xx > 0 & xx <= sort(xx)[3])
    }
  }
  if (is.null(fix)) {
    X <- all.traps
  }
  if (!is.null(fix)) {
    X <- rbind(all.traps, fix)
  }
  for (m in 1:ndesigns) {
    Qbest <- 10^10
    if (is.null(clusters)) {
      X.current <- sample(1:nrow(all.traps), ntraps)
      if (is.null(fix)) {
        X <- all.traps[X.current, ]
      }
      else {
        X <- rbind(all.traps[X.current, ], fix)
      }
    }
    else {
      X.current <- sample(1:nrow(clusters), ntraps)
      which.sites <- apply(clusters[X.current, ], 2, sum) > 
        0
      if (is.null(fix)) {
        X <- all.traps[which.sites, ]
      }
      else {
        X <- rbind(all.traps[which.sites, ], fix)
      }
    }
    # new Qfn (does En, Er, and outputs everything)
    Qinit <- Qfn(X = X, S = statespace, N = N, sigma = sigma, 
                 lambda0 = lambda0, 
                 occasions = occasions, detector = detector, D_per_mask_cell = D_per_mask_cell, 
                 transitions = trans, divby = divby, addback = addback, g_n234 = g_n234, pen_wt = pen_wt)
    Qinit <- c(Qinit, max(Qinit[c(4,5)]))
    Q <- Qinit[crit] #max(Qinit[c(4,5)]) # max(-En,-Er)
    Qhistory <- c(Qhistory, Q)
    Qdesign <- c(Qdesign, m)
    cat("Initial Q: ", Qinit, fill = TRUE)
    if (is.nan(Q)) {
      Dlist[[m]] <- list(Q = NA, X = X, X.current = X.current)
      next
    }
    repeat {
      for (i in 1:ntraps) {
        chk <- NN[X.current[i], ]
        chk[X.current] <- 0
        x.consider <- (1:ncol(NN))[chk == 1]
        qtest <- rep(10^10, length(x.consider))
        if (length(x.consider) > 0) {
          for (j in 1:length(x.consider)) {
            Xtest.current <- X.current
            Xtest.current[i] <- x.consider[j]
            if (!is.null(clusters)) {
              which.test <- apply(clusters[Xtest.current, 
                                           ], 2, sum) > 0
              if (is.null(fix)) {
                Xtest <- all.traps[which.test, ]
              }
              else {
                Xtest <- rbind(all.traps[which.test, 
                                         ], fix)
              }
            }
            else {
              if (is.null(fix)) {
                Xtest <- all.traps[Xtest.current, ]
              }
              else {
                Xtest <- rbind(all.traps[Xtest.current, 
                                         ], fix)
              }
            }
            # new Qfn (does En, Er, and outputs everything)
            Qinit <- Qfn(X = Xtest, S = statespace, N = N, 
                         sigma = sigma, lambda0 = lambda0, 
                         occasions = occasions, detector = detector, D_per_mask_cell = D_per_mask_cell, 
                         transitions = trans, divby = divby, addback = addback, g_n234 = g_n234, pen_wt = pen_wt)
            Qinit <- c(Qinit, max(Qinit[c(4,5)]))
            qtest[j] <- Qinit[crit] #max(Qinit[c(4,5)]) # max(-En,-Er)
          }
        }
        else {
          qtest <- NaN
        }
        if (any(is.nan(qtest))) {
          Dlist[[m]] <- list(Q = NA, X = X, X.current = X.current)
          next
        }
        if (min(qtest) < Q) {
          Q <- min(qtest)
          kp <- qtest == min(qtest)
          X.current[i] <- x.consider[kp][1]
          if (is.null(clusters)) {
            X <- all.traps[X.current, ]
          }
          else {
            which.sites <- apply(clusters[X.current, 
                                          ], 2, sum) > 0
            X <- all.traps[which.sites, ]
          }
        }
      }
      if (Qbest == Q) {
        break
      }
      if (Q < Qbest) 
        Qbest <- Q
      if (Q > Qbest) 
        cat("ERROR", fill = TRUE)
      Qhistory <- c(Qhistory, Q)
      Qdesign <- c(Qdesign, m)
      cat("Q: ", Qinit, fill = TRUE)
    }
    Dlist[[m]] <- list(Q = Qbest, X = X, X.current = X.current)
    m <- m + 1
  }
  Qvec <- rep(NA, length(Dlist))
  Xid <- matrix(NA, nrow = ntraps, ncol = length(Dlist))
  Xlst <- list()
  for (i in 1:length(Dlist)) {
    Qvec[i] <- Dlist[[i]]$Q
    Xid[, i] <- Dlist[[i]]$X.current
    Xlst[[i]] <- Dlist[[i]]$X
  }
  design.rank <- order(Qvec)
  tmp.xl <- list()
  tmp.Qh <- NULL
  tmp.Qi <- NULL
  for (i in 1:ndesigns) {
    tmp.xl[[i]] <- Xlst[[design.rank[i]]]
    tmp.Qh <- c(tmp.Qh, Qhistory[Qdesign == which(design.rank == 
                                                    i)])
    tmp.Qi <- c(tmp.Qi, Qdesign[Qdesign == which(design.rank == 
                                                   i)])
  }
  Qvec <- Qvec[design.rank]
  Xid <- Xid[, design.rank]
  Xlst <- tmp.xl
  Qdata <- data.frame(score = tmp.Qh, design = tmp.Qi)
  output <- list(Qvec = Qvec, Xid = Xid, Xlst = Xlst, Qdata = Qdata, 
                 all.traps = all.traps, statespace = statespace)
  return(output)
}