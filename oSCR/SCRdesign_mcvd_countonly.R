SCRdesign = function (statespace = NULL, all.traps = NULL, clusters = NULL, 
                            fix = NULL, clust.mids = clust.mids, ntraps = 9, ndesigns = 10, 
                            nn = 19, beta0 = -0.6, sigma = 2, N = 100, 
                            D_per_mask_cell = 1) 
{
  if (is.null(statespace))
    stop("Must supply a 'statespace' object (coordinates of the study area)")
  if (is.null(all.traps)) 
    stop("Must supply a 'statespace' object (coordinates of the study area)")
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
    Qinit <- Qfn(X = X, S = statespace, N = N, sigma = sigma, beta0 = beta0, 
                 D_per_mask_cell = D_per_mask_cell)
    Q <- max(Qinit[c(4,5)]) # max(-En,-Er))
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
            Qinit <- Qfn(X = Xtest, S = statespace, N = N, sigma = sigma, beta0 = beta0, 
                         D_per_mask_cell = D_per_mask_cell)
            qtest[j] <- max(Qinit[c(4,5)]) # max(-En,-Er)
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
