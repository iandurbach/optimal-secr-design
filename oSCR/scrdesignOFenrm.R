scrdesignOFenrm = function (v, alltraps, statespace, pop, N = 100, sigma, lambda0, 
                            occasions, detector, D, crit = 3, g_n234, pen_wt, sum_wt, ##
                            transitions = NULL, divby = 1, addback = 0) {
  
  # penalty for too clustered
  td <- e2dist(alltraps[v,], alltraps[v,]) 
  n2ish <- sum((td <= 2.499 * sigma) & (td > 1.499 * sigma))
  n3ish <- sum((td <= 3.499 * sigma) & (td > 2.499 * sigma))
  n4ish <- sum((td <= 4.449 * sigma) & (td > 3.499 * sigma))
  pen <- pen_wt * (max(0, g_n234[1] - n2ish) + max(0, g_n234[2] - n3ish) + max(0, g_n234[3] - n4ish))
  
  traps <- read.traps(data = data.frame(alltraps[v,]), detector = detector)
  
  if(length(lambda0) > 1){
    lambda0_mat <- matrix(lambda0[v], nrow = 1, ncol = length(v))
    enrm <- EnrmL(D = D, traps = traps, mask = statespace, noccasions = occasions, detectpar = list(lambda0 = lambda0_mat, sigma = sigma))
  } else {
    enrm <- Enrm(D = D, traps = traps, mask = statespace, noccasions = occasions, detectpar = list(lambda0 = lambda0[1], sigma = sigma))
  }
  
  c(-enrm[1], -enrm[3], pen-(min(enrm[1],enrm[2])))[crit]    
  
}


