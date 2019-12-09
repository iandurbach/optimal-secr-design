Qfn <- function (X, S, N = 100, sigma, beta0, 
                 occasions, detector, D_per_mask_cell, ##
                 transitions = NULL, divby = 1, addback = 0)  ##
{
  gr <- S
  traps <- X
  ntraps <- nrow(X)
  pp <- e2dist(traps, gr) 
  # start noneuc bit
  if (!is.null(transitions)) {  ##
    lcd <- costDistance(transitions, X, S)  ##
    pp <- lcd / divby + addback ##
  } #
  # end noneuc bit
  pph <- c() ##
  for(i in 1:occasions){ ##
    p0 <- exp(beta0[i]) ##
    pph_s <- p0 * exp(-pp * pp/(2 * sigma * sigma)) ##
    pph <- rbind(pph, pph_s) ##
  } ##
  # mean detection probability
  pp <- 1 - exp(-pph) 
  notpp <- 1 - pp 
  notcap <- exp(colSums(log(notpp))) 
  pbar <- mean(1 - notcap)
  p0times <- 1 - pbar
  cap.trap.j <- 1 - (notpp)
  bb <- cap.trap.j/((1 - pp))
  bling <- matrix(notcap, ncol = length(notcap), nrow = ntraps * occasions, ##
                  byrow = TRUE) * bb
  p1time <- colSums(bling)
  p2bar <- mean(1 - p0times - p1time)
  # Murray's criterion: first captures and recaptures (all new below, except last line)
  firstcaps <- sum((1 - exp(-apply(pph, 2, sum))) * D_per_mask_cell)
  if(detector == "count"){
    recaps <- sum(apply(pph,2,sum) * D_per_mask_cell) - firstcaps
  } else if (detector == "multi") {
    recaps <- 0
    for(i in 1:occasions){
      start_row <- ntraps * (i-1) + 1
      end_row <- ntraps * i
      p_s <- 1 - exp(-apply(pph[start_row:end_row, ], 2, sum))
      recaps <- recaps + sum(p_s * D_per_mask_cell)
    }
    recaps <- recaps - firstcaps
  } else if (detector == "proximity") {
    recaps <- sum(t(t(pp) * D_per_mask_cell)) - firstcaps
  }
  c(1 - pbar, 1 - p2bar,  1 - (pbar + p2bar), -firstcaps, -recaps)
}
