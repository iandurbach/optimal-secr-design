Qfn <- function (X, S, N = 100, sigma, beta0, 
                 D_per_mask_cell, ##
                 transitions = NULL, divby = 1, addback = 0)  ##
{
  p0 <- exp(beta0)
  gr <- S
  traps <- X
  ntraps <- nrow(X)
  pp <- e2dist(traps, gr) 
  pph <- p0 * exp(-pp * pp/(2 * sigma * sigma)) # previously pp then overwritten, but need for MC
  pp <- 1 - exp(-pph) 
  notpp <- 1 - pp 
  notcap <- exp(colSums(log(notpp))) 
  pbar <- mean(1 - notcap)
  p0times <- 1 - pbar
  cap.trap.j <- 1 - (notpp)
  bb <- cap.trap.j/((1 - pp))
  bling <- matrix(notcap, ncol = length(notcap), nrow = ntraps, 
                  byrow = TRUE) * bb
  p1time <- colSums(bling)
  p2bar <- mean(1 - p0times - p1time)
  # Murray's criterion: first captures and recaptures
  firstcaps <- sum((1 - exp(-apply(pph, 2, sum))) * D_per_mask_cell) ##
  recaps <- sum(apply(pph,2,sum) * D_per_mask_cell) - firstcaps ##
  c(1 - pbar, 1 - p2bar,  1 - (pbar + p2bar), -firstcaps, -recaps)
}
