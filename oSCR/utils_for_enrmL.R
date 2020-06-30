############################################################################################
# difference is significant only for large g0 
dfcast <- function (detectfn = 'HN', detectpar=list(g0 = 0.2, sigma = 25, z = NULL, w = NULL), matchsigma = 1, warning = TRUE) {
  
  if (!(detectfn %in% 14:18) ) {
    lambda0 <- -log(1- detectpar$g0)
    cast <- function (sigma2) {
      if (detectfn == 0)
        (detectpar$g0 * exp(-0.5 * matchsigma^2)) - 
        (1 - exp(- (lambda0 * exp(-0.5 * (matchsigma*detectpar$sigma)^2 / sigma2^2)))) 
      else if (detectfn == 1)
        (detectpar$g0 * (1 - exp(-matchsigma^-detectpar$z))) - 
        (1 - exp(- (lambda0 * (1 - exp(- (matchsigma*detectpar$sigma/sigma2)^-detectpar$z))))) 
      else if (detectfn == 2)
        (detectpar$g0 * exp(-matchsigma)) - 
        (1 - exp(- (lambda0 * exp(- (matchsigma*detectpar$sigma) / sigma2)))) 
      else stop ("invalid detectfn for dfcast")
    }
    detectpar <- list(lambda0 = lambda0, 
                      sigma = uniroot(cast, interval=c(0, detectpar$sigma))$root, 
                      z = detectpar$z,
                      w = detectpar$w)
    detectfn <- detectfn + 14   ## HN -> HHN, HR -> HHR, EX -> HEX
    if (warning)
      warning (call. = FALSE, "approximating detection function ",  secr:::.localstuff$DFN[detectfn+1], 
               paste0(" lambda0 = ", round(detectpar$lambda0,4), 
                      ", sigma = ", round(detectpar$sigma,1)))
  }
  return(list(detectfn = detectfn, detectpar = detectpar)) 
}

############################################################################################
defaultmodel <- function (CL, detectfn) {
  if (detectfn %in% c(0:8))
    model <- list(g0 = ~ 1, sigma = ~ 1)
  else if (detectfn %in% c(9))
    model <- list(b0 = ~ 1, b1 = ~ 1)
  else if (detectfn %in% c(10:11))
    model <- list(beta0 = ~ 1, beta1 = ~ 1)
  else ## detectfn %in% c(14:18))
    model <- list(lambda0 = ~ 1, sigma = ~ 1)
  if (!CL) model <- c(list(D = ~1), model)
  model
}
############################################################################################

replacedefaults <- function (default, user) replace(default, names(user), user)

##############################################################################

.local <- new.env()
.local$packageType <- ""  ## " pre-release"