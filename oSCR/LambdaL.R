##############################################################################
## package 'secrdesign'
## Lambda.R
## 2017-05-04
## 2017-07-15 getRSE renamed minnrRSE; optimalSpacing moved to optimalSpacing.R
## 2017-07-21 getdetectpar moved to getdetectpar.R
##############################################################################

LambdaL <- function (traps, mask, detectpar, noccasions, detectfn = 
                      c('HHN', 'HHR', 'HEX','HAN','HCG', 'HN', 'HR', 'EX')) {
  if (ms(traps) | ms(mask))
    stop ("Lambda does not accept multi-session traps or mask")
  if (is.character(detectfn))
    detectfn <- match.arg(detectfn)
  detectfn <- secr:::valid.detectfn(detectfn, valid = c(0,1,2,14:18))
  truncate <- ifelse(is.null(detectpar$truncate), 1e+10, detectpar$truncate)
  dfc <- dfcast (detectfn, detectpar)  # transforms detectfn 0 to 14, 2 to 16
  detectfn <- dfc$detectfn
  detectpar <- dfc$detectpar
  detectpars <- unlist(detectpar[secr:::parnames(detectfn)])
  
  temp <- .C("myLambdaC", 
             as.double(detectpars),
             as.integer(nrow(traps)),
             as.integer(nrow(mask)),
             as.double(unlist(traps)),
             as.double(unlist(mask)),
             as.integer(detectfn),
             L = double(nrow(mask)),
             resultcode = integer(1))
  if (temp$resultcode != 0)
    stop ("error in external function 'LambdaC'")
  L <- temp$L * noccasions
  covariates(mask) <- data.frame(Lambda = temp$L * noccasions)
  
  dettype <- secr:::detectorcode(traps, noccasions = noccasions)
  temp <- .C("mysumpkC", 
             as.integer(dettype),
             as.double(detectpars),
             as.integer(nrow(traps)),
             as.integer(nrow(mask)),
             as.double(unlist(traps)),
             as.double(unlist(mask)),
             as.integer(detectfn),
             sumpk = double(nrow(mask)),
             sumq2 = double(nrow(mask)),
             resultcode = integer(1))
  if (temp$resultcode != 0)
    stop ("error in external function 'sumpkc'")
  covariates(mask)$sumpk <- temp$sumpk * noccasions
  covariates(mask)$sumq2 <- temp$sumq2 
  return(mask)
}

EnrmL <- function (D, ...) {
  L <- LambdaL(...)
  allargs <- list(...)
  tr <- which(sapply(allargs, inherits, "traps"))
  detect <- detector(allargs[[tr]])[1]
  if (!detect %in% c("multi","proximity","count"))
    stop("only for 'multi', 'proximity' and 'count' detectors")
  cellarea <- attr(L, 'area')
  Lam <- covariates(L)$Lambda
  En <- sum( D * cellarea * (1 - exp(-Lam))) # expected n; multi, proximity, count
  if (detect %in% "count") { 
    Ec <- sum(D * cellarea * Lam )
    Em <- sum(D * cellarea * (covariates(L)$Lambda - (1-exp(-covariates(L)$Lambda))) * 
                (1-covariates(L)$sumq2))
  }
  else {   ## detect %in% c("multi", "proximity")
    Ec <- sum(D * cellarea * covariates(L)$sumpk )
    Em <- sum(D * cellarea * ((covariates(L)$sumpk - (1-exp(-covariates(L)$Lambda))) * 
                                (1-covariates(L)$sumq2)))
  }
  Er <- Ec - En        # expected r; count
  c(En = En, Er = Er, Em = Em)
}