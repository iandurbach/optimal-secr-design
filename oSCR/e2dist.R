# from oSCR package

e2dist <- function (x, y) {
  if (!is.matrix(x)) 
    x <- as.matrix(x)
  if (!is.matrix(y)) 
    y <- as.matrix(y)
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}