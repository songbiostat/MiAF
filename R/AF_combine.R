#' AF Combination of P-values
#'
#'@description This function combines p-values using adaptive Fisher method.
#'
#' @param p P-values to be combined. A matrix with dimenstions K by B.
#'          If an object of p-values from perm.mar.cov() is used, K is the number
#'          of OTUs, and B is the number of permutataion plus 1.
#'          The first column is observed marginal p-values,
#'          and the rest columns are p-values for permutation.
#' @param weight Weights given to the p-values.
#' @param log Indicator of whether p-values are on the log scale.
#'
#' @return A list containing the following elements:
#'         \item{pvalue}{P-value of AF statistic for observed sample (the first element)
#'                       and all permuted samples.}
#'         \item{which.selected}{The index of features to be selected.}
#'
#' @export
#'
AF_combine <- function(p, weight = 1, log = TRUE) {
  if(!log) {
    p <- log(p)
  }
  B <- ncol(p)
  wp <- p * weight
  s <- apply(apply(wp, 2, sort), 2, cumsum)
  ps.all <- t(apply(s, 1, rank, ties.method = "max"))/B
  rm(s)
  gc()
  AF.all <- apply(ps.all, 2, min)
  which.select <- order(wp[, 1])[1:which.min(ps.all[, 1])]
  p.all <- rank(AF.all, ties.method = "max")/B
  AF.list <- list(
    pvalue = p.all,
    which.selected = which.select
  )
  return(AF.list)
}
