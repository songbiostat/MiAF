#' AF Combination of P-values
#'
#'@description This function combines p-values using adaptive Fisher method.
#'
#' @param p P-values to be combined.
#' @param p.perm A permutation set of p-values.
#'               A matrix with dimensions length(p) by B (the number of permutations).
#' @param weight Weights given to the p-values.
#' @param log Indicator of whether p-values are on the log scale.
#'
#' @return A list containing the following elements:
#'         \item{pvalue}{P-value of AF statistic.}
#'         \item{stat}{AF statistic.}
#'         \item{p.perm}{P-values of AF statistics for all permuted samples.}
#'         \item{stat.perm}{AF statistics for all permuted samples.}
#'         \item{n.selected}{The number of features to be selected.}
#'         \item{which.selected}{The index of features to be selected.}
#'
#' @export
#'
AF_combine <- function(p, p.perm, weight = 1, log = TRUE) {
  if(!log) {
    p <- log(p)
    p.perm <- log(p.perm)
  }
  B <- ncol(p.perm)
  wp <- p * weight
  wp.perm <- p.perm * weight
  wp.sort <- sort(wp)
  wp.perm.sort <- apply(wp.perm, 2, sort)
  s <- cumsum(wp.sort)
  s.perm <- apply(wp.perm.sort, 2, cumsum)
  s.all <- cbind(s, s.perm)
  ps.all <- t(apply(s.all, 1, rank, ties.method = "max"))/(B + 1)
  AF.all <- apply(ps.all, 2, min)
  AF.perm <- AF.all[2:(B + 1)]
  AF <- AF.all[1]
  n.select <- which.min(ps.all[, 1])
  which.select <- order(wp)[1:n.select]
  p.all <- rank(AF.all, ties.method = "max")/(B + 1)
  p.AF <- p.all[1]
  p.AF.perm <- p.all[2:(B + 1)]
  AF.list <- list(
    pvalue = p.AF,
    stat = AF,
    p.perm = p.AF.perm,
    stat.perm = AF.perm,
    n.selected = n.select,
    which.selected = which.select
  )
  return(AF.list)
}
