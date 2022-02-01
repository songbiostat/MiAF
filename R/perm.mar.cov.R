#' Score Statistics Permutation
#'
#' @param Y A numeric vector for continuous or binary responses
#' @param X The matrix of transformed taxon proportion.
#' @param binary The indicator of whether the response is binary.
#' @param cov The matrix/data frame of covariates.
#' @param nperm The nnumber of permutations.
#' @param seed Set seed.
#'
#' @return A list containing the following elements:
#'         \item{pval.1side1}{Lower-tail p-values of one-sided score tests.}
#'         \item{weight}{Standard deviation of transformed microbial abundance,
#'                      which will be used to calculate the weight.}
#'
#' @import statmod stats
#'
#' @export
#'
perm.mar.cov <- function(Y, X, binary = FALSE, cov = NULL, nperm, seed) {
  n <- nrow(X)
  if(is.null(cov)) cov <- rep(1, n)

  ## fit the null model
  model <- ifelse(binary, "binomial", "gaussian")
  null.model <- glm(trait ~ ., family = model, data = data.frame(trait = Y, cov))

  ## regress each OTU on covariates
  Xres <- resid(lm(X ~ ., data = data.frame(cov)))
  weight <- apply(X, 2, sd)
  weight[apply(Xres, 2, sd) == 0] <- 0

  ## calculate score test statistics U
  ## permute X residuals and calculate permutation Us
  ## calculate p-values in return
  U <- glm.scoretest(null.model, Xres)
  set.seed(seed)
  perm.id <- replicate(nperm, sample(1:n))
  U.perm <- apply(perm.id, 2, function(x) glm.scoretest(null.model, Xres[x, ]))
  Up <- cbind(U, U.perm)
  rm(U, U.perm, Xres, X)
  gc()

  return(list(pval.1side = pnorm(Up), weight = weight))
  ## pval.1side are the one-sided p-values; can be further used by AF_combine functions
}
