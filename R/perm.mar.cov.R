#' Score Statistics Permutation
#'
#' @param Y A numeric vector for continuous or binary responses
#' @param X1 The matrix of unweighted taxon proportion for all non-rooted nodes.
#' @param X2 The matrix of weighted taxon proportion for all non-rooted nodes.
#' @param X3 The matrix of square-root transformed taxon proportion for all non-rooted nodes.
#' @param X4 The matrix of taxon proportion for tip nodes only.
#' @param binary The indicator of whether the response is binary.
#' @param cov The matrix/data frame of covariates.
#' @param nperm The nnumber of permutations.
#'
#' @return A list containing the following elements:
#'         \item{pval.1side1}{Lower-tail p-values of one-sided score tests using unweighted taxon proportion.}
#'         \item{weight1}{Weights for all non-rooted nodes in the test using unweighted taxon proportion.}
#'         \item{pval.1side2}{Lower-tail p-values of one-sided score tests using weighted taxon proportion.}
#'         \item{weight2}{Weights for all non-rooted nodes in the test using weighted taxon proportion.}
#'         \item{pval.1side3}{Lower-tail p-values of one-sided score tests using square-root transformed taxon proportion.}
#'         \item{weight3}{Weights for all non-rooted nodes in the test using square-root transformed taxon proportion.}
#'         \item{pval.1side4}{Lower-tail p-values of one-sided score tests using taxon proportion for tip nodes only.}
#'         \item{weight4}{Weights for all non-rooted nodes in the test using taxon proportion for tip nodes only.}
#'
#' @import statmod stats
#'
#' @export
#'
perm.mar.cov <- function(Y, X1, X2, X3, X4, binary = FALSE, cov = NULL, nperm) {
  n <- nrow(X1)
  if(is.null(cov)) cov <- rep(1, length(Y))

  ## fit the null model
  model <- ifelse(binary, "binomial", "gaussian")
  data1 <- data.frame(trait = Y, cov)
  null.model <- glm(trait ~ ., family = model, data = data1)

  ## regress each OTU on covariates
  K1 <- ncol(X1)
  Xres1 <- matrix(NA, n, K1)
  for (k in 1:K1) {
    data2 <- data.frame(response = X1[, k], cov)
    X.null <- glm(response ~ ., data = data2)
    Xhat <- fitted.values(X.null)
    Xres1[, k] <- X1[, k] - Xhat
  }
  weight1 <- apply(X1, 2, sd)
  weight1[apply(Xres1, 2, sd) == 0] <- 0

  K2 <- ncol(X2)
  Xres2 <- matrix(NA, n, K2)
  for (k in 1:K2) {
    data2 <- data.frame(response = X2[, k], cov)
    X.null <- glm(response ~ ., data = data2)
    Xhat <- fitted.values(X.null)
    Xres2[, k] <- X2[, k] - Xhat
  }
  weight2 <- apply(X2, 2, sd)
  weight2[apply(Xres2, 2, sd) == 0] <- 0

  K3 <- ncol(X3)
  Xres3 <- matrix(NA, n, K3)
  for (k in 1:K3) {
    data2 <- data.frame(response = X3[, k], cov)
    X.null <- glm(response ~ ., data = data2)
    Xhat <- fitted.values(X.null)
    Xres3[, k] <- X3[, k] - Xhat
  }
  weight3 <- apply(X3, 2, sd)
  weight3[apply(Xres3, 2, sd) == 0] <- 0

  K4 <- ncol(X4)
  Xres4 <- matrix(NA, n, K4)
  for (k in 1:K4) {
    data2 <- data.frame(response = X4[, k], cov)
    X.null <- glm(response ~ ., data = data2)
    Xhat <- fitted.values(X.null)
    Xres4[, k] <- X4[, k] - Xhat
  }
  weight4 <- apply(X4, 2, sd)
  weight4[apply(Xres4, 2, sd) == 0] <- 0

  ## calculate score test statistics U
  ## permute X residuals and calculate permutation Us
  ## calculate p-values in return
  order.perm <- replicate(nperm, sample(1:length(Y)))
  U1 <- glm.scoretest(null.model, Xres1)
  U.perm1 <- apply(order.perm, 2, function(x) glm.scoretest(null.model, Xres1[x, ]))
  Up1 <- cbind(U1, U.perm1)
  rm(U1, U.perm1, Xres1, data1, data2)
  gc()

  U2 <- glm.scoretest(null.model, Xres2)
  U.perm2 <- apply(order.perm, 2, function(x) glm.scoretest(null.model, Xres2[x, ]))
  Up2 <- cbind(U2, U.perm2)
  rm(U2, U.perm2, Xres2)
  gc()

  U3 <- glm.scoretest(null.model, Xres3)
  U.perm3 <- apply(order.perm, 2, function(x) glm.scoretest(null.model, Xres3[x, ]))
  Up3 <- cbind(U3, U.perm3)
  rm(U3, U.perm3, Xres3)
  gc()

  U4 <- glm.scoretest(null.model, Xres4)
  U.perm4 <- apply(order.perm, 2, function(x) glm.scoretest(null.model, Xres4[x, ]))
  Up4 <- cbind(U4, U.perm4)
  rm(U4, U.perm4, Xres4, order.perm, null.model)
  gc()

  return(list(pval.1side1 = pnorm(Up1), weight1 = weight1,
              pval.1side2 = pnorm(Up2), weight2 = weight2,
              pval.1side3 = pnorm(Up3), weight3 = weight3,
              pval.1side4 = pnorm(Up4), weight4 = weight4))
  ## pval.1side are the one-sided p-values; can be further used by AF_combine functions
}
