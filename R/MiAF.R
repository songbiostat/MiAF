#' Microbiome Adaptive Fisher Test
#'
#' @param Y A numeric vector for continuous or binary responses.
#' @param X A matrix of the OTU table.
#'          Row - samples, column - OTUs.
#' @param tree A rooted phylogenetic tree of R class "phylo".
#' @param cov A data frame/matrix for covariate adjustment.
#'            Row - samples, column - covariates.
#' @param model "gaussian" is for the linear regression model (continuous);
#'              "binomial" is for the logistic regression model (binary).
#' @param n.perm The number of permutations.
#' @param seed set seed.
#'
#' @return A list containing the following elements:
#'         \item{UniFrac}{Unweighted UniFrac-like test p-value.}
#'         \item{wUniFrac}{Weighted UniFrac-like test p-value.}
#'         \item{UniFrac5}{Generalized UniFrac-like test p-value.}
#'         \item{tip.abun}{Leaf-nodes-only test p-value.}
#'         \item{com}{The p-value of MiAF test that combines p-values of the above four tests.}
#'         \item{select}{A list of selected associated taxa.
#'                       lower: under-presented taxa based on lower-tail p-values.
#'                       upper: over-presented taxa based on upper-tail p-values.}
#'
#' @export
#'
MiAF <- function (Y, X, tree, cov = NULL, model = c("gaussian","binomial"), n.perm, seed)
{
  model <- match.arg(model)
  ext.prop <- cum_prop(X, tree)
  cum <- ext.prop$cum
  cum.labels <- ext.prop$cum.labels
  br.len <- ext.prop$br.len

  ## weighted test
  if(is.null(seed)) seed <- sample(10000, 1)
  p.w <- perm.mar.cov(Y, t(cum), binary = (model == "binomial"), cov = cov, nperm = n.perm, seed = seed)
  weight.w <- c(br.len) * p.w$weight
  weight.w.nz <- weight.w[weight.w != 0]
  p1.w <- p.w$pval.1side[weight.w != 0, ]
  cur.wUniFrac.l <- AF_combine(p1.w, weight = weight.w.nz, log = FALSE)
  cur.wUniFrac.u <- AF_combine(1 - p1.w, weight = weight.w.nz, log = FALSE)
  rm(p.w, ext.prop, p1.w, weight.w.nz)
  gc()
  cur.wUniFrac <- AF_combine(rbind(cur.wUniFrac.l$pvalue, cur.wUniFrac.u$pvalue), weight = 1, log = FALSE)
  sel.w <- select.l.u(cur.wUniFrac.l$which.selected, cur.wUniFrac.u$which.selected, com = cur.wUniFrac$which.selected, cum.labels[weight.w != 0])
  rm(weight.w, cur.wUniFrac.l, cur.wUniFrac.u)
  gc()

  ## unweighted test
  cum.u <- cum
  cum.u[cum.u != 0] <- 1
  p.u <- perm.mar.cov(Y, t(cum.u), binary = (model == "binomial"), cov = cov, nperm = n.perm, seed = seed)
  weight.u <- c(br.len) * p.u$weight
  weight.u.nz <- weight.u[weight.u != 0]
  p1.u <- p.u$pval.1side[weight.u != 0, ]
  cur.UniFrac.l <- AF_combine(p1.u, weight = weight.u.nz, log = FALSE)
  cur.UniFrac.u <- AF_combine(1 - p1.u, weight = weight.u.nz, log = FALSE)
  rm(p.u, p1.u, weight.u.nz)
  gc()
  cur.UniFrac <- AF_combine(rbind(cur.UniFrac.l$pvalue, cur.UniFrac.u$pvalue), weight = 1, log = FALSE)
  sel.u <- select.l.u(cur.UniFrac.l$which.selected, cur.UniFrac.u$which.selected, com = cur.UniFrac$which.selected, cum.labels[weight.u != 0])
  rm(cur.UniFrac.l, cur.UniFrac.u, weight.u)
  gc()

  ## square-root transformation test
  cum.5 <- sqrt(cum)
  p.5 <- perm.mar.cov(Y, t(cum.5), binary = (model == "binomial"), cov = cov, nperm = n.perm, seed = seed)
  weight.5 <- c(br.len) * p.5$weight
  weight.5.nz <- weight.5[weight.5 != 0]
  p1.5 <- p.5$pval.1side[weight.5 != 0, ]
  cur.UniFrac5.l <- AF_combine(p1.5, weight = weight.5.nz, log = FALSE)
  cur.UniFrac5.u <- AF_combine(1 - p1.5, weight = weight.5.nz, log = FALSE)
  rm(p.5, cum, p1.5, weight.5.nz)
  gc()
  cur.UniFrac5 <- AF_combine(rbind(cur.UniFrac5.l$pvalue, cur.UniFrac5.u$pvalue), weight = 1, log = FALSE)
  sel.5 <- select.l.u(cur.UniFrac5.l$which.selected, cur.UniFrac5.u$which.selected, com = cur.UniFrac5$which.selected, cum.labels[weight.5 != 0])
  rm(cur.UniFrac5.l, cur.UniFrac5.u, weight.5)
  gc()

  # leaf-nodes-only test
  tip.cum <- X/rowSums(X) #only tip
  p.tip <- perm.mar.cov(Y, tip.cum, binary = (model == "binomial"), cov = cov, nperm = n.perm, seed = seed)
  weight.tip <- p.tip$weight
  weight.tip.nz <- weight.tip[weight.tip != 0]
  p1.tip <- p.tip$pval.1side[weight.tip != 0, ]
  cur.tip.l <- AF_combine(p1.tip, weight = weight.tip.nz, log = FALSE)
  cur.tip.u <- AF_combine(1 - p1.tip, weight = weight.tip.nz, log = FALSE)
  rm(p.tip, p1.tip, weight.tip.nz)
  gc()
  cur.tip <- AF_combine(rbind(cur.tip.l$pvalue, cur.tip.u$pvalue), weight = 1, log = FALSE)
  sel.tip <- select.l.u(cur.tip.l$which.selected, cur.tip.u$which.selected, com = cur.tip$which.selected, colnames(X))
  rm(cur.tip.l, cur.tip.u, weight.tip)
  gc()


  ## ultimate combining
  AFcom <- AF_combine(rbind(cur.UniFrac$pvalue, cur.wUniFrac$pvalue, cur.UniFrac5$pvalue, cur.tip$pvalue), weight = 1, log = FALSE)
  sel <- select.com(sel.u, sel.w, sel.5, sel.tip, AFcom$which.selected)

  return(list(UniFrac = cur.UniFrac$pvalue[1],
              wUniFrac = cur.wUniFrac$pvalue[1],
              UniFrac5 = cur.UniFrac5$pvalue[1],
              tip.abun = cur.tip$pvalue[1],
              com = AFcom$pvalue[1],
              select = sel))
}
