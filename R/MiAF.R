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
MiAF <- function (Y, X, tree, cov = NULL, model = c("gaussian","binomial"), n.perm)
{
  model <- match.arg(model)
  ext.prop <- cum_prop(X, tree)
  cum <- ext.prop$cum
  cum.labels <- ext.prop$cum.labels
  br.len <- ext.prop$br.len
  br.len <- as.matrix(br.len)
  cum.u <- cum
  cum.u[cum.u != 0] <- 1 #unweighted
  cum.5 <- sqrt(cum) # 0.5-power
  lib.size <- rowSums(X)
  tip.cum <- X/lib.size #only tip

  p.u.w <- perm.mar.cov(Y, t(cum.u), t(cum), t(cum.5), tip.cum, binary = (model == "binomial"), cov = cov, nperm = n.perm)
  rm(ext.prop, cum, cum.u, cum.5, tip.cum)
  gc()

  weight.u <- c(br.len) * p.u.w$weight1
  weight.u.nz <- weight.u[weight.u != 0]
  label.u <- cum.labels[weight.u != 0]
  p1.u <- p.u.w$pval.1side1[, 1][weight.u != 0]
  p1.u.perm <- p.u.w$pval.1side1[, -1][weight.u != 0, ]
  ## unweighted test
  cur.UniFrac.l <- AF_combine(p1.u, p1.u.perm, weight = weight.u.nz, log = FALSE)
  sel.u.l <- cur.UniFrac.l$which.selected
  cur.UniFrac.u <- AF_combine(1 - p1.u, 1 - p1.u.perm, weight = weight.u.nz, log = FALSE)
  sel.u.u <- cur.UniFrac.u$which.selected
  rm(p1.u, p1.u.perm, weight.u, weight.u.nz)
  gc()
  cur.UniFrac <- AF_combine(c(cur.UniFrac.l$pvalue, cur.UniFrac.u$pvalue), rbind(cur.UniFrac.l$p.perm, cur.UniFrac.u$p.perm), weight = 1, log = FALSE)
  outcome.u <- cur.UniFrac$which.selected
  rm(cur.UniFrac.l, cur.UniFrac.u)
  gc()
  AF.UniFrac <- cur.UniFrac$pvalue
  sel.u <- select.l.u(sel.u.l, sel.u.u, com = outcome.u, label.u)
  rm(sel.u.l, sel.u.u, outcome.u)
  gc()

  weight.w <- c(br.len) * p.u.w$weight2
  weight.w.nz <- weight.w[weight.w != 0]
  label.w <- cum.labels[weight.w != 0]
  p1.w <- p.u.w$pval.1side2[, 1][weight.w != 0]
  p1.w.perm <- p.u.w$pval.1side2[, -1][weight.w != 0, ]
  ## weighted test
  cur.wUniFrac.l <- AF_combine(p1.w, p1.w.perm, weight = weight.w.nz, log = FALSE)
  sel.w.l <- cur.wUniFrac.l$which.selected
  cur.wUniFrac.u <- AF_combine(1 - p1.w, 1 - p1.w.perm, weight = weight.w.nz, log = FALSE)
  sel.w.u <- cur.wUniFrac.u$which.selected
  rm(p1.w, p1.w.perm, weight.w, weight.w.nz)
  gc()
  cur.wUniFrac <- AF_combine(c(cur.wUniFrac.l$pvalue, cur.wUniFrac.u$pvalue), rbind(cur.wUniFrac.l$p.perm, cur.wUniFrac.u$p.perm), weight = 1, log = FALSE)
  outcome.w <- cur.wUniFrac$which.selected
  rm(cur.wUniFrac.l, cur.wUniFrac.u)
  gc()
  AF.wUniFrac <- cur.wUniFrac$pvalue
  sel.w <- select.l.u(sel.w.l, sel.w.u, com = outcome.w, label.w)
  rm(sel.w.l, sel.w.u, outcome.w)
  gc()

  weight.5 <- c(br.len) * p.u.w$weight3
  weight.5.nz <- weight.5[weight.5 != 0]
  label.5 <- cum.labels[weight.5 != 0]
  p1.5 <- p.u.w$pval.1side3[, 1][weight.5 != 0]
  p1.5.perm <- p.u.w$pval.1side3[, -1][weight.5 != 0, ]
  ## 0.5 - transformation test
  cur.UniFrac5.l <- AF_combine(p1.5, p1.5.perm, weight = weight.5.nz, log = FALSE)
  sel.5.l <- cur.UniFrac5.l$which.selected
  cur.UniFrac5.u <- AF_combine(1 - p1.5, 1 - p1.5.perm, weight = weight.5.nz, log = FALSE)
  sel.5.u <- cur.UniFrac5.u$which.selected
  rm(p1.5, p1.5.perm, weight.5, weight.5.nz)
  gc()
  cur.UniFrac5 <- AF_combine(c(cur.UniFrac5.l$pvalue, cur.UniFrac5.u$pvalue), rbind(cur.UniFrac5.l$p.perm, cur.UniFrac5.u$p.perm), weight = 1, log = FALSE)
  outcome.5 <- cur.UniFrac5$which.selected
  rm(cur.UniFrac5.l, cur.UniFrac5.u)
  gc()
  AF.UniFrac5 <- cur.UniFrac5$pvalue
  sel.5 <- select.l.u(sel.5.l, sel.5.u, com = outcome.5, label.5)
  rm(sel.5.l, sel.5.u, outcome.5)
  gc()

  weight.tip <- p.u.w$weight4
  weight.tip.nz <- weight.tip[weight.tip != 0]
  label.tip <- colnames(X)
  p1.tip <- p.u.w$pval.1side4[, 1][weight.tip != 0]
  p1.tip.perm <- p.u.w$pval.1side4[, -1][weight.tip != 0, ]
  ## only-tip-node test
  cur.tip.l <- AF_combine(p1.tip, p1.tip.perm, weight = weight.tip.nz, log = FALSE)
  sel.tip.l <- cur.tip.l$which.selected
  cur.tip.u <- AF_combine(1 - p1.tip, 1 - p1.tip.perm, weight = weight.tip.nz, log = FALSE)
  sel.tip.u <- cur.tip.u$which.selected
  rm(p1.tip, p1.tip.perm, weight.tip, weight.tip.nz)
  gc()
  cur.tip <- AF_combine(c(cur.tip.l$pvalue, cur.tip.u$pvalue), rbind(cur.tip.l$p.perm, cur.tip.u$p.perm), weight = 1, log = FALSE)
  outcome.tip <- cur.tip$which.selected
  rm(cur.tip.l, cur.tip.u)
  gc()
  AF.tip <- cur.tip$pvalue
  sel.tip <- select.l.u(sel.tip.l, sel.tip.u, com = outcome.tip, label.tip)
  rm(sel.tip.l, sel.tip.u, outcome.tip)
  gc()

  ## ultimate combining
  AFcom <- AF_combine(c(AF.UniFrac, AF.wUniFrac, AF.UniFrac5, AF.tip), rbind(cur.UniFrac$p.perm, cur.wUniFrac$p.perm, cur.UniFrac5$p.perm, cur.tip$p.perm), weight = 1, log = FALSE)
  p.com <- AFcom$pvalue
  outcome.com <- AFcom$which.selected
  sel <- select.com(sel.u, sel.w, sel.5, sel.tip, outcome.com)

  return(list(UniFrac = AF.UniFrac,
              wUniFrac = AF.wUniFrac,
              UniFrac5 = AF.UniFrac5,
              tip.abun = AF.tip,
              com = p.com,
              select = sel))
}
