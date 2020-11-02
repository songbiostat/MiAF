library(statmod)

## Correct the typo in simulateData function in MiSPU package
## The original return command: return(list(informative.OTU = OTU,whole.OTU = OTU))
simulateData.new <- function (nSam = 100, s = 12, ncluster = 20, mu = 1000, size = 25)
{
  data(throat.tree, envir = environment())
  ## The estimate of Dirichlet-multinomial distribution.
  data(dd, envir = environment())
  tree <- get("throat.tree", envir = environment())
  dd1 = get("dd", envir = environment())
  tree.dist <- cophenetic(tree)
  ## Partitioning Around Medoids
  obj <- pam(as.dist(tree.dist), ncluster, diss = TRUE)
  clustering <- obj$clustering
  otu.ids <- tree$tip.label
  p.est = dd1$pi
  names(p.est) <- names(dd1$pi)
  theta <- dd1$theta
  gplus <- (1 - theta)/theta
  p.est <- p.est[otu.ids]
  g.est <- p.est * gplus
  p.clus <- sort(tapply(p.est, clustering, sum), decreasing = T)
  scale2 = function(x) as.numeric(scale(x))
  comm <- matrix(0, nSam, length(g.est))
  rownames(comm) <- 1:nrow(comm)
  colnames(comm) <- names(g.est)
  comm.p <- comm
  nSeq <- rnbinom(nSam, mu = mu, size = size)
  for (i in 1:nSam) {
    comm.p[i, ] <- rdirichlet(1, g.est)[1, ]
    comm[i, ] <- rmultinom(1, nSeq[i], prob = comm.p[i, ])[,1]
  }
  otu.ids <- names(which(clustering == s))
  OTU = comm[, otu.ids]
  return(list(informative.OTU = OTU, whole.OTU = comm))
}

## Score Test
perm.mar.cov <- function(Y, X, binary = FALSE, cov = NULL, nperm = 1000, seed = NULL, ...) {
  K <- ncol(X)
  n <- nrow(X)
  if(is.null(cov)) cov <- rep(1, length(Y))
  
  ## fit the null model
  model <- ifelse(binary, "binomial", "gaussian")
  data1 <- data.frame(trait = Y, cov)
  null.model <- glm(trait ~ ., family = model, data = data1)
  
  ## regress each OTU on covariates
  Xres <- matrix(NA, n, K)
  for (k in 1:K) {
    data2 <- data.frame(response = X[, k], cov)
    X.null <- glm(response ~ ., data = data2)
    Xhat <- fitted.values(X.null)
    Xres[, k] <- X[, k] - Xhat
  }
 
  weight <- apply(X, 2, sd)
  weight[apply(Xres, 2, sd) == 0] <- 0

  ## calculate score test statistics U
  U <- glm.scoretest(null.model, Xres)

  ## permute X residuals and calculate permutation Us
  order.perm <- replicate(nperm, sample(1:length(Y)))
  U.perm <- apply(order.perm, 2, function(x) glm.scoretest(null.model, Xres[x, ]))

  ## Calculate two-sided and one-sided p-values
  Up <- cbind(U, U.perm)
  p1 <- pnorm(Up)
  
  return(list(U = Up, pval.1side = p1, weight = weight))
  ## U are the score statistics
  ## pval.1side is the lower 1-sided p-values; can be further used by minp and AF.perm... functions
}


## Compute AF P-values
AF.wr <- function(p, p.perm, weight = 1, log = TRUE) {
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
  AF.all <- -log(apply(ps.all, 2, min))
  AF.perm <- AF.all[2:(B + 1)]
  AF <- AF.all[1]
  n.select <- which.min(ps.all[, 1])
  which.select <- (1:length(p))[order(wp)][1:n.select]
  p.all <- rank(-AF.all, ties.method = "max")/(B + 1)
  p.AF <- p.all[1]
  p.AF.perm <- p.all[2:(B + 1)]
  AF.list <- list(
    pvalue = p.AF,
    stat = AF,
    stat.perm = AF.perm,
    p.perm = p.AF.perm,
    n.selected = n.select,
    which.selected = which.select
  )
  return(AF.list)
}

## Comebine P-values using AF method or minP method
combine.method <- function(tests, tests.perm, method = c("AF", "minP")) {
  method <- match.arg(method)
  if (method == "AF") {
    res <- AF.wr(tests, tests.perm, weight = 1, log = FALSE)
  } else if (method == "minP") {
    tests.all <- cbind(tests, tests.perm)
    pmin.all <- apply(tests.all, 2, min)
    pval.all <- rank(pmin.all, ties.method = "max")/(length(pmin.all))
    res <- list(
      pvalue = pval.all[1],
      stat = pmin.all[1],
      stat.perm = pmin.all[-1],
      p.perm = pval.all[-1]
    )
  }
  return(res)
}

## Y is the outcome, X is the OTU table, tree is the phylogenetic tree
## cov is the covariates which will be adjusted for
## model: gaussian - linear model, binomial - logistic regression
## returns P-values from MiAFu (unweighted), MiAFw (weighted) and MiAF (combined)
MiAF <- function (Y, X, tree, cov = NULL, model = c("gaussian","binomial"), n.perm = 1000)
{
  model <- match.arg(model)
  GuniF.cum <- GUniFrac_cum(X, tree)
  cum <- GuniF.cum$cum     # "extended" OTU table 
  br.len <- GuniF.cum$br.len    # branch length
  br.len <- as.matrix(br.len)
  tmp.cum <- cum  
  tmp.cum[tmp.cum != 0] <- 1    # unweighted version
  ncum <- length(c(br.len))

  xx <- t(rbind(tmp.cum, cum))
  p.u.w <- perm.mar.cov(Y, xx, binary = (model == "binomial"), cov = cov, res = TRUE, n.perm = n.perm, seed = NULL)
  w2 <- p.u.w$weight

  w2.u <- w2[1:ncum]
  weight.u <- c(br.len) * w2.u    # weight for unweighted version
  weight.u.nz <- weight.u[weight.u != 0]    # remove those OTUs whose weight is zero

  w2.w <- w2[-(1:ncum)]
  weight.w <- c(br.len) * w2.w    # weight for weighted version
  weight.w.nz <- weight.w[weight.w != 0]    # remove those OTUs whose weight is zero

  p1.u <- p.u.w$pval.1side[1:ncum, 1][weight.u != 0]     # 1-sided P-values for unweighted version
  p1.w <- p.u.w$pval.1side[-(1:ncum), 1][weight.w != 0]     # 1-sided P-values for weighted version

  p1.u.perm <- p.u.w$pval.1side[1:ncum, -1][weight.u != 0, ]    # permutation of 1-sided P-values for unweighted version
  p1.w.perm <- p.u.w$pval.1side[-(1:ncum), -1][weight.w != 0, ]     # permutation of 1-sided P-values for weighted version

  cur.UniFrac.l <- AF.wr(p1.u, p1.u.perm, weight = weight.u.nz, log = FALSE)     # lower 1-sided AF P-value for unweighted version
  cur.UniFrac.u <- AF.wr(1 - p1.u, 1 - p1.u.perm, weight = weight.u.nz, log = FALSE)     # upper 1-sided AF P-value for unweighted version
  cur.UniFrac <- combine.method(c(cur.UniFrac.l$pvalue, cur.UniFrac.u$pvalue), rbind(cur.UniFrac.l$p.perm, cur.UniFrac.u$p.perm), method = "AF")     # combine lower and upper 1-sided P-values for unweighted version
  AF.UniFrac <- cur.UniFrac$pvalue

  cur.wUniFrac.l <- AF.wr(p1.w, p1.w.perm, weight = weight.w.nz, log = FALSE)     # lower 1-sided AF P-value for weighted version
  cur.wUniFrac.u <- AF.wr(1 - p1.w, 1 - p1.w.perm, weight = weight.w.nz, log = FALSE)     # upper 1-sided AF P-value for weighted version  
  cur.wUniFrac <- combine.method(c(cur.wUniFrac.l$pvalue, cur.wUniFrac.u$pvalue), rbind(cur.wUniFrac.l$p.perm, cur.wUniFrac.u$p.perm), method = "AF")     # combine lower and upper 1-sided P-values for weighted version
  AF.wUniFrac <- cur.wUniFrac$pvalue

  p.com <- combine.method(c(cur.UniFrac$pvalue, cur.wUniFrac$pvalue), rbind(cur.UniFrac$p.perm, cur.wUniFrac$p.perm), method = "AF")$pvalue     # combine unweighted and weighted versions of P-values

  return(list(UniFrac = AF.UniFrac, wUniFrac = AF.wUniFrac, com = p.com))
}
