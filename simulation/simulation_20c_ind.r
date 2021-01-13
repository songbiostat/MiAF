## This is a simulation example for the independent case of scenario 1.
## Under scenario 1 where OTUs are divided into 20 clusters, and the outcome is truly associated with the 1st cluster.
## We will assume the outcome is truly associated with each cluster in turn.
## In this independent case, covariates Z2 is independent of OTUs.

library(MiSPU)
library(boot)
library(MiRKAT)
library(vegan)
library(GUniFrac)
source("GUniFrac.R")
source("function.R")

data(throat.otu.tab)
data(throat.tree)
data(throat.meta)
tree <- throat.tree

n <- 100     ## number of samples
clu <- 20     ## OTUs are divided into 20 clusters under scenario 1
## The effect size is set to be 0, 0.6, 0.8, 1.2, 1.6 and 2 in the independent case
beta <- 1.6     ## effect size
rep.pow <- 1000     ## number of replicates 
n.perm <- 1000     ## number of permutation

p.MiRKATw <- c()
p.MiRKATu <- c()
p.MiRKAT5 <- c()
p.MiRKATopt <- c()
p.aMiSPUw <- c()
p.aMiSPUu <- c()
p.aMiSPU <- c()
p.AFw <- c()
p.AFu <- c()
p.AFcom <- c()

for (i in 1:rep.pow){
  Z1 <- rbinom(n, 1, 0.5)
  ## assume the outcome is truly associated with the 1st (s) cluster (s = 1)
  OTU <- simulateData.new(nSam = n, s = 1, ncluster = clu, mu = 1000, size = 25)
  X <- OTU$whole.OTU
  X.inf <- OTU$informative.OTU
  Z2 <- rnorm(n, 0, 1)
  Z <- cbind (Z1, Z2)
  P <- inv.logit(0.5 * scale(Z1 + Z2, center = TRUE, scale = TRUE) + beta * scale(apply(X.inf, 1, sum), center = TRUE, scale = TRUE))
  Y <- rbinom(n, 1, P)

  unifrac <- GUniFrac::GUniFrac(X, tree)$unifracs
  dw <- unifrac[,,"d_1"] # Weighted UniFrac
  du <- unifrac[,,"d_UW"] # Unweighted UniFrac
  d5 <- unifrac[,,"d_0.5"]   # GUniFrac with alpha 0.5
  dbc <- as.matrix(vegdist(X , method="bray"))
  
  du = D2K(du)
  dw = D2K(dw)
  d5 = D2K(d5)
  dbc = D2K(dbc)
  ks = list(du, dw, d5, dbc)
  
  ## MiRKAT method
  RKAT <- MiRKAT(Y, Z, Ks = ks, out_type = "D", nperm = 1000, method = "davies")
  p.MiRKATu[i] <- RKAT$indivP[1]
  p.MiRKATw[i] <- RKAT$indivP[2]
  p.MiRKAT5[i] <- RKAT$indivP[3]
  p.MiRKATopt[i] <- RKAT$omnibus_p
  
  ## aMiSPU method
  SPU <- MiSPU(Y, X, throat.tree, Z, model = "binomial", pow = c(2:8, Inf), n.perm = n.perm)
  p.aMiSPUw[i] <- SPU$Weighted$pvs[9]
  p.aMiSPUu[i] <- SPU$Unweighted$pvs[9]
  p.aMiSPU[i] <- SPU$aMiSPU$pvalue
  
  ## AF method
  AF <- MiAF (Y, X, throat.tree, Z, model = "binomial", n.perm = n.perm)
  p.AFw[i] <- AF$wUniFrac
  p.AFu[i] <- AF$UniFrac
  p.AFcom[i] <- AF$com
}

power.MiRKATw <- mean(p.MiRKATw <= 0.05)
power.MiRKATu <- mean(p.MiRKATu <= 0.05)
power.MiRKAT5 <- mean(p.MiRKAT5 <= 0.05)
power.MiRKATopt <- mean(p.MiRKATopt <= 0.05)
power.aMiSPUw <- mean(p.aMiSPUw <= 0.05)
power.aMiSPUu <- mean(p.aMiSPUu <= 0.05)
power.aMiSPU <- mean(p.aMiSPU <= 0.05)
power.AFw <- mean(p.AFw <= 0.05)
power.AFu <- mean(p.AFu <= 0.05)
power.AFcom <- mean(p.AFcom <= 0.05)
