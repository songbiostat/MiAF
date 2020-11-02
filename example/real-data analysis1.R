library(MiSPU)
library(MiRKAT)
library(ggplot2)
library(ggtree)
source("function_select.R")
source("GUniFrac_select.R")

data(throat.otu.tab)
data(throat.tree)
data(throat.meta)

set.seed(201908)
otu.tab <- throat.otu.tab
tree <- throat.tree

id <- rownames(throat.otu.tab)
meta <- throat.meta[id,]
y <- as.numeric(meta$SmokingStatus == 'Smoker')
gender <- as.numeric(meta$Sex == 'Female')
antibio <- as.numeric(meta$AntibioticUsePast3Months_TimeFromAntibioticUsage)
zcov <- cbind(gender, antibio)

if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
  stop("The OTU table contains unknown OTUs! OTU names
       in the OTU table and the tree should match!" )
}


# Get the subtree if tree contains more OTUs
absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
if (length(absent) != 0) {
  tree <- drop.tip(tree, absent)
  warning("The tree has more OTU than the OTU table!")
}

otu.tab.r <- Rarefy(otu.tab)$otu.tab.rff

unifracs.r= GUniFrac::GUniFrac(otu.tab.r, tree)$unifracs
D.weighted = unifracs.r[,,"d_1"]
D.unweighted = unifracs.r[,,"d_UW"]
D.BC= as.matrix(vegdist(otu.tab.r, method="bray"))

duw = D2K(D.unweighted)
dw = D2K(D.weighted)
dbc = D2K(D.BC)
ks = list(duw, dw, dbc)

RKATr <- MiRKAT(y, X = zcov, Ks = ks, out_type = "D", nperm = 9999, method= "permutation")
SPUr <- MiSPU::MiSPU(y, otu.tab.r, throat.tree, zcov, model = "binomial", pow = c(2:8, Inf), n.perm = 9999)
AFr <- MiAF(y, otu.tab.r, throat.tree, zcov, model = "binomial", n.perm = 9999)


## plot associated taxon selected by AF method
if(is.null(tree$node.label)) tree$node.label <- paste("Node", 1:tree$Nnode, sep = "")
all.labels <- c(tree$tip.label, tree$node.label)
col <- rep(3, length(all.labels))
col[which(all.labels %in% AFr$select$lower)] <- 1
col[which(all.labels %in% AFr$select$upper)] <- 2
si <- rep(3, length(all.labels))
si[which(all.labels %in% AFr$select$lower)] <- 1
si[which(all.labels %in% AFr$select$upper)] <- 1
cols <- c("1" = "red3", "3" = "grey50", "2" = "blue3")
siz <- c("1" = 0.3,  "3" = 0.2)

pdf("select.pdf", width = 5, height = 10)
ggtree(tree, aes(color=as.factor(col), size = as.factor(si))) + 
  scale_colour_manual(values = cols) + scale_size_manual(values = siz) +
  geom_text2(aes(label=all.labels, subset = isTip &(all.labels %in% AFr$select$lower)), hjust = -.5) +
  geom_text2(aes(label=all.labels, subset = isTip &(all.labels %in% AFr$select$upper)), hjust = -.5) +
  geom_text2(aes(label=all.labels, subset = isTip &(! all.labels %in% AFr$select$lower) & (! all.labels %in% AFr$select$upper)), hjust = -.5)
dev.off()
