library(MiSPU)
library(MiRKAT)
library(GUniFrac)
library(phytools)
source("function_select.R")
source("GUniFrac_select.R")

## The processed otu table and meta data are available at openly available 
## in MicrobiomeHD at http://doi.org/10.5281/zenodo.1146764. 
otu <- read.table("RDP/hiv_lozupone.otu_table.100.denovo.rdp_assigned")
meta <- read.table("hiv_lozupone.metadata.txt", sep = "\t", header = TRUE)

tree <- read.newick("tree.nwk")
taxa_name <- rownames(otu)
rownames(otu) <- gsub("(.*)_", "", taxa_name) 
meta <- meta[-(1:6),]
id <- colnames(otu)
meta <- meta[meta$new_sample_name %in% id, ]
otu <- as.matrix(otu)
meta$new_sample_name <- as.character(meta$new_sample_name)
otu <- otu[, meta$new_sample_name]
otu <- t(otu)
meta <- meta[meta$time_point == 1,]
otu <- otu[meta$new_sample_name, ]
y <- ifelse(meta$DiseaseState == "H", 0, 1)
age <- meta$age
bodymass <- meta$body_mass_index
zcov <- cbind(bodymass, age)

if (sum(!(colnames(otu) %in% tree$tip.label)) != 0) {
  stop("The OTU table contains unknown OTUs! OTU names
       in the OTU table and the tree should match!" )
}

# Get the subtree if tree contains more OTUs
absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu))]
if (length(absent) != 0) {
  tree <- drop.tip(tree, absent)
  warning("The tree has more OTU than the OTU table!")
}

set.seed(12345)

otu.tab.r <- Rarefy(otu)$otu.tab.rff

unifracs.r= GUniFrac::GUniFrac(otu.tab.r, tree)$unifracs
D.weighted = unifracs.r[,,"d_1"]
D.unweighted = unifracs.r[,,"d_UW"]
D.5 = unifracs.r[,,"d_0.5"]
D.BC= as.matrix(vegdist(otu.tab.r, method="bray"))

duw = D2K(D.unweighted)
dw = D2K(D.weighted)
dbc = D2K(D.BC)
ks = list(duw, dw, dbc)

RKATr <- MiRKAT(y, X = NULL, Ks = ks, out_type = "D", nperm = 9999, method= "permutation")
SPUr <- MiSPU::MiSPU(y, otu.tab.r, tree, cov = NULL, model = "binomial", pow = c(2:8, Inf), n.perm = 9999)
AFr <- MiAF (y, otu.tab.r, tree, cov = NULL, model = "binomial", n.perm = 9999)
