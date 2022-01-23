#' Taxon Proportion for All Non-rooted Nodes in a Phylogenetic Tree
#'
#' @description Calculate the taxon proportion for all non-rooted nodes
#'              in a phylogenetic tree, as well as return the labels for those nodes.
#'              The code is adapted from GUniFrac package contributed by
#'              Jun Chen (Jun 12, 2015).
#'
#' @param otu.tab OTU count table.
#'                Row - n samples, column - OTUs
#' @param tree A rooted phylogenetic tree of R class "phylo".
#'
#' @return A list containing the following elements:
#'         \item{cum}{Taxon proportion for all non-rooted nodes.
#'                    The nodes are arranged in the same order as the second column of tree$edge.}
#'         \item{br.len}{Branch length for all non-rooted nodes.}
#'         \item{cum.labels}{Labels for all non-rooted nodes.
#'                           If node.label is not provided by the tree of R class "phylo", we assign "Node k" to the kth node. }
#'
#' @export
#'
#' @import ape
#'
cum_prop <- function (otu.tab, tree) {
	if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!")

	# Convert into proportions
	otu.tab <- as.matrix(otu.tab)
	row.sum <- rowSums(otu.tab)
	otu.tab <- otu.tab / row.sum
	n <- nrow(otu.tab)

	# Construct the returning array
	if (is.null(rownames(otu.tab))) {
		rownames(otu.tab) <- paste("comm", 1:n, sep="_")
	}

	# Check OTU name consistency
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

	# Reorder the otu.tab matrix if the OTU orders are different
	tip.label <- tree$tip.label
	otu.tab <- otu.tab[, tip.label]


	ntip <- length(tip.label)
	nbr <- nrow(tree$edge)
	edge <- tree$edge
	edge2 <- edge[, 2]
  br.len <- tree$edge.length  # branch length


    #  Accumulate OTU proportions up the tree
	cum <- matrix(0, nbr, n)							# Branch abundance matrix
	for (i in 1:ntip) {
		tip.loc <- which(edge2 == i)
		cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]
		node <- edge[tip.loc, 1]						# Assume the direction of edge
		node.loc <- which(edge2 == node)
		while (length(node.loc)) {
			cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]
			node <- edge[node.loc, 1]
			node.loc <- which(edge2 == node)
		}
	}

	if(is.null(tree$node.label)) tree$node.label <- paste("Node", 1:tree$Nnode, sep = "")
	all.labels <- c(tree$tip.label, tree$node.label)
	cum.labels <- all.labels[tree$edge[, 2]]

	out = list(cum = cum, br.len = br.len, cum.labels = cum.labels)
	return(out)
}
