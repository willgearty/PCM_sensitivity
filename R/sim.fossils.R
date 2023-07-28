# Written by Will Gearty 2/10/2023

library(TreeSim)
library(FossilSim)
library(geiger)

# This uses modified source code from FossilSim to pick the
# correct proportion of fossils and the correct number of extant tips.
# The chance of sampling an extinct species is dependent on its branch length.
# Non-uniform fossil sampling is imposed by temporarily rescaling the tree
# (options can be found in geiger::rescale).
# n: number of total taxa in the final tree
# prop_extinct: proportion of the final taxa that are extinct
# numbsim: number of simulations
# lambda: birth/speciation rate
# mu: death/extinction rate
# complete: if TRUE, unsampled ancestors are included in the final tree
#           (meaning it will have more than n tips)
# model: a transformation model passed to geiger::rescale used only for fossil sampling
#        (for temporal sampling biases/trends)
# ...: other arguments passed to geiger::rescale (e.g., `a` for the EB `model`)
sim.fbd.taxa.prop <- function(n, prop_extinct, numbsim, lambda, mu,
                              complete = FALSE, model = NULL, progress = TRUE, ...)
{
  trees <- sim.bd.taxa(n, numbsim, lambda, mu, frac = 1, complete = TRUE)

  if (progress) pb <- txtProgressBar(max = numbsim, style = 3)
  for(i in 1:length(trees))
  {
    t <- trees[[i]]
    f <- sim.fossils(prop_extinct * n, tree = t, model = model, ...)

    tree <- SAtree.from.fossils(t,f)$tree

    node.ages <- FossilSim:::n.ages(tree)

    origin <- max(node.ages) + tree$root.edge

    if( !complete ) {
      tree <- FossilSim:::drop.unsampled(tree, frac = 1 - prop_extinct, n = -1)
      node.ages <- FossilSim:::n.ages(tree)
    }

    trees[[i]] <- tree
    trees[[i]]$root.edge <- origin - max(node.ages)

    trees[[i]] <- SAtree(trees[[i]], complete)
    trees[[i]]$tip.label <- paste("t", 1:Ntip(trees[[i]]), sep = "")
    if (progress) setTxtProgressBar(pb, i)
  }
  if (progress) close(pb)
  class(trees) <- c("multiPhylo", "list")
  return(trees)
}

# Sample fossil species across a phylogeny
# The sampling process can be biased temporally by specifying a `model` that is
# used to temporarily transform the tree using geiger::rescale.
# Note that a given extinct species can only be sampled at most once.
# n: number of fossil tips to sample
# tree: non-ultrametric phylogeny
# model: a transformation model passed to geiger::rescale used only for fossil sampling
#        (for temporal sampling biases/trends)
# ...: other arguments passed to geiger::rescale (e.g., `a` for the EB `model`)
sim.fossils <- function(n, tree = NULL, model = NULL, ...) {
  # scale tree however the user wants
  scaled_tree <- tree
  if (!is.null(model)) {
    scaled_tree <- rescale(scaled_tree, model, ...)
  }
  # get branch lengths of rescaled tree
  rescaled_taxonomy <- sim.taxonomy(scaled_tree, beta = 1)
  # sample species using the rescaled branch lengths
  sps <- sample.int(n = nrow(rescaled_taxonomy), size = n,
                        prob = sapply(rescaled_taxonomy$start - rescaled_taxonomy$end,
                                      function(x) max(0.00001, x))) # each branch always has a small chance
  # sample the fossils from the sampled species
  fdf <- fossils()
  taxonomy <- sim.taxonomy(tree, beta = 1)
  for (sp in sps){
    start <- max(taxonomy$start[which(taxonomy$sp == sp)])
    end <- min(taxonomy$end[which(taxonomy$sp == sp)])
    h <- runif(1, min = end, max = start)
    fdf <- rbind(fdf, data.frame(sp = sp, edge = sp, hmin = h, hmax = h, stringsAsFactors = F))
  }
  fdf <- fossils(fdf)
  return(fdf)
}
