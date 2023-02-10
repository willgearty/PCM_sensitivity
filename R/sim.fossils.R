library(TreeSim)
library(FossilSim)
library(geiger)

# this uses modified source code from FossilSim to pick the
# correct proportion of fossils and the correct number of extant tips
# the chance of sampling an extinct species is dependent its the branch length
# non-uniform fossil sampling is imposed by rescaling the tree
# options can be found in geiger::rescale
sim.fbd.taxa.prop <- function(n, prop_extinct, numbsim, lambda, mu,
                              complete = FALSE, model = NULL, ...)
{
  trees <- TreeSim::sim.bd.taxa(n, numbsim, lambda, mu, frac = 1, complete = TRUE)

  pb <- txtProgressBar(max = numbsim, style = 3)
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
    setTxtProgressBar(pb, i)
  }
  close(pb)
  class(trees) <- c("multiPhylo", "list")
  return(trees)
}

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
                        prob = rescaled_taxonomy$start - rescaled_taxonomy$end)
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
