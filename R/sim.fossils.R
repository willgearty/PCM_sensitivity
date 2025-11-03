# Written by Will Gearty 2/10/2023
# Last updated 6/27/2025

library(TreeSim)
library(FossilSim)
library(geiger)
library(dplyr)

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
                              complete = FALSE, model = function(t) 1, progress = TRUE, ...)
{
  trees <- sim.bd.taxa(n, numbsim, lambda, mu, frac = 1, complete = TRUE)

  if (progress) pb <- txtProgressBar(max = numbsim, style = 3)
  for(i in 1:length(trees))
  {
    t <- trees[[i]]
    f <- sim.fossil.tips(prop_extinct * n, tree = t, model = model, ...)

    tree <- SAtree.from.fossils(t, f)$tree

    node.ages <- FossilSim:::n.ages(tree)

    origin <- max(node.ages) + tree$root.edge

    if ( !complete ) {
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

# TODO: document args
# determine a sampling "model" (i.e., a temporal bias for fossil sampling)
sim.fossil.tips <- function(n, tree, model = function(t) 1, ...) {
  while (TRUE) {
    # simulate a large number of fossil occurrences with sim.fossils.poisson
    # until at least n extinct edges have fossils
    foss <- sim.fossils.poisson(50, tree, root.edge = FALSE)
    # filter to only those fossils on extinct branches (exclude extant tips and internal branches)
    foss_sub <- subset(foss, edge <= Ntip(tree))
    if (length(unique(foss_sub$edge)) >= n) break
  }
  # based on the relative age of the fossils and the sampling "model",
  # assign each fossil a "recovery potential"
  max_age <- max(FossilSim:::n.ages(tree))
  foss_sub$rel_age <- (max_age - foss_sub$hmin) / max_age
  foss_sub$recovery_potential <- model(foss_sub$rel_age)
  # for each branch, sum the recovery potentials for all occurrences on that branch
  branch_recovery <- foss_sub %>%
    group_by(edge) %>%
    summarise(recovery_potential = sum(recovery_potential, na.rm = TRUE)) %>%
    ungroup()
  # sample the desired # of branches based on the recovery potentials using sample(x, n, prob)
  # x is the branch ids
  # n is the number of desired extinct tips
  # prob is the recovery potentials
  sampled_branches <- sample(branch_recovery$edge, n,
                             prob = branch_recovery$recovery_potential)
  # return a fossil occurrence for each of those branches
  foss_sub %>%
    filter(edge %in% sampled_branches) %>%
    group_by(edge) %>%
    slice(1) %>%
    ungroup() %>%
    select(-rel_age, -recovery_potential) %>%
    as.fossils()
}

