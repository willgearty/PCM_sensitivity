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
  # assumes that 0.9 is the maximum relative extinction rate
  n_extant_tips <- n * ((lambda / mu) / (1 / 0.9))
  n_extinct_tips <- prop_extinct * n
  trees <- sim.bd.taxa(n_extant_tips, numbsim,
                       lambda, mu, frac = 1, complete = TRUE)
  # make sure all trees have enough extinct tips
  while (TRUE) {
    n_extinct <- sapply(trees, function(tr) length(is.extinct(tr)))
    # ensure there is a surplus of extinct tips to sample from later
    enough_extinct <- n_extinct > (n_extinct_tips * 1.2)
    if (any(!enough_extinct)) {
      cat("more")
      more_trees <- sim.bd.taxa(n_extant_tips, sum(!enough_extinct),
                                lambda, mu, frac = 1, complete = TRUE)
      trees[!enough_extinct] <- more_trees
    } else {
      break
    }
  }

  if (progress) pb <- txtProgressBar(max = numbsim, style = 3)
  for (i in 1:length(trees))
  {
    t <- trees[[i]]
    f <- sim.fossil.tips(n_extinct_tips, tree = t, model = model, ...)

    tree <- SAtree.from.fossils(t, f)$tree

    node.ages <- FossilSim:::n.ages(tree)

    origin <- max(node.ages) + tree$root.edge

    if ( !complete ) {
      tree <- FossilSim:::drop.unsampled(tree, frac = (n - n_extinct_tips) / n_extant_tips, n = -1)
      node.ages <- FossilSim:::n.ages(tree)
    }

    trees[[i]] <- tree
    trees[[i]]$root.edge <- origin - max(node.ages)

    trees[[i]] <- SAtree(trees[[i]], complete)
    trees[[i]]$tip.label <- paste("t", 1:Ntip(trees[[i]]), sep = "")
    if (progress) setTxtProgressBar(pb, i)
  }
  if (progress) close(pb)
  # class(trees) <- c("multiPhylo", "list")
  return(trees)
}

# This function simulates fossil occurrences along the extinct branches of a
# tree using sim.fossils.poisson, and then samples from the branches based on
# aggregate weighting under a time-varying model. This allows the user to impose
# a temporal bias for fossil sampling.
# n: number of extinct tips to sample
# tree: a phylo object representing the tree to sample from
# model: a function that takes a relative age (between 0 and 1) and returns a
#        recovery potential (a non-negative number)
sim.fossil.tips <- function(n, tree, model = function(t) 1, ...) {
  while (TRUE) {
    # simulate a large number of fossil occurrences with sim.fossils.poisson
    # until at least n extinct edges have fossils
    foss <- sim.fossils.poisson.data.table(50, tree, root.edge = FALSE)
    # filter to only those fossils on extinct branches (exclude extant tips and internal branches)
    foss_sub <- subset(foss, edge %in% match(is.extinct(tree), tree$tip.label))
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
  # note: sample() uses floor(n)
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

