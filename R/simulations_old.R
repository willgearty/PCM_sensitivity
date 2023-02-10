# this represents Will's stream of consciousness investigating how to
# simulate the trees we want. Enjoy!

library(ape)
library(phytools)
library(paleotree)

# Cooper et al 2016 https://academic.oup.com/biolinnean/article/118/1/64/2440254
# https://github.com/nhcooper123/OhYou/blob/master/R/OhYou_simulations.R
# lambda = 1
# mu = 0, 0.25, 0.5, 0.75
# n_tips = 25, 50, 100, 150, 200, 500, 1000


# we want to get a set of simulation trees with a specific number of tips
# and a specific proportion of fossils

# get proportions of extinct and extant tips
get_props <- function(tree) {
  n_extant <- length(getExtant(tree))
  n_extinct <- length(getExtinct(tree))
  n_total <- n_extant + n_extinct
  data.frame(raw = c(n_extant, n_extinct, n_total),
             prop = c(n_extant, n_extinct, n_total) / n_total,
             row.names = c("extant", "extinct", "total"))
}

library(TreeSim)
# sim.bd.taxa.age doesn't include fossils

# sim.bd.age is not reliable to get the number of tips we want
# it also often fails
trees <- sim.bd.age(2, numbsim = 1, lambda = 2, mu = 0.5, frac = 0.6)
lapply(trees, get_props)

# sim.bd.taxa is much faster, but we can't specify the age
# but we could scale the trees
# can specify the number of extant tips, but can not specify the number of
# extinct tips
trees <- sim.bd.taxa(50, numbsim = 10, lambda = 2, mu = 0.5, frac = .6)
lapply(trees, get_props)

# you could brute force the algorithm to just keep going until it
# finds a simulation with the number of extinct tips that you want
# the problem with this is that it is hard to vary the fossil proportion
# it will require lots of trial and error to figure out what lambda and mu work
sim.bd.taxa.brute <- function(n, n_extinct, numbsim, lambda, mu, frac = 1,
                              complete = TRUE, stochsampling = FALSE) {
  lst <- list()
  pb <- txtProgressBar(max = numbsim)
  for (i in 1:numbsim) {
    while(TRUE) {
      tree <- sim.bd.taxa(n, 1, lambda = lambda, mu = mu, frac = frac,
                          complete = complete, stochsampling = stochsampling)
      tree <- tree[[1]]
      if (length(getExtinct(tree)) == n_extinct) break
    }
    lst[[i]] <- tree
    setTxtProgressBar(pb, i)
  }
  close(pb)
  lst
}
trees <- sim.bd.taxa.brute(50, 20, numbsim = 10, lambda = 2, mu = 0.5, frac = 1)
lapply(trees, get_props)
# lambda = 2, mu = 1.9 appears to give a fossil prop of ~.95 (with lots of noise)
# lambda = 2, mu = 1 appears to give a fossil prop of ~.5
# lambda = 2, mu = .2 appears to give a fossil prop of ~.1

# this is tougher to brute force, but still works:
# 1. a lot of the tries result in failures (everything goes extinct)
# 2. there's a lot of noise with high rates, resulting in a wide range of outcomes
trees <- sim.bd.taxa.brute(5, 95, numbsim = 10, lambda = 2, mu = 1.8, frac = 1)
lapply(trees, get_props)

# but all of this requires that we tweak the rates to get the proportions we want

# instead, we probably want to have set rates, then just have higher fossil/extant sampling for diff sims

library(FossilSim)

trees <- sim.fbd.taxa(75, numbsim = 3, lambda = 1, mu = .5, psi = 0.1, frac = .6)
lapply(trees, get_props)
# this seems promising...can exactly control the # of extant taxa
# and can control the fossil sampling rate
# we'll still need to brute force it though to get the exact fossil prop we want

sim.fbd.taxa.brute1 <- function(n, n_extinct, numbsim, lambda, mu, psi, frac = 1,
                                complete = FALSE) {
  lst <- list()
  pb <- txtProgressBar(max = numbsim)
  for (i in 1:numbsim) {
    while(TRUE) { # guarantee we get the number of extinct taxa that we want
      tree <- sim.fbd.taxa(n, 1, lambda = lambda, mu = mu, psi = psi, frac = frac,
                           complete = complete)
      tree <- tree[[1]]
      if (length(getExtinct(tree)) == n_extinct) break
    }
    lst[[i]] <- tree
    setTxtProgressBar(pb, i)
  }
  close(pb)
  lst
}

trees <- sim.fbd.taxa.brute1(75, 25, numbsim = 10, lambda = 1, mu = .5, psi = 0.1, frac = .6)
lapply(trees, get_props)
# lambda = 2, mu = 1.9 appears to give a fossil prop of ~.95 (with lots of noise)
# lambda = 2, mu = 1 appears to give a fossil prop of ~.5
# lambda = 2, mu = .2 appears to give a fossil prop of ~.1

# this one uses the source code of FossilSim to significantly speed up the brute forcing
# it will also estimate the fossil sampling rate for you
sim.fbd.taxa.brute2 <- function(n, n_extinct, numbsim, lambda, mu, complete = FALSE)
{
  trees <- TreeSim::sim.bd.taxa(n, numbsim, lambda, mu, frac = 1, complete = TRUE)
  # simulate the trait(s) here!!

  pb <- txtProgressBar(max = numbsim)
  for(i in 1:length(trees))
  {
    t <- trees[[i]]
    # estimate the best sampling rate for the desired number of fossils
    lambda <- n_extinct / (sum(t$edge.length) + t$root.edge)

    while(TRUE) { # guarantee we get the number of extinct taxa that we want
      f <- sim.fossils.poisson(tree = t, rate = lambda)
      if (nrow(f) == n_extinct) break
    }

    tree <- SAtree.from.fossils(t,f)$tree

    node.ages <- FossilSim:::n.ages(tree)

    origin <- max(node.ages) + tree$root.edge

    if( !complete ) {
      tree <- FossilSim:::drop.unsampled(tree, frac = 1, n = n)
      node.ages <- FossilSim:::n.ages(tree)
    }

    trees[[i]] <- tree
    trees[[i]]$root.edge <- origin - max(node.ages)

    trees[[i]] <- SAtree(trees[[i]], complete)
    trees[[i]]$tip.label <- paste("t", 1:Ntip(trees[[i]]), sep = "")
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(trees)
}
trees <- sim.fbd.taxa.brute2(75, 25, numbsim = 10, lambda = 1, mu = .5)
lapply(trees, get_props)

trees <- sim.fbd.taxa.brute2(25, 75, numbsim = 10, lambda = 1, mu = .5)
lapply(trees, get_props)

# sim.fossils.intervals is how we could get more fossils towards the root/tips
