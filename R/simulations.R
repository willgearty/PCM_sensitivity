library(ape)
library(phytools)
library(paleotree)
library(TreeSim)
library(FossilSim)
library(geiger)

# Cooper et al 2016 https://academic.oup.com/biolinnean/article/118/1/64/2440254
# https://github.com/nhcooper123/OhYou/blob/master/R/OhYou_simulations.R
# lambda = 1
# mu = 0, 0.25, 0.5, 0.75
# n_tips = 25, 50, 100, 150, 200, 500, 1000

# we want to get a set of simulation trees with a specific number of tips
# and a specific proportion of fossils

# get proportions of extinct and extant tips
# use this to make sure the code works
get_props <- function(tree) {
  n_extant <- length(getExtant(tree))
  n_extinct <- length(getExtinct(tree))
  n_total <- n_extant + n_extinct
  data.frame(raw = c(n_extant, n_extinct, n_total),
             prop = c(n_extant, n_extinct, n_total) / n_total,
             row.names = c("extant", "extinct", "total"))
}

# load in simulation functions
source("sim.fossils.R")

trees <- sim.fbd.taxa.prop(100, .25, numbsim = 10, lambda = 1, mu = .25)
# check that we get the correct proportions
lapply(trees, get_props)

# Now we want to vary the preservation rate through time...
trees1 <- sim.fbd.taxa.prop(100, .1, numbsim = 10, lambda = 1, mu = .25,
                            model = "EB", a = 2)
trees2 <- sim.fbd.taxa.prop(100, .1, numbsim = 10, lambda = 1, mu = .25,
                            model = "EB", a = -2)
# check that the rescaling shifts which extinct species we are sampling
par(mfrow=c(2, 1))
plot(trees1[[1]]) # should have the fossils near the tips
plot(trees2[[1]]) # should have the fossils near the root
layout(1)
