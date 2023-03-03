#Bethany      3rd March 2023
#Simulate continuous trait evolution across a phylogeny using different models

library(ape)
library(geiger)
library(phytools)

tree <- read.tree("data/testtree.txt")
plot(tree)

#Brownian motion in geiger::sim.char()
BM_traits <- sim.char(phy = tree, par = 0.02, nsim = 1, model = "BM", root = 1)
names(BM_traits) <- c(paste0('t', seq(1, 100, 1)))
contMap(tree, BM_traits)

