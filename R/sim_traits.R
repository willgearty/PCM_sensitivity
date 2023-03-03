#Bethany      3rd March 2023
#Simulate continuous trait evolution across a phylogeny using different models

library(ape)
library(geiger)
library(phytools)
library(mvMORPH)

tree <- read.tree("data/testtree.txt")
plot(tree)

#BM in geiger::sim.char()
BM_geiger <- sim.char(phy = tree, par = 0.02, nsim = 1, model = "BM", root = 1)
names(BM_geiger) <- c(paste0('t', seq(1, 100, 1)))
contMap(tree, BM_geiger)

#BM, BM with trend, bounded and OU in phytools::fastBM()
#Note that internal nodes can also be estimated
BM_phytools <- fastBM(tree, a = 1, mu = 0, sig2 = 1, bounds = c(-Inf, Inf),
                      internal = F, nsim = 1)
names(BM_phytools) <- c(paste0('t', seq(1, 100, 1)))
contMap(tree, BM_phytools)

OU_phytools <- fastBM(tree, a = 1, sig2 = 1, alpha = 1, theta = 1, internal = F,
                      nsim = 1)
names(OU_phytools) <- c(paste0('t', seq(1, 100, 1)))
contMap(tree, OU_phytools)

#BM, OU and ACDC / Early Burst in mvMORPH::mvSIM()
BM_mvMORPH <- mvSIM(tree = tree, nsim = 1, model = "BM1",
                    param = list(alpha = 1, theta = 1, sigma = 0.5))
names(BM_mvMORPH) <- c(paste0('t', seq(1, 100, 1)))
contMap(tree, BM_mvMORPH)

OU_mvMORPH <- mvSIM(tree = tree, nsim = 1, model = "OU1",
                    param = list(alpha = 1, theta = 1, sigma = 0.5))
names(OU_mvMORPH) <- c(paste0('t', seq(1, 100, 1)))
contMap(tree, OU_mvMORPH)

EB_mvMORPH <- mvSIM(tree = tree, nsim = 1, model = "BM1",
                    param = list(alpha = 1, theta = 1, beta = 1, sigma = 0.5))
names(EB_mvMORPH) <- c(paste0('t', seq(1, 100, 1)))
contMap(tree, EB_mvMORPH)
