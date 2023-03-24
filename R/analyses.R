# Written by Will Gearty, Bethany Allen, and Pedro Godoy

# Load libraries
library(ape)
library(phytools)
library(TreeSim)
library(FossilSim)
library(geiger)
library(mvMORPH)
library(pbapply)

# Load functions
source("R/sim.fossils.R")

# Simulate trees of various sizes, etc
n_tips <- c(50, 100, 200, 500, 1000)
fossil_props <- c(0.05, 0.1, 0.25, 0.5, 0.95)
lambdas <- 1
mus <- c(0.25, .9)
n_sim <- 100

pblapply(n_tips, function(n_tip) {
  lapply(fossil_props, function(fossil_prop) {
    lapply(lambdas, function(lambda) {
      lapply(mus, function(mu) {
        trees <- sim.fbd.taxa.prop(n_tip, fossil_prop, numbsim = n_sim, lambda = lambda, mu = mu)
        lapply(trees, function(tree) {
          tree$n_tip <- n_tip
          tree$fossil_prop <- fossil_prop
          tree$lambda <- lambda
          tree$mu <- mu
        })
      })
    })
  })
})


# TODO: Simulate traits with different models
