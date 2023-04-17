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

tree_list <- lapply(n_tips, function(n_tip) {
  lapply(fossil_props, function(fossil_prop) {
    lapply(lambdas, function(lambda) {
      lapply(mus, function(mu) {
        trees <- sim.fbd.taxa.prop(n_tip, fossil_prop, numbsim = n_sim, lambda = lambda, mu = mu)
        lapply(trees, function(tree) {
          tree$n_tip <- n_tip
          tree$fossil_prop <- fossil_prop
          tree$lambda <- lambda
          tree$mu <- mu
          tree
        })
      })
    })
  })
})

saveRDS(tree_list, "./data/tree_simulations.RDS")

# Simulate traits -------------------------------------------------------

tree_list <- readRDS("./data/tree_simulations.RDS")

traits <- tree_list

for (i in 1:length(n_tips)) {
  for (j in 1:length(fossil_props)) {
    for (k in 1:length(lambdas)) {
      for (l in 1:length(mus)) {
        for (m in 1:n_sim) {
          tree <- tree_list[[i]][[j]][[k]][[l]][[m]]

          #Simulate weak BM trait evolution
          wBM_trait <- mvSIM(tree = tree, nsim = 1, model = "BM1",
                             param = list(trend = FALSE,
                                          theta = 0, #ancestral state
                                          sigma = 0.1 #strength of drift
                             ))

          #Simulate strong BM trait evolution
          sBM_trait <- mvSIM(tree = tree, nsim = 1, model = "BM1",
                             param = list(trend = FALSE,
                                          theta = 0, #ancestral state
                                          sigma = 0.5 #strength of drift
                             ))

          #Simulate weak trended BM trait evolution
          wtrend_trait <- mvSIM(tree = tree, nsim = 1, model = "BM1",
                                param = list(trend = 0.1, #strength of trend
                                             theta = 0, #ancestral state
                                             sigma = 0.1 #strength of drift
                                ))

          #Simulate strong trended BM trait evolution
          strend_trait <- mvSIM(tree = tree, nsim = 1, model = "BM1",
                                param = list(trend = 0.3, #strength of trend
                                             theta = 0, #ancestral state
                                             sigma = 0.1 #strength of drift
                                ))

          #Simulate weak OU ("SSP") trait evolution, centred on ancestral state
          wOUc_trait <- mvSIM(tree = tree, nsim = 1, model = "OU1",
                              param = list(alpha = 0.01, #strength of selection
                                           theta = 0, #ancestral state
                                           sigma = 0.1 #strength of drift
                              ))

          #Simulate strong OU ("SSP") trait evolution, centred on ancestral state
          sOUc_trait <- mvSIM(tree = tree, nsim = 1, model = "OU1",
                              param = list(alpha = 10, #strength of selection
                                           theta = 0, #ancestral state
                                           sigma = 0.1 #strength of drift
                              ))

          #Simulate weak OU ("SSP") trait evolution, with shifted optimum
          wOUs_trait <- mvSIM(tree = tree, nsim = 1, model = "OU1",
                              param = list(root = TRUE,
                                           alpha = 0.01, #strength of selection
                                           theta = c(0, 1), #ancestral state,
                                           #optimum
                                           sigma = 0.1 #strength of drift
                              ))

          #Simulate strong OU ("SSP") trait evolution, with shifted optimum
          sOUs_trait <- mvSIM(tree = tree, nsim = 1, model = "OU1",
                              param = list(root = TRUE,
                                           alpha = 10, #strength of selection
                                           theta = c(0, 1), #ancestral state,
                                           #optimum
                                           sigma = 0.1 #strength of drift
                              ))

          #Simulate weak AC trait evolution
          wAC_trait <- mvSIM(tree = tree, nsim = 1, model = "EB",
                             param = list(theta = 0, #ancestral state
                                          beta = 0.1, #exponential rate
                                          sigma = 0.001 #strength of drift
                             ))

          #Simulate strong AC trait evolution
          sAC_trait <- mvSIM(tree = tree, nsim = 1, model = "EB",
                             param = list(theta = 0, #ancestral state
                                          beta = 0.3, #exponential rate
                                          sigma = 0.001 #strength of drift
                             ))

          #Simulate weak DC ("Early Burst") trait evolution
          wDC_trait <- mvSIM(tree = tree, nsim = 1, model = "EB",
                             param = list(theta = 0, #ancestral state
                                          beta = -0.1, #exponential rate
                                          sigma = 0.001 #strength of drift
                             ))

          #Simulate strong DC ("Early Burst") trait evolution
          sDC_trait <- mvSIM(tree = tree, nsim = 1, model = "EB",
                             param = list(theta = 1, #ancestral state
                                          beta = -0.3, #exponential
                                          sigma = 0.001 #strength of drift
                             ))

          #Paste together trait simulations
          trait_list <- list(wBM_trait, sBM_trait, wtrend_trait,
                             strend_trait, wOUc_trait, sOUc_trait,
                             wOUs_trait, sOUs_trait, wAC_trait,
                             sAC_trait, wDC_trait, sDC_trait)
          names(trait_list) <- c("wBM", "sBM", "wtrend", "strend",
                                 "wOUc", "sOUc", "wOUs", "sOUs", "wAC",
                                 "sAC", "wDC", "sDC")

          #Replace tree with trait list
          traits[[i]][[j]][[k]][[l]][[m]] <- trait_list
        }
      }
    }
  }
}

saveRDS(traits, "./data/trait_simulations.RDS")
