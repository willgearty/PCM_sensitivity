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

# Settings -----------------------------------------------------------------
n_tips <- c(50, 100, 200, 500, 1000)
fossil_props <- c(0.05, 0.1, 0.25, 0.5, 0.95)
lambdas <- 1
mus <- c(0.25, .9)
n_sim <- 100

# Simulate trees -----------------------------------------------------------
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

## getting rid of possible zero-length branches
# by adding 0.00001 to zero-length branches

tree_list <- lapply(tree_list, function(ntip) {
  lapply(ntip, function(fossil_props) {
    lapply(fossil_props, function(lambdas) {
      lapply(lambdas, function(mus) {
        lapply(mus, function(tree) {
          zero <- tree$edge.length == 0
          tree$edge.length[zero] <- tree$edge.length[zero] + 0.00001
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

pb <- txtProgressBar(max = length(n_tips) * length(fossil_props) * length(lambdas) *
                       length(mus) * n_sim,
                     style = 3)
n <- 0
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
          setTxtProgressBar(pb, n)
          n <- n + 1
        }
      }
    }
  }
}

saveRDS(traits, "./data/trait_simulations.RDS")

# A hack to split the results into separate files
mods <- c("wBM", "sBM", "wtrend", "strend",
          "wOUc", "sOUc", "wOUs", "sOUs", "wAC",
          "sAC", "wDC", "sDC")
for (mod in mods) assign(mod, tree_list)

for (mod in mods) {
  for (i in 1:length(n_tips)) {
    for (j in 1:length(fossil_props)) {
      for (k in 1:length(lambdas)) {
        for (l in 1:length(mus)) {
          for (m in 1:n_sim) {
              tmp <- get(mod)
              tmp[[i]][[j]][[k]][[l]][[m]] <- traits[[i]][[j]][[k]][[l]][[m]][[mod]]
              assign(mod, tmp)
          }
        }
      }
    }
  }
}

for (mod in mods) saveRDS(get(mod), paste0("./data/", mod, ".RDS"))


# Model fitting -------------------------------------------------------


# loading trees
tree_list <- readRDS("./data/tree_simulations.RDS")


## getting rid of possible zero-length branches
# by adding 0.00001 to zero-length branches

tree_list <- lapply(tree_list, function(ntip) {
  lapply(ntip, function(fossil_props) {
    lapply(fossil_props, function(lambdas) {
      lapply(lambdas, function(mus) {
        lapply(mus, function(tree) {
          zero <- tree$edge.length == 0
          tree$edge.length[zero] <- tree$edge.length[zero] + 0.00001
          tree
        })
      })
    })
  })
})


# loading traits
wBM_trait  <- readRDS("./data/wBM.RDS")
sBM_trait  <- readRDS("./data/sBM.RDS")

wtrend_trait <- readRDS("./data/wtrend.RDS")
strend_trait  <- readRDS("./data/strend.RDS")

wOUc_trait  <- readRDS("./data/wOUc.RDS")
sOUc_trait  <- readRDS("./data/sOUc.RDS")

wOUs_trait  <- readRDS("./data/wOUs.RDS")
sOUs_trait  <- readRDS("./data/sOUs.RDS")

wAC_trait  <- readRDS("./data/wAC.RDS")
sAC_trait  <- readRDS("./data/sAC.RDS")

wDC_trait  <- readRDS("./data/wDC.RDS")
sDC_trait  <- readRDS("./data/sDC.RDS")


# create empty list to store all results
model_fitting_results <- list()

for (mod in mods) {

# create lists that will be used in the loop
assign(paste0("model_fitting_",mod), tree_list)
simulated_traits <- get(paste0(mod,"_trait"))


## use mvMORPH to fit 6 models (for each of the 12 sets of data simulated with the 12 models)

for (i in 1:length(n_tips)) {
  for (j in 1:length(fossil_props)) {
     for (k in 1:length(lambdas)) {
       for (l in 1:length(mus)) {
         for (m in 1:n_sim) {
          tree <- tree_list[[i]][[j]][[k]][[l]][[m]]
          data <- simulated_traits[[i]][[j]][[k]][[l]][[m]]

          # BM
          assign(paste0("fit_", mod, "_BM"),
                 mvBM(tree = tree, data = data,
                    model = "BM1", method = "rpf"))

          # trend
          assign(paste0("fit_", mod, "_trend"),
                 mvBM(tree = tree, data = data,
                    model = "BM1", method = "rpf",
                    param = list(trend = TRUE)))

          # OU 1 theta
          assign(paste0("fit_", mod, "_OU1"),
                 mvOU(tree = tree, data = data,
                    model="OU1", param=list(root=FALSE)))

          # OU 2 theta
          assign(paste0("fit_", mod, "_OU2"),
                 mvOU(tree = tree, data = data,
                    model="OU1", param=list(root=TRUE)))

          # AC
          assign(paste0("fit_", mod, "_AC"),
                 mvEB(tree = tree, data = data,
                    param=list(up=1)))

          # DC
          assign(paste0("fit_", mod, "_DC"),
                 mvEB(tree = tree, data = data))

# create a list with the results for the current dataset/model
obj_list <- list(get(paste0("fit_", mod, "_BM")),
                 get(paste0("fit_", mod, "_trend")),
                 get(paste0("fit_", mod, "_OU1")),
                 get(paste0("fit_", mod, "_OU2")),
                 get(paste0("fit_", mod, "_AC")),
                 get(paste0("fit_", mod, "_DC")))

# assign names to the list elements
names(obj_list) <- c("BM", "trend", "OU1", "OU2", "AC", "DC")

#Replace tree with trait list
eval(parse(text = paste0("model_fitting_", mod, "[[", i, "]]",
                           "[[", j, "]]", "[[", k, "]]",
                           "[[", l, "]]", "[[", m, "]] <- obj_list")))

        }
      }
    }
  }
}

# add the list to the results list
model_fitting_results[[mod]] <-  get(paste0("model_fitting_", mod))
}


