#Model fitting with AIC value extraction for PCM sensitivity study
setwd("C:\\Users\\frisk\\Dropbox\\Work in progress\\PCM_sensitivity\\PCM_sensitivity")

#Written by Will Gearty, Bethany Allen, and Pedro Godoy, edited by Alfio Alessandro Chiarenza

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

# Model fitting -------------------------------------------------------


# loading trees
tree_list <- readRDS("./data/tree_simulations.RDS")


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


# information needed for the loops:
mods <- c("wBM", "sBM", "wtrend", "strend",
          "wOUc", "sOUc", "wOUs", "sOUs", "wAC",
          "sAC", "wDC", "sDC")
n_tips <- c(50, 100, 200, 500, 1000)
fossil_props <- c(0.05, 0.1, 0.25, 0.5, 0.95)
lambdas <- 1
mus <- c(0.25, .9)
n_sim <- 100


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

            #IF THERE ARE ISSUES TRY WITH: # set the model name
            #mod <- paste0("n", n_tips[i], "_fp", fossil_props[j], "_l", lambdas[k], "_m", mus[l], "_s", m)

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

            # create a list with the AIC values for the current dataset/model
            obj_list <- list(AIC(get(paste0("fit_", mod, "_BM"))),
                             AIC(get(paste0("fit_", mod, "_trend"))),
                             AIC(get(paste0("fit_", mod, "_OU1"))),
                             AIC(get(paste0("fit_", mod, "_OU2"))),
                             AIC(get(paste0("fit_", mod, "_AC"))),
                             AIC(get(paste0("fit_", mod, "_DC"))))

            # assign names to the list elements
            names(obj_list) <- c("BM", "trend", "OU1", "OU2", "AC", "DC")

            ##Replace tree with trait list and store the AIC values in the model_fitting_results object
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

# save all results in a single file (but we might want to separate by model?)
saveRDS(model_fitting_results, "./data/model_fitting_results.RDS")

