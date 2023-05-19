#Model fitting with AIC value extraction for PCM sensitivity study
setwd("C:\\Users\\frisk\\Dropbox\\Work in progress\\PCM_sensitivity\\PCM_sensitivity")

#Written by Will Gearty, Bethany Allen, and Pedro Godoy, edited by Alfio Alessandro Chiarenza

# Load libraries
library(ape)
library(phytools)
library(geiger)
library(mvMORPH)
library(pbapply)

# Load functions
source("R/sim.fossils.R")

# Model fitting -------------------------------------------------------


# loading trees
tree_list <- readRDS("./data/tree_simulations.RDS")


# loading traits
wBM_trait  <- readRDS("./data/simulated_traits/wBM.RDS")
sBM_trait  <- readRDS("./data/simulated_traits/sBM.RDS")

wtrend_trait <- readRDS("./data/simulated_traits/wtrend.RDS")
strend_trait  <- readRDS("./data/simulated_traits/strend.RDS")

wOUc_trait  <- readRDS("./data/simulated_traits/wOUc.RDS")
sOUc_trait  <- readRDS("./data/simulated_traits/sOUc.RDS")

wOUs_trait  <- readRDS("./data/simulated_traits/wOUs.RDS")
sOUs_trait  <- readRDS("./data/simulated_traits/sOUs.RDS")

wAC_trait  <- readRDS("./data/simulated_traits/wAC.RDS")
sAC_trait  <- readRDS("./data/simulated_traits/sAC.RDS")

wDC_trait  <- readRDS("./data/simulated_traits/wDC.RDS")
sDC_trait  <- readRDS("./data/simulated_traits/sDC.RDS")


# information needed for the loops:
mods <- c("wBM", "sBM", "wtrend", "strend",
          "wOUc", "sOUc", "wOUs", "sOUs", "wAC",
          "sAC", "wDC", "sDC")
n_tips <- c(50, 100, 200, 500, 1000)
fossil_props <- c(0.05, 0.1, 0.25, 0.5, 0.95)
lambdas <- 1
mus <- c(0.25, .9)
n_sim <- 100

# for each of the trait evolution models used in the simulations
for (mod in mods) {

  # pull simulated trait values
  simulated_traits <- get(paste0(mod,"_trait"))

  # create objects to store outputs
  model_fitting_results <- tree_list
  theta_estimates <- tree_list; sigma_estimates <- tree_list

  ## use mvMORPH to fit 6 models (for each of the 12 sets of data simulated with the 12 models)

  for (i in 1:length(n_tips)) {
    for (j in 1:length(fossil_props)) {
      for (k in 1:length(lambdas)) {
        for (l in 1:length(mus)) {
          for (m in 1:n_sim) {

            # for testing
            # i <- 1; j <- 1; k <- 1; l <- 1; m <- 1

            tree <- tree_list[[i]][[j]][[k]][[l]][[m]]
            data <- simulated_traits[[i]][[j]][[k]][[l]][[m]]

            # BM
            fit_BM <- mvBM(tree = tree, data = data,
                        model = "BM1", method = "rpf",
                        diagnostic = FALSE, echo = FALSE)

            # trend
            fit_trend <- mvBM(tree = tree, data = data,
                        model = "BM1", method = "rpf",
                        param = list(trend = TRUE),
                        diagnostic = FALSE, echo = FALSE)

            # OU 1 theta
            fit_OU1 <- mvOU(tree = tree, data = data,
                        model="OU1", param=list(root=FALSE),
                        diagnostic = FALSE, echo = FALSE)

            # OU 2 theta
            fit_OU2 <- mvOU(tree = tree, data = data,
                        model="OU1", param=list(root=TRUE),
                        diagnostic = FALSE, echo = FALSE)

            # AC
            fit_AC <- mvEB(tree = tree, data = data,
                        param=list(up=1),
                        diagnostic = FALSE, echo = FALSE)

            # DC
            fit_DC <- mvEB(tree = tree, data = data,
                           diagnostic = FALSE, echo = FALSE)


            # create a list with the AIC values
            AIC_list <- list(fit_BM$AIC, fit_trend$AIC, fit_OU1$AIC,
                             fit_OU2$AIC, fit_AC$AIC, fit_DC$AIC)

            # create a list with the theta values
            theta_list <- list(fit_BM$theta, fit_trend$theta, fit_OU1$theta,
                             fit_OU2$theta, fit_AC$theta, fit_DC$theta)

            # create a list with the sigma values
            sigma_list <- list(fit_BM$sigma, fit_trend$sigma, fit_OU1$sigma,
                             fit_OU2$sigma, fit_AC$sigma, fit_DC$sigma)


            # assign names to the list elements
            names(AIC_list) <- c("BM", "trend", "OU1", "OU2", "AC", "DC")
            names(theta_list) <- c("BM", "trend", "OU1", "OU2", "AC", "DC")
            names(sigma_list) <- c("BM", "trend", "OU1", "OU2", "AC", "DC")

            ##Replace tree with AIC or estimated parameters
            model_fitting_results[[i]][[j]][[k]][[l]][[m]] <- AIC_list
            theta_estimates[[i]][[j]][[k]][[l]][[m]] <- theta_list
            sigma_estimates[[i]][[j]][[k]][[l]][[m]] <- sigma_list

          }
        }
      }
    }
  }

  # save results in file named with input model
  saveRDS(model_fitting_results,
          paste0("./data/model_fitting/model_fitting_", mod, ".RDS"))
  saveRDS(theta_estimates,
          paste0("./data/model_fitting/theta_estimates_", mod, ".RDS"))
  saveRDS(sigma_estimates,
          paste0("./data/model_fitting/sigma_estimates_", mod, ".RDS"))
}
