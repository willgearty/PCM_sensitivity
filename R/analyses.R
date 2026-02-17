# Written by Will Gearty, Bethany Allen, Pedro Godoy, and Alfio Alessandro Chiarenza

# Load libraries
#devtools::install_github("willgearty/pcmtools")
library(ape); library(phytools); library(geiger); library(TreeSim)
library(FossilSim); library(mvMORPH); library(pbapply); library(dplyr)
library(tibble); library(tidyr); library(ggplot2); library(pcmtools)
library(deeptime); library(future); library(future.apply)
library(forcats); library(ggh4x); library(data.table)

# Load functions
source("R/sim.fossils.poisson.R") # uses data.table::rbindlist for speed
source("R/sim.fossils.R")

# Download trees if missing ----------------------------------------------------
# Normal trees
# https://1drv.ms/u/s!ArhYkoKadYP1gP91QcgXJtFRpCnYmA?e=hg8ylU

# Scaled trees
# https://1drv.ms/u/s!ArhYkoKadYP1gP92kMLrgsNOjzGKBQ?e=eK3eFQ

# Settings -----------------------------------------------------------------
n_tips <- c(50, 100, 200, 500, 1000)
fossil_props <- c(0, 0.1, 0.25, 0.5, 0.95)
lambdas <- 1
mus <- c(0.5, 0.9)
models <- list("root" = function(x) 10 ^ -(x - 1) - 1,
               "random" = function(x) 1,
               "recent" = function(x) 100 ^ x - 1)
n_sim <- 100

# get all unique combinations of parameters
settings <- expand.grid(n_tip = n_tips, fossil_prop = fossil_props,
                        lambda = lambdas, mu = mus, model = names(models))

plan(multisession, workers = 8)
# Simulate trees -----------------------------------------------------------
set.seed(1234)
tree_df <- future_mapply(function(n_tip, fossil_prop, lambda, mu, model) {
  trees <- sim.fbd.taxa.prop(n_tip, fossil_prop, numbsim = n_sim,
                             lambda = lambda, mu = mu,
                             model = models[[model]], progress = FALSE)
  tmp <- data.frame(n_tip = rep(n_tip, n_sim),
                    fossil_prop = rep(fossil_prop, n_sim),
                    lambda = rep(lambda, n_sim),
                    mu = rep(mu, n_sim),
                    beta = model,
                    sim = seq_len(n_sim))
  tmp$tree <- trees
  tmp
}, n_tip = settings$n_tip, fossil_prop = settings$fossil_prop,
lambda = settings$lambda, mu = settings$mu, model = settings$model,
SIMPLIFY = FALSE, future.seed = TRUE) %>% bind_rows()

## getting rid of possible zero-length branches
# by adding 0.00001 to zero-length branches
tree_df$tree <- lapply(tree_df$tree, function(tree) {
  zero <- tree$edge.length == 0
  tree$edge.length[zero] <- tree$edge.length[zero] + 0.00001
  tree
})
saveRDS(tree_df, "./data/tree_simulations.RDS")

# rescale trees to height of 1
# not sure why, but geiger::rescale wasn't working here
rescaleTree <- function(tree, scale) {
  tree$edge.length <- tree$edge.length / max(nodeHeights(tree)[, 2]) * scale
  return(tree)
}
tree_df_scaled <- tree_df
tree_df_scaled$tree <- lapply(tree_df$tree, function(tree) {
  rescaleTree(tree, 1)
})

saveRDS(tree_df_scaled, "./data/tree_simulations_scaled.RDS")

# Simulate traits -------------------------------------------------------
tree_df <- readRDS("./data/tree_simulations.RDS")

set.seed(1234)
wBM_trait <- pblapply(tree_df$tree, function(tree) {
  #Simulate weak BM trait evolution
  mvSIM(tree = tree, nsim = 1, model = "BM1",
        param = list(trend = FALSE,
                     theta = 0, #ancestral state
                     sigma = 0.1 #strength of drift
        ))
})

sBM_trait <- pblapply(tree_df$tree, function(tree) {
  #Simulate strong BM trait evolution
  mvSIM(tree = tree, nsim = 1, model = "BM1",
        param = list(trend = FALSE,
                     theta = 0, #ancestral state
                     sigma = 0.5 #strength of drift
        ))
})

wtrend_trait <- pblapply(tree_df$tree, function(tree) {
  #Simulate weak trended BM trait evolution
  mvSIM(tree = tree, nsim = 1, model = "BM1",
        param = list(trend = 0.1, #strength of trend
                     theta = 0, #ancestral state
                     sigma = 0.1 #strength of drift
        ))
})

strend_trait <- pblapply(tree_df$tree, function(tree) {
  #Simulate strong trended BM trait evolution
  mvSIM(tree = tree, nsim = 1, model = "BM1",
        param = list(trend = 0.3, #strength of trend
                     theta = 0, #ancestral state
                     sigma = 0.1 #strength of drift
        ))
})

wOUc_trait <- pblapply(tree_df$tree, function(tree) {
  max_height <- max(nodeHeights(tree))
  #Simulate weak OU trait evolution, centered on ancestral state ("SSP")
  mvSIM(tree = tree, nsim = 1, model = "OU1",
        param = list(alpha = log(2) / max_height, #strength of selection
                     theta = 0, #ancestral state
                     sigma = 0.1 #strength of drift
        ))
})

sOUc_trait <- pblapply(tree_df$tree, function(tree) {
  max_height <- max(nodeHeights(tree))
  #Simulate strong OU ("SSP") trait evolution, centered on ancestral state ("SSP")
  mvSIM(tree = tree, nsim = 1, model = "OU1",
        param = list(alpha = log(2) / (max_height / 5), #strength of selection
                     theta = 0, #ancestral state
                     sigma = 0.1 #strength of drift
        ))
})

wOUs_trait <- pblapply(tree_df$tree, function(tree) {
  max_height <- max(nodeHeights(tree))
  #Simulate weak OU trait evolution, with shifted optimum
  mvSIM(tree = tree, nsim = 1, model = "OU1",
        param = list(root = TRUE,
                     alpha = log(2) / max_height, #strength of selection
                     theta = c(0, 2), #ancestral state, optimum
                     sigma = 0.1 #strength of drift
        ))
})

sOUs_trait <- pblapply(tree_df$tree, function(tree) {
  max_height <- max(nodeHeights(tree))
  #Simulate strong OU trait evolution, with shifted optimum
  mvSIM(tree = tree, nsim = 1, model = "OU1",
        param = list(root = TRUE,
                     alpha = log(2) / (max_height / 5), #strength of selection
                     theta = c(0, 2), #ancestral state, optimum
                     sigma = 0.1 #strength of drift
        ))
})

wAC_trait <- pblapply(tree_df$tree, function(tree) {
  #Simulate weak AC trait evolution
  mvSIM(tree = tree, nsim = 1, model = "EB",
        param = list(theta = 0, #ancestral state
                     beta = 0.1, #exponential rate
                     sigma = 0.001 #strength of drift
        ))
})

sAC_trait <- pblapply(tree_df$tree, function(tree) {
  #Simulate strong AC trait evolution
  mvSIM(tree = tree, nsim = 1, model = "EB",
        param = list(theta = 0, #ancestral state
                     beta = 0.3, #exponential rate
                     sigma = 0.001 #strength of drift
        ))
})

wDC_trait <- pblapply(tree_df$tree, function(tree) {
  #Simulate weak DC ("Early Burst") trait evolution
  mvSIM(tree = tree, nsim = 1, model = "EB",
        param = list(theta = 0, #ancestral state
                     beta = -0.1, #exponential rate
                     sigma = 0.001 #strength of drift
        ))
})

sDC_trait <- pblapply(tree_df$tree, function(tree) {
  #Simulate strong DC ("Early Burst") trait evolution
  mvSIM(tree = tree, nsim = 1, model = "EB",
        param = list(theta = 1, #ancestral state
                     beta = -0.3, #exponential
                     sigma = 0.001 #strength of drift
        ))
})

# Save trait simulation results
mods <- c("wBM", "sBM", "wtrend", "strend",
          "wOUc", "sOUc", "wOUs", "sOUs",
          "wAC", "sAC", "wDC", "sDC")
for (mod in mods) saveRDS(get(paste0(mod, "_trait")),
                          paste0("./data/simulated_traits/", mod, ".RDS"))


# Model fitting -------------------------------------------------------
# loading trees
tree_df <- readRDS("./data/tree_simulations.RDS")

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

# information needed for the loop:
mods <- c("wBM", "sBM", "wtrend", "strend",
          "wOUc", "sOUc", "wOUs", "sOUs",
          "wAC", "sAC", "wDC", "sDC")

plan(multisession, workers = 8) # set up parallel processing
# for each of the trait evolution models used in the simulations
for (mod in mods) {
  print(paste("Fitting models to", mod, "simulations"))
  # pull simulated trait values
  simulated_traits <- get(paste0(mod,"_trait"))

  model_fitting_results <- pblapply(seq_len(nrow(tree_df)), function(tree_i) {
    tryMV <- function(expr) {
      lst <- tryCatch(expr,
                      error = function(e) {
                        cat(paste0("Error: ", conditionMessage(e)), "\n")
                        cat(paste0("tree index: ", tree_i, "\n"))
                        cat(paste0("original model: ", mod, "\n"))
                        list(AIC = NA, AICc = NA, theta = NA, alpha = NA, beta = NA,
                             sigma = NA, trend = NA, error = e)
                      }
      )
      # free up some memory by removing some things
      lst$LogLik <- NULL
      lst$hess.values <- NULL
      lst$param <- NULL
      lst$llik <- NULL # this function especially takes up so much unnecessary memory
      unclass(lst)
    }
    sim_trait <- simulated_traits[[tree_i]]
    tree <- tree_df$tree[[tree_i]]
    # BM
    fit_BM <- tryMV(mvBM(tree = tree, data = sim_trait,
                         model = "BM1", method = "rpf",
                         diagnostic = FALSE, echo = FALSE))

    # trend
    fit_trend <- tryMV(mvBM(tree = tree, data = sim_trait,
                            model = "BM1", method = "rpf",
                            param = list(trend = TRUE),
                            diagnostic = FALSE, echo = FALSE))

    # OU 1 theta
    fit_OU1 <- tryMV(mvOU(tree = tree, data = sim_trait,
                          model="OU1", param=list(root=FALSE),
                          diagnostic = FALSE, echo = FALSE))

    # OU 2 theta
    fit_OU2 <- tryMV(mvOU(tree = tree, data = sim_trait,
                          model="OU1", param=list(root=TRUE),
                          diagnostic = FALSE, echo = FALSE))

    # ACDC
    fit_ACDC <- tryMV(mvEB(tree = tree, data = sim_trait,
                           param=list(up=1),
                           diagnostic = FALSE, echo = FALSE))

    list(BM = fit_BM, trend = fit_trend, OU1 = fit_OU1,
         OU2 = fit_OU2, ACDC = fit_ACDC)
  }, cl = "future", future.seed = TRUE)
  saveRDS(model_fitting_results,
          paste0("./data/model_fitting/model_fitting_", mod, ".RDS"))
}
plan(sequential)

# Analyze results ---------------------------------------------------
model_results <- list()
mods <- c("wBM", "sBM", "wtrend", "strend",
          "wOUc", "sOUc", "wOUs", "sOUs",
          "wAC", "sAC", "wDC", "sDC")
for (mod in mods) {
  model_results[[mod]] <- readRDS(paste0("./data/model_fitting/model_fitting_", mod, ".RDS"))
}

# read in tree_df
tree_df <- readRDS("./data/tree_simulations.RDS")

# extract to data.frames
model_fits_df <- lapply(model_results, \(mod) {
  lapply(mod, \(tree) {
    sapply(tree, \(fit_model) {
      fit_model$AICc
    })
  }) %>% bind_rows() %>% cbind(tree_df %>% select(-tree), .)
}) %>% bind_rows(.id = "model")
colnames(model_fits_df)[10:11] <- c("OUc", "OUs")

# remove the simulations that we don't want
# 1. OUs and trend models without fossils
# 2. duplicate fossil-less simulations

model_fits_df_filt <- model_fits_df %>%
  filter(!(model %in% c("wOUs", "sOUs", "wtrend", "strend") & fossil_prop == 0)) %>%
  filter(!(fossil_prop == 0 & beta != "random"))

fit_models <- c("BM", "trend", "OUc", "OUs", "ACDC")
model_fits_df_long <- model_fits_df_filt %>%
  pivot_longer(cols = all_of(fit_models), names_to = "fit_model", values_to = "aicc") %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta, sim) %>%
  mutate(aicc_w = AICweights(aicc), aicc_d = aicc - min(aicc)) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = mods),
         fit_model = factor(fit_model, levels = fit_models),
         across(c(n_tip, fossil_prop, lambda, mu, sim, aicc_w), ~ as.numeric(.x))) %>%
  mutate(across(c(n_tip, fossil_prop, lambda, mu, beta), ~ as.factor(.x)))

param_estimates_df <- lapply(model_results, \(mod) {
  lapply(mod, \(tree) {
    lapply(tree, \(fit_model) {
      thetas <- rep_len(fit_model$theta, 2)
      bind_cols(alpha = fit_model$alpha[[1]], exp_rate = fit_model$beta[[1]],
                sigma = fit_model$sigma[[1]], trend = fit_model$trend[[1]],
                theta_0 = thetas[1], theta_1 = thetas[2],
                convergence = fit_model$convergence)
    }) %>% bind_rows(.id = "fit_model")
  }) %>% bind_rows() %>% cbind(tree_df %>% select(-tree) %>% slice(rep(1:n(), each = 5)), .)
}) %>% bind_rows(.id = "model")

param_estimates_df_clean <- param_estimates_df %>%
  mutate(fit_model = fct_recode(fit_model, OUc = "OU1", OUs = "OU2")) %>%
  mutate(model = factor(model, levels = mods),
         fit_model = factor(fit_model, levels = fit_models),
         across(c(n_tip, fossil_prop, lambda, mu, sim), ~ as.numeric(.x))) %>%
  mutate(across(c(n_tip, fossil_prop, lambda, mu, beta), ~ as.factor(.x))) %>%
  mutate(correct_model = case_when(
    model %in% c("wBM", "sBM") ~ "BM",
    model %in% c("wtrend", "strend") ~ "trend",
    model %in% c("wOUc", "sOUc") ~ "OUc",
    model %in% c("wOUs", "sOUs") ~ "OUs",
    model %in% c("wAC", "sAC", "wDC", "sDC") ~ "ACDC"
  )) %>%
  filter(fit_model == correct_model) %>%
  select(-correct_model) %>%
  mutate(rel_hl = (log(2) / alpha) / sapply(tree_df$tree, function(tree) max(nodeHeights(tree)))) %>%
  filter(!(model %in% c("wOUs", "sOUs", "wtrend", "strend") & fossil_prop == 0)) %>%
  filter(!(fossil_prop == 0 & beta != "random"))

# need to do some pivoting to get the two different theta values for each simulation
theta_estimates_df_long <- param_estimates_df_clean %>%
  pivot_longer(cols = c(theta_0, theta_1), names_to = "theta", values_to = "theta_val",
               names_prefix = "theta_") %>%
  filter(grepl("OU", model)) %>%
  mutate(model = factor(model, levels = c("wOUc", "sOUc", "wOUs", "sOUs")))

# clean up the old object
remove(model_results)

## AIC plot ------------------------------------------------------
gg1 <- ggplot(model_fits_df_long %>% filter(mu == "0.5")) +
  geom_violin(aes(x = factor(fossil_prop), y = aicc_w, color = fit_model)) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(rows = vars(n_tip), cols = vars(model)) +
  theme_bw(base_size = 20)
ggsave("./figures/AIC.pdf", gg1, width = 40, height = 20)

## best model plots -----------------------------------------------------
# proportion of best fitting models that are the correct model
model_fits_df_summ <- model_fits_df_long %>%
  mutate(correct_model = case_when(
    model %in% c("wBM", "sBM") ~ "BM",
    model %in% c("wtrend", "strend") ~ "trend",
    model %in% c("wOUc", "sOUc") ~ "OUc",
    model %in% c("wOUs", "sOUs") ~ "OUs",
    model %in% c("wAC", "sAC", "wDC", "sDC") ~ "ACDC"
  )) %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta, sim) %>%
  summarise(correct = ifelse(all(is.na(aicc_w)), NA,
                             any(unique(correct_model) ==
                                   fit_model[which(aicc_w == max(aicc_w, na.rm = TRUE))])),
            .groups = "drop") %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta) %>%
  summarise(prop_true = sum(correct)/n(), .groups = "drop") %>%
  mutate(beta = fct_recode(as.factor(beta), `root-biased` = "root", `random` = "random", `recent-biased` = "recent"))

gg2a <- ggplot(model_fits_df_summ %>% filter(mu == 0.5)) +
  geom_line(aes(x = fossil_prop, y = prop_true, color = n_tip,
                linetype = beta, group = interaction(n_tip, beta))) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations with Correct Best Fit Model", limits = c(0, 1)) +
  scale_color_brewer("# of tips", palette = "Dark2") +
  scale_linetype_discrete("Fossil Distribution") +
  theme_bw(base_size = 20) +
  facet_wrap(~model)
gg2b <- ggplot(model_fits_df_summ %>% filter(mu == 0.9)) +
  geom_line(aes(x = fossil_prop, y = prop_true, color = n_tip,
                linetype = beta, group = interaction(n_tip, beta))) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations with Correct Best Fit Model", limits = c(0, 1)) +
  scale_color_brewer("# of tips", palette = "Dark2") +
  scale_linetype_discrete("Fossil Distribution") +
  theme_bw(base_size = 20) +
  facet_wrap(~model)
gg2 <- ggarrange2(gg2a, gg2b, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Prop_Best.pdf", gg2, width = 16, height = 20)
ggsave("./figures/Prop_Best_50.pdf", gg2a, width = 16, height = 10)
ggsave("./figures/Prop_Best_90.pdf", gg2b, width = 16, height = 10)

# generating models for best fitting models
model_fits_df_summ2 <- model_fits_df_long %>%
  mutate(correct_model = case_when(
    model %in% c("wBM", "sBM") ~ "BM",
    model %in% c("wtrend", "strend") ~ "trend",
    model %in% c("wOUc", "sOUc") ~ "OUc",
    model %in% c("wOUs", "sOUs") ~ "OUs",
    model %in% c("wAC", "sAC", "wDC", "sDC") ~ "ACDC"
  )) %>%
  mutate(correct_model = factor(correct_model, levels = fit_models)) %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta, sim) %>%
  slice_max(aicc_w) %>%
  ungroup() %>%
  group_by(n_tip, fossil_prop, lambda, mu, beta, fit_model) %>%
  count(correct_model) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(beta = fct_recode(as.factor(beta), `root-biased` = "root", `random` = "random", `recent-biased` = "recent"))

gg2c <- ggplot(model_fits_df_summ2 %>% filter(mu == 0.5)) +
  geom_line(aes(x = fossil_prop, y = prop, color = correct_model,
                linetype = beta, group = interaction(correct_model, beta))) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations", limits = c(0, 1)) +
  scale_color_brewer("Simulated Model", palette = "Dark2") +
  scale_linetype_discrete("Fossil Distribution") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(fit_model), rows = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg2d <- ggplot(model_fits_df_summ2 %>% filter(mu == 0.9)) +
  geom_line(aes(x = fossil_prop, y = prop, color = correct_model,
                linetype = beta, group = interaction(correct_model, beta))) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations", limits = c(0, 1)) +
  scale_color_brewer("Simulated Model", palette = "Dark2") +
  scale_linetype_discrete("Fossil Distribution") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(fit_model), rows = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg2_b <- ggarrange2(gg2c, gg2d, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Prop_Sim.pdf", gg2_b, width = 16, height = 25)
ggsave("./figures/Prop_Sim_50.pdf", gg2c, width = 16, height = 12.5)
ggsave("./figures/Prop_Sim_90.pdf", gg2d, width = 16, height = 12.5)

gg2c_n <- ggplot(model_fits_df_summ2 %>% filter(mu == 0.5)) +
  geom_line(aes(x = fossil_prop, y = n, color = correct_model,
                linetype = beta, group = interaction(correct_model, beta))) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Number of Simulations") +
  scale_color_brewer("Simulated Model", palette = "Dark2") +
  scale_linetype_discrete("Fossil Distribution") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(fit_model), rows = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg2d_n <- ggplot(model_fits_df_summ2 %>% filter(mu == 0.9)) +
  geom_line(aes(x = fossil_prop, y = n, color = correct_model,
                linetype = beta, group = interaction(correct_model, beta))) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Number of Simulations") +
  scale_color_brewer("Simulated Model", palette = "Dark2") +
  scale_linetype_discrete("Fossil Distribution") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(fit_model), rows = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg2_b_n <- ggarrange2(gg2c_n, gg2d_n, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/N_Sim.pdf", gg2_b_n, width = 16, height = 25)
ggsave("./figures/N_Sim_50.pdf", gg2c_n, width = 16, height = 12.5)
ggsave("./figures/N_Sim_90.pdf", gg2d_n, width = 16, height = 12.5)

# when wrong answer, what is it?
model_fits_df_summ3 <- model_fits_df_long %>%
  mutate(correct_model = case_when(
    model %in% c("wBM", "sBM") ~ "BM",
    model %in% c("wtrend", "strend") ~ "trend",
    model %in% c("wOUc", "sOUc") ~ "OUc",
    model %in% c("wOUs", "sOUs") ~ "OUs",
    model %in% c("wAC", "sAC", "wDC", "sDC") ~ "ACDC"
  )) %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta, sim) %>%
  mutate(correct = ifelse(all(is.na(aicc_w)), NA,
                          any(unique(correct_model) ==
                                fit_model[which(aicc_w == max(aicc_w, na.rm = TRUE))]))) %>%
  filter(!correct) %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta, sim) %>%
  arrange(-aicc_w) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta) %>%
  count(fit_model) %>%
  mutate(prop_wrong = n / sum(n), prop_all = n / 100) %>%
  ungroup() %>%
  mutate(beta = fct_recode(as.factor(beta), `root-biased` = "root", `random` = "random", `recent-biased` = "recent"))

gg2e <- ggplot(model_fits_df_summ3 %>% filter(mu == 0.5)) +
  geom_col(data = . %>% filter(beta == "root-biased"),
           aes(x = as.numeric(fossil_prop) - .25, y = n, fill = fit_model, color = "blue"),
           position = "stack", width = .2) +
  geom_col(data = . %>% filter(beta == "random"),
           aes(x = as.numeric(fossil_prop), y = n, fill = fit_model, color = "green"),
           position = "stack", width = .2) +
  geom_col(data = . %>% filter(beta == "recent-biased"),
           aes(x = as.numeric(fossil_prop) + .25, y = n, fill = fit_model, color = "red"),
           position = "stack", width = .2) +
  scale_x_continuous("Proportion of Fossils in Tree", breaks = 1:5, labels = levels(model_fits_df_summ3$fossil_prop)) +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Fit Model", palette = "Dark2") +
  scale_color_identity("Fossil Distribution", guide = guide_legend(), labels = c("root-biased", "random", "recent-biased")) +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(model), rows = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))

gg2f <- ggplot(model_fits_df_summ3 %>% filter(mu == 0.9)) +
  geom_col(data = . %>% filter(beta == "root-biased"),
           aes(x = as.numeric(fossil_prop) - .25, y = n, fill = fit_model, color = "blue"),
           position = "stack", width = .2) +
  geom_col(data = . %>% filter(beta == "random"),
           aes(x = as.numeric(fossil_prop), y = n, fill = fit_model, color = "green"),
           position = "stack", width = .2) +
  geom_col(data = . %>% filter(beta == "recent-biased"),
           aes(x = as.numeric(fossil_prop) + .25, y = n, fill = fit_model, color = "red"),
           position = "stack", width = .2) +
  scale_x_continuous("Proportion of Fossils in Tree", breaks = 1:5, labels = levels(model_fits_df_summ3$fossil_prop)) +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Fit Model", palette = "Dark2") +
  scale_color_identity("Fossil Distribution", guide = guide_legend(), labels = c("root-biased", "random", "recent-biased")) +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(model), rows = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))

gg2_c <- ggarrange2(gg2e, gg2f, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Wrong_Best2.pdf", gg2_c, width = 40, height = 40)
ggsave("./figures/Wrong_Best2_50.pdf", gg2e, width = 40, height = 20)
ggsave("./figures/Wrong_Best2_90.pdf", gg2f, width = 40, height = 20)

gg2e2 <- ggplot(model_fits_df_summ3 %>% filter(mu == 0.5)) +
  geom_col(aes(x = fossil_prop, y = n, fill = fit_model), position = "stack") +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Best Fit Model", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_nested(beta + n_tip ~ model,
               labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg2f2 <- ggplot(model_fits_df_summ3 %>% filter(mu == 0.9)) +
  geom_col(aes(x = fossil_prop, y = n, fill = fit_model), position = "stack") +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Best Fit Model", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_nested(beta + n_tip ~ model,
               labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg2_c2 <- ggarrange2(gg2e2, gg2f2, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Wrong_Best.pdf", gg2_c2, width = 40, height = 40)
ggsave("./figures/Wrong_Best_50.pdf", gg2e2, width = 40, height = 20)
ggsave("./figures/Wrong_Best_90.pdf", gg2f2, width = 40, height = 20)

gg2e3 <- ggplot(model_fits_df_summ3 %>% filter(mu == 0.5, model %in% c("wBM", "sBM", "wtrend", "strend"))) +
  geom_col(aes(x = fossil_prop, y = n, fill = fit_model), position = "stack") +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Best Fit Model", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_nested(beta + n_tip ~ model,
               labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg2f3 <- ggplot(model_fits_df_summ3 %>% filter(mu == 0.9, model %in% c("wBM", "sBM", "wtrend", "strend"))) +
  geom_col(aes(x = fossil_prop, y = n, fill = fit_model), position = "stack") +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Best Fit Model", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_nested(beta + n_tip ~ model,
               labeller = labeller(n_tip = function(x) paste(x, "tips")))
ggsave("./figures/Wrong_Best_BM_50.pdf", gg2e3, width = 16, height = 20)
ggsave("./figures/Wrong_Best_BM_90.pdf", gg2f3, width = 16, height = 20)

gg2e4 <- ggplot(model_fits_df_summ3 %>% filter(mu == 0.5, model %in% c("wOUc", "sOUc", "wOUs", "sOUs"))) +
  geom_col(aes(x = fossil_prop, y = n, fill = fit_model), position = "stack") +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Best Fit Model", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_nested(beta + n_tip ~ model,
               labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg2f4 <- ggplot(model_fits_df_summ3 %>% filter(mu == 0.9, model %in% c("wOUc", "sOUc", "wOUs", "sOUs"))) +
  geom_col(aes(x = fossil_prop, y = n, fill = fit_model), position = "stack") +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Best Fit Model", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_nested(beta + n_tip ~ model,
               labeller = labeller(n_tip = function(x) paste(x, "tips")))
ggsave("./figures/Wrong_Best_OU_50.pdf", gg2e4, width = 16, height = 20)
ggsave("./figures/Wrong_Best_OU_90.pdf", gg2f4, width = 16, height = 20)

gg2e5 <- ggplot(model_fits_df_summ3 %>% filter(mu == 0.5, model %in% c("wAC", "sAC", "wDC", "sDC"))) +
  geom_col(aes(x = fossil_prop, y = n, fill = fit_model), position = "stack") +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Best Fit Model", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_nested(beta + n_tip ~ model,
               labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg2f5 <- ggplot(model_fits_df_summ3 %>% filter(mu == 0.9, model %in% c("wAC", "sAC", "wDC", "sDC"))) +
  geom_col(aes(x = fossil_prop, y = n, fill = fit_model), position = "stack") +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Best Fit Model", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_nested(beta + n_tip ~ model,
               labeller = labeller(n_tip = function(x) paste(x, "tips")))
ggsave("./figures/Wrong_Best_ACDC_50.pdf", gg2e5, width = 16, height = 20)
ggsave("./figures/Wrong_Best_ACDC_90.pdf", gg2f5, width = 16, height = 20)

ggplot(model_fits_df_summ3 %>% filter(mu == 0.9)) +
  geom_line(aes(x = fossil_prop, y = n, color = fit_model, linetype = beta,
                group = interaction(fit_model, beta))) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_color_brewer("Fit Model", palette = "Dark2") +
  scale_linetype_discrete("Fossil Distribution") +
  theme_bw(base_size = 20) +
  facet_grid(n_tip ~ model,
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
ggsave("./figures/test2.pdf", width = 40, height = 20)

# summarize proportion of correct models
gg2g <- ggplot(model_fits_df_summ %>% filter(mu == 0.5)) +
  geom_col(aes(x = fossil_prop, y = prop_true / 12, fill = model)) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations Correctly Identified",
                     limits = c(0, 1)) +
  scale_fill_brewer("Simulated Model", palette = "Paired") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(n_tip), rows = vars(beta),
             labeller = labeller(n_tip = function(x) paste(x, "tips"),
                                 beta = c("root" = "root-biased",
                                          "random" = "random",
                                          "recent" = "recent-biased")))

gg2h <- ggplot(model_fits_df_summ %>% filter(mu == 0.9)) +
  geom_col(aes(x = fossil_prop, y = prop_true / 12, fill = model)) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations Correctly Identified",
                     limits = c(0, 1)) +
  scale_fill_brewer("Simulated Model", palette = "Paired") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(n_tip), rows = vars(beta),
             labeller = labeller(n_tip = function(x) paste(x, "tips"),
                                 beta = c("root" = "root-biased",
                                          "random" = "random",
                                          "recent" = "recent-biased")))

gg2_d <- ggarrange2(gg2g, gg2h, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Prop_Best_Stacked.pdf", gg2_d, width = 18.53, height = 20)
ggsave("./figures/Prop_Best_Stacked_50.pdf", gg2g, width = 18.53, height = 10)
ggsave("./figures/Prop_Best_Stacked_90.pdf", gg2h, width = 18.53, height = 10)

gg2i <- ggplot(model_fits_df_summ %>% filter(mu == 0.5)) +
  geom_col(aes(x = fossil_prop, y = prop_true / 12)) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations Correctly Identified",
                     limits = c(0, 1)) +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(n_tip), rows = vars(beta),
             labeller = labeller(n_tip = function(x) paste(x, "tips"),
                                 beta = c("root" = "root-biased",
                                          "random" = "random",
                                          "recent" = "recent-biased")))

gg2j <- ggplot(model_fits_df_summ %>% filter(mu == 0.9)) +
  geom_col(aes(x = fossil_prop, y = prop_true / 12)) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations Correctly Identified",
                     limits = c(0, 1)) +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(n_tip), rows = vars(beta),
             labeller = labeller(n_tip = function(x) paste(x, "tips"),
                                 beta = c("root" = "root-biased",
                                          "random" = "random",
                                          "recent" = "recent-biased")))

gg2_e <- ggarrange2(gg2i, gg2j, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Prop_Best_Combined.pdf", gg2_e, width = 16, height = 20)
ggsave("./figures/Prop_Best_Combined_50.pdf", gg2i, width = 16, height = 10)
ggsave("./figures/Prop_Best_Combined_90.pdf", gg2j, width = 16, height = 10)

# proportions of simulations with clear best model
model_fits_df_summ4 <- model_fits_df_long %>%
  mutate(correct = case_when(
    model %in% c("wBM", "sBM") ~ "BM",
    model %in% c("wtrend", "strend") ~ "trend",
    model %in% c("wOUc", "sOUc") ~ "OUc",
    model %in% c("wOUs", "sOUs") ~ "OUs",
    model %in% c("wAC", "sAC", "wDC", "sDC") ~ "ACDC"
  ) == fit_model) %>%
  mutate(best = ifelse(is.na(aicc_d), FALSE, aicc_d == 0)) %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta, sim) %>%
  summarise(best_model = fit_model[which(aicc_d == 0)],
            correct = ifelse(any(correct & best), "correct", "incorrect"),
            clear = ifelse(all(aicc_d[!is.na(aicc_d) & !best] > 2), "clear", "unclear"),
            .groups = "drop") %>%
  mutate(cor_clear = factor(interaction(correct, clear, sep = " & "),
                            levels = c("incorrect & clear", "incorrect & unclear",
                                       "correct & unclear", "correct & clear")),
         model_summ = factor(case_when(
           model %in% c("wBM", "sBM") ~ "BM",
           model %in% c("wtrend", "strend") ~ "trend",
           model %in% c("wOUc", "sOUc") ~ "OUc",
           model %in% c("wOUs", "sOUs") ~ "OUs",
           model %in% c("wAC", "sAC", "wDC", "sDC") ~ "ACDC"
         ), levels = c("BM", "trend", "OUc", "OUs", "ACDC")))

model_fits_df_summ4 %>%
  summarise(t(table(correct)/n()))

gg2k <- ggplot(model_fits_df_summ4 %>% filter(mu == 0.5)) +
  geom_bar(data = . %>% filter(beta == "root-biased"),
           aes(x = as.numeric(fossil_prop) - .25,
               fill = cor_clear, color = "blue"),
           position = "stack", width = .2) +
  geom_bar(data = . %>% filter(beta == "random"),
           aes(x = as.numeric(fossil_prop),
               fill = cor_clear, color = "green"),
           position = "stack", width = .2) +
  geom_bar(data = . %>% filter(beta == "recent-biased"),
           aes(x = as.numeric(fossil_prop) + .25,
               fill = cor_clear, color = "red"),
           position = "stack", width = .2) +
  scale_x_continuous("Proportion of Fossils in Tree", breaks = 1:5, labels = levels(model_fits_df_summ3$fossil_prop)) +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Fit Status", palette = "Dark2") +
  scale_color_identity("Fossil Distribution", labels = c("root-biased", "random", "recent-biased"),
                       guide = guide_legend(direction = "horizontal",
                                            title.position = "top",
                                            label.position = "bottom",
                                            label.hjust = 1, label.vjust = .5,
                                            label.theme = element_text(angle = 90, size = 18),
                                            keyheight = grid::unit(7, "lines"),
                                            override.aes = list(linewidth = 1))) +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(model), rows = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips"))) +
  theme(legend.spacing.x = unit(.2, 'lines'))

gg2l <- ggplot(model_fits_df_summ4 %>% filter(mu == 0.9)) +
  geom_bar(data = . %>% filter(beta == "root-biased"),
           aes(x = as.numeric(fossil_prop) - .25,
               fill = cor_clear, color = "blue"),
           position = "stack", width = .2) +
  geom_bar(data = . %>% filter(beta == "random"),
           aes(x = as.numeric(fossil_prop),
               fill = cor_clear, color = "green"),
           position = "stack", width = .2) +
  geom_bar(data = . %>% filter(beta == "recent-biased"),
           aes(x = as.numeric(fossil_prop) + .25,
               fill = cor_clear, color = "red"),
           position = "stack", width = .2) +
  scale_x_continuous("Proportion of Fossils in Tree", breaks = 1:5, labels = levels(model_fits_df_summ3$fossil_prop)) +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Fit Status", palette = "Dark2") +
  scale_color_identity("Fossil Distribution", labels = c("root-biased", "random", "recent-biased"),
                       guide = guide_legend(direction = "horizontal",
                                            title.position = "top",
                                            label.position = "bottom",
                                            label.hjust = 1, label.vjust = .5,
                                            label.theme = element_text(angle = 90, size = 18),
                                            keyheight = grid::unit(7, "lines"),
                                            override.aes = list(linewidth = 1))) +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(model), rows = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips"))) +
  theme(legend.spacing.x = unit(.2, 'lines'))

gg2_f <- ggarrange2(gg2k, gg2l, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Prop_Correct_Clear.pdf", gg2_f, width = 40, height = 40)
ggsave("./figures/Prop_Correct_Clear_50.pdf", gg2k, width = 40, height = 20)
ggsave("./figures/Prop_Correct_Clear_90.pdf", gg2l, width = 40, height = 20)

# and the same but combined across simulated models
gg2m <- ggplot(model_fits_df_summ4 %>% filter(mu == 0.5)) +
  geom_bar(aes(x = fossil_prop, y = after_stat(count), fill = cor_clear),
           position = "fill") +
  geom_bar(data = model_fits_df_summ4 %>% filter(mu == 0.5, beta == "random", fossil_prop == 0),
           aes(x = fossil_prop, y = after_stat(count), fill = cor_clear), alpha = 0.5,
           position = "fill", layout = "fixed_cols") +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Proportion of Simulations", limits = c(0, 1)) +
  scale_fill_brewer("Fit Status", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(n_tip), rows = vars(beta),
             labeller = labeller(n_tip = function(x) paste(x, "tips"),
                                 beta = c("root" = "root-biased",
                                          "random" = "random",
                                          "recent" = "recent-biased")))

gg2n <- ggplot(model_fits_df_summ4 %>% filter(mu == 0.9, !(beta != "random" & fossil_prop == 0))) +
  geom_bar(aes(x = fossil_prop, y = after_stat(count), fill = cor_clear),
           position = "fill") +
  geom_bar(data = model_fits_df_summ4 %>% filter(mu == 0.9, beta == "random", fossil_prop == 0),
           aes(x = fossil_prop, y = after_stat(count), fill = cor_clear), alpha = 0.5,
           position = "fill", layout = "fixed_cols") +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Proportion of Simulations", limits = c(0, 1)) +
  scale_fill_brewer("Fit Status", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(n_tip), rows = vars(beta),
             labeller = labeller(n_tip = function(x) paste(x, "tips"),
                                 beta = c("root" = "root-biased",
                                          "random" = "random",
                                          "recent" = "recent-biased")))

gg2_g <- ggarrange2(gg2m, gg2n, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Prop_Correct_Clear_Combined.pdf", gg2_g, width = 18.53, height = 20)
ggsave("./figures/Prop_Correct_Clear_Combined_50.pdf", gg2m, width = 18.53, height = 10)
ggsave("./figures/Prop_Correct_Clear_Combined_90.pdf", gg2n, width = 18.53, height = 10)

# same but split out by simulated model, not phylogeny size
gg2o <- ggplot(model_fits_df_summ4 %>% filter(mu == 0.5)) +
  geom_bar(aes(x = fossil_prop, y = after_stat(count), fill = cor_clear),
           position = "fill") +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Proportion of Simulations", limits = c(0, 1)) +
  scale_fill_brewer("Fit Status", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(beta), rows = vars(model),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))

gg2p <- ggplot(model_fits_df_summ4 %>% filter(mu == 0.9)) +
  geom_bar(aes(x = fossil_prop, y = after_stat(count), fill = cor_clear),
           position = "fill") +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Proportion of Simulations", limits = c(0, 1)) +
  scale_fill_brewer("Fit Status", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(beta), rows = vars(model),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))

gg2_h <- ggarrange2(gg2o, gg2p, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Prop_Correct_Clear_Combined2.pdf", gg2_h, width = 10, height = 40)
ggsave("./figures/Prop_Correct_Clear_Combined2_50.pdf", gg2o, width = 10, height = 20)
ggsave("./figures/Prop_Correct_Clear_Combined2_90.pdf", gg2p, width = 10, height = 20)

# same but collapsed by model groupings
model_fits_df_summ4_summ <- model_fits_df_summ4 %>%
  group_by(model_summ, fossil_prop, lambda, mu, beta, cor_clear) %>%
  count() %>%
  ungroup() %>%
  group_by(model_summ, fossil_prop, lambda, mu, beta) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

gg2q <- ggplot(model_fits_df_summ4_summ %>% filter(mu == 0.5)) +
  geom_col(aes(x = fossil_prop, y = prop, fill = cor_clear)) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Proportion of Simulations", limits = c(0, 1)) +
  scale_fill_brewer("Fit Status", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(model_summ), rows = vars(beta))

gg2r <- ggplot(model_fits_df_summ4_summ %>% filter(mu == 0.9)) +
  geom_col(aes(x = fossil_prop, y = prop, fill = cor_clear)) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Proportion of Simulations", limits = c(0, 1)) +
  scale_fill_brewer("Fit Status", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(model_summ), rows = vars(beta))

gg2_i <- ggarrange2(gg2q, gg2r, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Prop_Correct_Clear_Combined3.pdf", gg2_i, width = 18.53, height = 20)
ggsave("./figures/Prop_Correct_Clear_Combined3_50.pdf", gg2q, width = 18.53, height = 10)
ggsave("./figures/Prop_Correct_Clear_Combined3_90.pdf", gg2r, width = 18.53, height = 10)

### heatmap --------------------------------------------------------
model_fits_df_summ5 <- model_fits_df_summ4 %>%
  mutate(fossil_treatment = factor(case_when(
    fossil_prop == 0 ~ "none",
    .default = beta
  ), levels = c("none", "recent", "random", "root"))) %>%
  group_by(model, lambda, mu, fossil_treatment, correct) %>%
  count(best_model) %>%
  ungroup(fossil_treatment, correct) %>%
  complete(fossil_treatment, best_model, fill = list(n = 0)) %>%
  group_by(model, lambda, mu, fossil_treatment) %>%
  mutate(prop = if(sum(n) == 0) NA else n / sum(n)) %>%
  mutate(perc = round(prop, 2) * 100) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(correct = ifelse(!is.na(prop) &&
                            grepl(substr(model,
                                         nchar(as.character(model)) - 1,
                                         nchar(as.character(model))),
                                  best_model), "correct", "incorrect")) %>%
  ungroup()

viridis_palette <- viridis::viridis_pal()(101)
gg2s <- ggplot(model_fits_df_summ5 %>% filter(mu == 0.5)) +
  geom_tile(aes(x = model, y = best_model, fill = prop), color = NA) +
  geom_tile(data = model_fits_df_summ5 %>% filter(mu == 0.5, correct == "correct"),
            aes(x = model, y = best_model, color = correct), fill = NA, linewidth = 1) +
  scale_x_discrete("Simulated Model", expand = expansion()) +
  scale_y_discrete("Best Fit Model", expand = expansion()) +
  coord_cartesian(clip = "off") +
  scale_fill_viridis_c("Proportion of Simulations", limits = c(0, 1),
                       guide = guide_colorbar(
                         theme = theme(legend.title = element_text(vjust = 1),
                                       legend.key.width = unit(15, "lines")),
                         order = 1
                       )) +
  scale_color_manual(NULL, values = c("correct" = "red"),
                     labels = c("correct" = "Correct Model"),
                     guide = guide_legend(
                       theme = theme(legend.text = element_text(size = 20, vjust = 1),
                                     legend.key.justification = "top")
                     )) +
  theme_classic(base_size = 20, ink = "black") +
  facet_grid(rows = vars(fossil_treatment),
             labeller = labeller(fossil_treatment = c("none" = "no fossils",
                                                      "root" = "root-biased\nfossils",
                                                      "random" = "random\nfossils",
                                                      "recent" = "recent-biased\nfossils"))) +
  theme(legend.position = "bottom",,
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

gg2t <- ggplot(model_fits_df_summ5 %>% filter(mu == 0.9)) +
  geom_tile(aes(x = model, y = best_model, fill = prop), color = NA) +
  geom_tile(data = model_fits_df_summ5 %>% filter(mu == 0.9, correct == "correct"),
            aes(x = model, y = best_model, color = correct), fill = NA, linewidth = 1) +
  geom_text(aes(x = model, y = best_model, label = perc / 100),
            color = deeptime:::white_or_black(viridis_palette[
              model_fits_df_summ5 %>% filter(mu == 0.9) %>% pull(perc) + 1]),
            size = 4) +
  scale_x_discrete("Simulated Model", expand = expansion()) +
  scale_y_discrete("Best Fit Model", expand = expansion()) +
  coord_cartesian(clip = "off") +
  scale_fill_viridis_c("Proportion of Simulations", limits = c(0, 1),
                       guide = guide_colorbar(
                         theme = theme(legend.title = element_text(vjust = 1),
                                       legend.key.width = unit(15, "lines")),
                         order = 1
                       )) +
  scale_color_manual(NULL, values = c("correct" = "red"),
                     labels = c("correct" = "Correct Model"),
                     guide = guide_legend(
                       theme = theme(legend.text = element_text(size = 20, vjust = 1),
                                     legend.key.justification = "top")
                     )) +
  theme_classic(base_size = 20, ink = "black") +
  facet_grid(rows = vars(fossil_treatment),
             labeller = labeller(fossil_treatment = c("none" = "no fossils",
                                                      "root" = "root-biased\nfossils",
                                                      "random" = "random\nfossils",
                                                      "recent" = "recent-biased\nfossils"))) +
  theme(legend.position = "bottom",,
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("./figures/Best_Fit_Heatmap_50.pdf", gg2s, width = 10, height = 12)
ggsave("./figures/Best_Fit_Heatmap_90.pdf", gg2t, width = 10, height = 12)

## sigma plot ------------------------------------------------------
correct_sigmas <- data.frame(model = factor(c("wBM", "sBM"), levels = c("wBM", "sBM")),
                             sigma = c(0.1, 0.5))

gg3a <- ggplot(param_estimates_df_clean %>% filter(model %in% c("wBM", "sBM"), mu == 0.5)) +
  geom_hline(data = correct_sigmas, aes(yintercept = sigma), linewidth = 1.25) +
  geom_violin(aes(x = beta, y = sigma, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Sigma") +
  scale_x_discrete("Fossil Sampling Bias", labels = c("root-\nbiased", "random", "recent-\nbiased")) +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_grid(rows = vars(model), cols = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips"))) +
  coord_cartesian(ylim = c(0, 1))
gg3b <- ggplot(param_estimates_df_clean %>% filter(model %in% c("wBM", "sBM"), mu == 0.9)) +
  geom_hline(data = correct_sigmas, aes(yintercept = sigma), linewidth = 1.25) +
  geom_violin(aes(x = beta, y = sigma, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Sigma") +
  scale_x_discrete("Fossil Sampling Bias", labels = c("root-\nbiased", "random", "recent-\nbiased")) +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_grid(rows = vars(model), cols = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips"))) +
  coord_cartesian(ylim = c(0, 1))
gg3 <- ggarrange2(gg3a, gg3b, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Sigmas.pdf", gg3, width = 16, height = 14)
ggsave("./figures/Sigmas_50.pdf", gg3a, width = 16, height = 7)
ggsave("./figures/Sigmas_90.pdf", gg3b, width = 16, height = 7)

## theta plot ------------------------------------------------------
boxplot.stats(theta_estimates_df_long$theta_val)
theta_estimates_df_long_clean <- theta_estimates_df_long %>%
  group_by(mu, model, theta) %>%
  filter(!theta_val %in% boxplot.stats(theta_val)$out) %>% # remove outliers
  mutate(theta_text = ifelse((theta == 0 & model %in% c("wOUc", "sOUc")) |
                               (theta == 1 & model %in% c("wOUs", "sOUs")), "θ", "root")) %>%
  mutate(model_theta = factor(paste(model, theta, sep = "_"),
                              levels = c("wOUc_0", "sOUc_0", "wOUs_0", "wOUs_1", "sOUs_0", "sOUs_1"))) %>%
  filter(!is.na(model_theta)) %>%
  ungroup()

theta_range_stats <- theta_estimates_df_long_clean %>%
  group_by(beta, fossil_prop, model, theta, n_tip, mu) %>%
  summarise(range = list(setNames(range(theta_val), c("range_min", "range_max"))),
            n = n(),
            .groups = "drop") %>%
  unnest_wider(range) %>%
  complete(beta, fossil_prop, model, theta, n_tip, mu, fill = list(n = 0)) %>%
  group_by(mu, model, theta) %>%
  mutate(num_y = max(range_max, na.rm = TRUE) + (max(range_max, na.rm = TRUE) - min(range_min, na.rm = TRUE)) * .075) %>%
  ungroup() %>%
  mutate(num_x = c("-1" = 1, "0" = 2, "1" = 3)[beta] +
           c("0" = -.375, "0.1" = -.1875, "0.25" = 0, "0.5" = .1875, "0.95" = .375)[fossil_prop]) %>%
  mutate(model_theta = factor(paste(model, theta, sep = "_"),
                              levels = c("wOUc_0", "sOUc_0", "wOUs_0", "wOUs_1", "sOUs_0", "sOUs_1"))) %>%
  filter(!is.na(model_theta)) %>%
  mutate(theta_text = ifelse((theta == 0 & model %in% c("wOUc", "sOUc")) |
                               (theta == 1 & model %in% c("wOUs", "sOUs")), "θ", "root"))

correct_thetas <- data.frame(model = factor(c("wOUc", "sOUc", "wOUs", "wOUs", "sOUs", "sOUs"),
                                            levels = c("wOUc", "sOUc", "wOUs", "sOUs")),
                             theta_text = c("θ", "θ", "root", "θ", "root", "θ"),
                             theta_val = c(0, 0, 0, 2, 0, 2))
gg4a <- ggplot(theta_estimates_df_long_clean %>% filter(mu == 0.5)) +
  geom_hline(data = correct_thetas, aes(yintercept = theta_val), linewidth = 1.25, color = "grey30") +
  #geom_violin(aes(x = beta, y = theta_val, fill = fossil_prop), color = "yellow",
  #            scale = "width", drop = FALSE, position = position_dodge(preserve = "single")) +
  geom_violin(data = ~ .x %>% mutate(x_num = c("-1" = 1, "0" = 2, "1" = 3)[beta] +
                                       c("0" = -.375, "0.1" = -.1875, "0.25" = 0, "0.5" = .1875, "0.95" = .375)[fossil_prop]),
              aes(x = x_num, y = theta_val, fill = fossil_prop, group = as.factor(x_num), color = fossil_prop),
              scale = "width", drop = FALSE, position = "identity") +
  geom_text(data = theta_range_stats %>% filter(mu == 0.5, n > 0),
            aes(x = num_x, y = num_y, label = n, color = fossil_prop), size = 1.9) +
  scale_y_continuous("Estimated Theta", expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete("Fossil Sampling Bias", limits = factor(1:3), labels = c("root-\nbiased", "random", "recent-\nbiased")) +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_nested(model + theta_text ~ n_tip, scales = "free_y",
               labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg4b <- ggplot(theta_estimates_df_long_clean %>% filter(mu == 0.9)) +
  geom_hline(data = correct_thetas, aes(yintercept = theta_val), linewidth = 1.25, color = "grey30") +
  #geom_violin(aes(x = beta, y = theta_val, fill = fossil_prop), color = "yellow",
  #            scale = "width", drop = FALSE, position = position_dodge(preserve = "single")) +
  geom_violin(data = ~ .x %>% mutate(x_num = c("-1" = 1, "0" = 2, "1" = 3)[beta] +
                                       c("0" = -.375, "0.1" = -.1875, "0.25" = 0, "0.5" = .1875, "0.95" = .375)[fossil_prop]),
              aes(x = x_num, y = theta_val, fill = fossil_prop, group = as.factor(x_num), color = fossil_prop),
              scale = "width", drop = FALSE, position = "identity") +
  geom_text(data = theta_range_stats %>% filter(mu == 0.9, n > 0),
            aes(x = num_x, y = num_y, label = n, color = fossil_prop), size = 1.9) +
  scale_y_continuous("Estimated Theta", expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete("Fossil Sampling Bias", limits = factor(1:3), labels = c("root-\nbiased", "random", "recent-\nbiased")) +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_nested(model + theta_text ~ n_tip, scales = "free_y",
               labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg4 <- ggarrange2(gg4a, gg4b, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Thetas.pdf", gg4, width = 16, height = 22, device = cairo_pdf)
ggsave("./figures/Thetas_50.pdf", gg4a, width = 16, height = 11, device = cairo_pdf)
ggsave("./figures/Thetas_90.pdf", gg4b, width = 16, height = 11, device = cairo_pdf)

## half-life plot ----------------------------------------------
# half-lives relative to max tree height
correct_rhls <- data.frame(model = factor(c("wOUc", "sOUc", "wOUs", "sOUs"),
                                          levels = c("wOUc", "sOUc", "wOUs", "sOUs")),
                           rel_hl = c(1, 1/5, 1, 1/5))

gg5a <- ggplot(param_estimates_df_clean %>%
                 filter(model %in% c("wOUc", "wOUs", "sOUc", "sOUs"), mu == 0.5)) +
  geom_hline(data = correct_rhls, aes(yintercept = rel_hl)) +
  geom_violin(aes(x = beta, y = rel_hl, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Relative Phylogenetic Half-life") +
  scale_x_discrete("Fossil Sampling Bias", labels = c("root-\nbiased", "random", "recent-\nbiased")) +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_grid(rows = vars(model), cols = vars(n_tip), scales = "free_y",
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg5b <- ggplot(param_estimates_df_clean %>%
                 filter(model %in% c("wOUc", "sOUc", "wOUs", "sOUs"), mu == 0.9)) +
  geom_hline(data = correct_rhls, aes(yintercept = rel_hl)) +
  geom_violin(aes(x = beta, y = rel_hl, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Relative Phylogenetic Half-life") +
  scale_x_discrete("Fossil Sampling Bias", labels = c("root-\nbiased", "random", "recent-\nbiased")) +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_grid(rows = vars(model), cols = vars(n_tip), scales = "free_y",
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg5 <- ggarrange2(gg5a, gg5b, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Halflives.pdf", gg5, width = 22, height = 22)
ggsave("./figures/Halflives_50.pdf", gg5a, width = 22, height = 11)
ggsave("./figures/Halflives_90.pdf", gg5b, width = 22, height = 11)

## trend plot ----------------------------------------------
correct_trends <- data.frame(model = factor(c("wtrend", "strend"), levels = c("wtrend", "strend")),
                             trend = c(0.1, 0.3))

gg6a <- ggplot(param_estimates_df_clean %>% filter(model %in% c("wtrend", "strend"), mu == 0.5)) +
  geom_hline(data = correct_trends, aes(yintercept = trend)) +
  geom_violin(aes(x = beta, y = trend, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Trend") +
  scale_x_discrete("Fossil Sampling Bias", labels = c("root-\nbiased", "random", "recent-\nbiased")) +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2", breaks = c(0.1, 0.25, 0.5, 0.95), limits = factor(c(0, 0.1, 0.25, 0.5, 0.95))) +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2", breaks = c(0.1, 0.25, 0.5, 0.95), limits = factor(c(0, 0.1, 0.25, 0.5, 0.95))) +
  theme_bw(base_size = 20) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_grid(rows = vars(model), cols = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips"))) +
  coord_cartesian(ylim = c(-0.1, 0.5))
gg6b <- ggplot(param_estimates_df_clean %>% filter(model %in% c("wtrend", "strend"), mu == 0.9)) +
  geom_hline(data = correct_trends, aes(yintercept = trend)) +
  geom_violin(aes(x = beta, y = trend, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Trend") +
  scale_x_discrete("Fossil Sampling Bias", labels = c("root-\nbiased", "random", "recent-\nbiased")) +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2", breaks = c(0.1, 0.25, 0.5, 0.95), limits = factor(c(0, 0.1, 0.25, 0.5, 0.95))) +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2", breaks = c(0.1, 0.25, 0.5, 0.95), limits = factor(c(0, 0.1, 0.25, 0.5, 0.95))) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_grid(rows = vars(model), cols = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips"))) +
  coord_cartesian(ylim = c(-0.1, 0.5))
gg6 <- ggarrange2(gg6a, gg6b, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Trends.pdf", gg6, width = 16, height = 14)
ggsave("./figures/Trends_50.pdf", gg6a, width = 16, height = 7)
ggsave("./figures/Trends_90.pdf", gg6b, width = 16, height = 7)

## beta plot -----------------------------------------------
correct_betas <- data.frame(model = factor(c("wAC", "sAC", "wDC", "sDC"), levels = c("wAC", "sAC", "wDC", "sDC")),
                            exp_rate = c(0.1, 0.3, -0.1, -0.3))

gg7a <- ggplot(param_estimates_df_clean %>% filter(model %in% c("wAC", "sAC", "wDC", "sDC"), mu == 0.5)) +
  geom_hline(data = correct_betas, aes(yintercept = exp_rate)) +
  geom_violin(aes(x = beta, y = exp_rate, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Beta") +
  scale_x_discrete("Fossil Sampling Bias", labels = c("root-\nbiased", "random", "recent-\nbiased")) +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_grid(rows = vars(model), cols = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg7b <- ggplot(param_estimates_df_clean %>% filter(model %in% c("wAC", "sAC", "wDC", "sDC"), mu == 0.9)) +
  geom_hline(data = correct_betas, aes(yintercept = exp_rate)) +
  geom_violin(aes(x = beta, y = exp_rate, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Beta") +
  scale_x_discrete("Fossil Sampling Bias", labels = c("root-\nbiased", "random", "recent-\nbiased")) +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_grid(rows = vars(model), cols = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg7 <- ggarrange2(gg7a, gg7b, nrow = 2, draw = FALSE, labels = c("mu = 0.5", "mu = 0.9"))
ggsave("./figures/Betas.pdf", gg7, width = 16, height = 14)
ggsave("./figures/Betas_50.pdf", gg7a, width = 16, height = 7)
ggsave("./figures/Betas_90.pdf", gg7b, width = 16, height = 7)

## heights of fossils --------------------------------------
tree_df_scaled <- readRDS("./data/tree_simulations_scaled.RDS")
fossil_heights <- lapply(seq_along(tree_df_scaled$tree), FUN = function(i) {
  tree <- tree_df_scaled$tree[[i]]
  tree$edge.length <- tree$edge.length / max(nodeHeights(tree)[, 2])
  n_heights <- nodeHeights(tree)[which(tree$edge <= Ntip(tree))]
  tmp <- n_heights[n_heights < (1 - 0.000001)] # because computer math
  if(length(tmp > 0)) {
    cbind.data.frame(height = tmp, tree_df_scaled[i, 1:ncol(tree_df_scaled) - 1])
  }
}) %>%
  do.call(rbind, .) %>%
  mutate(beta = fct_recode(as.factor(beta), `root-biased` = "root", `random` = "random", `recent-biased` = "recent"))

gg5 <- ggplot(fossil_heights) +
  geom_histogram(aes(x = height, fill = factor(beta)), position = "dodge",
                 binwidth = 0.05, boundary = 1) +
  scale_y_continuous("# of Fossils") +
  scale_x_continuous("Relative Height in Phylogeny") +
  scale_fill_brewer(NULL, palette = "Set1") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(mu), rows = vars(n_tip), scales = "free_y",
             labeller = labeller(n_tip = function(x) paste(x, "tips"),
                                 mu = function(x) paste("\u03BC = ", x))) +
  theme(legend.position = "top", panel.spacing.x = unit(1.5, "lines"))
ggsave("./figures/Fossil_heights.pdf", gg5, width = 12, height = 8)

fossil_heights_bins <- fossil_heights %>%
  mutate(height_bin = cut(height, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)) %>%
  group_by(mu, n_tip, fossil_prop, beta, height_bin) %>%
  count() %>%
  ungroup() %>%
  mutate(perc_foss = n / (n_tip * fossil_prop * 100) * 100,
         height_bin_cont = as.numeric(height_bin) / 20 - 0.025)

gg6 <- ggplot(fossil_heights_bins %>% filter(mu == 0.5)) +
  geom_col(aes(x = height_bin_cont, y = perc_foss, fill = factor(beta)),
           position = "dodge") +
  scale_x_continuous("Relative Height in Phylogeny") +
  scale_y_continuous("% of Fossils in Phylogenies with Treatment") +
  scale_fill_brewer(NULL, palette = "Set1") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(fossil_prop), rows = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips"),
                                 fossil_prop = function(x) {
                                   paste0(as.numeric(x) * 100, "% fossils")
                                 })) +
  theme(legend.position = "top", panel.spacing.x = unit(1.5, "lines"))
gg7 <- ggplot(fossil_heights_bins %>% filter(mu == 0.90)) +
  geom_col(aes(x = height_bin_cont, y = perc_foss, fill = factor(beta)),
           position = "dodge") +
  scale_x_continuous("Relative Height in Phylogeny") +
  scale_y_continuous("% of Fossils in Phylogenies with Treatment") +
  scale_fill_brewer(NULL, palette = "Set1") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(fossil_prop), rows = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips"),
                                 fossil_prop = function(x) {
                                   paste0(as.numeric(x) * 100, "% fossils")
                                 })) +
  theme(legend.position = "top", panel.spacing.x = unit(1.5, "lines"))
ggsave("./figures/Fossil_heights_50.pdf", gg6, width = 12, height = 12)
ggsave("./figures/Fossil_heights_90.pdf", gg7, width = 12, height = 12)

# plots for schematic ####
# TODO: finish schematic components and combine them (and make it pretty)
## intial full trees ####
par(mfrow = c(5, 1), mar = c(0, 0, 0, 0))
set.seed(1234)
ex_trees <- sapply(n_tips, FUN = function(n) {
  sim.bd.taxa(n, 1, 1, 0.9, frac = 1, complete = TRUE)[[1]]
}, USE.NAMES = TRUE, simplify = FALSE)

lapply(ex_trees, FUN = function(tree) {
  plot.phylo(tree, show.tip.label = FALSE, no.margin = TRUE)
})

tree_50 <- ladderize(ex_trees[[1]], right = FALSE)
tree_200 <- ladderize(ex_trees[[3]], right = FALSE)

layout(1)
plot.phylo(tree_50, show.tip.label = FALSE, no.margin = TRUE)
dev.print(svg, filename = "./figures/Example_Full_Tree_50tips.svg", width = 8, height = 6)
dev.print(png, filename = "./figures/Example_Full_Tree_50tips.png", width = 8, height = 6,
          units = "in", res = 300)

## generate some example fossils ####
# an example for just showing fossils (otherwise there are too many)
set.seed(1234)
ex_foss_low <- sim.fossils.poisson(1, tree_50, root.edge = FALSE)

layout(1)
# do this instead of using FossilSim so we can remove the plot margins
xx <- node.depth.edgelength(tree_50)
yy <- node.height(tree_50)
ex_foss_low$x <- max(xx) - ex_foss_low$hmax
ex_foss_low$y <- yy[ex_foss_low$edge]

plot.phylo(tree_50, show.tip.label = FALSE, no.margin = TRUE)
points(ex_foss_low$x, ex_foss_low$y, col = "red", pch = 18, cex = 0.75)
dev.print(svg, filename = "./figures/Example_Fossils_LowRate_50tips.svg", width = 8, height = 6)
dev.print(png, filename = "./figures/Example_Fossils_LowRate_50tips.png", width = 8, height = 6,
          units = "in", res = 300)

# now examples with the real rate to generate recovery potentials
set.seed(1234)
ex_foss_high <- sim.fossils.poisson(50, tree_50, root.edge = FALSE)

# calculate recovery potentials under different models
foss_sub <- subset(ex_foss_high, edge <= Ntip(tree_50))
ex_recov <- lapply(models, function(model) {
  max_age <- max(FossilSim:::n.ages(tree_50))
  foss_sub$rel_age <- (max_age - foss_sub$hmin) / max_age
  foss_sub$recovery_potential <- model(foss_sub$rel_age)
  # for each branch, sum the recovery potentials for all occurrences on that branch
  branch_recovery <- foss_sub %>%
    group_by(edge) %>%
    summarise(recovery_potential = sum(recovery_potential, na.rm = TRUE)) %>%
    ungroup()
})

par(mfrow = c(3, 1), mar = c(0, 0, 0, 0))
for (i in 1:3) {
  cols <- rep("black", nrow(tree_50$edge))
  cols[match(ex_recov[[i]]$edge, tree_50$edge[, 2])] <- viridisLite::magma(1000)[cut(log10(ex_recov[[i]]$recovery_potential), 1000)]
  edge_wdt <- rep(.5, nrow(tree_50$edge))
  edge_wdt[match(ex_recov[[i]]$edge, tree_50$edge[, 2])] <- 3
  plot.phylo(tree_50, edge.color = cols, show.tip.label = FALSE, edge.width = edge_wdt)
  if (i == 1) {
    fields::image.plot(legend.only = TRUE, zlim = c(0,1), col = viridisLite::magma(1000),
                       horizontal = TRUE, smallplot = c(0.05, 0.4, 0.7, 0.8),
                       legend.lab = "Relative Recovery Potential", legend.line = -4)
  }
}
dev.print(svg, filename = "./figures/Example_Recovery_Potentials_50tips.svg", width = 6, height = 8)
dev.print(png, filename = "./figures/Example_Recovery_Potentials_50tips.png", width = 6, height = 8,
          units = "in", res = 300)

## resulting trees ####
tree_sub <- tree_df %>%
  filter(n_tip == 50, mu == 0.9, sim == 1, fossil_prop == 0.5) %>%
  mutate(beta = fct_recode(as.factor(beta), `root-biased` = "root", `random` = "random", `recent-biased` = "recent"))

par(mfrow = c(1, 3), mar = c(7.5, 0, 0, 0))
for (i in seq_len(nrow(tree_sub))) {
  plot.phylo(ladderize(tree_sub$tree[[i]], right = FALSE), show.tip.label = FALSE)
  palaeoverse::axis_geo_phylo(side = 1, lab_size = 1.5, cex.axis = 1.5,
                              lwd = 1.5, height = 0.1,
                              skip = c("Oligocene", "Holocene"),
                              title = tree_sub$beta[i])
}
dev.print(svg, filename = "./figures/Example_Resulting_Trees_50tips.svg", width = 9, height = 4.5)
dev.print(png, filename = "./figures/Example_Resulting_Trees_50tips.png", width = 9, height = 4.5,
          units = "in", res = 300)

## example trait simulations ####
ind <- which(with(tree_df, n_tip == 200 & fossil_prop == 0.5 & mu == 0.9 & sim == 90 & beta == "random"))

par(mfrow = c(2, 2), mar = c(2, 0, 2, 2))
phenogram(tree_df$tree[[ind]], sBM_trait[[ind]][,1], fsize = FALSE, spread.labels = FALSE)
phenogram(tree_df$tree[[ind]], strend_trait[[ind]][,1], fsize = FALSE, spread.labels = FALSE)
phenogram(tree_df$tree[[ind]], sDC_trait[[ind]][,1], fsize = FALSE, spread.labels = FALSE)
phenogram(tree_df$tree[[ind]], sOUs_trait[[ind]][,1], fsize = FALSE, spread.labels = FALSE)

## testing theta values for OU model ####
max_height <- max(nodeHeights(tree_df$tree[[ind]]))
theta_vals <- c(1, 2, 3, 4, 5)
halflives <- max_height / c(10, 5, 2, 1)
ou_vals <- expand.grid(theta = theta_vals, halflife = halflives)
ou_vals$alpha <- log(2) / ou_vals$halflife

set.seed(1234)
test_ou <- mapply(function(alpha, theta) mvSIM(tree_df$tree[[ind]], nsim = 2, model = "OU1",
                        param = list(root = TRUE,
                                     alpha = alpha, #strength of selection
                                     theta = c(0, theta), #ancestral state, optimum
                                     sigma = 0.1 #strength of drift
                        )),
                  alpha = ou_vals$alpha, theta = ou_vals$theta, SIMPLIFY = FALSE)

par(mfrow = c(4, 5), mar = c(2, 4, 2, 0))
#phenogram(tree_df$tree[[ind]], sBM_trait[[ind]][,1], fsize = FALSE, spread.labels = FALSE)
#title("BM")
for (i in seq_along(test_ou)) {
  phenogram(tree_df$tree[[ind]], test_ou[[i]][,1], fsize = FALSE, spread.labels = FALSE, color = "blue")
  for (j in 2:ncol(test_ou[[i]])) {
    phenogram(tree_df$tree[[ind]], test_ou[[i]][, j], fsize = FALSE, spread.labels = FALSE, add = TRUE, color = "red")
  }
  title(paste("theta = ", ou_vals$theta[i], ", halflife = ", round(ou_vals$halflife[i], 1), sep = ""))
}

test_ou_ape <- mapply(function(alpha, theta) rTraitCont(tree_df$tree[[ind]], model = "OU",
                                                        alpha = alpha, theta = theta, root.value = 0),
                  alpha = ou_vals$alpha, theta = ou_vals$theta, SIMPLIFY = FALSE)
par(mfrow = c(4, 5), mar = c(2, 4, 2, 0))
for (i in seq_along(test_ou_ape)) {
  phenogram(tree_df$tree[[ind]], test_ou_ape[[i]], fsize = FALSE, spread.labels = FALSE)
  title(paste("theta = ", ou_vals$theta[i], ", halflife = ", round(ou_vals$halflife[i], 3), sep = ""))
}
