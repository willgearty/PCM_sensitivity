# Written by Will Gearty, Bethany Allen, and Pedro Godoy

# Load libraries
#devtools::install_github("willgearty/pcmtools")
library(ape); library(phytools); library(geiger); library(TreeSim)
library(FossilSim); library(mvMORPH); library(pbapply); library(dplyr)
library(tibble); library(tidyr); library(ggplot2); library(pcmtools)
library(deeptime); library(future); library(forcats)

# Load functions
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
mus <- c(0.25, 0.9)
betas <- c(-3, 0, 3)
n_sim <- 100

# get all unique combinations of parameters
settings <- expand.grid(n_tip = n_tips, fossil_prop = fossil_props,
                        lambda = lambdas, mu = mus, beta = betas)

# Simulate trees -----------------------------------------------------------
set.seed(1234)
tree_df <- pbmapply(function(n_tip, fossil_prop, lambda, mu, beta) {
  trees <- sim.fbd.taxa.prop(n_tip, fossil_prop, numbsim = n_sim,
                             lambda = lambda, mu = mu,
                             model = "EB", a = beta, progress = FALSE)
  tmp <- data.frame(n_tip = rep(n_tip, n_sim),
                    fossil_prop = rep(fossil_prop, n_sim),
                    lambda = rep(lambda, n_sim),
                    mu = rep(mu, n_sim),
                    beta = rep(beta, n_sim),
                    sim = seq_len(n_sim))
  tmp$tree <- trees
  tmp
}, n_tip = settings$n_tip, fossil_prop = settings$fossil_prop,
lambda = settings$lambda, mu = settings$mu, beta = settings$beta,
SIMPLIFY = FALSE) %>% do.call(rbind, .)

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
  #Simulate weak OU ("SSP") trait evolution, centred on ancestral state
  mvSIM(tree = tree, nsim = 1, model = "OU1",
        param = list(alpha = log(2) / max_height, #strength of selection
                     theta = 0, #ancestral state
                     sigma = 0.1 #strength of drift
        ))
})

sOUc_trait <- pblapply(tree_df$tree, function(tree) {
  max_height <- max(nodeHeights(tree))
  #Simulate strong OU ("SSP") trait evolution, centred on ancestral state
  mvSIM(tree = tree, nsim = 1, model = "OU1",
        param = list(alpha = log(2) / (max_height / 10), #strength of selection
                     theta = 0, #ancestral state
                     sigma = 0.1 #strength of drift
        ))
})

wOUs_trait <- pblapply(tree_df$tree, function(tree) {
  max_height <- max(nodeHeights(tree))
  #Simulate weak OU ("SSP") trait evolution, with shifted optimum
  mvSIM(tree = tree, nsim = 1, model = "OU1",
        param = list(root = TRUE,
                     alpha = log(2) / max_height, #strength of selection
                     theta = c(0, 1), #ancestral state, optimum
                     sigma = 0.1 #strength of drift
        ))
})

sOUs_trait <- pblapply(tree_df$tree, function(tree) {
  max_height <- max(nodeHeights(tree))
  #Simulate strong OU ("SSP") trait evolution, with shifted optimum
  mvSIM(tree = tree, nsim = 1, model = "OU1",
        param = list(root = TRUE,
                     alpha = log(2) / (max_height / 10), #strength of selection
                     theta = c(0, 1), #ancestral state, optimum
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

plan(multisession)
# for each of the trait evolution models used in the simulations
for (mod in mods) {
  if (mod == "wBM") next
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
  }, cl = "future")
  saveRDS(model_fitting_results,
          paste0("./data/model_fitting/model_fitting_", mod, ".RDS"))
}
plan(sequential)

# Analyze results ---------------------------------------------------
# ('settings' must be created)
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

fit_models <- c("BM", "trend", "OU1", "OU2", "ACDC")
model_fits_df_long <- model_fits_df %>%
  pivot_longer(cols = all_of(fit_models), names_to = "fit_model", values_to = "aicc") %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta, sim) %>%
  mutate(aicc_w = AICweights(aicc), aicc_d = aicc - min(aicc)) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = mods),
         fit_model = factor(fit_model, levels = fit_models),
         across(c(n_tip, fossil_prop, lambda, mu, beta, sim, aicc_w), ~ as.numeric(.x))) %>%
  mutate(across(c(n_tip, fossil_prop, lambda, mu, beta), ~ as.factor(.x))) %>%
  mutate(beta = fct_recode(beta, old = "-3", even = "0", young = "3"))

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

param_estimates_df <- param_estimates_df %>%
  mutate(model = factor(model, levels = mods),
         fit_model = factor(fit_model, levels = fit_models),
         across(c(n_tip, fossil_prop, lambda, mu, beta, sim), ~ as.numeric(.x))) %>%
  mutate(across(c(n_tip, fossil_prop, lambda, mu, beta), ~ as.factor(.x))) %>%
  mutate(correct_model = case_when(
    model %in% c("wBM", "sBM") ~ "BM",
    model %in% c("wtrend", "strend") ~ "trend",
    model %in% c("wOUc", "sOUc") ~ "OU1",
    model %in% c("wOUs", "sOUs") ~ "OU2",
    model %in% c("wAC", "sAC", "wDC", "sDC") ~ "ACDC"
  )) %>%
  filter(fit_model == correct_model) %>%
  select(-correct_model) %>%
  mutate(rel_hl = (log(2) / alpha) / sapply(tree_df$tree, function(tree) max(nodeHeights(tree)))) %>%
  mutate(beta = fct_recode(beta, old = "-3", even = "0", young = "3"))

# need to do some pivoting to get the two different theta values for each simulation
theta_estimates_df_long <- param_estimates_df %>%
  pivot_longer(cols = c(theta_0, theta_1), names_to = "theta", values_to = "theta_val",
               names_prefix = "theta_") %>%
  mutate(model_theta = paste(model, theta, sep = "_")) %>%
  filter(model_theta %in% c("wOUc_0", "sOUc_0", "wOUs_0", "sOUs_0", "wOUs_1", "sOUs_1")) %>%
  group_by(model_theta) %>%
  filter(!theta_val %in% boxplot.stats(theta_val)$out) %>% # remove outliers
  ungroup()

# clean up the old object
remove(model_results)

## AIC plot ------------------------------------------------------
gg1 <- ggplot(model_fits_df_long %>% filter(mu == "0.25")) +
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
    model %in% c("wOUc", "sOUc") ~ "OU1",
    model %in% c("wOUs", "sOUs") ~ "OU2",
    model %in% c("wAC", "sAC", "wDC", "sDC") ~ "ACDC"
  )) %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta, sim) %>%
  summarise(correct = ifelse(all(is.na(aicc_w)), NA,
                             any(unique(correct_model) ==
                                   fit_model[which(aicc_w == max(aicc_w, na.rm = TRUE))])),
            .groups = "drop") %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta) %>%
  summarise(prop_true = sum(correct)/n(), .groups = "drop")

gg2a <- ggplot(model_fits_df_summ %>% filter(mu == 0.25)) +
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
gg2 <- ggarrange2(gg2a, gg2b, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Prop_Best.pdf", gg2, width = 16, height = 20)

# generating models for best fitting models
model_fits_df_summ2 <- model_fits_df_long %>%
  mutate(correct_model = case_when(
    model %in% c("wBM", "sBM") ~ "BM",
    model %in% c("wtrend", "strend") ~ "trend",
    model %in% c("wOUc", "sOUc") ~ "OU1",
    model %in% c("wOUs", "sOUs") ~ "OU2",
    model %in% c("wAC", "sAC", "wDC", "sDC") ~ "ACDC"
  )) %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta, sim) %>%
  slice_max(aicc_w) %>%
  ungroup() %>%
  group_by(n_tip, fossil_prop, lambda, mu, beta, fit_model) %>%
  count(correct_model) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

gg2c <- ggplot(model_fits_df_summ2 %>% filter(mu == 0.25)) +
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
gg2_b <- ggarrange2(gg2c, gg2d, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Prop_Sim.pdf", gg2_b, width = 16, height = 25)

# when wrong answer, what is it?
model_fits_df_summ3 <- model_fits_df_long %>%
  mutate(correct_model = case_when(
    model %in% c("wBM", "sBM") ~ "BM",
    model %in% c("wtrend", "strend") ~ "trend",
    model %in% c("wOUc", "sOUc") ~ "OU1",
    model %in% c("wOUs", "sOUs") ~ "OU2",
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
  mutate(prop = n / sum(n)) %>%
  ungroup()

gg2e <- ggplot(model_fits_df_summ3 %>% filter(mu == 0.25)) +
  geom_col(data = . %>% filter(beta == "old"),
           aes(x = as.numeric(fossil_prop) - .25, y = n, fill = fit_model, color = "blue"),
           position = "stack", width = .2) +
  geom_col(data = . %>% filter(beta == "even"),
           aes(x = as.numeric(fossil_prop), y = n, fill = fit_model, color = "green"),
           position = "stack", width = .2) +
  geom_col(data = . %>% filter(beta == "young"),
           aes(x = as.numeric(fossil_prop) + .25, y = n, fill = fit_model, color = "red"),
           position = "stack", width = .2) +
  scale_x_continuous("Proportion of Fossils in Tree", breaks = 1:5, labels = levels(model_fits_df_summ3$fossil_prop)) +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Fit Model", palette = "Dark2") +
  scale_color_identity("Fossil Distribution", guide = guide_legend(), labels = c("old", "even", "young")) +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(model), rows = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))

gg2f <- ggplot(model_fits_df_summ3 %>% filter(mu == 0.9)) +
  geom_col(data = . %>% filter(beta == "old"),
           aes(x = as.numeric(fossil_prop) - .25, y = n, fill = fit_model, color = "blue"),
           position = "stack", width = .2) +
  geom_col(data = . %>% filter(beta == "even"),
           aes(x = as.numeric(fossil_prop), y = n, fill = fit_model, color = "green"),
           position = "stack", width = .2) +
  geom_col(data = . %>% filter(beta == "young"),
           aes(x = as.numeric(fossil_prop) + .25, y = n, fill = fit_model, color = "red"),
           position = "stack", width = .2) +
  scale_x_continuous("Proportion of Fossils in Tree", breaks = 1:5, labels = levels(model_fits_df_summ3$fossil_prop)) +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Fit Model", palette = "Dark2") +
  scale_color_identity("Fossil Distribution", guide = guide_legend(), labels = c("old", "even", "young")) +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(model), rows = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))

gg2_c <- ggarrange2(gg2e, gg2f, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Wrong_Best.pdf", gg2_c, width = 40, height = 40)

# summarize proportion of correct models
gg2g <- ggplot(model_fits_df_summ %>% filter(mu == 0.25)) +
  geom_col(aes(x = fossil_prop, y = prop_true / 12, fill = model)) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations Correctly Identified",
                     limits = c(0, 1)) +
  scale_fill_brewer("Simulated Model", palette = "Paired") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(n_tip), rows = vars(beta),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))

gg2h <- ggplot(model_fits_df_summ %>% filter(mu == 0.9)) +
  geom_col(aes(x = fossil_prop, y = prop_true / 12, fill = model)) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations Correctly Identified",
                     limits = c(0, 1)) +
  scale_fill_brewer("Simulated Model", palette = "Paired") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(n_tip), rows = vars(beta),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))

gg2_d <- ggarrange2(gg2g, gg2h, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Prop_Best_Stacked.pdf", gg2_d, width = 18.53, height = 20)

gg2i <- ggplot(model_fits_df_summ %>% filter(mu == 0.25)) +
  geom_col(aes(x = fossil_prop, y = prop_true / 12)) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations Correctly Identified",
                     limits = c(0, 1)) +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(n_tip), rows = vars(beta),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))

gg2j <- ggplot(model_fits_df_summ %>% filter(mu == 0.9)) +
  geom_col(aes(x = fossil_prop, y = prop_true / 12)) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations Correctly Identified",
                     limits = c(0, 1)) +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(n_tip), rows = vars(beta),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))

gg2_e <- ggarrange2(gg2i, gg2j, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Prop_Best_Combined.pdf", gg2_e, width = 16, height = 20)

# proportions of simulations with clear best model
model_fits_df_summ4 <- model_fits_df_long %>%
  mutate(correct = case_when(
    model %in% c("wBM", "sBM") ~ "BM",
    model %in% c("wtrend", "strend") ~ "trend",
    model %in% c("wOUc", "sOUc") ~ "OU1",
    model %in% c("wOUs", "sOUs") ~ "OU2",
    model %in% c("wAC", "sAC", "wDC", "sDC") ~ "ACDC"
  ) == fit_model) %>%
  mutate(best = ifelse(is.na(aicc_d), FALSE, aicc_d == 0)) %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, beta, sim) %>%
  summarise(correct = ifelse(any(correct & best), "correct", "incorrect"),
            clear = ifelse(all(aicc_d[!is.na(aicc_d) & !best] > 2), "clear", "unclear"),
            .groups = "drop") %>%
  mutate(cor_clear = factor(interaction(correct, clear, sep = " & "),
                            levels = c("correct & clear", "correct & unclear",
                                       "incorrect & clear", "incorrect & unclear")))

gg2k <- ggplot(model_fits_df_summ4 %>% filter(mu == 0.25)) +
  geom_bar(data = . %>% filter(beta == "old"),
           aes(x = as.numeric(fossil_prop) - .25,
               fill = cor_clear, color = "blue"),
           position = "stack", width = .2) +
  geom_bar(data = . %>% filter(beta == "even"),
           aes(x = as.numeric(fossil_prop),
               fill = cor_clear, color = "green"),
           position = "stack", width = .2) +
  geom_bar(data = . %>% filter(beta == "young"),
           aes(x = as.numeric(fossil_prop) + .25,
               fill = cor_clear, color = "red"),
           position = "stack", width = .2) +
  scale_x_continuous("Proportion of Fossils in Tree", breaks = 1:5, labels = levels(model_fits_df_summ3$fossil_prop)) +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Fit Status", palette = "Dark2") +
  scale_color_identity("Fossil Distribution", labels = c("old", "even", "young"),
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
  geom_bar(data = . %>% filter(beta == "old"),
           aes(x = as.numeric(fossil_prop) - .25,
               fill = cor_clear, color = "blue"),
           position = "stack", width = .2) +
  geom_bar(data = . %>% filter(beta == "even"),
           aes(x = as.numeric(fossil_prop),
               fill = cor_clear, color = "green"),
           position = "stack", width = .2) +
  geom_bar(data = . %>% filter(beta == "young"),
           aes(x = as.numeric(fossil_prop) + .25,
               fill = cor_clear, color = "red"),
           position = "stack", width = .2) +
  scale_x_continuous("Proportion of Fossils in Tree", breaks = 1:5, labels = levels(model_fits_df_summ3$fossil_prop)) +
  scale_y_continuous("Number of Simulations", limits = c(0, 100)) +
  scale_fill_brewer("Fit Status", palette = "Dark2") +
  scale_color_identity("Fossil Distribution", labels = c("old", "even", "young"),
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

gg2_f <- ggarrange2(gg2k, gg2l, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Prop_Correct_Clear.pdf", gg2_f, width = 40, height = 40)

# and the same but combined across simulated models
gg2m <- ggplot(model_fits_df_summ4 %>% filter(mu == 0.25)) +
  geom_bar(aes(x = fossil_prop, y = after_stat(count) / 1200, fill = cor_clear)) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Proportion of Simulations", limits = c(0, 1)) +
  scale_fill_brewer("Fit Status", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(n_tip), rows = vars(beta),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))

gg2n <- ggplot(model_fits_df_summ4 %>% filter(mu == 0.9)) +
  geom_bar(aes(x = fossil_prop, y = after_stat(count) / 1200, fill = cor_clear)) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Proportion of Simulations", limits = c(0, 1)) +
  scale_fill_brewer("Fit Status", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(n_tip), rows = vars(beta),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))

gg2_g <- ggarrange2(gg2m, gg2n, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Prop_Correct_Clear_Combined.pdf", gg2_g, width = 18.53, height = 20)

## sigma plot ------------------------------------------------------
correct_sigmas <- data.frame(model = factor(c("wBM", "sBM"), levels = c("wBM", "sBM")),
                             sigma = c(0.1, 0.5))

gg3a <- ggplot(param_estimates_df %>% filter(model %in% c("wBM", "sBM"), mu == 0.25)) +
  geom_hline(data = correct_sigmas, aes(yintercept = sigma), linewidth = 1.25) +
  geom_violin(aes(x = beta, y = sigma, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Sigma") +
  scale_x_discrete("Fossil Sampling Bias") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(rows = vars(model), cols = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips"))) +
  coord_cartesian(ylim = c(0, 1))
gg3b <- ggplot(param_estimates_df %>% filter(model %in% c("wBM", "sBM"), mu == 0.9)) +
  geom_hline(data = correct_sigmas, aes(yintercept = sigma), linewidth = 1.25) +
  geom_violin(aes(x = beta, y = sigma, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Sigma") +
  scale_x_discrete("Fossil Sampling Bias") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(rows = vars(model), cols = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips"))) +
  coord_cartesian(ylim = c(0, 1))
gg3 <- ggarrange2(gg3a, gg3b, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Sigmas.pdf", gg3, width = 16, height = 12)

## theta plot ------------------------------------------------------
theta_estimates_df_long <- theta_estimates_df_long %>%
  mutate(model_theta = factor(model_theta, levels = c("wOUc_0", "sOUc_0", "wOUs_0", "sOUs_0", "wOUs_1", "sOUs_1"))) %>%
  complete(beta, fossil_prop, model_theta, n_tip, mu) %>%
  mutate(transparent = ifelse(is.na(theta_val), "transparent", "black")) %>%
  mutate(theta_val = ifelse(is.na(theta_val), 0, theta_val))
theta_boxplot_stats <- theta_estimates_df_long %>%
  group_by(beta, fossil_prop, model_theta, n_tip, mu) %>%
  summarise(boxplot = list(setNames(boxplot.stats(theta_val)$stats,
                                    c('lower_whisker','lower_hinge','median','upper_hinge','upper_whisker'))),
            .groups = "drop") %>%
  unnest_wider(boxplot) %>%
  complete(beta, fossil_prop, model_theta, n_tip, mu)

correct_thetas <- data.frame(model_theta = factor(c("wOUc_0", "sOUc_0", "wOUs_0", "sOUs_0", "wOUs_1", "sOUs_1"),
                                                  levels = c("wOUc_0", "sOUc_0", "wOUs_0", "sOUs_0", "wOUs_1", "sOUs_1")),
                             theta_val = c(rep(0, 4), 1, 1))
gg4a <- ggplot(theta_estimates_df_long %>% filter(mu == 0.25)) +
  geom_hline(data = correct_thetas, aes(yintercept = theta_val), linewidth = 1.25, color = "grey30") +
  geom_boxplot(aes(x = beta, y = theta_val, fill = fossil_prop, color = transparent),
               position = position_dodge2(padding = .2, preserve = "single")) +
  scale_y_continuous("Estimated Theta") +
  scale_x_discrete("Fossil Sampling Bias") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_identity(NULL, ) +
  theme_bw(base_size = 20) +
  facet_grid(rows = vars(model_theta), cols = vars(n_tip), scales = "free_y",
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg4b <- ggplot(theta_estimates_df_long %>% filter(mu == 0.9)) +
  geom_hline(data = correct_thetas, aes(yintercept = theta_val), linewidth = 1.25, color = "grey30") +
  geom_boxplot(aes(x = beta, y = theta_val, fill = fossil_prop, color = transparent),
               position = position_dodge2(padding = .2, preserve = "single")) +
  scale_y_continuous("Estimated Theta") +
  scale_x_discrete("Fossil Sampling Bias") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_identity(NULL, ) +
  theme_bw(base_size = 20) +
  facet_grid(rows = vars(model_theta), cols = vars(n_tip), scales = "free_y",
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg4 <- ggarrange2(gg4a, gg4b, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Thetas.pdf", gg4, width = 16, height = 20)

## half-life plot ----------------------------------------------
# half-lives relative to max tree height
correct_rhls <- data.frame(model = factor(c("wOUc", "sOUc", "wOUs", "sOUs"),
                                          levels = c("wOUc", "sOUc", "wOUs", "sOUs")),
                           rel_hl = c(1, .1, 1, .1))

gg5a <- ggplot(param_estimates_df %>%
                 filter(model %in% c("wOUc", "wOUs", "sOUc", "sOUs"), mu == 0.25)) +
  geom_hline(data = correct_rhls, aes(yintercept = rel_hl)) +
  geom_violin(aes(x = beta, y = rel_hl, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Relative Phylogenetic Half-life") +
  scale_x_discrete("Fossil Sampling Bias") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(rows = vars(model), cols = vars(n_tip), scales = "free_y",
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg5b <- ggplot(param_estimates_df %>%
                 filter(model %in% c("wOUc", "sOUc", "wOUs", "sOUs"), mu == 0.9)) +
  geom_hline(data = correct_rhls, aes(yintercept = rel_hl)) +
  geom_violin(aes(x = beta, y = rel_hl, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Relative Phylogenetic Half-life") +
  scale_x_discrete("Fossil Sampling Bias") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(rows = vars(model), cols = vars(n_tip), scales = "free_y",
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg5 <- ggarrange2(gg5a, gg5b, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Halflives.pdf", gg5, width = 22, height = 20)

## trend plot ----------------------------------------------
correct_trends <- data.frame(model = factor(c("wtrend", "strend"), levels = c("wtrend", "strend")),
                             trend = c(0.1, 0.3))

gg6a <- ggplot(param_estimates_df %>% filter(model %in% c("wtrend", "strend"), mu == 0.25)) +
  geom_hline(data = correct_trends, aes(yintercept = trend)) +
  geom_violin(aes(x = beta, y = trend, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Trend") +
  scale_x_discrete("Fossil Sampling Bias") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(rows = vars(model), cols = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips"))) +
  coord_cartesian(ylim = c(-0.1, 0.5))
gg6b <- ggplot(param_estimates_df %>% filter(model %in% c("wtrend", "strend"), mu == 0.9)) +
  geom_hline(data = correct_trends, aes(yintercept = trend)) +
  geom_violin(aes(x = beta, y = trend, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Trend") +
  scale_x_discrete("Fossil Sampling Bias") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(rows = vars(model), cols = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips"))) +
  coord_cartesian(ylim = c(-0.1, 0.5))
gg6 <- ggarrange2(gg6a, gg6b, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Trends.pdf", gg6, width = 16, height = 12)

## beta plot -----------------------------------------------
correct_betas <- data.frame(model = factor(c("wAC", "sAC", "wDC", "sDC"), levels = c("wAC", "sAC", "wDC", "sDC")),
                            exp_rate = c(0.1, 0.3, -0.1, -0.3))

gg7a <- ggplot(param_estimates_df %>% filter(model %in% c("wAC", "sAC", "wDC", "sDC"), mu == 0.25)) +
  geom_hline(data = correct_betas, aes(yintercept = exp_rate)) +
  geom_violin(aes(x = beta, y = exp_rate, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Beta") +
  scale_x_discrete("Fossil Sampling Bias") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(rows = vars(model), cols = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg7b <- ggplot(param_estimates_df %>% filter(model %in% c("wAC", "sAC", "wDC", "sDC"), mu == 0.9)) +
  geom_hline(data = correct_betas, aes(yintercept = exp_rate)) +
  geom_violin(aes(x = beta, y = exp_rate, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Beta") +
  scale_x_discrete("Fossil Sampling Bias") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_grid(rows = vars(model), cols = vars(n_tip),
             labeller = labeller(n_tip = function(x) paste(x, "tips")))
gg7 <- ggarrange2(gg7a, gg7b, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Betas.pdf", gg7, width = 16, height = 12)

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
  mutate(beta = fct_recode(as.factor(beta), old = "-3", even = "0", young = "3"))

gg5 <- ggplot(fossil_heights) +
  geom_histogram(aes(x = height, fill = factor(mu)), position = "dodge",
                 binwidth = 0.05, boundary = 1) +
  scale_y_continuous("# of Fossils") +
  scale_x_continuous("Relative Height in Phylogeny") +
  scale_fill_discrete("mu") +
  theme_bw(base_size = 20) +
  facet_grid(cols = vars(beta), rows = vars(n_tip), scales = "free_y",
             labeller = labeller(n_tip = function(x) paste(x, "tips"))) +
  theme(legend.position = "top")
ggsave("./figures/Fossil_heights.pdf", gg5, width = 12, height = 8)

