# Written by Will Gearty, Bethany Allen, and Pedro Godoy

# Load libraries
library(ape)
library(phytools)
library(TreeSim)
library(FossilSim)
library(geiger)
library(mvMORPH)
library(pbapply)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
#devtools::install_github("willgearty/pcmtools")
library(pcmtools)
library(deeptime)
library(future)

# Load functions
source("R/sim.fossils.R")

# Load trees if missing ----------------------------------------------------
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
model_results <- list()
for (mod in mods) {
  model_results[[mod]] <- readRDS(paste0("./data/model_fitting/model_fitting_", mod, ".RDS"))
}

# extract  to data.frames
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
  group_by(model, n_tip, fossil_prop, lambda, mu, sim) %>%
  mutate(aicc_w = AICweights(aicc)) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = mods),
         fit_model = factor(fit_model, levels = fit_models),
         across(c(n_tip, fossil_prop, lambda, mu, sim, aicc_w), ~ as.numeric(.x))) %>%
  mutate(across(c(n_tip, fossil_prop, lambda, mu), ~ as.factor(.x)))

param_estimates_df <- lapply(model_results, \(mod) {
  lapply(mod, \(tree) {
    lapply(tree, \(fit_model) {
      thetas <- rep_len(fit_model$theta, 2)
      bind_cols(alpha = fit_model$alpha[[1]], exp_rate = fit_model$beta[[1]],
                sigma = fit_model$sigma[[1]],
                theta_0 = thetas[1], theta_1 = thetas[2],
                convergence = fit_model$convergence)
    }) %>% bind_rows(.id = "fit_model")
  }) %>% bind_rows() %>% cbind(tree_df %>% select(-tree) %>% slice(rep(1:n(), each = 5)), .)
}) %>% bind_rows(.id = "model")

param_estimates_df <- param_estimates_df %>%
  mutate(model = factor(model, levels = mods),
         fit_model = factor(fit_model, levels = fit_models),
         across(c(n_tip, fossil_prop, lambda, mu, sim), ~ as.numeric(.x))) %>%
  mutate(across(c(n_tip, fossil_prop, lambda, mu), ~ as.factor(.x))) %>%
  mutate(correct_model = case_when(
    model %in% c("wBM", "sBM") ~ "BM",
    model %in% c("wtrend", "strend") ~ "trend",
    model %in% c("wOUc", "sOUc") ~ "OU1",
    model %in% c("wOUs", "sOUs") ~ "OU2",
    model %in% c("wAC", "sAC", "wDC", "sDC") ~ "ACDC"
  )) %>%
  filter(fit_model == correct_model) %>%
  select(-correct_model)

theta_estimates_df_long <- param_estimates_df %>%
  pivot_longer(cols = c(theta_0, theta_1), names_to = "theta", values_to = "theta_val",
               names_prefix = "theta_") %>%
  mutate(model_theta = paste(model, theta, sep = "_")) %>%
  filter(model_theta %in% c("wOUc_0", "sOUc_0", "wOUs_0", "sOUs_0", "wOUs_1", "sOUs_1")) %>%
  group_by(model_theta) %>%
  filter(!theta_val %in% boxplot.stats(theta_val)$out) %>% # remove outliers
  ungroup()

# clean up the old object
remove(model_fits)

# make some plots!

gg1 <- ggplot(model_fits_df_long %>% filter(mu == "0.25")) +
  geom_violin(aes(x = factor(fossil_prop), y = aicc_w, color = fit_model)) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(rows = vars(n_tip), cols = vars(model)) +
  theme_bw(base_size = 20)
ggsave("./figures/AIC.pdf", gg1, width = 40, height = 20)

model_fits_df_summ <- model_fits_df_long %>%
  mutate(correct_model = case_when(
    model %in% c("wBM", "sBM") ~ "BM",
    model %in% c("wtrend", "strend") ~ "trend",
    model %in% c("wOUc", "sOUc") ~ "OU1",
    model %in% c("wOUs", "sOUs") ~ "OU2",
    model %in% c("wAC", "sAC", "wDC", "sDC") ~ "ACDC"
  )) %>%
  group_by(model, n_tip, fossil_prop, lambda, mu, sim) %>%
  summarise(correct = ifelse(all(is.na(aicc_w)), NA,
                             any(unique(correct_model) ==
                                   fit_model[which(aicc_w == max(aicc_w, na.rm = TRUE))])),
            .groups = "drop") %>%
  group_by(model, n_tip, fossil_prop, lambda, mu) %>%
  summarise(prop_true = sum(correct)/n(), .groups = "drop")

gg2 <- ggplot(model_fits_df_summ) +
  geom_line(aes(x = fossil_prop, y = prop_true, color = n_tip,
                linetype = mu, group = interaction(n_tip, mu))) +
  scale_x_discrete("Proportion of Fossils in Tree") +
  scale_y_continuous("Prop. of Simulations with Correct Best Fit Model", limits = c(0, 1)) +
  scale_color_brewer("# of tips", palette = "Dark2") +
  scale_linetype_discrete("Rel. Ext. Rate") +
  theme_bw(base_size = 20) +
  facet_wrap(~model)
ggsave("./figures/Prop_Best.pdf", gg2, width = 16, height = 12)

correct_sigmas <- data.frame(model = factor(mods, levels = mods),
                             sigma = c(0.1, 0.5, rep(0.1, 6), rep(0.001, 4)))
gg3a <- ggplot(param_estimates_df %>%
                filter(! model %in% c("wAC", "sAC", "wDC", "sDC"), mu == 0.25)) +
  geom_hline(data = correct_sigmas %>%
               filter(! model %in% c("wAC", "sAC", "wDC", "sDC")),
             aes(yintercept = sigma)) +
  geom_violin(aes(x = n_tip, y = sigma, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Sigma") +
  scale_x_discrete("Number of Tips") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_wrap(~model, nrow = 2) +
  coord_cartesian(ylim = c(0, 1))
gg3b <- ggplot(param_estimates_df %>%
         filter(! model %in% c("wAC", "sAC", "wDC", "sDC"), mu == 0.9)) +
  geom_hline(data = correct_sigmas %>%
               filter(! model %in% c("wAC", "sAC", "wDC", "sDC")),
             aes(yintercept = sigma)) +
  geom_violin(aes(x = n_tip, y = sigma, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Sigma") +
  scale_x_discrete("Number of Tips") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_wrap(~model, nrow = 2) +
  coord_cartesian(ylim = c(0, 1))
gg3 <- ggarrange2(gg3a, gg3b, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Sigmas.pdf", gg3, width = 22, height = 20)

correct_thetas <- data.frame(model_theta = c("wOUc_0", "sOUc_0", "wOUs_0", "sOUs_0", "wOUs_1", "sOUs_1"),
                             theta_val = c(rep(0, 4), 1, 1))
gg4a <- ggplot(theta_estimates_df_long %>% filter(mu == 0.25)) +
  geom_hline(data = correct_thetas, aes(yintercept = theta_val)) +
  geom_violin(aes(x = n_tip, y = theta_val, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Theta") +
  scale_x_discrete("Number of Tips") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_wrap(~model_theta, nrow = 2, scales = "free_y")
gg4b <- ggplot(theta_estimates_df_long %>% filter(mu == 0.9)) +
  geom_hline(data = correct_thetas, aes(yintercept = theta_val)) +
  geom_violin(aes(x = n_tip, y = theta_val, fill = fossil_prop, color = fossil_prop)) +
  scale_y_continuous("Estimated Theta") +
  scale_x_discrete("Number of Tips") +
  scale_fill_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  scale_color_brewer("prop. of tips\nthat are\nfossils", palette = "Dark2") +
  theme_bw(base_size = 20) +
  facet_wrap(~model_theta, nrow = 2, scales = "free_y")
gg4 <- ggarrange2(gg4a, gg4b, nrow = 2, draw = FALSE, labels = c("mu = 0.25", "mu = 0.9"))
ggsave("./figures/Thetas.pdf", gg4, width = 22, height = 20)

# proportional tip heights of fossils
tree_df_scaled <- readRDS("./data/tree_simulations_scaled.RDS")
fossil_heights <- function(tree_list) {
  lapply(tree_list, function(tree) {
    tree$edge.length <- tree$edge.length / max(nodeHeights(tree)[, 2])
    n_heights <- nodeHeights(tree)[which(tree$edge <= Ntip(tree))]
    n_heights[n_heights < (1 - 0.000001)] # because computer math
  }) %>% do.call(c, .)
}
hist(fossil_heights(tree_df_scaled$tree))
