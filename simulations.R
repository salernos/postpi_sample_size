#===============================================================================
#
#  PROGRAM: simulations.R
#
#  AUTHORS: Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: R script to reproduce the main simulation results for the 
#           manuscript entitled:
#
#           Sample Size Considerations for Post-Prediction Inference
#
#  UPDATED: 2024-08-03
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- LOAD NECESSARY PACKAGES ---------------------------------------------------

devtools::install_github("ipd-tools/ipd")

library(pacman)

p_load(ipd, gam, ggridges, grid, gtable, tidyverse, update = F)

#=== HELPER FUNCTIONS ==========================================================

#--- DATA GENERATING MECHANISM -------------------------------------------------

dgm_main <- function(true_beta, x1, x2, x3, x4, err) {

  return(true_beta * x1 + 0.5 * x2 + 3 * x3^3 + 4 * x4^2 + err)
}

#--- SIMULATION FUNCTION -------------------------------------------------------

sim_func <- function(nsims = 100, n_t = 100, n_l = 100, n_u = 100, nboot = 100,

  dgm, true_beta) {

  n_tot <- n_t + n_l + n_u

  est_oracle <- est_naive <- est_classic <- est_postpi1 <- est_postpi2 <-

    est_ppi <- est_plusplus <- est_popinf <- rep(NA, nsims)

  err_oracle <- err_naive <- err_classic <- err_postpi1 <- err_postpi2 <-

    err_ppi <- err_plusplus <- err_popinf <- rep(NA, nsims)

  for(sim in 1:nsims) {

    if (sim %% 100 == 0) cat("sim:", sim, "of", nsims, "\n")

    #--- SIMULATE DATA ---------------------------------------------------------

    #-- GENERATE COVARIATES

    x1 <- rnorm(n_tot)
    x2 <- rnorm(n_tot)
    x3 <- rnorm(n_tot)
    x4 <- rnorm(n_tot)

    err <- rnorm(n_tot, sd = 4)

    #-- GENERATE OUTCOME FROM THE DATA GENERATING FUNCTION

    y <- dgm(true_beta, x1, x2, x3, x4, err)

    #-- CREATE INDICATOR OF TRAINING, LABELED, AND UNLABELED DATA

    set_ind <- c(rep("train", n_t), rep("lab", n_l), rep("unlab", n_u))

    #-- CONSTRUCT FULL DATASET

    dat <- data.frame(y, x1, x2, x3, x4, set_ind)

    #-- TRAINING DATA

    dat_t <- dat |>

      filter(set_ind == "train") |>

      select(y, x1, x2, x3, x4)

    #-- FIT PREDICTION MODEL

    fit_pred <- gam(

      y ~ s(x1, df = 3) + s(x2, df = 3) + s(x3, df = 3) + s(x4, df = 3),

      data = dat_t)

    #-- ADD PREDICTIONS TO DATA

    dat$pred <- predict(fit_pred, dat)

    #-- LABELED DATA

    dat_l = dat |>

      filter(set_ind == "lab")

    #-- UNLABELED DATA

    dat_u = dat |>

      filter(set_ind == "unlab")

    #--- FIT VARIOUS METHODS ---------------------------------------------------

    #-- 0. ORACLE ESTIMATE

    fit_oracle <- lm(y ~ x1, data = dat_u)

    est_oracle[sim] <- fit_oracle$coefficients[2]

    err_oracle[sim] <- summary(fit_oracle)$coefficients[2, 2]

    #-- 1. NAIVE ESTIMATE

    fit_naive <- lm(pred ~ x1, data = dat_u)

    est_naive[sim] <- fit_naive$coefficients[2]

    err_naive[sim] <- summary(fit_naive)$coefficients[2, 2]

    #-- 2. CLASSIC ESTIMATE

    fit_classic <- lm(y ~ x1, data = dat_l)

    est_classic[sim] <- fit_classic$coefficients[2]

    err_classic[sim] <- summary(fit_classic)$coefficients[2, 2]

    #-- 3. ORIGINAL POSTPI ESTIMATE

    fit_postpi1 <- ipd(y - pred ~ x1, data = dat_l, unlabeled_data = dat_u,

      method = "postpi_boot", model = "ols", rel_func = "gam", scale_se = F,

      nboot = nboot)

    est_postpi1[sim] <- fit_postpi1$coefficients[2]

    err_postpi1[sim] <- fit_postpi1$se[2]

    #-- 4. SCALED POSTPI ESTIMATE

    fit_postpi2 <- ipd(y - pred ~ x1, data = dat_l, unlabeled_data = dat_u,

      method = "postpi_boot", model = "ols", rel_func = "gam", n_t = n_t,

      nboot = nboot)

    est_postpi2[sim] <- fit_postpi2$coefficients[2]

    err_postpi2[sim] <- fit_postpi2$se[2]

    #-- 5. PPI

    fit_ppi <- ipd(y - pred ~ x1, data = dat_l, unlabeled_data = dat_u,

      method = "ppi", model = "ols")

    est_ppi[sim] <- fit_ppi$coefficients[2]

    err_ppi[sim] <- fit_ppi$se[2]

    #-- 6. PPI++

    fit_plusplus <- ipd(y - pred ~ x1, data = dat_l,

      unlabeled_data = dat_u, method = "ppi_plusplus", model = "ols")

    est_plusplus[sim] <- fit_plusplus$coefficients[2]

    err_plusplus[sim] <- fit_plusplus$se[2]

    #-- 7. POPInf

    fit_popinf <- ipd(y - pred ~ x1, data = dat_l,

      unlabeled_data = dat_u, method = "popinf", model = "ols")

    est_popinf[sim] <- fit_popinf$coefficients[2]

    err_popinf[sim] <- fit_popinf$se[2]
  }

  results <- data.frame(

    est_oracle,   err_oracle,
    est_naive,    err_naive,
    est_classic,  err_classic,
    est_postpi1,  err_postpi1,
    est_postpi2,  err_postpi2,
    est_ppi,      err_ppi,
    est_plusplus, err_plusplus,
    est_popinf,   err_popinf)

  results <- results |>

    mutate(

      #- 0. Oracle

      upr_oracle  = est_oracle + 1.96 * err_oracle,
      lwr_oracle  = est_oracle - 1.96 * err_oracle,
      cov_oracle  = ((upr_oracle > true_beta) & (lwr_oracle < true_beta)),
      tval_oracle = est_oracle / err_oracle,
      pval_oracle = 2 * pt(-abs(tval_oracle), df = (n_u - 2)),

      #- 1. Naive

      upr_naive  = est_naive + 1.96 * err_naive,
      lwr_naive  = est_naive - 1.96 * err_naive,
      cov_naive  = ((upr_naive > true_beta) & (lwr_naive < true_beta)),
      tval_naive = est_naive / err_naive,
      pval_naive = 2 * pt(-abs(tval_naive), df = (n_u - 2)),

      #- 2. Classic

      upr_classic  = est_classic + 1.96 * err_classic,
      lwr_classic  = est_classic - 1.96 * err_classic,
      cov_classic  = ((upr_classic > true_beta) & (lwr_classic < true_beta)),
      tval_classic = est_classic / err_classic,
      pval_classic = 2 * pt(-abs(tval_classic), df = (n_l - 2)),

      #- 3. Original PostPI

      upr_postpi1  = est_postpi1 + 1.96 * err_postpi1,
      lwr_postpi1  = est_postpi1 - 1.96 * err_postpi1,
      cov_postpi1  = ((upr_postpi1 > true_beta) & (lwr_postpi1 < true_beta)),
      tval_postpi1 = est_postpi1 / err_postpi1,
      pval_postpi1 = 2 * pt(-abs(tval_postpi1), df = (n_u - 2)),

      #- 4. Scaled PostPI

      upr_postpi2  = est_postpi2 + 1.96 * err_postpi2,
      lwr_postpi2  = est_postpi2 - 1.96 * err_postpi2,
      cov_postpi2  = ((upr_postpi2 > true_beta) & (lwr_postpi2 < true_beta)),
      tval_postpi2 = est_postpi2 / err_postpi2,
      pval_postpi2 = 2 * pt(-abs(tval_postpi2), df = (n_u - 2)),

      #- 5. PPI

      upr_ppi  = est_ppi + 1.96 * err_ppi,
      lwr_ppi  = est_ppi - 1.96 * err_ppi,
      cov_ppi  = ((upr_ppi > true_beta) & (lwr_ppi < true_beta)),
      tval_ppi = est_ppi / err_ppi,
      pval_ppi = 2 * pt(-abs(tval_ppi), df = (n_l + n_u - 2)),

      #- 6. PPI++

      upr_plusplus  = est_plusplus + 1.96 * err_plusplus,
      lwr_plusplus  = est_plusplus - 1.96 * err_plusplus,
      cov_plusplus  = ((upr_plusplus > true_beta) & (lwr_plusplus < true_beta)),
      tval_plusplus = est_plusplus / err_plusplus,
      pval_plusplus = 2 * pt(-abs(tval_plusplus), df = (n_l + n_u - 2)),

      #- 7. POPInf

      upr_popinf  = est_popinf + 1.96 * err_popinf,
      lwr_popinf  = est_popinf - 1.96 * err_popinf,
      cov_popinf  = ((upr_popinf > true_beta) & (lwr_popinf < true_beta)),
      tval_popinf = est_popinf / err_popinf,
      pval_popinf = 2 * pt(-abs(tval_popinf), df = (n_l + n_u - 2))
    )

  return(results)
}

#=== SIMULATIONS ===============================================================

#--- SETTINGS ------------------------------------------------------------------

#-- Number of Simulated Replicates

nsims <- 1000

#-- Number of Bootstrap Replicates

nboot <- 100

#-- Effect Sizes

true_beta_opts <- 0:1

#-- Sample Sizes

n_opts <- 1:3

dgm_opts <- 1

settings <- expand.grid(true_beta_opts, n_opts, dgm_opts)

results <- vector("list", nrow(settings))

#--- LOOP OVER SIMULATION SETTINGS ---------------------------------------------

tic <- proc.time()

for(i in 1:length(results)) {

  cat("\nSetting", i, "of", length(results), "---\n\n")

  #-- Get Setting Parameters

  #- Effect Size

  true_beta <- settings[i, 1]

  #- Sample Sizes

  if (settings[i, 2] == 1) {

    n_t <- 500; n_l <- 500; n_u <- 500

  } else if (settings[i, 2] == 2) {

    n_t <- 100; n_l <- 500; n_u <- 1000

  } else if (settings[i, 2] == 3) {

    n_t <- 1000; n_l <- 100; n_u <- 500
  }

  #- Data generating Mechanism

  if (settings[i, 3] == 1) {

    dgm <- dgm_main

  } else if (settings[i, 3] == 2) {

    dgm <- dgm_wang

  } else if (settings[i, 3] == 3) {

    dgm <- dgm_miao

  } else if (settings[i, 3] == 4) {

    dgm <- dgm_angelopoulos
  }

  #-- Run Simulations

  results_sim <- sim_func(nsims, n_t, n_l, n_u, nboot, dgm, true_beta)

  #-- Store Results

  save(results_sim, file = paste0("sim_nt", n_t, "_nl", n_l, "_nu", n_u,

    "_beta", true_beta, "_dgm", deparse(substitute(dgm)), ".RData"))

  results[[i]] <- results_sim
}

toc <- proc.time()

cat("Total Time:", toc - tic, "\n")

#=== SUMMARIZE RESULTS =========================================================

#--- LABELS --------------------------------------------------------------------

cbp <- c(

  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

methods <- c(

  "oracle" = "Oracle", "naive" = "Naive", "classic" = "Classical",
  "ppi" = "PPI", "plusplus" = "PPI++", "popinf" = "POP-Inf",
  "postpi1" = "Original PostPI", "postpi2" = "Proposed Update")

n_vec <- c(

  "n[train] == 500*','~n == 500*','~N == 500",
  "n[train] == 100*','~n == 500*','~N == 1000",
  "n[train] == 1000*','~n == 100*','~N == 500")

#--- FIGURE 1: P-VALUES UNDER THE NULL -----------------------------------------

results_pval <- bind_rows(

  lapply(1:length(results), function(i) {

    results[[i]] |>

      mutate(

        beta1 = settings[i, 1],

        n = n_vec[settings[i, 2]] |>

          factor(levels = n_vec, labels = n_vec)) |>

      select(n, beta1, starts_with("pval_")) |>

      pivot_longer(starts_with("pval_"))})) |>

  filter(beta1 == 0) |>

  mutate(

    name = gsub("pval_", "", name) |>

      factor(levels = names(methods), labels = methods))

fig_pval <- results_pval |>

  ggplot(aes(sample = value, group = name, color = name, fill = name)) +

    theme_bw() + coord_equal() +

    facet_wrap(~ n, labeller = label_parsed) +

    geom_abline(slope = 1, intercept = 0, linetype = 2) +

    stat_qq(distribution = stats::qunif) +

    scale_color_manual(values = cbp) +

    scale_fill_manual(values = cbp) +

    xlim(0, 1) + ylim(0, 1) +

    labs(

      x = "\nTheoretical Quantiles",

      y = "Empirical Quantiles\n",

      color = "Method:",

      fill = "Method:") +

    theme(

      legend.text = element_text(size = 14),

      legend.title = element_text(face = "bold", size = 16),

      strip.background = element_rect(fill = cbp[6]),

      strip.text = element_text(face = "bold", color = "#FFFFFF", size = 14),

      axis.text = element_text(size = 14),

      axis.title = element_text(face = "bold", size = 16)) +

    guides(

      color = guide_legend(override.aes = list(size = 5)))

fig_pval

ggsave("fig_pval.png", fig_pval, width = 16, height = 8)

#--- FIGURE 2: COVERAGE UNDER THE ALTERNATIVE ----------------------------------

results_cov <- bind_rows(

  lapply(1:length(results), function(i) {

    results[[i]] |>

      mutate(

        beta1 = settings[i, 1],

        n = n_vec[settings[i, 2]] |>

          factor(levels = n_vec, labels = n_vec)) |>

      select(n, beta1, starts_with("cov_")) |>

      pivot_longer(starts_with("cov_"))})) |>

  filter(beta1 != 0) |>

  mutate(

    name = gsub("cov_", "", name) |>

      factor(levels = names(methods), labels = methods)) |>

  group_by(n, beta1, name) |>

  summarise(cov = mean(value))

fig_cov <- results_cov |>

  ggplot(aes(x = cov, y = name, group = name, color = name, fill = name)) +

    theme_bw() +

    facet_wrap(~ n, labeller = label_parsed) +

    geom_vline(xintercept = 0.95, linetype = 2) +

    geom_point(size = 5) +

    scale_y_discrete(limit = rev) +

    scale_color_manual(values = cbp) +

    scale_fill_manual(values = cbp) +

    xlim(0, 1) +

    labs(

      x = "\nCoverage Probability",

      y = "Method\n",

      color = "Method:",

      fill = "Method:") +

  theme(

    legend.text = element_text(size = 14),

    legend.title = element_text(face = "bold", size = 16),

    strip.background = element_rect(fill = cbp[6]),

    strip.text = element_text(face = "bold", color = "#FFFFFF", size = 14),

    axis.text = element_text(size = 14),

    axis.title = element_text(face = "bold", size = 16)) +

  guides(

    color = guide_legend(override.aes = list(size = 5)))

fig_cov

ggsave("fig_cov.png", fig_cov, width = 16, height = 8)

#--- FIGURE 3: POWER UNDER THE ALTERNATIVE -------------------------------------

results_power <- bind_rows(

  lapply(1:length(results), function(i) {

    results[[i]] |>

      mutate(

        beta1 = settings[i, 1],

        n = n_vec[settings[i, 2]] |>

          factor(levels = n_vec, labels = n_vec)) |>

      select(n, beta1, starts_with("pval_")) |>

      pivot_longer(starts_with("pval_"))})) |>

  filter(beta1 != 0) |>

  mutate(

    name = gsub("pval_", "", name) |>

      factor(levels = names(methods), labels = methods))

powers <- results_power |>

  mutate(samps = str_extract_all(n, "\\d+")) |>

  unnest_wider(samps, names_sep = "_") |>

  rename(nt = samps_1, nl = samps_2, nu = samps_3) |>

  mutate(across(nt:nu, as.numeric)) |>

  group_by(n, beta1, name) |>

  summarize(power = mean(value < 0.05),

    nt = mean(nt), nl = mean(nl), nu = mean(nu)) |>

  mutate(

    sampsize = case_when(

      name %in% c("Oracle", "Naive", "Original PostPI", "Proposed Update") ~ nu,

      name == "Classical" ~ nl,

      name %in% c("PPI", "PPI++", "POP-Inf") ~ nl + nu),

    label = paste0("\nPower: ", power*100, "%; Samples Used: ",

      prettyNum(sampsize, big.mark=",")))

fig_power <- results_power |>

  mutate(

    sampused = case_when(

      name %in% c("Oracle", "Naive", "Original PostPI", "Proposed Update") ~

        "n[u]~Samples",

      name == "Classical" ~ "n[l]~Samples",

      name %in% c("PPI", "PPI++", "POP-Inf") ~ "n[l] + n[u]~Samples") |>

      factor(levels = c("n[l]~Samples", "n[u]~Samples",

        "n[l] + n[u]~Samples"))) |>

  left_join(powers, by = c("n", "name")) |>

  ggplot(aes(x = value, y = name, group = name, color = name, fill = name)) +

  theme_bw() +

  facet_grid(vars(sampused), vars(n),

    labeller = label_parsed, scales = "free_y") +

  stat_density_ridges(from = 0, to = 1,

    alpha = 0.5, quantile_lines = T, linewidth = 1.1) +

  geom_text(#data = powers,

            color = "black",

    mapping = aes(x = 0.5, y = name, label = label), size = 4) +

  scale_y_discrete(limits = rev) +

  scale_color_manual(values = cbp) +

  scale_fill_manual(values = cbp) +

  xlim(0, 1) +

  labs(

    x = "\nP-Values",

    y = "Method\n",

    color = "Method:",

    fill = "Method:") +

  theme(

    legend.text = element_text(size = 14),

    legend.title = element_text(face = "bold", size = 16),

    strip.background = element_rect(fill = cbp[6]),

    strip.text = element_text(face = "bold", color = "#FFFFFF", size = 14),

    axis.text = element_text(size = 14),

    axis.title = element_text(face = "bold", size = 16)) +

  guides(

    color = guide_legend(override.aes = list(size = 5)))

fig_power

fig_power_gt <- ggplot_gtable(ggplot_build(fig_power))

gtable_show_layout(fig_power_gt)

fig_power_gt$heights[10] <- (1/4) * fig_power_gt$heights[10]
fig_power_gt$heights[14] <- (3/4) * fig_power_gt$heights[14]

grid.draw(fig_power_gt)

ggsave("fig_power.png", fig_power, width = 16, height = 8)

fig_power_simplified <- results_power |>

  mutate(

    sampused = case_when(

      name %in% c("Oracle", "Naive", "Original PostPI", "Proposed Update") ~

        "n[u]~Samples",

      name == "Classical" ~ "n[l]~Samples",

      name %in% c("PPI", "PPI++", "POP-Inf") ~ "n[l] + n[u]~Samples") |>

      factor(levels = c("n[l]~Samples", "n[u]~Samples",

                        "n[l] + n[u]~Samples"))) |>

  left_join(powers, by = c("n", "name")) |>

  filter(name %in% c("Classical", "Oracle", "Naive",

    "Original PostPI", "Proposed Update")) |>

  ggplot(aes(x = value, y = name, group = name, color = name, fill = name)) +

  theme_bw() +

  facet_grid(vars(sampused), vars(n),

             labeller = label_parsed, scales = "free_y") +

  stat_density_ridges(from = 0, to = 1,

                      alpha = 0.5, quantile_lines = T, linewidth = 1.1) +

  geom_text(#data = powers,

    color = "black",

    mapping = aes(x = 0.5, y = name, label = label), size = 4) +

  scale_y_discrete(limits = rev) +

  scale_color_manual(values = cbp) +

  scale_fill_manual(values = cbp) +

  xlim(0, 1) +

  labs(

    x = "\nP-Values",

    y = "Method\n",

    color = "Method:",

    fill = "Method:") +

  theme(

    legend.text = element_text(size = 14),

    legend.title = element_text(face = "bold", size = 16),

    strip.background = element_rect(fill = cbp[6]),

    strip.text = element_text(face = "bold", color = "#FFFFFF", size = 14),

    axis.text = element_text(size = 14),

    axis.title = element_text(face = "bold", size = 16)) +

  guides(

    color = guide_legend(override.aes = list(size = 5)))

fig_power_simplified

fig_power_simplified_gt <- ggplot_gtable(ggplot_build(fig_power_simplified))

gtable_show_layout(fig_power_simplified_gt)

fig_power_simplified_gt$heights[10] <- (1/4) * fig_power_simplified_gt$heights[10]

grid.draw(fig_power_simplified_gt)

ggsave("fig_power_simplified.png", fig_power_simplified, width = 16, height = 8)

#--- FIGURE 4: POINT ESTIMATES UNDER THE ALTERNATIVE ---------------------------

results_est <- bind_rows(

  lapply(1:length(results), function(i) {

    results[[i]] |>

      mutate(

        beta1 = settings[i, 1],

        n = n_vec[settings[i, 2]] |>

          factor(levels = n_vec, labels = n_vec)) |>

      select(n, beta1, starts_with("est_")) |>

      pivot_longer(starts_with("est_"))})) |>

  filter(beta1 != 0) |>

  mutate(

    name = gsub("est_", "", name) |>

      factor(levels = names(methods), labels = methods))

fig_est <- results_est |>

  ggplot(aes(y = name, x = value, group = name, fill = name)) +

    theme_bw() +

    facet_wrap(~ n, labeller = label_parsed) +

    geom_vline(xintercept = 1, linetype = 2) +

    geom_boxplot() +

    scale_y_discrete(limit = rev) +

    scale_color_manual(values = cbp) +

    scale_fill_manual(values = cbp) +

    labs(

      x = "\nEstimated Effect Size",

      y = "Method\n",

      color = "Method:",

      fill = "Method:") +

  theme(

    legend.text = element_text(size = 14),

    legend.title = element_text(face = "bold", size = 16),

    strip.background = element_rect(fill = cbp[6]),

    strip.text = element_text(face = "bold", color = "#FFFFFF", size = 14),

    axis.text = element_text(size = 14),

    axis.title = element_text(face = "bold", size = 16)) +

  guides(

    color = guide_legend(override.aes = list(size = 5)))

fig_est

ggsave("fig_est.png", fig_est, width = 16, height = 8)

#--- FIGURE 5: STANDARD ERRORS UNDER THE ALTERNATIVE ---------------------------

results_err <- bind_rows(

  lapply(1:length(results), function(i) {

    results[[i]] |>

      mutate(

        beta1 = settings[i, 1],

        n = n_vec[settings[i, 2]] |>

          factor(levels = n_vec, labels = n_vec)) |>

      select(n, beta1, starts_with("err_")) |>

      pivot_longer(starts_with("err_"))})) |>

  filter(beta1 != 0) |>

  mutate(

    name = gsub("err_", "", name) |>

      factor(levels = names(methods), labels = methods))

results_err_ref <- results_err |>

  filter(name == "Oracle") |>

  group_by(n) |>

  summarise(oracle_value = mean(value))

fig_err <- results_err |>

  left_join(results_err_ref, "n") |>

  ggplot(aes(y = name, x = value, group = name, fill = name)) +

    theme_bw() +

    facet_wrap(~ n, labeller = label_parsed) +

    geom_vline(aes(xintercept = oracle_value), linetype = 2) +

    geom_boxplot() +

    scale_y_discrete(limit = rev) +

    scale_color_manual(values = cbp) +

    scale_fill_manual(values = cbp) +

    labs(

      x = "\nEstimated Standard Error",

      y = "Method\n",

      color = "Method:",

      fill = "Method:") +

    theme(

      legend.text = element_text(size = 14),

      legend.title = element_text(face = "bold", size = 16),

      strip.background = element_rect(fill = cbp[6]),

      strip.text = element_text(face = "bold", color = "#FFFFFF", size = 14),

      axis.text = element_text(size = 14),

      axis.title = element_text(face = "bold", size = 16)) +

    guides(

      color = guide_legend(override.aes = list(size = 5)))

fig_err

ggsave("fig_err.png", fig_err, width = 16, height = 8)

#--- TABLE 1: FULL SIMULATION RESULTS ------------------------------------------

results_full <- bind_rows(

  lapply(1:length(results), function(i) {

    results[[i]] |>

      mutate(

        beta1 = settings[i, 1],

        n = n_vec[settings[i, 2]] |>

          factor(levels = n_vec, labels = n_vec))})) |>

  pivot_longer(cols = -c(n, beta1),

    names_to = c(".value", "method"), names_sep = "_")

results_summary <- results_full |>

  group_by(n, beta1, method) |>

  summarise(

    Bias = mean(est - beta1),

    MSE = mean((est - beta1)^2),

    `Average CI Width` = mean(upr - lwr),

    Coverage = mean(cov),

    `Type I Error/Power` = mean(pval < 0.05)) |>

  separate(n, into = c("n_t", "n_l", "n_u"), sep = "\\*','~") |>

  mutate(

    across(c(n_t, n_l, n_u), ~ as.numeric(gsub("[^0-9]", "", .))),

    across(c(Bias, MSE, `Average CI Width`), ~ sprintf("%.3f", .)),

    method = factor(method, levels = names(methods), labels = methods))

View(results_summary)

results_summary |>

  gt::gt() |>

  gt::as_latex() |>

  cat()

#=== END =======================================================================
