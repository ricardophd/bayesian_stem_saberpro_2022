# Reproduce Bayesian STEM selection analysis from the packaged processed data.
# This script is designed for the GitHub bundle and does not require the raw
# Saber 11, Saber Pro, crosswalk, or SNIES files.

find_repo_root <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
    return(normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE))
  }

  sourced_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(sourced_file)) {
    script_dir <- dirname(normalizePath(sourced_file, winslash = "/"))
    return(normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE))
  }

  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(wd, "data", "processed", "analytic_model_data.csv"))) {
    return(wd)
  }
  if (basename(wd) == "scripts") {
    return(normalizePath(file.path(wd, ".."), winslash = "/", mustWork = FALSE))
  }
  stop("Could not infer repository root. Run from the repository root or use Rscript scripts/02_reproduce_from_processed.R.")
}

ROOT <- find_repo_root()
setwd(ROOT)

LOCAL_LIB <- file.path(ROOT, "r_libs")
if (dir.exists(LOCAL_LIB)) {
  .libPaths(c(normalizePath(LOCAL_LIB, winslash = "/"), .libPaths()))
}

required_packages <- c("data.table", "ggplot2", "scales", "MASS")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "), "\n",
    "Run scripts/00_install_packages.R first, or install them manually."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

set.seed(as.integer(Sys.getenv("BAYES_STEM_SEED", "20260616")))
POSTERIOR_DRAWS <- as.integer(Sys.getenv("BAYES_STEM_DRAWS", "8000"))
ROPE <- as.numeric(Sys.getenv("BAYES_STEM_ROPE", "0.01"))

DATA_FILE <- file.path(ROOT, "data", "processed", "analytic_model_data.csv")
OUT_DIR <- file.path(ROOT, "outputs")
FIG_DIR <- file.path(OUT_DIR, "figures")
TABLE_DIR <- file.path(OUT_DIR, "tables")
MODEL_DIR <- file.path(OUT_DIR, "models")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)

quintile_levels <- paste0("Q", 1:5)
subject_specs <- data.table(
  subject = c("math", "science", "language"),
  model_col = c("math_q", "science_q", "language_q"),
  label = c("Mathematics", "Natural sciences", "Language")
)

theme_publication <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

posterior_interval <- function(draws) {
  c(
    median = stats::median(draws),
    lower = stats::quantile(draws, 0.025, names = FALSE),
    upper = stats::quantile(draws, 0.975, names = FALSE)
  )
}

order_subject_facets <- function(dt) {
  dt <- copy(dt)
  dt[, subject_label := factor(subject_label, levels = subject_specs$label)]
  dt
}

summarise_beta_binomial <- function(dt) {
  long_counts <- rbindlist(lapply(seq_len(nrow(subject_specs)), function(i) {
    q_col <- subject_specs$model_col[[i]]
    out <- dt[, .(
      n = .N,
      stem = sum(STEM_AREA)
    ), by = .(gender, quintile = get(q_col))]
    out[, subject := subject_specs$subject[[i]]]
    out[, subject_label := subject_specs$label[[i]]]
    out[]
  }), use.names = TRUE)

  long_counts[, non_stem := n - stem]
  long_counts[, observed_rate := stem / n]

  summary <- rbindlist(lapply(seq_len(nrow(long_counts)), function(i) {
    d <- rbeta(POSTERIOR_DRAWS, long_counts$stem[[i]] + 1, long_counts$non_stem[[i]] + 1)
    ints <- posterior_interval(d)
    data.table(
      subject = long_counts$subject[[i]],
      subject_label = long_counts$subject_label[[i]],
      gender = long_counts$gender[[i]],
      quintile = long_counts$quintile[[i]],
      n = long_counts$n[[i]],
      stem = long_counts$stem[[i]],
      observed_rate = long_counts$observed_rate[[i]],
      posterior_median = ints[["median"]],
      lower_95 = ints[["lower"]],
      upper_95 = ints[["upper"]]
    )
  }))

  contrasts <- rbindlist(lapply(split(summary, by = c("subject", "quintile"), keep.by = TRUE), function(cell) {
    if (!all(c("Men", "Women") %in% as.character(cell$gender))) {
      return(NULL)
    }
    men <- cell[gender == "Men"]
    women <- cell[gender == "Women"]
    men_draws <- rbeta(POSTERIOR_DRAWS, men$stem + 1, men$n - men$stem + 1)
    women_draws <- rbeta(POSTERIOR_DRAWS, women$stem + 1, women$n - women$stem + 1)
    diff <- women_draws - men_draws
    ints <- posterior_interval(diff)
    data.table(
      subject = women$subject,
      subject_label = women$subject_label,
      quintile = women$quintile,
      diff_median = ints[["median"]],
      lower_95 = ints[["lower"]],
      upper_95 = ints[["upper"]],
      pr_women_lower = mean(diff < 0),
      pr_women_higher = mean(diff > 0),
      pr_equivalent_1pp = mean(abs(diff) < ROPE)
    )
  }), fill = TRUE)

  list(summary = summary, contrasts = contrasts)
}

make_positive_definite <- function(vcov_mat) {
  eig <- eigen(vcov_mat, symmetric = TRUE)
  values <- pmax(eig$values, 1e-8)
  eig$vectors %*% diag(values, nrow = length(values)) %*% t(eig$vectors)
}

draw_mvn <- function(mu, vcov_mat, ndraws) {
  vcov_pd <- make_positive_definite(vcov_mat)
  z <- matrix(rnorm(ndraws * length(mu)), nrow = ndraws, ncol = length(mu))
  draws <- sweep(z %*% chol(vcov_pd), 2, mu, "+")
  colnames(draws) <- names(mu)
  draws
}

fit_bayesian_logistic_laplace <- function(dt) {
  cells <- dt[, .(
    stem = sum(STEM_AREA),
    n = .N
  ), by = .(gender, math_q, science_q, language_q, nse_model, sb11_period)]
  cells[, non_stem := n - stem]
  cells <- droplevels(cells)

  rhs_formula <- ~ gender * math_q + gender * science_q + gender * language_q + nse_model + sb11_period
  X_full <- model.matrix(rhs_formula, data = cells)
  keep_cols <- colSums(abs(X_full)) > 0
  keep_cols["(Intercept)"] <- TRUE
  X <- X_full[, keep_cols, drop = FALSE]
  y <- cells$stem
  n <- cells$n

  prior_sd <- rep(1.5, ncol(X))
  names(prior_sd) <- colnames(X)
  prior_sd["(Intercept)"] <- 2.5
  prior_sd[grepl(":", names(prior_sd), fixed = TRUE)] <- 1.0

  log1pexp <- function(z) ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))

  neg_logpost <- function(beta) {
    eta <- as.vector(X %*% beta)
    if (any(!is.finite(eta))) return(Inf)
    log_lik <- sum(y * eta - n * log1pexp(eta))
    log_prior <- sum(stats::dnorm(beta, mean = 0, sd = prior_sd, log = TRUE))
    -(log_lik + log_prior)
  }

  neg_grad <- function(beta) {
    eta <- as.vector(X %*% beta)
    p <- stats::plogis(eta)
    grad_logpost <- as.vector(crossprod(X, y - n * p)) - beta / (prior_sd^2)
    -grad_logpost
  }

  start <- rep(0, ncol(X))
  names(start) <- colnames(X)
  glm_fit <- try(
    glm(cbind(stem, non_stem) ~ gender * math_q + gender * science_q + gender * language_q + nse_model + sb11_period,
        data = cells,
        family = binomial()),
    silent = TRUE
  )
  if (!inherits(glm_fit, "try-error")) {
    glm_coef <- coef(glm_fit)
    common <- intersect(names(start), names(glm_coef))
    start[common] <- ifelse(is.finite(glm_coef[common]), glm_coef[common], 0)
  }

  opt <- optim(
    par = start,
    fn = neg_logpost,
    gr = neg_grad,
    method = "BFGS",
    control = list(maxit = 1000, reltol = 1e-9)
  )
  if (opt$convergence != 0L) {
    warning("optim did not fully converge. Code: ", opt$convergence, "; message: ", opt$message)
  }

  beta_hat <- opt$par
  eta_hat <- as.vector(X %*% beta_hat)
  p_hat <- stats::plogis(eta_hat)
  weights <- n * p_hat * (1 - p_hat)
  hessian <- crossprod(X * weights, X) + diag(1 / (prior_sd^2), nrow = length(prior_sd))
  vcov_mat <- tryCatch(solve(hessian), error = function(e) MASS::ginv(hessian))
  beta_draws <- draw_mvn(beta_hat, vcov_mat, POSTERIOR_DRAWS)

  coef_summary <- data.table(
    term = names(beta_hat),
    estimate = as.numeric(beta_hat),
    sd = sqrt(diag(vcov_mat)),
    lower_95 = beta_hat + stats::qnorm(0.025) * sqrt(diag(vcov_mat)),
    upper_95 = beta_hat + stats::qnorm(0.975) * sqrt(diag(vcov_mat)),
    prior_sd = prior_sd
  )

  list(
    formula = rhs_formula,
    cells = cells,
    coef_summary = coef_summary,
    beta_mode = beta_hat,
    vcov = vcov_mat,
    beta_draws = beta_draws,
    coef_names = names(beta_hat),
    opt = opt
  )
}

newdata_for_subject <- function(model_dt, subject_id) {
  row <- subject_specs[subject == subject_id]
  if (nrow(row) != 1L) stop("Unknown subject: ", subject_id)

  out <- CJ(
    gender = factor(c("Men", "Women"), levels = levels(model_dt$gender)),
    quintile = quintile_levels,
    sorted = FALSE
  )
  out[, math_q := factor("Q3", levels = quintile_levels)]
  out[, science_q := factor("Q3", levels = quintile_levels)]
  out[, language_q := factor("Q3", levels = quintile_levels)]
  out[, (row$model_col) := factor(quintile, levels = quintile_levels)]
  out[, nse_model := factor(names(which.max(table(model_dt$nse_model))), levels = levels(model_dt$nse_model))]
  out[, sb11_period := factor(names(which.max(table(model_dt$sb11_period))), levels = levels(model_dt$sb11_period))]
  out[, subject := row$subject]
  out[, subject_label := row$label]
  out[]
}

predict_posterior <- function(fit, newdata) {
  X_new_full <- model.matrix(fit$formula, data = newdata)
  missing_cols <- setdiff(fit$coef_names, colnames(X_new_full))
  if (length(missing_cols) > 0L) {
    X_new_full <- cbind(
      X_new_full,
      matrix(0, nrow = nrow(X_new_full), ncol = length(missing_cols), dimnames = list(NULL, missing_cols))
    )
  }
  X_new <- X_new_full[, fit$coef_names, drop = FALSE]
  stats::plogis(X_new %*% t(fit$beta_draws))
}

summarise_predictions <- function(fit, model_dt) {
  pieces <- lapply(subject_specs$subject, function(s) {
    nd <- newdata_for_subject(model_dt, s)
    draws <- predict_posterior(fit, nd)
    pred_summary <- rbindlist(lapply(seq_len(nrow(nd)), function(i) {
      ints <- posterior_interval(draws[i, ])
      data.table(
        subject = nd$subject[[i]],
        subject_label = nd$subject_label[[i]],
        gender = nd$gender[[i]],
        quintile = nd$quintile[[i]],
        predicted_median = ints[["median"]],
        lower_95 = ints[["lower"]],
        upper_95 = ints[["upper"]]
      )
    }))

    contrasts <- rbindlist(lapply(quintile_levels, function(q) {
      men_idx <- which(nd$gender == "Men" & nd$quintile == q)
      women_idx <- which(nd$gender == "Women" & nd$quintile == q)
      diff <- draws[women_idx, ] - draws[men_idx, ]
      ints <- posterior_interval(diff)
      data.table(
        subject = nd$subject[[women_idx]],
        subject_label = nd$subject_label[[women_idx]],
        quintile = q,
        diff_median = ints[["median"]],
        lower_95 = ints[["lower"]],
        upper_95 = ints[["upper"]],
        pr_women_lower = mean(diff < 0),
        pr_women_higher = mean(diff > 0),
        pr_equivalent_1pp = mean(abs(diff) < ROPE)
      )
    }))

    offsets <- rbindlist(lapply(1:4, function(k) {
      men_idx <- which(nd$gender == "Men" & nd$quintile == quintile_levels[[k]])
      women_idx <- which(nd$gender == "Women" & nd$quintile == quintile_levels[[k + 1L]])
      diff <- draws[women_idx, ] - draws[men_idx, ]
      ints <- posterior_interval(diff)
      data.table(
        subject = nd$subject[[women_idx]],
        subject_label = nd$subject_label[[women_idx]],
        comparison = paste0("Women ", quintile_levels[[k + 1L]], " - Men ", quintile_levels[[k]]),
        diff_median = ints[["median"]],
        lower_95 = ints[["lower"]],
        upper_95 = ints[["upper"]],
        pr_women_higher = mean(diff > 0)
      )
    }))

    list(summary = pred_summary, contrasts = contrasts, offsets = offsets)
  })

  list(
    summary = rbindlist(lapply(pieces, `[[`, "summary")),
    contrasts = rbindlist(lapply(pieces, `[[`, "contrasts")),
    offsets = rbindlist(lapply(pieces, `[[`, "offsets"))
  )
}

save_plot <- function(plot, filename, width = 8, height = 5.5) {
  ggsave(file.path(FIG_DIR, paste0(filename, ".png")), plot, width = width, height = height, dpi = 600)
  ggsave(file.path(FIG_DIR, paste0(filename, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
}

plot_observed <- function(obs_summary) {
  obs_summary <- order_subject_facets(obs_summary)
  ggplot(obs_summary, aes(x = quintile, y = posterior_median, color = gender, group = gender)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.12, linewidth = 0.45) +
    facet_wrap(~ subject_label, nrow = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
    scale_color_manual(values = c("Men" = "#2D5F8B", "Women" = "#C85A3E")) +
    labs(
      x = "Saber 11 quintile",
      y = "Posterior STEM rate"
    ) +
    theme_publication()
}

plot_predictions <- function(pred_summary) {
  pred_summary <- order_subject_facets(pred_summary)
  ggplot(pred_summary, aes(x = quintile, y = predicted_median, color = gender, group = gender)) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 2.2) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = gender), alpha = 0.16, color = NA) +
    facet_wrap(~ subject_label, nrow = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
    scale_color_manual(values = c("Men" = "#2D5F8B", "Women" = "#C85A3E")) +
    scale_fill_manual(values = c("Men" = "#2D5F8B", "Women" = "#C85A3E")) +
    labs(
      x = "Saber 11 quintile",
      y = "Predicted probability"
    ) +
    theme_publication()
}

plot_differences <- function(contrast_summary) {
  contrast_summary <- order_subject_facets(contrast_summary)
  ggplot(contrast_summary, aes(x = quintile, y = diff_median)) +
    geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
    geom_pointrange(aes(ymin = lower_95, ymax = upper_95), color = "#6A3D2A", linewidth = 0.55) +
    facet_wrap(~ subject_label, nrow = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      x = "Saber 11 quintile",
      y = "Probability difference"
    ) +
    theme_publication()
}

plot_offset_contrasts <- function(offset_summary) {
  offset_summary <- order_subject_facets(offset_summary)
  ggplot(offset_summary, aes(x = comparison, y = diff_median)) +
    geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
    geom_pointrange(aes(ymin = lower_95, ymax = upper_95), color = "#4D4D4D", linewidth = 0.5) +
    facet_wrap(~ subject_label, nrow = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      x = NULL,
      y = "Probability difference"
    ) +
    coord_flip() +
    theme_publication()
}

message("Reading processed data: ", DATA_FILE)
analysis_dt <- fread(DATA_FILE)
analysis_dt[, gender := factor(gender, levels = c("Men", "Women"))]
analysis_dt[, nse_model := factor(nse_model)]
analysis_dt[, sb11_period := factor(sb11_period)]
analysis_dt[, math_q := factor(math_q, levels = quintile_levels)]
analysis_dt[, science_q := factor(science_q, levels = quintile_levels)]
analysis_dt[, language_q := factor(language_q, levels = quintile_levels)]
analysis_dt[, STEM_AREA := as.integer(STEM_AREA)]

fwrite(
  analysis_dt[, .(n = .N, stem = sum(STEM_AREA), stem_rate = mean(STEM_AREA)), by = .(gender)],
  file.path(TABLE_DIR, "03_sample_summary_by_gender.csv")
)

message("Analytic sample size: ", nrow(analysis_dt))
message("STEM rate: ", percent(mean(analysis_dt$STEM_AREA), accuracy = 0.1))

message("Summarising descriptive beta-binomial posteriors ...")
beta_results <- summarise_beta_binomial(analysis_dt)
fwrite(beta_results$summary, file.path(TABLE_DIR, "04_observed_beta_binomial_stem_rates.csv"))
fwrite(beta_results$contrasts, file.path(TABLE_DIR, "05_observed_beta_binomial_gender_contrasts.csv"))
save_plot(plot_observed(beta_results$summary), "01_observed_stem_rate_by_quintile_gender", width = 9, height = 4.8)

message("Fitting Bayesian logistic model ...")
bayes_fit <- fit_bayesian_logistic_laplace(analysis_dt)
fwrite(bayes_fit$coef_summary, file.path(TABLE_DIR, "06_bayesian_logistic_coefficients.csv"))
saveRDS(bayes_fit, file.path(MODEL_DIR, "bayesian_logistic_laplace_fit.rds"))

message("Generating posterior predictions and contrasts ...")
pred_results <- summarise_predictions(bayes_fit, analysis_dt)
fwrite(pred_results$summary, file.path(TABLE_DIR, "07_bayesian_predicted_probabilities.csv"))
fwrite(pred_results$contrasts, file.path(TABLE_DIR, "08_bayesian_gender_contrasts_same_quintile.csv"))
fwrite(pred_results$offsets, file.path(TABLE_DIR, "09_bayesian_one_quintile_offset_contrasts.csv"))

save_plot(plot_predictions(pred_results$summary), "02_bayesian_predicted_probability_by_quintile_gender", width = 9, height = 4.8)
save_plot(plot_differences(pred_results$contrasts), "03_bayesian_gender_difference_same_quintile", width = 9, height = 4.8)
save_plot(plot_offset_contrasts(pred_results$offsets), "04_one_quintile_offset_contrast", width = 9.5, height = 5.2)

run_info <- data.table(
  item = c("posterior_draws", "rope_probability_points", "analytic_n", "stem_rate", "input_data", "output_folder"),
  value = c(
    as.character(POSTERIOR_DRAWS),
    as.character(ROPE),
    as.character(nrow(analysis_dt)),
    as.character(mean(analysis_dt$STEM_AREA)),
    DATA_FILE,
    OUT_DIR
  )
)
fwrite(run_info, file.path(TABLE_DIR, "10_run_info.csv"))

message("Done.")
