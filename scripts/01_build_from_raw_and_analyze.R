# Bayesian analysis: Saber 11 performance and STEM selection in Saber Pro 2022
# This script rebuilds the linked Saber 11 -> Saber Pro 2022 analytic sample,
# creates cohort-specific Saber 11 performance quintiles, and estimates the
# probability of choosing a STEM field by gender and performance quintile.

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

ROOT <- get_script_dir()
setwd(ROOT)

LOCAL_LIB <- file.path(ROOT, "r_libs_bayes_stem")
if (dir.exists(LOCAL_LIB)) {
  .libPaths(c(normalizePath(LOCAL_LIB, winslash = "/"), .libPaths()))
}

required_packages <- c("data.table", "readr", "openxlsx", "ggplot2", "scales", "MASS")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "), "\n",
    "Install them in the project library with:\n",
    "dir.create('r_libs_bayes_stem', showWarnings = FALSE)\n",
    ".libPaths(c(normalizePath('r_libs_bayes_stem'), .libPaths()))\n",
    "install.packages(c('data.table','readr','openxlsx','ggplot2','scales'), ",
    "lib=normalizePath('r_libs_bayes_stem'), repos='https://cloud.r-project.org')"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(openxlsx)
  library(ggplot2)
  library(scales)
})

set.seed(as.integer(Sys.getenv("BAYES_STEM_SEED", "20260616")))

QUICK_RUN <- identical(Sys.getenv("BAYES_STEM_QUICK", "0"), "1")
POSTERIOR_DRAWS <- as.integer(Sys.getenv("BAYES_STEM_DRAWS", ifelse(QUICK_RUN, "1000", "8000")))
QUICK_N_MAX <- as.integer(Sys.getenv("BAYES_STEM_NMAX", "50000"))
QUINTILE_REFERENCE <- Sys.getenv("BAYES_STEM_QUINTILE_REFERENCE", "saber11_cohort")
ROPE <- as.numeric(Sys.getenv("BAYES_STEM_ROPE", "0.01"))

OUT_DIR <- file.path(ROOT, "outputs", "bayesian_stem_saberpro_2022")
FIG_DIR <- file.path(OUT_DIR, "figures")
TABLE_DIR <- file.path(OUT_DIR, "tables")
MODEL_DIR <- file.path(OUT_DIR, "models")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)

message("Output folder: ", OUT_DIR)
if (QUICK_RUN) {
  message("BAYES_STEM_QUICK=1: reading at most ", QUICK_N_MAX, " rows per raw file.")
}

icfes_delimiter <- intToUtf8(172)
quintile_levels <- paste0("Q", 1:5)

subject_specs <- data.table(
  subject = c("math", "science", "language"),
  quintile_col = c(
    "PUNT_MATEMATICAS_QUINTILE",
    "PUNT_C_NATURALES_QUINTILE",
    "PUNT_LENGUAJE_QUINTILE"
  ),
  model_col = c("math_q", "science_q", "language_q"),
  label = c("Matematicas", "Ciencias naturales", "Lenguaje")
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

normalize_name <- function(x) {
  y <- enc2utf8(as.character(x))
  y <- iconv(y, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
  y <- toupper(gsub("[^A-Z0-9]+", "_", y))
  gsub("^_|_$", "", y)
}

normalize_label <- function(x) {
  y <- enc2utf8(as.character(x))
  y <- iconv(y, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
  y <- toupper(trimws(y))
  gsub("\\s+", " ", y)
}

find_col <- function(nms, normalized_target) {
  idx <- match(normalized_target, normalize_name(nms))
  if (is.na(idx)) {
    stop("Column not found: ", normalized_target)
  }
  nms[[idx]]
}

clean_code <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.0+$", "", x)
  x[x == ""] <- NA_character_
  x
}

clean_gender <- function(x) {
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[y %in% c("M", "MASCULINO", "HOMBRE", "MALE")] <- "Men"
  out[y %in% c("F", "FEMENINO", "MUJER", "FEMALE")] <- "Women"
  out
}

as_number <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = TRUE)))
}

detect_delimiter <- function(path) {
  line <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  candidates <- list(
    pipe = "|",
    notsign = icfes_delimiter,
    comma = ",",
    semicolon = ";",
    tab = "\t"
  )
  counts <- vapply(candidates, function(delim) {
    matches <- gregexpr(delim, line, fixed = TRUE)[[1]]
    if (identical(matches, -1L)) 0L else length(matches)
  }, integer(1))
  if (max(counts) == 0L) {
    stop("Could not detect delimiter for: ", path)
  }
  candidates[[which.max(counts)]]
}

read_delim_selected <- function(path, wanted_cols, n_max = Inf) {
  delim <- detect_delimiter(path)
  header <- readr::read_delim(
    path,
    delim = delim,
    n_max = 0,
    show_col_types = FALSE,
    locale = readr::locale(encoding = "UTF-8"),
    progress = FALSE
  )
  original_names <- names(header)
  upper_names <- toupper(original_names)
  wanted_upper <- toupper(wanted_cols)
  present_upper <- intersect(wanted_upper, upper_names)
  if (length(present_upper) == 0L) {
    stop("None of the requested columns were found in: ", path)
  }
  selected_original <- original_names[match(present_upper, upper_names)]

  if (is.finite(n_max)) {
    dt <- readr::read_delim(
      path,
      delim = delim,
      col_select = tidyselect::all_of(selected_original),
      col_types = readr::cols(.default = readr::col_character()),
      n_max = n_max,
      show_col_types = FALSE,
      locale = readr::locale(encoding = "UTF-8"),
      progress = FALSE,
      lazy = FALSE
    )
  } else {
    dt <- readr::read_delim(
      path,
      delim = delim,
      col_select = tidyselect::all_of(selected_original),
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE,
      locale = readr::locale(encoding = "UTF-8"),
      progress = FALSE,
      lazy = FALSE
    )
  }

  setDT(dt)
  setnames(dt, toupper(names(dt)))
  dt
}

make_quintile <- function(x) {
  x_num <- as_number(x)
  out <- rep(NA_character_, length(x_num))
  ok <- is.finite(x_num)
  if (sum(ok) == 0L) {
    return(factor(out, levels = quintile_levels, ordered = TRUE))
  }

  breaks <- unique(as.numeric(stats::quantile(
    x_num[ok],
    probs = seq(0, 1, 0.2),
    na.rm = TRUE,
    type = 7
  )))

  if (length(breaks) == 6L) {
    out[ok] <- as.character(cut(
      x_num[ok],
      breaks = breaks,
      labels = quintile_levels,
      include.lowest = TRUE,
      right = TRUE
    ))
  } else {
    ranks <- data.table::frank(x_num, ties.method = "average", na.last = "keep")
    r_breaks <- unique(as.numeric(stats::quantile(
      ranks[ok],
      probs = seq(0, 1, 0.2),
      na.rm = TRUE,
      type = 7
    )))
    out[ok] <- as.character(cut(
      ranks[ok],
      breaks = r_breaks,
      labels = quintile_levels[seq_len(length(r_breaks) - 1L)],
      include.lowest = TRUE,
      right = TRUE
    ))
  }
  factor(out, levels = quintile_levels, ordered = TRUE)
}

add_quintiles <- function(dt, score_cols) {
  for (score_col in score_cols) {
    if (score_col %in% names(dt)) {
      dt[, (score_col) := as_number(get(score_col))]
      dt[, paste0(score_col, "_QUINTILE") := make_quintile(get(score_col))]
    }
  }
  dt
}

read_crosswalk <- function() {
  key_2019 <- fread(
    file.path(ROOT, "cruces-SB11-SBPRO", "llaveSB112006_SBPro2019.txt"),
    sep = ",",
    colClasses = "character"
  )
  key_2023 <- fread(
    file.path(ROOT, "cruces-SB11-SBPRO", "Llave_Saber11_SaberPro.txt"),
    sep = ";",
    colClasses = "character"
  )

  setnames(key_2019, "estu_consecutivo_11", "estu_consecutivo_sb11", skip_absent = TRUE)
  setnames(key_2023, "estu_consecutivo_11", "estu_consecutivo_sb11", skip_absent = TRUE)
  setnames(key_2023, "estu_consecutivo_sbpro", "estu_consecutivo_PRO", skip_absent = TRUE)

  keep <- c("estu_consecutivo_sb11", "estu_consecutivo_PRO")
  key <- rbindlist(
    list(key_2019[, ..keep], key_2023[, ..keep]),
    use.names = TRUE,
    fill = TRUE
  )
  key <- key[!is.na(estu_consecutivo_sb11) & !is.na(estu_consecutivo_PRO)]
  key <- unique(key, by = c("estu_consecutivo_sb11", "estu_consecutivo_PRO"))

  key[, YEAR_SB11 := sub("^(?:SABER11|SB11)([0-9]{4}).*", "\\1", estu_consecutivo_sb11)]
  key[, SB11_PERIOD := sub("^(?:SABER11|SB11)([0-9]{5}).*", "\\1", estu_consecutivo_sb11)]
  key[, YEAR_SBPRO := sub("^EK([0-9]{4}).*", "\\1", estu_consecutivo_PRO)]
  key[YEAR_SBPRO == "2022"]
}

read_saberpro_2022 <- function(key_analysis) {
  wanted <- c(
    "PERIODO",
    "ESTU_CONSECUTIVO",
    "ESTU_GENERO",
    "ESTU_SNIES_PRGMACADEMICO",
    "ESTU_PRGM_ACADEMICO",
    "GRUPOREFERENCIA",
    "ESTU_NUCLEO_PREGRADO",
    "ESTU_NIVEL_PRGM_ACADEMICO"
  )
  n_max <- if (QUICK_RUN) QUICK_N_MAX else Inf
  files <- file.path(
    ROOT,
    "dataSaberPro",
    c("SaberPro_Genericas_20221.TXT", "SaberPro_Genericas_20222.TXT")
  )
  sbpro <- rbindlist(lapply(files, read_delim_selected, wanted_cols = wanted, n_max = n_max), fill = TRUE)
  setnames(sbpro, c("PERIODO", "ESTU_GENERO"), c("PERIODO_SBPRO", "ESTU_GENERO_SBPRO"), skip_absent = TRUE)
  sbpro[, ESTU_SNIES_PRGMACADEMICO := clean_code(ESTU_SNIES_PRGMACADEMICO)]

  sbpro <- merge(
    sbpro,
    key_analysis[, .(estu_consecutivo_PRO, estu_consecutivo_sb11, SB11_PERIOD, YEAR_SB11)],
    by.x = "ESTU_CONSECUTIVO",
    by.y = "estu_consecutivo_PRO",
    all = FALSE
  )

  snies <- read.xlsx(file.path(ROOT, "dataSaberPro", "programas_SNIES.xlsx"), sheet = "Programas")
  setDT(snies)
  code_col <- find_col(names(snies), "CODIGO_SNIES_DEL_PROGRAMA")
  area_col <- find_col(names(snies), "AREA_DE_CONOCIMIENTO")
  cine_broad_col <- find_col(names(snies), "CINE_F_2013_AC_CAMPO_AMPLIO")
  cine_specific_col <- find_col(names(snies), "CINE_F_2013_AC_CAMPO_ESPECIFIC")
  nucleo_col <- find_col(names(snies), "NUCLEO_BASICO_DEL_CONOCIMIENTO")

  snies_small <- snies[, .(
    snies_code = clean_code(get(code_col)),
    area_conocimiento = as.character(get(area_col)),
    cine_campo_amplio = as.character(get(cine_broad_col)),
    cine_campo_especifico = as.character(get(cine_specific_col)),
    nucleo_basico = as.character(get(nucleo_col))
  )]
  snies_small <- unique(snies_small[!is.na(snies_code)], by = "snies_code")

  sbpro <- merge(
    sbpro,
    snies_small,
    by.x = "ESTU_SNIES_PRGMACADEMICO",
    by.y = "snies_code",
    all.x = TRUE
  )

  stem_areas <- normalize_label(c(
    "Ingenieria, arquitectura, urbanismo y afines",
    "Matematicas y ciencias naturales",
    "Agronomia, veterinaria y afines",
    "Ciencias de la salud"
  ))
  sbpro[, area_key := normalize_label(area_conocimiento)]
  sbpro[, STEM_AREA := as.integer(area_key %in% stem_areas)]
  sbpro[is.na(area_conocimiento) | area_key == "SIN INFORMACION", STEM_AREA := NA_integer_]
  sbpro[]
}

read_saber11_panel <- function(key_analysis) {
  specs <- data.table(
    sb11_period = c("20142", "20151", "20152", "20161", "20162", "20171", "20172", "20181", "20182"),
    year = c(2014L, 2015L, 2015L, 2016L, 2016L, 2017L, 2017L, 2018L, 2018L),
    file = c(
      "SB11_20142.txt",
      "SB11_20151.txt",
      "SB11_20152.txt",
      "SB11_20161.txt",
      "SB11_20162.txt",
      "SB11_20171.TXT",
      "SB11_20172.TXT",
      "SB11_20181.TXT",
      "SB11_20182.TXT"
    )
  )
  specs <- specs[sb11_period %in% key_analysis$SB11_PERIOD]
  if (nrow(specs) == 0L) {
    stop("No Saber 11 files in the 2014-2 to 2018-2 panel match the 2022 Saber Pro key.")
  }

  wanted <- c(
    "PERIODO",
    "ESTU_CONSECUTIVO",
    "ESTU_GENERO",
    "ESTU_NSE_INDIVIDUAL",
    "COLE_CALENDARIO",
    "COLE_GENERO",
    "COLE_NATURALEZA",
    "COLE_BILINGUE",
    "COLE_JORNADA",
    "COLE_CARACTER",
    "PUNT_LECTURA_CRITICA",
    "PUNT_LENGUAJE",
    "PUNT_MATEMATICAS",
    "PUNT_C_NATURALES",
    "PUNT_SOCIALES_CIUDADANAS",
    "PUNT_INGLES",
    "PUNT_GLOBAL"
  )
  score_cols <- c(
    "PUNT_LENGUAJE",
    "PUNT_MATEMATICAS",
    "PUNT_C_NATURALES",
    "PUNT_SOCIALES_CIUDADANAS",
    "PUNT_INGLES",
    "PUNT_GLOBAL"
  )
  n_max <- if (QUICK_RUN) QUICK_N_MAX else Inf

  pieces <- vector("list", nrow(specs))
  for (i in seq_len(nrow(specs))) {
    message("Reading Saber 11 cohort ", specs$sb11_period[[i]], " ...")
    dt <- read_delim_selected(file.path(ROOT, "dataSaber11", specs$file[[i]]), wanted, n_max = n_max)

    if ("PUNT_LECTURA_CRITICA" %in% names(dt) && !"PUNT_LENGUAJE" %in% names(dt)) {
      setnames(dt, "PUNT_LECTURA_CRITICA", "PUNT_LENGUAJE")
    }

    dt[, SB11_PERIOD_SOURCE := specs$sb11_period[[i]]]
    dt[, YEAR_SB11_SOURCE := specs$year[[i]]]

    if (identical(QUINTILE_REFERENCE, "saber11_cohort")) {
      dt <- add_quintiles(dt, score_cols)
    }

    dt <- merge(
      dt,
      key_analysis[, .(estu_consecutivo_sb11, estu_consecutivo_PRO)],
      by.x = "ESTU_CONSECUTIVO",
      by.y = "estu_consecutivo_sb11",
      all = FALSE
    )

    if (identical(QUINTILE_REFERENCE, "linked_saberpro")) {
      dt <- add_quintiles(dt, score_cols)
    }

    setnames(dt, c("PERIODO", "ESTU_GENERO"), c("PERIODO_SB11", "ESTU_GENERO_SB11"), skip_absent = TRUE)
    pieces[[i]] <- dt
  }

  panel <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  panel[panel == ""] <- NA
  panel[]
}

prepare_analysis_data <- function(panel_sb11, sbpro_2022) {
  keep_sbpro <- c(
    "ESTU_CONSECUTIVO",
    "PERIODO_SBPRO",
    "ESTU_GENERO_SBPRO",
    "ESTU_SNIES_PRGMACADEMICO",
    "ESTU_PRGM_ACADEMICO",
    "GRUPOREFERENCIA",
    "ESTU_NUCLEO_PREGRADO",
    "ESTU_NIVEL_PRGM_ACADEMICO",
    "area_conocimiento",
    "cine_campo_amplio",
    "cine_campo_especifico",
    "nucleo_basico",
    "STEM_AREA"
  )
  dt <- merge(
    panel_sb11,
    sbpro_2022[, ..keep_sbpro],
    by.x = "estu_consecutivo_PRO",
    by.y = "ESTU_CONSECUTIVO",
    all = FALSE
  )

  dt[, gender := factor(clean_gender(ESTU_GENERO_SB11), levels = c("Men", "Women"))]
  dt[, nse_model := toupper(trimws(as.character(ESTU_NSE_INDIVIDUAL)))]
  dt[is.na(nse_model) | nse_model == "", nse_model := "NO INFO"]
  dt[, nse_model := factor(nse_model)]
  dt[, sb11_period := factor(SB11_PERIOD_SOURCE)]

  for (j in seq_len(nrow(subject_specs))) {
    q_col <- subject_specs$quintile_col[[j]]
    m_col <- subject_specs$model_col[[j]]
    dt[, (m_col) := factor(as.character(get(q_col)), levels = quintile_levels, ordered = FALSE)]
  }

  dt <- dt[
    !is.na(STEM_AREA) &
      !is.na(gender) &
      !is.na(math_q) &
      !is.na(science_q) &
      !is.na(language_q)
  ]
  dt[]
}

posterior_interval <- function(draws) {
  c(
    median = stats::median(draws),
    lower = stats::quantile(draws, 0.025, names = FALSE),
    upper = stats::quantile(draws, 0.975, names = FALSE)
  )
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

  rows <- lapply(seq_len(nrow(long_counts)), function(i) {
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
  })
  summary <- rbindlist(rows)

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

  log1pexp <- function(z) {
    ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))
  }

  neg_logpost <- function(beta) {
    eta <- as.vector(X %*% beta)
    if (any(!is.finite(eta))) {
      return(Inf)
    }
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

  message("Fitting Bayesian logistic model with ", ncol(X), " coefficients and ", nrow(cells), " aggregated cells ...")
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
    modal_nse = names(which.max(table(dt$nse_model))),
    modal_period = names(which.max(table(dt$sb11_period))),
    opt = opt
  )
}

newdata_for_subject <- function(model_dt, subject_id) {
  row <- subject_specs[subject == subject_id]
  if (nrow(row) != 1L) {
    stop("Unknown subject: ", subject_id)
  }
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

order_subject_facets <- function(dt) {
  dt <- copy(dt)
  dt[, subject_label := factor(subject_label, levels = subject_specs$label)]
  dt
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
      title = "Observed STEM selection by Saber 11 performance quintile",
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
      title = "Bayesian predicted probability of STEM selection",
      subtitle = "Other Saber 11 performance areas fixed at Q3; NSE and cohort fixed at their modal values",
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
      title = "Posterior gender difference at the same Saber 11 quintile",
      subtitle = "Women minus men; values below zero indicate lower STEM probability for women",
      x = "Saber 11 quintile",
      y = "Probability difference"
    ) +
    theme_publication()
}

plot_probability_women_lower <- function(contrast_summary) {
  contrast_summary <- order_subject_facets(contrast_summary)
  ggplot(contrast_summary, aes(x = quintile, y = pr_women_lower)) +
    geom_hline(yintercept = 0.5, color = "grey35", linewidth = 0.35) +
    geom_col(fill = "#587C6D", width = 0.72) +
    facet_wrap(~ subject_label, nrow = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      title = "Posterior probability that women have lower STEM probability than men",
      x = "Saber 11 quintile",
      y = "Pr(Women < Men)"
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
      title = "One-quintile offset contrast",
      subtitle = "Women one quintile higher minus men one quintile lower",
      x = NULL,
      y = "Probability difference"
    ) +
    coord_flip() +
    theme_publication()
}

message("Reading crosswalk ...")
key_analysis <- read_crosswalk()
fwrite(key_analysis[, .N, by = .(YEAR_SB11, SB11_PERIOD)][order(YEAR_SB11, SB11_PERIOD)],
       file.path(TABLE_DIR, "00_key_saberpro2022_by_saber11_period.csv"))

message("Reading Saber Pro 2022 ...")
sbpro_2022 <- read_saberpro_2022(key_analysis)
fwrite(sbpro_2022[, .N, by = .(area_conocimiento, STEM_AREA)][order(-N)],
       file.path(TABLE_DIR, "01_saberpro2022_area_counts.csv"))

message("Reading Saber 11 panel and creating quintiles ...")
panel_sb11 <- read_saber11_panel(key_analysis)

message("Preparing linked analytic data ...")
analysis_dt <- prepare_analysis_data(panel_sb11, sbpro_2022)
fwrite(analysis_dt[, .(
  estu_consecutivo_sb11 = ESTU_CONSECUTIVO,
  estu_consecutivo_PRO,
  gender,
  nse_model,
  sb11_period,
  math_q,
  science_q,
  language_q,
  STEM_AREA,
  area_conocimiento,
  ESTU_PRGM_ACADEMICO
)], file.path(TABLE_DIR, "02_analysis_data_minimal.csv"))

sample_summary <- analysis_dt[, .(
  n = .N,
  stem = sum(STEM_AREA),
  stem_rate = mean(STEM_AREA)
), by = .(gender)]
fwrite(sample_summary, file.path(TABLE_DIR, "03_sample_summary_by_gender.csv"))

message("Analytic sample size: ", nrow(analysis_dt))
message("STEM rate: ", percent(mean(analysis_dt$STEM_AREA), accuracy = 0.1))

message("Summarising beta-binomial descriptive posteriors ...")
beta_results <- summarise_beta_binomial(analysis_dt)
fwrite(beta_results$summary, file.path(TABLE_DIR, "04_observed_beta_binomial_stem_rates.csv"))
fwrite(beta_results$contrasts, file.path(TABLE_DIR, "05_observed_beta_binomial_gender_contrasts.csv"))

p_observed <- plot_observed(beta_results$summary)
save_plot(p_observed, "01_observed_stem_rate_by_quintile_gender", width = 9, height = 4.8)

message("Fitting Bayesian logistic model ...")
bayes_fit <- fit_bayesian_logistic_laplace(analysis_dt)
fwrite(bayes_fit$coef_summary, file.path(TABLE_DIR, "06_bayesian_logistic_coefficients.csv"))
saveRDS(bayes_fit, file.path(MODEL_DIR, "bayesian_logistic_laplace_fit.rds"))

message("Generating posterior predictions and contrasts ...")
pred_results <- summarise_predictions(bayes_fit, analysis_dt)
fwrite(pred_results$summary, file.path(TABLE_DIR, "07_bayesian_predicted_probabilities.csv"))
fwrite(pred_results$contrasts, file.path(TABLE_DIR, "08_bayesian_gender_contrasts_same_quintile.csv"))
fwrite(pred_results$offsets, file.path(TABLE_DIR, "09_bayesian_one_quintile_offset_contrasts.csv"))

p_pred <- plot_predictions(pred_results$summary)
p_diff <- plot_differences(pred_results$contrasts)
p_lower <- plot_probability_women_lower(pred_results$contrasts)
p_offset <- plot_offset_contrasts(pred_results$offsets)

save_plot(p_pred, "02_bayesian_predicted_probability_by_quintile_gender", width = 9, height = 4.8)
save_plot(p_diff, "03_bayesian_gender_difference_same_quintile", width = 9, height = 4.8)
save_plot(p_lower, "04_probability_women_lower_than_men_same_quintile", width = 9, height = 4.8)
save_plot(p_offset, "05_one_quintile_offset_contrast", width = 9.5, height = 5.2)

run_info <- data.table(
  item = c(
    "quick_run",
    "posterior_draws",
    "quintile_reference",
    "rope_probability_points",
    "analytic_n",
    "stem_rate",
    "output_folder"
  ),
  value = c(
    as.character(QUICK_RUN),
    as.character(POSTERIOR_DRAWS),
    QUINTILE_REFERENCE,
    as.character(ROPE),
    as.character(nrow(analysis_dt)),
    as.character(mean(analysis_dt$STEM_AREA)),
    OUT_DIR
  )
)
fwrite(run_info, file.path(TABLE_DIR, "10_run_info.csv"))

message("Done. Key outputs:")
message(" - ", file.path(TABLE_DIR, "08_bayesian_gender_contrasts_same_quintile.csv"))
message(" - ", file.path(FIG_DIR, "02_bayesian_predicted_probability_by_quintile_gender.png"))
message(" - ", file.path(FIG_DIR, "03_bayesian_gender_difference_same_quintile.png"))
