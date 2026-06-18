# Install the R packages needed to reproduce the processed-data analysis.
# Packages are installed into a local project library, not into the system library.

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
  if (file.exists(file.path(wd, "README.md")) && dir.exists(file.path(wd, "scripts"))) {
    return(wd)
  }
  if (basename(wd) == "scripts") {
    return(normalizePath(file.path(wd, ".."), winslash = "/", mustWork = FALSE))
  }
  wd
}

root <- find_repo_root()

lib <- file.path(root, "r_libs")
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(normalizePath(lib, winslash = "/"), .libPaths()))

packages <- c("data.table", "ggplot2", "scales")
missing <- setdiff(packages, rownames(installed.packages()))

if (length(missing) > 0) {
  install.packages(missing, lib = normalizePath(lib, winslash = "/"), repos = "https://cloud.r-project.org")
}

message("Library path: ", normalizePath(lib, winslash = "/"))
message("All required packages are available.")
