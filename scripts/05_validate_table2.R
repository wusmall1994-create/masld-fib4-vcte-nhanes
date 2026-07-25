project_dir <- normalizePath(".", winslash = "/")
suppressPackageStartupMessages(library(survey))

d <- read.csv(file.path(project_dir, "data", "derived", "core_analytic_cohort.csv"), check.names = FALSE)
logical_columns <- c("adult", "valid_vcte", "complete_fib4", "masld_formal",
                     "fib4_low_age")
for (column in logical_columns) {
  if (column %in% names(d)) d[[column]] <- d[[column]] == "True"
}
d <- d[d$adult & d$valid_vcte & d$complete_fib4, ]
d$lsm_ge_8 <- d$LSM >= 8.0
d$lsm_ge_8_6 <- d$LSM >= 8.6
d$lsm_ge_10 <- d$LSM >= 10.0

estimate_domain <- function(design, domain_expr, variable, df_design) {
  dom <- subset(design, eval(domain_expr))
  fit <- svymean(as.formula(paste0("~", variable)), dom, na.rm = TRUE)
  estimate <- as.numeric(coef(fit)[1])
  standard_error <- as.numeric(SE(fit)[1])
  critical <- qt(0.975, df = df_design)
  c(estimate = estimate, standard_error = standard_error,
    ci_lower = estimate - critical * standard_error,
    ci_upper = estimate + critical * standard_error)
}

rows <- list()
index <- 1
for (cycle_name in c("2017-2020", "2021-2023")) {
  x <- d[d$cycle == cycle_name, ]
  x$low_fib4 <- as.numeric(x$fib4_low_age)
  design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight,
                      data = x, nest = TRUE)
  design_df <- degf(design)

  for (cutoff in c(8.0, 8.6, 10.0)) {
    suffix <- if (cutoff == 8.0) "8" else if (cutoff == 8.6) "8_6" else "10"
    elevated_var <- paste0("lsm_ge_", suffix)
    x_elevated <- as.numeric(x[[elevated_var]])
    x$elevated <- x_elevated
    x$discordant <- as.numeric(x$fib4_low_age & x_elevated == 1)
    design <- update(design, elevated = x$elevated, discordant = x$discordant,
                     low_fib4 = x$low_fib4)

    definitions <- list(
      lsm_elevated = list(quote(masld_formal), "elevated"),
      discordant_in_masld = list(quote(masld_formal), "discordant"),
      lsm_elevated_among_low_fib4 = list(quote(masld_formal & fib4_low_age), "elevated"),
      low_fib4_among_lsm = list(substitute(masld_formal & VAR == 1,
                                           list(VAR = as.name(elevated_var))), "low_fib4")
    )

    for (metric_name in names(definitions)) {
      result <- estimate_domain(design, definitions[[metric_name]][[1]],
                                definitions[[metric_name]][[2]], design_df)
      rows[[index]] <- data.frame(
        cycle = cycle_name, lsm_cutoff = cutoff, metric = metric_name,
        design_df = design_df,
        estimate_percent = 100 * result["estimate"],
        standard_error_percent = 100 * result["standard_error"],
        ci_lower_percent = 100 * result["ci_lower"],
        ci_upper_percent = 100 * result["ci_upper"]
      )
      index <- index + 1
    }
  }
}

out <- do.call(rbind, rows)
write.csv(out, file.path(project_dir, "results", "tables", "r_survey_table2_validation.csv"), row.names = FALSE)
python_result <- read.csv(file.path(project_dir, "results", "tables", "revision_table2_bidirectional.csv"))
metric_stems <- c(
  lsm_elevated = "lsm_elevated",
  discordant_in_masld = "discordant_in_masld",
  lsm_elevated_among_low_fib4 = "lsm_elevated_among_low_fib4",
  low_fib4_among_lsm = "low_fib4_among_lsm"
)
differences <- c()
for (i in seq_len(nrow(out))) {
  p <- python_result[python_result$cycle == out$cycle[i] &
                     python_result$lsm_cutoff == out$lsm_cutoff[i], ]
  stem <- metric_stems[[out$metric[i]]]
  differences <- c(
    differences,
    abs(out$estimate_percent[i] - p[[paste0(stem, "_percent")]]),
    abs(out$ci_lower_percent[i] - p[[paste0(stem, "_ci_lower")]]),
    abs(out$ci_upper_percent[i] - p[[paste0(stem, "_ci_upper")]])
  )
}
max_difference <- max(differences)
cat(sprintf("Maximum absolute Python-R difference: %.3e percentage points\n", max_difference))
if (max_difference > 1e-8) stop("Independent R validation exceeded tolerance")
cat(sprintf("R %s; survey %s\n", getRversion(), packageVersion("survey")))
print(out, row.names = FALSE)
