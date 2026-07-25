project_dir <- normalizePath(".", winslash = "/")
suppressPackageStartupMessages(library(survey))

d <- read.csv(file.path(project_dir, "data", "derived", "core_analytic_cohort.csv"), check.names = FALSE)
logical_columns <- c("adult", "valid_vcte", "complete_fib4", "masld_formal", "fib4_low_age")
for (column in logical_columns) d[[column]] <- d[[column]] == "True"
d$female <- as.numeric(d$RIAGENDR == 2)
d$diabetes <- as.numeric(
  d$DIQ010 %in% 1 | d$DIQ050 %in% 1 | d$DIQ070 %in% 1 |
    (!is.na(d$LBXGH) & d$LBXGH >= 6.5)
)
d$lsm_ge_8 <- as.numeric(d$LSM >= 8.0)
d$age_10 <- d$RIDAGEYR / 10
d$bmi_5 <- d$BMXBMI / 5
d$tg_50 <- d$LBXSTR / 50
d$hdl_10 <- d$LBDHDD / 10
d$sbp_10 <- d$mean_sbp / 10
d$cap_10 <- d$CAP / 10
d$race_ethnicity <- factor(
  d$RIDRETH3,
  levels = c(3, 1, 2, 4, 6, 7),
  labels = c("Non-Hispanic White", "Mexican American", "Other Hispanic",
             "Non-Hispanic Black", "Non-Hispanic Asian", "Other or multiracial")
)

labels <- c(
  age_10 = "Age, per 10 years", female = "Female vs male",
  bmi_5 = "BMI, per 5 kg/m2", diabetes = "Diabetes",
  tg_50 = "Triglycerides, per 50 mg/dL",
  hdl_10 = "HDL-C, per 10 mg/dL", sbp_10 = "Systolic BP, per 10 mmHg",
  cap_10 = "CAP, per 10 dB/m"
)
covariates <- names(labels)
model_formula <- as.formula(paste("lsm_ge_8 ~", paste(covariates, collapse = " + ")))
race_adjusted_formula <- as.formula(paste(
  "lsm_ge_8 ~", paste(c(covariates, "race_ethnicity"), collapse = " + ")
))
rows <- list()
race_adjusted_rows <- list()
diagnostics <- list()

for (cycle_name in c("2017-2020", "2021-2023")) {
  x <- d[d$cycle == cycle_name & d$adult & !is.na(d$weight) & d$weight > 0, ]
  eligible_domain <- x$valid_vcte & x$complete_fib4 & x$masld_formal & x$fib4_low_age
  model_variables <- c("lsm_ge_8", covariates, "race_ethnicity")
  complete_covariates <- complete.cases(x[, model_variables])
  x$model_domain <- eligible_domain & complete_covariates
  design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight,
                      data = x, nest = TRUE)
  model_design <- subset(design, model_domain)
  fit <- svyglm(model_formula, design = model_design, family = quasibinomial())
  race_adjusted_fit <- svyglm(race_adjusted_formula, design = model_design, family = quasibinomial())
  design_df <- degf(design)
  critical <- qt(0.975, df = design_df)
  b <- coef(fit)[covariates]
  se <- sqrt(diag(vcov(fit)))[covariates]
  for (variable in covariates) {
    rows[[length(rows) + 1]] <- data.frame(
      cycle = cycle_name, variable = variable, label = labels[[variable]],
      adjusted_or = exp(b[[variable]]),
      ci_lower = exp(b[[variable]] - critical * se[[variable]]),
      ci_upper = exp(b[[variable]] + critical * se[[variable]]),
      p_value = 2 * pt(abs(b[[variable]] / se[[variable]]), df = design_df, lower.tail = FALSE),
      design_df = design_df
    )
    race_b <- coef(race_adjusted_fit)[[variable]]
    race_se <- sqrt(diag(vcov(race_adjusted_fit)))[[variable]]
    race_adjusted_rows[[length(race_adjusted_rows) + 1]] <- data.frame(
      cycle = cycle_name, variable = variable, label = labels[[variable]],
      adjusted_or = exp(race_b),
      ci_lower = exp(race_b - critical * race_se),
      ci_upper = exp(race_b + critical * race_se),
      p_value = 2 * pt(abs(race_b / race_se), df = design_df, lower.tail = FALSE),
      design_df = design_df
    )
  }
  model_data <- x[x$model_domain, covariates]
  correlation <- cor(model_data, use = "pairwise.complete.obs")
  diag(correlation) <- NA
  diagnostics[[length(diagnostics) + 1]] <- data.frame(
    cycle = cycle_name, eligible_low_fib4_n = sum(eligible_domain),
    model_n = sum(x$model_domain),
    excluded_for_missing_covariates_n = sum(eligible_domain & !complete_covariates),
    excluded_for_missing_covariates_percent =
      100 * sum(eligible_domain & !complete_covariates) / sum(eligible_domain),
    event_n = sum(x$lsm_ge_8[x$model_domain]), design_df = design_df,
    max_absolute_predictor_correlation = max(abs(correlation), na.rm = TRUE),
    convergence = fit$converged
  )
}

write.csv(do.call(rbind, rows), file.path(project_dir, "results", "tables", "revision_table6_logistic_primary.csv"), row.names = FALSE)
write.csv(do.call(rbind, race_adjusted_rows), file.path(project_dir, "results", "tables", "revision_table9_race_adjusted_sensitivity.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostics), file.path(project_dir, "results", "tables", "revision_table6_logistic_primary_diagnostics.csv"), row.names = FALSE)
print(do.call(rbind, diagnostics), row.names = FALSE)
print(do.call(rbind, rows), row.names = FALSE)
