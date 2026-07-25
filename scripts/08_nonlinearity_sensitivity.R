project_dir <- normalizePath(".", winslash = "/")
suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(splines))

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

continuous <- c("age_10", "bmi_5", "tg_50", "hdl_10", "sbp_10", "cap_10")
all_covariates <- c("age_10", "female", "bmi_5", "diabetes", "tg_50", "hdl_10", "sbp_10", "cap_10")
rows <- list()

for (cycle_name in c("2017-2020", "2021-2023")) {
  x <- d[d$cycle == cycle_name & d$adult & !is.na(d$weight) & d$weight > 0, ]
  complete_covariates <- complete.cases(x[, c("lsm_ge_8", all_covariates)])
  x$model_domain <- x$valid_vcte & x$complete_fib4 & x$masld_formal &
    x$fib4_low_age & complete_covariates
  design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight,
                      data = x, nest = TRUE)
  model_design <- subset(design, model_domain)
  linear_formula <- as.formula(paste("lsm_ge_8 ~", paste(all_covariates, collapse = " + ")))
  linear_fit <- svyglm(linear_formula, design = model_design, family = quasibinomial())

  for (variable in continuous) {
    other_terms <- setdiff(all_covariates, variable)
    spline_formula <- as.formula(paste(
      "lsm_ge_8 ~ ns(", variable, ", df = 3) +", paste(other_terms, collapse = " + ")
    ))
    spline_fit <- svyglm(spline_formula, design = model_design, family = quasibinomial())
    comparison <- anova(linear_fit, spline_fit, method = "Wald")
    p_value <- comparison$p
    if (is.null(p_value)) p_value <- comparison[["p"]]
    rows[[length(rows) + 1]] <- data.frame(
      cycle = cycle_name,
      variable = variable,
      model_n = sum(x$model_domain),
      event_n = sum(x$lsm_ge_8[x$model_domain]),
      design_df = degf(design),
      spline_df = 3,
      wald_p_for_nonlinearity = as.numeric(p_value)
    )
  }
}

out <- do.call(rbind, rows)
write.csv(out, file.path(project_dir, "results", "tables", "revision_table8_nonlinearity_sensitivity.csv"), row.names = FALSE)
print(out, row.names = FALSE)
