from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "data" / "derived"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)


def design_df(base):
    psus = base[["SDMVSTRA", "SDMVPSU"]].drop_duplicates().shape[0]
    strata = base["SDMVSTRA"].nunique()
    return max(1, psus - strata)


def domain_mean_ci(base, domain, variable):
    """Taylor-linearized ratio estimate for a survey domain.

    The full cycle-specific analysis sample is retained. Observations outside
    the domain contribute a zero linearized value, preserving strata and PSU
    structure instead of treating the domain as a newly sampled dataset.
    """
    if isinstance(variable, str):
        values = base[variable]
    else:
        values = pd.Series(variable, index=base.index)
    d = base[["weight", "SDMVSTRA", "SDMVPSU"]].copy()
    d["value"] = values
    dom = pd.Series(domain, index=base.index).fillna(False).astype(bool)
    valid = dom & d["value"].notna() & d["weight"].gt(0)
    denominator = d.loc[valid, "weight"].sum()
    estimate = np.sum(d.loc[valid, "weight"] * d.loc[valid, "value"]) / denominator
    d["linearized"] = 0.0
    d.loc[valid, "linearized"] = d.loc[valid, "weight"] * (d.loc[valid, "value"] - estimate)
    clusters = d.groupby(["SDMVSTRA", "SDMVPSU"], as_index=False)["linearized"].sum()
    variance_total = 0.0
    single_psu_strata = 0
    for _, stratum in clusters.groupby("SDMVSTRA"):
        n_psu = len(stratum)
        if n_psu < 2:
            single_psu_strata += 1
            continue
        variance_total += n_psu / (n_psu - 1) * np.square(
            stratum["linearized"] - stratum["linearized"].mean()
        ).sum()
    standard_error = np.sqrt(variance_total) / denominator
    critical = t.ppf(0.975, design_df(base))
    return {
        "estimate": estimate,
        "lower": estimate - critical * standard_error,
        "upper": estimate + critical * standard_error,
        "se": standard_error,
        "design_df": design_df(base),
        "single_psu_strata": single_psu_strata,
        "domain_n": int(valid.sum()),
    }


def weighted_quantile(values, weights, probabilities):
    values = np.asarray(values)
    weights = np.asarray(weights)
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[keep], weights[keep]
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cumulative = np.cumsum(weights) - 0.5 * weights
    cumulative /= weights.sum()
    return np.interp(probabilities, cumulative, values)


def prepare():
    d = pd.read_csv(DERIVED / "core_analytic_cohort.csv")
    d = d[d["adult"] & d["valid_vcte"] & d["complete_fib4"]].copy()
    d["masld_domain"] = d["masld_formal"].astype(bool)
    d["female"] = (d["RIAGENDR"] == 2).astype(int)
    race_labels = {
        1: "Mexican American",
        2: "Other Hispanic",
        3: "Non-Hispanic White",
        4: "Non-Hispanic Black",
        6: "Non-Hispanic Asian",
        7: "Other or multiracial",
    }
    d["race_ethnicity"] = d["RIDRETH3"].map(race_labels)
    for code, label in race_labels.items():
        d[f"race_{code}"] = (d["RIDRETH3"] == code).astype(int)
    d["obesity"] = (d["BMXBMI"] >= 30).astype(int)
    diabetes_treatment = d["DIQ050"].eq(1) | d["DIQ070"].eq(1)
    d["diabetes"] = (
        d["DIQ010"].eq(1) | diabetes_treatment | d["LBXGH"].ge(6.5)
    ).astype(int)
    d["hypertension_cmrf"] = d["cmrf_bp"].astype(int)
    d["dyslipidemia_cmrf"] = (d["cmrf_tg"].astype(bool) | d["cmrf_hdl"].astype(bool)).astype(int)
    additional_cmrf = (
        d["cmrf_glucose"].astype(int)
        + d["cmrf_bp"].astype(int)
        + d["cmrf_tg"].astype(int)
        + d["cmrf_hdl"].astype(int)
    )
    d["metabolic_cluster"] = (
        d["diabetes"].eq(1) | ((d["BMXBMI"] >= 30) & (additional_cmrf >= 2))
    ).astype(int)
    return d


def estimate_row(base, domain, variable, scale=100):
    result = domain_mean_ci(base, domain, variable)
    return {
        "estimate": scale * result["estimate"],
        "ci_lower": scale * result["lower"],
        "ci_upper": scale * result["upper"],
        "domain_n": result["domain_n"],
        "design_df": result["design_df"],
        "single_psu_strata": result["single_psu_strata"],
    }


def baseline_table(d):
    rows = []
    variables = [
        ("Age, years", "RIDAGEYR", "mean"),
        ("Female, %", "female", "percent"),
        ("Mexican American, %", "race_1", "percent"),
        ("Other Hispanic, %", "race_2", "percent"),
        ("Non-Hispanic White, %", "race_3", "percent"),
        ("Non-Hispanic Black, %", "race_4", "percent"),
        ("Non-Hispanic Asian, %", "race_6", "percent"),
        ("Other or multiracial, %", "race_7", "percent"),
        ("BMI, kg/m²", "BMXBMI", "mean"),
        ("Obesity, %", "obesity", "percent"),
        ("Diabetes, %", "diabetes", "percent"),
        ("Hypertension criterion, %", "hypertension_cmrf", "percent"),
        ("Hypertriglyceridemia criterion, %", "cmrf_tg", "percent"),
        ("Low HDL-C criterion, %", "cmrf_hdl", "percent"),
        ("Glucose criterion, %", "cmrf_glucose", "percent"),
        ("Low FIB-4, %", "fib4_low_age", "percent"),
        ("CAP, dB/m", "CAP", "mean"),
        ("LSM, kPa", "LSM", "quantile"),
        ("FIB-4", "fib4", "quantile"),
    ]
    for cycle, base in d.groupby("cycle"):
        domain = base["masld_domain"]
        for label, variable, kind in variables:
            if kind == "quantile":
                z = base.loc[domain & base[variable].notna(), [variable, "weight"]]
                q1, median, q3 = weighted_quantile(z[variable], z["weight"], [0.25, 0.5, 0.75])
                rows.append({"cycle": cycle, "characteristic": label, "kind": kind,
                             "estimate": median, "ci_lower": np.nan, "ci_upper": np.nan,
                             "q1": q1, "q3": q3, "unweighted_n": len(z)})
            else:
                scale = 100 if kind == "percent" else 1
                result = estimate_row(base, domain, variable, scale=scale)
                rows.append({"cycle": cycle, "characteristic": label, "kind": kind,
                             "estimate": result["estimate"], "ci_lower": result["ci_lower"],
                             "ci_upper": result["ci_upper"], "q1": np.nan, "q3": np.nan,
                             "unweighted_n": result["domain_n"]})
    return pd.DataFrame(rows)


def main_table(d):
    rows = []
    for cycle, base in d.groupby("cycle"):
        masld = base["masld_domain"]
        for cutoff in (8.0, 8.6, 10.0):
            elevated = (base["LSM"] >= cutoff).astype(int)
            low = base["fib4_low_age"].astype(bool)
            discordant = (low & elevated.eq(1)).astype(int)
            a = estimate_row(base, masld, elevated)
            b = estimate_row(base, masld, discordant)
            c = estimate_row(base, masld & low, elevated)
            d_rev = estimate_row(base, masld & elevated.eq(1), low.astype(int))
            rows.append({
                "cycle": cycle, "lsm_cutoff": cutoff,
                "masld_n": int(masld.sum()), "low_fib4_n": int((masld & low).sum()),
                "lsm_elevated_n": int((masld & elevated.eq(1)).sum()),
                "discordant_n": int((masld & low & elevated.eq(1)).sum()),
                "lsm_elevated_percent": a["estimate"], "lsm_elevated_ci_lower": a["ci_lower"], "lsm_elevated_ci_upper": a["ci_upper"],
                "discordant_in_masld_percent": b["estimate"], "discordant_in_masld_ci_lower": b["ci_lower"], "discordant_in_masld_ci_upper": b["ci_upper"],
                "lsm_elevated_among_low_fib4_percent": c["estimate"], "lsm_elevated_among_low_fib4_ci_lower": c["ci_lower"], "lsm_elevated_among_low_fib4_ci_upper": c["ci_upper"],
                "low_fib4_among_lsm_percent": d_rev["estimate"], "low_fib4_among_lsm_ci_lower": d_rev["ci_lower"], "low_fib4_among_lsm_ci_upper": d_rev["ci_upper"],
                "design_df": a["design_df"], "single_psu_strata": a["single_psu_strata"],
            })
    return pd.DataFrame(rows)


def subgroup_table(d):
    rows = []
    for cycle, base in d.groupby("cycle"):
        masld = base["masld_domain"]
        low = base["fib4_low_age"].astype(bool)
        domains = {
            "All operational MASLD": masld,
            "Diabetes": masld & base["diabetes"].eq(1),
            "Exploratory metabolic-risk clustering": masld & base["metabolic_cluster"].eq(1),
        }
        for cutoff in (8.0, 8.6, 10.0):
            elevated = (base["LSM"] >= cutoff).astype(int)
            for label, domain in domains.items():
                overall = estimate_row(base, domain, (low & elevated.eq(1)).astype(int))
                among_low = estimate_row(base, domain & low, elevated)
                rows.append({
                    "cycle": cycle, "lsm_cutoff": cutoff, "subgroup": label,
                    "subgroup_n": int(domain.sum()), "low_fib4_n": int((domain & low).sum()),
                    "discordant_n": int((domain & low & elevated.eq(1)).sum()),
                    "discordant_in_subgroup_percent": overall["estimate"],
                    "discordant_in_subgroup_ci_lower": overall["ci_lower"],
                    "discordant_in_subgroup_ci_upper": overall["ci_upper"],
                    "lsm_elevated_among_low_fib4_percent": among_low["estimate"],
                    "lsm_elevated_among_low_fib4_ci_lower": among_low["ci_lower"],
                    "lsm_elevated_among_low_fib4_ci_upper": among_low["ci_upper"],
                })
    return pd.DataFrame(rows)


def age_rule_table(d):
    rows = []
    for cycle, base in d.groupby("cycle"):
        masld = base["masld_domain"]
        elevated = (base["LSM"] >= 8.0).astype(int)
        age_groups = {"≤65 years": base["RIDAGEYR"] <= 65, ">65 years": base["RIDAGEYR"] > 65}
        rules = {"Age-adjusted": base["fib4_low_age"].astype(bool),
                 "Fixed FIB-4 <1.3": base["fib4_low_fixed"].astype(bool)}
        for age_label, age_domain in age_groups.items():
            for rule_label, low in rules.items():
                domain = masld & age_domain
                low_prevalence = estimate_row(base, domain, low.astype(int))
                among_low = estimate_row(base, domain & low, elevated)
                rows.append({
                    "cycle": cycle, "age_group": age_label, "fib4_rule": rule_label,
                    "subgroup_n": int(domain.sum()), "low_fib4_n": int((domain & low).sum()),
                    "discordant_n": int((domain & low & elevated.eq(1)).sum()),
                    "low_fib4_percent": low_prevalence["estimate"],
                    "low_fib4_ci_lower": low_prevalence["ci_lower"],
                    "low_fib4_ci_upper": low_prevalence["ci_upper"],
                    "lsm_elevated_among_low_fib4_percent": among_low["estimate"],
                    "lsm_elevated_among_low_fib4_ci_lower": among_low["ci_lower"],
                    "lsm_elevated_among_low_fib4_ci_upper": among_low["ci_upper"],
                })
    return pd.DataFrame(rows)


def etiology_complete_sensitivity_table(d):
    """Repeat the main conditional estimate where alcohol/HBV/HCV are observed."""
    rows = []
    for cycle, base in d.groupby("cycle"):
        masld = base["masld_domain"]
        known = base["etiology_known"].astype(bool)
        low = base["fib4_low_age"].astype(bool)
        missingness = estimate_row(base, masld, (~known).astype(int))
        for cutoff in (8.0, 8.6, 10.0):
            elevated = (base["LSM"] >= cutoff).astype(int)
            domain = masld & known & low
            estimate = estimate_row(base, domain, elevated)
            rows.append({
                "cycle": cycle,
                "lsm_cutoff": cutoff,
                "masld_n": int(masld.sum()),
                "etiology_complete_masld_n": int((masld & known).sum()),
                "etiology_complete_low_fib4_n": int(domain.sum()),
                "discordant_n": int((domain & elevated.eq(1)).sum()),
                "etiology_incomplete_percent": missingness["estimate"],
                "etiology_incomplete_ci_lower": missingness["ci_lower"],
                "etiology_incomplete_ci_upper": missingness["ci_upper"],
                "lsm_elevated_among_low_fib4_percent": estimate["estimate"],
                "lsm_elevated_among_low_fib4_ci_lower": estimate["ci_lower"],
                "lsm_elevated_among_low_fib4_ci_upper": estimate["ci_upper"],
            })
    return pd.DataFrame(rows)


if __name__ == "__main__":
    data = prepare()
    outputs = {
        "revision_table1_full_baseline.csv": baseline_table(data),
        "revision_table2_bidirectional.csv": main_table(data),
        "revision_table3_high_risk.csv": subgroup_table(data),
        "revision_table5_age_rule.csv": age_rule_table(data),
        "revision_table7_etiology_complete_sensitivity.csv": etiology_complete_sensitivity_table(data),
    }
    for filename, frame in outputs.items():
        frame.to_csv(TABLES / filename, index=False, encoding="utf-8-sig")
        print(filename, len(frame))
