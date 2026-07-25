from pathlib import Path

import numpy as np
import pandas as pd

from analysis import domain_mean_ci

ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "data" / "derived"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)


def weighted_mean(base, domain, variable):
    return domain_mean_ci(base, domain, variable)


def weighted_sd(base, domain, variable, mean):
    values = base[variable]
    valid = domain & values.notna() & base["weight"].gt(0)
    weights = base.loc[valid, "weight"]
    return np.sqrt(np.sum(weights * (values[valid] - mean) ** 2) / weights.sum())


def continuous_smd(mean_a, sd_a, mean_b, sd_b):
    pooled = np.sqrt((sd_a ** 2 + sd_b ** 2) / 2)
    return (mean_a - mean_b) / pooled if pooled > 0 else np.nan


def binary_smd(p_a, p_b):
    pooled = np.sqrt((p_a * (1 - p_a) + p_b * (1 - p_b)) / 2)
    return (p_a - p_b) / pooled if pooled > 0 else np.nan


def main():
    data = pd.read_csv(DERIVED / "core_analytic_cohort.csv")
    # Selection-bias comparisons target the full MEC-examined adult sample.
    # Use MEC weights here; phlebotomy weights are zero for nonrespondents and
    # therefore are not suitable for comparing complete and excluded adults.
    data["weight"] = data["selection_weight"]
    data = data[data["adult"] & data["weight"].gt(0) & data["SDMVSTRA"].notna() & data["SDMVPSU"].notna()].copy()
    data["complete_analysis_data"] = data["valid_vcte"] & data["complete_fib4"]
    data["female"] = (data["RIAGENDR"] == 2).astype(float)
    data["obesity"] = np.where(data["BMXBMI"].notna(), (data["BMXBMI"] >= 30).astype(float), np.nan)
    diabetes_treatment = data["DIQ050"].eq(1) | data["DIQ070"].eq(1)
    diabetes_observed = (
        data["DIQ010"].notna()
        | data["DIQ050"].notna()
        | data["DIQ070"].notna()
        | data["LBXGH"].notna()
    )
    data["diabetes"] = np.where(
        diabetes_observed,
        (data["DIQ010"].eq(1) | diabetes_treatment | data["LBXGH"].ge(6.5)).astype(float),
        np.nan,
    )
    variables = [
        ("Age, years", "RIDAGEYR", "continuous"),
        ("Female, %", "female", "binary"),
        ("BMI, kg/m²", "BMXBMI", "continuous"),
        ("Obesity, %", "obesity", "binary"),
        ("Diabetes, %", "diabetes", "binary"),
    ]
    rows = []
    summary = []
    for cycle, base in data.groupby("cycle"):
        included = base["complete_analysis_data"]
        excluded = ~included
        inclusion = weighted_mean(base, pd.Series(True, index=base.index), "complete_analysis_data")
        summary.append({
            "cycle": cycle,
            "adult_n": len(base),
            "complete_n": int(included.sum()),
            "excluded_n": int(excluded.sum()),
            "weighted_complete_percent": 100 * inclusion["estimate"],
            "weighted_complete_ci_lower": 100 * inclusion["lower"],
            "weighted_complete_ci_upper": 100 * inclusion["upper"],
        })
        for label, variable, kind in variables:
            a = weighted_mean(base, included, variable)
            b = weighted_mean(base, excluded, variable)
            if kind == "continuous":
                sd_a = weighted_sd(base, included, variable, a["estimate"])
                sd_b = weighted_sd(base, excluded, variable, b["estimate"])
                smd = continuous_smd(a["estimate"], sd_a, b["estimate"], sd_b)
                scale = 1
            else:
                smd = binary_smd(a["estimate"], b["estimate"])
                scale = 100
            rows.append({
                "cycle": cycle,
                "characteristic": label,
                "complete_estimate": scale * a["estimate"],
                "complete_ci_lower": scale * a["lower"],
                "complete_ci_upper": scale * a["upper"],
                "complete_variable_n": a["domain_n"],
                "excluded_estimate": scale * b["estimate"],
                "excluded_ci_lower": scale * b["lower"],
                "excluded_ci_upper": scale * b["upper"],
                "excluded_variable_n": b["domain_n"],
                "standardized_mean_difference": smd,
                "absolute_smd": abs(smd),
            })
    pd.DataFrame(summary).to_csv(TABLES / "selection_bias_flow_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(rows).to_csv(TABLES / "selection_bias_comparison.csv", index=False, encoding="utf-8-sig")
    print(pd.DataFrame(summary).to_string(index=False))
    print(f"Wrote {len(rows)} comparison rows to {TABLES / 'selection_bias_comparison.csv'}")


if __name__ == "__main__":
    main()
