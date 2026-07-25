"""Supplementary analyses requested during pre-submission review.

A. Absolute weighted population totals (millions) for key domains, with
   Taylor-linearized standard errors on the full cycle design sample.
B. Weight-scheme sensitivity for August 2021-August 2023: repeat the main
   bidirectional estimates and the diabetes subgroup using the MEC
   examination weight (WTMEC2YR) instead of the phlebotomy weight (WTPH2YR).
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis import domain_mean_ci, design_df, prepare  # noqa: E402

TABLES = ROOT / "results" / "tables"


def weighted_total_ci(base, domain):
    """Estimated population total for a domain indicator, with Taylor SE."""
    d = base[["weight", "SDMVSTRA", "SDMVPSU"]].copy()
    dom = pd.Series(domain, index=base.index).fillna(False).astype(bool)
    valid = dom & d["weight"].gt(0)
    d["linearized"] = 0.0
    d.loc[valid, "linearized"] = d.loc[valid, "weight"]
    total = d["linearized"].sum()
    clusters = d.groupby(["SDMVSTRA", "SDMVPSU"], as_index=False)["linearized"].sum()
    variance_total = 0.0
    for _, stratum in clusters.groupby("SDMVSTRA"):
        n_psu = len(stratum)
        if n_psu < 2:
            continue
        variance_total += n_psu / (n_psu - 1) * np.square(
            stratum["linearized"] - stratum["linearized"].mean()
        ).sum()
    se = np.sqrt(variance_total)
    critical = t.ppf(0.975, design_df(base))
    return {
        "total_millions": total / 1e6,
        "ci_lower_millions": (total - critical * se) / 1e6,
        "ci_upper_millions": (total + critical * se) / 1e6,
        "domain_n": int(valid.sum()),
    }


def population_totals(d):
    rows = []
    for cycle, base in d.groupby("cycle"):
        masld = base["masld_domain"]
        low = base["fib4_low_age"].astype(bool)
        elevated8 = (base["LSM"] >= 8.0)
        diabetes = base["diabetes"].eq(1)
        domains = {
            "Operational MASLD": masld,
            "Operational MASLD + low FIB-4": masld & low,
            "Operational MASLD + low FIB-4 + LSM>=8.0 kPa": masld & low & elevated8,
            "Operational MASLD + diabetes": masld & diabetes,
            "Operational MASLD + diabetes + low FIB-4": masld & diabetes & low,
            "Operational MASLD + diabetes + low FIB-4 + LSM>=8.0 kPa":
                masld & diabetes & low & elevated8,
        }
        for label, domain in domains.items():
            r = weighted_total_ci(base, domain)
            rows.append({"cycle": cycle, "domain": label, **r})
    return pd.DataFrame(rows)


def mec_weight_sensitivity(d):
    """Repeat main + diabetes estimates for 2021-2023 with WTMEC2YR."""
    base = d[d["cycle"] == "2021-2023"].copy()
    base["weight"] = base["selection_weight"]
    base = base[base["weight"].gt(0)].copy()
    masld = base["masld_domain"]
    low = base["fib4_low_age"].astype(bool)
    rows = []
    for cutoff in (8.0, 8.6, 10.0):
        elevated = (base["LSM"] >= cutoff).astype(int)
        discordant = (low & elevated.eq(1)).astype(int)
        a = domain_mean_ci(base, masld, elevated)
        b = domain_mean_ci(base, masld, discordant)
        c = domain_mean_ci(base, masld & low, elevated)
        rev = domain_mean_ci(base, masld & elevated.eq(1), low.astype(int))
        diab = masld & base["diabetes"].eq(1)
        dc = domain_mean_ci(base, diab & low, elevated)
        rows.append({
            "weight": "WTMEC2YR (MEC)", "lsm_cutoff": cutoff,
            "masld_n": int(masld.sum()),
            "low_fib4_n": int((masld & low).sum()),
            "lsm_elevated_percent": 100 * a["estimate"],
            "lsm_elevated_ci": f'{100*a["lower"]:.1f}-{100*a["upper"]:.1f}',
            "discordant_in_masld_percent": 100 * b["estimate"],
            "discordant_in_masld_ci": f'{100*b["lower"]:.1f}-{100*b["upper"]:.1f}',
            "lsm_elevated_among_low_fib4_percent": 100 * c["estimate"],
            "lsm_elevated_among_low_fib4_ci": f'{100*c["lower"]:.1f}-{100*c["upper"]:.1f}',
            "low_fib4_among_lsm_percent": 100 * rev["estimate"],
            "low_fib4_among_lsm_ci": f'{100*rev["lower"]:.1f}-{100*rev["upper"]:.1f}',
            "diabetes_low_fib4_lsm_elevated_percent": 100 * dc["estimate"],
            "diabetes_low_fib4_lsm_elevated_ci": f'{100*dc["lower"]:.1f}-{100*dc["upper"]:.1f}',
            "design_df": a["design_df"],
        })
    return pd.DataFrame(rows)


if __name__ == "__main__":
    data = prepare()

    totals = population_totals(data)
    totals.to_csv(TABLES / "supplementary_population_totals.csv",
                  index=False, encoding="utf-8-sig")
    print("=== A. Weighted population totals (millions) ===")
    print(totals.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    sens = mec_weight_sensitivity(data)
    sens.to_csv(TABLES / "supplementary_mec_weight_sensitivity.csv",
                index=False, encoding="utf-8-sig")
    print("\n=== B. 2021-2023 with WTMEC2YR (vs WTPH2YR primary: 20.4/16.4/18.9/80.4, diabetes 28.8) ===")
    print(sens.to_string(index=False))
