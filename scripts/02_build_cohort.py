from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
TABLES = ROOT / "results" / "tables"
DERIVED.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)


def load(name: str) -> pd.DataFrame:
    path = RAW / name
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}. Run scripts/01_download_nhanes.py first.")
    return pd.read_sas(path, format="xport")


def build_cycle(prefix: str) -> pd.DataFrame:
    if prefix == "2017-2020":
        files = {
            "demo": "P_DEMO.xpt",
            "bmx": "P_BMX.xpt",
            "bio": "P_BIOPRO.xpt",
            "cbc": "P_CBC.xpt",
            "lux": "P_LUX.xpt",
            "bp": "P_BPXO.xpt",
            "ghb": "P_GHB.xpt",
            "hdl": "P_HDL.xpt",
            "tg": "P_TRIGLY.xpt",
            "diq": "P_DIQ.xpt",
            "bpq": "P_BPQ.xpt",
            "alq": "P_ALQ.xpt",
            "heq": "P_HEQ.xpt",
            "hepbd": "P_HEPBD.xpt",
            "hepc": "P_HEPC.xpt",
        }
        exam_weight = "WTMECPRP"
    else:
        files = {
            "demo": "DEMO_L.xpt",
            "bmx": "BMX_L.xpt",
            "bio": "BIOPRO_L.xpt",
            "cbc": "CBC_L.xpt",
            "lux": "LUX_L.xpt",
            "bp": "BPXO_L.xpt",
            "ghb": "GHB_L.xpt",
            "hdl": "HDL_L.xpt",
            "tg": "TRIGLY_L.xpt",
            "diq": "DIQ_L.xpt",
            "bpq": "BPQ_L.xpt",
            "alq": "ALQ_L.xpt",
            "heq": "HEQ_L.xpt",
            "hepbd": "HEPBD_L.xpt",
            "hepc": "HEPC_L.xpt",
        }
        exam_weight = "WTMEC2YR"

    parts = {key: load(value) for key, value in files.items()}
    df = parts["demo"][["SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3", exam_weight,
                         "SDMVSTRA", "SDMVPSU"]]
    df = df.merge(parts["bmx"][["SEQN", "BMXBMI", "BMXWAIST"]], on="SEQN", how="inner")
    bio_columns = ["SEQN", "LBXSATSI", "LBXSASSI", "LBXSTR"]
    if prefix == "2021-2023":
        bio_columns.append("WTPH2YR")
    df = df.merge(parts["bio"][bio_columns], on="SEQN", how="inner")
    df = df.merge(parts["cbc"][["SEQN", "LBXPLTSI"]], on="SEQN", how="inner")
    df = df.merge(parts["lux"][["SEQN", "LUAXSTAT", "LUANMVGP", "LUXSMED",
                                 "LUXSIQRM", "LUXCAPM", "LUAPNME"]], on="SEQN", how="inner")
    for key in ("bp", "ghb", "hdl", "tg", "diq", "bpq", "alq", "heq", "hepbd", "hepc"):
        keep = [c for c in parts[key].columns if c == "SEQN" or c not in df.columns]
        df = df.merge(parts[key][keep], on="SEQN", how="left")
    df = df.rename(columns={exam_weight: "selection_weight", "LBXSATSI": "ALT", "LBXSASSI": "AST",
                            "LBXPLTSI": "platelets", "LUXSMED": "LSM", "LUXCAPM": "CAP"})
    # The August 2021-August 2023 release introduced WTPH2YR to account for
    # phlebotomy nonresponse. FIB-4 requires blood analytes, so WTPH2YR is the
    # analysis weight for that cycle. The MEC weight is retained separately for
    # complete-case versus excluded-participant comparisons.
    if prefix == "2021-2023":
        df = df.rename(columns={"WTPH2YR": "weight"})
    else:
        df["weight"] = df["selection_weight"]
    df["cycle"] = prefix
    df["fib4"] = df["RIDAGEYR"] * df["AST"] / (df["platelets"] * np.sqrt(df["ALT"]))
    df["adult"] = df["RIDAGEYR"] >= 20
    df["valid_vcte"] = (
        (df["LUAXSTAT"] == 1)
        & (df["LUANMVGP"] >= 10)
        & df["LSM"].notna()
        & df["CAP"].notna()
        & ((df["LUXSIQRM"].isna()) | (df["LUXSIQRM"] < 30))
    )
    df["complete_fib4"] = (
        df[["RIDAGEYR", "ALT", "AST", "platelets"]].notna().all(axis=1)
        & (df["ALT"] > 0)
        & (df["platelets"] > 0)
    )
    df["steatosis274"] = df["CAP"] >= 274
    asian = df["RIDRETH3"] == 6
    bmi_criterion = np.where(asian, df["BMXBMI"] >= 23, df["BMXBMI"] >= 25)
    male_waist_cutoff = np.where(asian, 90, 94)
    waist_criterion = np.where(
        df["RIAGENDR"] == 1,
        df["BMXWAIST"] >= male_waist_cutoff,
        df["BMXWAIST"] >= 80,
    )
    df["overweight_cmf"] = bmi_criterion | waist_criterion
    # This is deliberately labelled a proxy: formal MASLD requires all five
    # cardiometabolic criteria, which are added only after the event audit passes.
    df["masld_proxy"] = df["steatosis274"] & df["overweight_cmf"]

    systolic = [c for c in ("BPXOSY1", "BPXOSY2", "BPXOSY3") if c in df]
    diastolic = [c for c in ("BPXODI1", "BPXODI2", "BPXODI3") if c in df]
    df["mean_sbp"] = df[systolic].mean(axis=1)
    df["mean_dbp"] = df[diastolic].mean(axis=1)
    diabetes_report = df["DIQ010"].isin([1, 3])
    diabetes_treatment = df["DIQ050"].eq(1) | df["DIQ070"].eq(1)
    bp_treatment_var = "BPQ050A" if prefix == "2017-2020" else "BPQ150"
    lipid_treatment_var = "BPQ100D" if prefix == "2017-2020" else "BPQ101D"
    bp_treatment = df[bp_treatment_var].eq(1)
    lipid_treatment = df[lipid_treatment_var].eq(1)
    # The MASLD glycaemic criterion includes self-reported prediabetes. DIQ010
    # captures diabetes/borderline diabetes and DIQ160 separately captures a
    # clinician diagnosis of prediabetes.
    prediabetes_report = df["DIQ160"].eq(1)
    df["cmrf_glucose"] = (
        diabetes_report | prediabetes_report | diabetes_treatment | (df["LBXGH"] >= 5.7)
    )
    df["cmrf_bp"] = (df["mean_sbp"] >= 130) | (df["mean_dbp"] >= 85) | bp_treatment
    df["cmrf_tg"] = (df["LBXSTR"] >= 150) | lipid_treatment
    low_hdl = np.where(df["RIAGENDR"] == 1, df["LBDHDD"] <= 40, df["LBDHDD"] <= 50)
    df["cmrf_hdl"] = low_hdl | lipid_treatment
    df["any_cmrf"] = df[["overweight_cmf", "cmrf_glucose", "cmrf_bp", "cmrf_tg", "cmrf_hdl"]].any(axis=1)

    annual_days = df["ALQ121"].map({0: 0, 1: 365, 2: 300, 3: 182, 4: 104, 5: 52,
                                     6: 30, 7: 12, 8: 9, 9: 4.5, 10: 1.5})
    avg_drinks = df["ALQ130"].where(df["ALQ130"].between(0, 30), np.nan)
    df["alcohol_g_day"] = annual_days * avg_drinks * 14 / 365
    df.loc[df["ALQ111"].eq(2) | df["ALQ121"].eq(0), "alcohol_g_day"] = 0
    alcohol_limit = np.where(df["RIAGENDR"] == 1, 30, 20)
    df["high_alcohol"] = df["alcohol_g_day"] >= alcohol_limit
    df["alcohol_known"] = df["alcohol_g_day"].notna()
    df["hbv_known"] = df["LBDHBG"].isin([1, 2, 3])
    df["hcv_known"] = df["LBXHCR"].isin([1, 2, 3])
    df["etiology_known"] = df["alcohol_known"] & df["hbv_known"] & df["hcv_known"]
    df["viral_hepatitis"] = df["LBDHBG"].eq(1) | df["LBXHCR"].eq(1)
    df["masld_formal"] = df["steatosis274"] & df["any_cmrf"] & ~df["high_alcohol"] & ~df["viral_hepatitis"]
    df["fib4_low_fixed"] = df["fib4"] < 1.3
    age_adjusted_cut = np.where(df["RIDAGEYR"] > 65, 2.0, 1.3)
    df["fib4_low_age"] = df["fib4"] < age_adjusted_cut
    return df


def weighted_pct(df: pd.DataFrame, mask: pd.Series) -> float:
    ok = df["weight"].notna() & (df["weight"] > 0)
    if not ok.any():
        return np.nan
    return 100 * df.loc[ok & mask, "weight"].sum() / df.loc[ok, "weight"].sum()


def audit(df: pd.DataFrame) -> list[dict]:
    rows = []
    base = df[df["adult"] & df["valid_vcte"] & df["complete_fib4"]].copy()
    populations = {
        "all_valid_adults": base,
        "masld_proxy": base[base["masld_proxy"]],
        "masld_formal": base[base["masld_formal"]],
    }
    for pop_name, pop in populations.items():
        for lsm_cut in (8.0, 8.6, 10.0):
            elevated = pop["LSM"] >= lsm_cut
            for fib4_name in ("fib4_low_fixed", "fib4_low_age"):
                missed = elevated & pop[fib4_name]
                rows.append({
                    "cycle": pop["cycle"].iloc[0],
                    "population": pop_name,
                    "lsm_cut": lsm_cut,
                    "fib4_rule": fib4_name,
                    "n_population": len(pop),
                    "n_lsm_elevated": int(elevated.sum()),
                    "n_low_fib4_lsm_elevated": int(missed.sum()),
                    "pct_of_lsm_elevated_low_fib4_unweighted": (
                        100 * missed.sum() / elevated.sum() if elevated.sum() else np.nan
                    ),
                    "pct_population_lsm_elevated_weighted": weighted_pct(pop, elevated),
                    "pct_population_low_fib4_lsm_elevated_weighted": weighted_pct(pop, missed),
                })
    return rows


def main() -> None:
    cycles = [build_cycle("2017-2020"), build_cycle("2021-2023")]
    combined = pd.concat(cycles, ignore_index=True)
    results = pd.DataFrame([row for cycle in cycles for row in audit(cycle)])
    combined.to_csv(DERIVED / "core_analytic_cohort.csv", index=False)
    results.to_csv(TABLES / "feasibility_counts.csv", index=False)
    display = results[
        (results["population"] == "masld_formal")
        & (results["fib4_rule"] == "fib4_low_age")
    ]
    print(display.to_string(index=False))
    print("\nFiles written:")
    print(DERIVED / "core_analytic_cohort.csv")
    print(TABLES / "feasibility_counts.csv")


if __name__ == "__main__":
    main()
