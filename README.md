# Low FIB-4 and elevated liver stiffness in operational MASLD

Reproducible code for a repeated cross-sectional analysis of NHANES
2017–March 2020 and August 2021–August 2023. The study quantifies the
classification discordance between age-adjusted low-risk FIB-4 and elevated
VCTE liver stiffness among adults meeting an operational MASLD
definition, evaluates diabetes and age subgroups, assesses item-nonresponse
selection, and explores correlates of LSM >=8.0 kPa among participants with
low FIB-4.

VCTE is treated as a non-invasive comparator, not a histological reference
standard. The outputs must not be interpreted as biopsy-confirmed false-negative
rates or evidence that additional VCTE testing improves clinical outcomes.

## Repository structure

```text
.
|-- data/
|   |-- raw/            # downloaded NHANES XPT files; ignored by Git
|   `-- derived/        # participant-level analytic data; ignored by Git
|-- results/
|   |-- tables/         # manuscript-ready CSV outputs
|   `-- figures/        # PNG and PDF figures
|-- scripts/
|   |-- 01_download_nhanes.py
|   |-- 02_build_cohort.py
|   |-- analysis.py
|   |-- 04_selection_bias.py
|   |-- 05_validate_table2.R
|   |-- 06_logistic_regression.R
|   |-- 07_make_figures.py
|   |-- 08_nonlinearity_sensitivity.R
|   `-- 09_supplementary_population_and_weight_sensitivity.py
|-- install_r_packages.R
|-- requirements.txt
`-- run_all.py
```

## Data source

All input data are public NHANES files distributed by the US National Center
for Health Statistics. The repository does not redistribute XPT files or the
participant-level derived cohort. `scripts/01_download_nhanes.py` downloads the
required files directly from the official CDC/NCHS server.

The two survey periods are analysed separately. They must not be interpreted as
a longitudinal cohort or pooled temporal trend.

## Main operational definitions

- Adults: age >=20 years.
- Valid VCTE: complete examination, at least 10 valid measurements, available
  CAP and LSM, and reported LSM IQR/median <30% where available.
- Steatosis: CAP >=274 dB/m.
- Operational MASLD: steatosis plus at least one available
  cardiometabolic risk factor, after excluding estimated alcohol intake
  >=30 g/day in men or >=20 g/day in women, HBsAg positivity, or detectable
  HCV RNA. Overweight or central obesity uses ethnicity-specific BMI and waist
  thresholds, including a 90-cm waist threshold for Asian men.
- Low FIB-4: <1.3 at age <=65 years and <2.0 above age 65 years.
- Elevated LSM: >=8.0 kPa in the primary analysis; 8.6 and 10.0 kPa in
  sensitivity analyses.
- Diabetes subgroup: self-reported physician-diagnosed diabetes, HbA1c >=6.5%,
  or current glucose-lowering treatment.
- Cardiometabolic glycaemic criterion: self-reported diabetes, borderline
  diabetes or prediabetes, HbA1c >=5.7%, or current glucose-lowering treatment.
- Low HDL-C criterion: <=40 mg/dL in men or <=50 mg/dL in women, or current
  lipid-lowering treatment.

These definitions follow the manuscript analysis plan but do not reproduce
every element of a clinical MASLD diagnosis from the public NHANES variables.

## Software

The final analysis was run with Python 3.13.14 and R 4.6.1. Python package
versions are pinned in `requirements.txt`. The R analysis uses `survey` 4.5.

### Install Python dependencies

```bash
python -m venv .venv
```

Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

macOS/Linux:

```bash
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### Install the R dependency

```bash
Rscript install_r_packages.R
```

Ensure `Rscript` is available on `PATH` before running the complete pipeline.

## Reproduce all results

From the repository root:

```bash
python run_all.py
```

The pipeline performs the following steps:

1. downloads official NHANES XPT files;
2. merges components and constructs the analytic variables;
3. generates baseline, bidirectional-discordance, subgroup, threshold and age
   sensitivity tables using Taylor-linearized domain estimates;
4. compares complete analytic observations with excluded adults;
5. independently reproduces Table 2 with R `survey`;
6. fits cycle-specific survey-weighted exploratory logistic models;
7. checks continuous-variable functional form with survey-weighted natural
   spline sensitivity models;
8. repeats the linear model with additional race/ethnicity adjustment;
9. exports the participant-flow and different-denominator figures;
10. computes weighted population totals and the 2021-2023 MEC-weight
    sensitivity analysis.

The R validation output can be compared with
`results/tables/revision_table2_bidirectional.csv`. In the final development
environment, the maximum absolute difference across unrounded point estimates
and confidence limits was approximately 1.3e-12 percentage points.

## Survey weights

The 2017-March 2020 analyses use the prepandemic MEC examination weight
`WTMECPRP`. Because FIB-4 requires blood analytes, the August 2021-August 2023
primary and regression analyses use the cycle-specific phlebotomy weight
`WTPH2YR`, as recommended by NCHS. The complete-case versus excluded-adult
comparison continues to use the MEC examination weight `WTMEC2YR`, because
phlebotomy nonrespondents have a zero phlebotomy weight and are part of the
population being compared in that analysis.

## Logistic regression

The exploratory models are restricted to participants with operational MASLD
and low FIB-4. The binary outcome is LSM >=8.0 kPa. Primary models adjust for age,
sex, BMI, diabetes, serum triglycerides, HDL-C, systolic blood pressure and CAP while
retaining survey weights, strata and PSUs. Continuous-variable odds ratios are
scaled as stated in the output table. The regression is associative and is not
a clinical prediction model or referral rule. A sensitivity model additionally
adjusts for race/ethnicity.

## Expected key checks

- Operational MASLD sample sizes: 2,785 and 1,757.
- Design degrees of freedom: 25 and 15.
- Low-FIB-4 participants with LSM >=8.0 kPa: approximately 15.4% and 18.9%.
- Logistic model sizes: 2,150 participants/331 events and 1,433/274 events.

If these checks differ materially, confirm the downloaded file versions,
software environment and code version before using the outputs.

## Privacy and repository contents

NHANES public data are de-identified, but participant-level files are excluded
from Git by default to keep the repository small and to follow data-minimisation
practice. Do not force-add `data/raw` or `data/derived` unless the target
repository and journal explicitly require redistribution.

## Citation and archival release

Before public release, replace the placeholders in `CITATION.cff` and `LICENSE`.
After pushing to GitHub, create a versioned GitHub release (for example `v1.0.0`)
and archive that release with Zenodo. Zenodo will assign a DOI; place the DOI
and the GitHub release URL in the manuscript's Code Availability statement.

Suggested statement after archiving:

> All code used for cohort construction, survey-weighted domain estimation,
> complete-case comparison, independent R survey validation, regression,
> tables and figures is available at [GitHub release URL] and archived at
> [Zenodo DOI], version 1.0.0.

## License

Code is released under the MIT License. NHANES data remain subject to the terms
and documentation of CDC/NCHS and are not covered by this software license.
