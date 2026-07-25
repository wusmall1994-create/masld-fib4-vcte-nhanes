"""Run the complete reproducibility pipeline from the repository root."""
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def run(command):
    print("\n>", " ".join(map(str, command)), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main():
    run([sys.executable, "scripts/01_download_nhanes.py"])
    run([sys.executable, "scripts/02_build_cohort.py"])
    run([sys.executable, "scripts/analysis.py"])
    run([sys.executable, "scripts/04_selection_bias.py"])
    rscript = shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript was not found on PATH. Install R and the survey package.")
    run([rscript, "scripts/05_validate_table2.R"])
    run([rscript, "scripts/06_logistic_regression.R"])
    run([rscript, "scripts/08_nonlinearity_sensitivity.R"])
    run([sys.executable, "scripts/07_make_figures.py"])
    print("\nComplete. See results/tables and results/figures.")


if __name__ == "__main__":
    main()
