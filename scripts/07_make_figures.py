"""Generate the two manuscript figures from derived data and final tables."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "data" / "derived"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)


def make_flowchart() -> None:
    data = pd.read_csv(DERIVED / "core_analytic_cohort.csv")
    flow = {}
    for cycle, x in data.groupby("cycle"):
        adults = x["adult"]
        valid = adults & x["valid_vcte"]
        complete = valid & x["complete_fib4"]
        masld = complete & x["masld_formal"]
        flow[cycle.replace("-", "–")] = [
            ("Adults", int(adults.sum())), ("Valid VCTE", int(valid.sum())),
            ("Complete FIB-4", int(complete.sum())),
            ("Operational MASLD", int(masld.sum())),
        ]
    fig, axes = plt.subplots(1, 2, figsize=(9, 5.2))
    for ax, (cycle, steps) in zip(axes, flow.items()):
        ax.axis("off")
        ax.set_title(cycle, fontsize=13, weight="bold")
        ys = np.linspace(0.86, 0.14, len(steps))
        for i, ((label, n), y) in enumerate(zip(steps, ys)):
            ax.text(0.5, y, f"{label}\nn={n:,}", ha="center", va="center", fontsize=10,
                    bbox=dict(boxstyle="round,pad=0.45", facecolor="#EEF3F8", edgecolor="#1F4D78"))
            if i < len(steps) - 1:
                ax.annotate("", xy=(0.5, ys[i + 1] + 0.07), xytext=(0.5, y - 0.07),
                            arrowprops=dict(arrowstyle="->", color="#555555"))
    fig.suptitle("Participant selection", fontsize=14, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIGURES / "figure1_participant_flow.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "figure1_participant_flow.pdf", bbox_inches="tight")
    plt.close(fig)


def make_bidirectional() -> None:
    table = pd.read_csv(TABLES / "revision_table2_bidirectional.csv")
    table = table[table["lsm_cutoff"] == 8.0]
    categories = ["All MASLD", "Among low FIB-4", "Among elevated LSM"]
    stems = ["discordant_in_masld_percent", "lsm_elevated_among_low_fib4_percent",
             "low_fib4_among_lsm_percent"]
    x = np.arange(len(categories))
    width = 0.34
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    for i, (cycle, color) in enumerate([("2017-2020", "#4C78A8"), ("2021-2023", "#F58518")]):
        row = table[table["cycle"] == cycle].iloc[0]
        values = [row[s] for s in stems]
        bars = ax.bar(x + (i - 0.5) * width, values, width,
                      label=cycle.replace("-", "–"), color=color)
        ax.bar_label(bars, labels=[f"{v:.1f}%" for v in values], padding=2, fontsize=8)
    ax.set_xticks(x, categories)
    ax.set_ylabel("Survey-weighted proportion (%)")
    ax.set_ylim(0, 92)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIGURES / "figure2_different_denominators.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "figure2_different_denominators.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    make_flowchart()
    make_bidirectional()
