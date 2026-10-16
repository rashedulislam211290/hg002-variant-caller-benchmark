#!/usr/bin/env python3
"""Step 6: build Figure 1 and a summary table from the hap.py output.

Usage (from the repository root):
    python 04_make_figure.py --results results --out figures
Reads results/bcftools.summary.csv and results/freebayes.summary.csv.
"""
import argparse
import os
from decimal import Decimal, ROUND_HALF_UP

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CALLERS = {"bcftools": "BCFtools", "freebayes": "FreeBayes"}
COLORS = {"BCFtools": "#2a78d6", "FreeBayes": "#eb6834"}
TYPES = {"SNP": "SNV", "INDEL": "Indel"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def load(results_dir):
    rows = []
    for key, name in CALLERS.items():
        df = pd.read_csv(os.path.join(results_dir, f"{key}.summary.csv"))
        df = df[df["Filter"] == "ALL"]
        for _, r in df.iterrows():
            fp, gt, al = int(r["QUERY.FP"]), int(r["FP.gt"]), int(r["FP.al"])
            rows.append({
                "Type": TYPES[r["Type"]], "Caller": name,
                "TP": int(r["TRUTH.TP"]), "FN": int(r["TRUTH.FN"]), "FP": fp,
                "FP.gt": gt, "FP.al": al, "FP.other": fp - gt - al,
                "Recall": r["METRIC.Recall"], "Precision": r["METRIC.Precision"],
                "F1": r["METRIC.F1_Score"],
            })
    return pd.DataFrame(rows).set_index(["Type", "Caller"])


def plot(d, out_dir):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.linewidth": 0.6})
    fig = plt.figure(figsize=(7.2, 6.4), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], hspace=0.55, wspace=0.42)

    # (A) recall / precision / F1
    ax = fig.add_subplot(gs[0, :])
    metrics = [("Recall", "Recall"), ("Precision", "Precision"), ("F1", "F1 score")]
    rows = [(t, m) for t in ("SNV", "Indel") for m in metrics]
    y = np.arange(len(rows))[::-1].astype(float)
    y[3:] -= 0.6
    lo_all = 1.0
    for yi, (t, (col, label)) in zip(y, rows):
        a, b = d.loc[(t, "BCFtools"), col], d.loc[(t, "FreeBayes"), col]
        lo_all = min(lo_all, a, b)
        ax.plot([a, b], [yi, yi], color="#c3c2b7", lw=1.2, zorder=1)
        for tool, v in (("BCFtools", a), ("FreeBayes", b)):
            ax.scatter(v, yi, s=42, color=COLORS[tool], edgecolor="white", linewidth=1.2, zorder=3)
        ax.text(max(a, b) + 0.0015, yi, f"{max(a, b):.4f}", va="center", ha="left", fontsize=6.5, color=MUTED)
        ax.text(min(a, b) - 0.0015, yi, f"{min(a, b):.4f}", va="center", ha="right", fontsize=6.5, color=MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{t}  {label}" for t, (_, label) in rows])
    ax.set_xlim(np.floor((lo_all - 0.01) * 200) / 200, 1.004)
    ax.set_xlabel("Metric value")
    ax.xaxis.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    for tool in COLORS:
        ax.scatter([], [], s=42, color=COLORS[tool], label=tool)
    ax.legend(loc="upper left", frameon=False, ncol=2, handletextpad=0.3, fontsize=7.5)
    ax.set_title("A", loc="left", fontweight="bold", fontsize=11, x=-0.16)

    # (B) FN and FP counts
    ax = fig.add_subplot(gs[1, 0])
    groups = [("SNV", "FN"), ("SNV", "FP"), ("Indel", "FN"), ("Indel", "FP")]
    x = np.array([0, 1, 2.4, 3.4])
    w = 0.36
    ymax = 0
    for j, tool in enumerate(COLORS):
        vals = [d.loc[(t, tool), e] for t, e in groups]
        ymax = max(ymax, max(vals))
        xs = x + (j - 0.5) * (w + 0.03)
        ax.bar(xs, vals, w, color=COLORS[tool], label=tool, linewidth=0)
        for xb, v in zip(xs, vals):
            ax.text(xb, v + ymax * 0.015, f"{v:,}", ha="center", va="bottom", fontsize=6, color=MUTED, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels(["FN", "FP", "FN", "FP"])
    for xc, lab in ((0.5, "SNV"), (2.9, "Indel")):
        ax.text(xc, -0.20, lab, transform=ax.get_xaxis_transform(), ha="center", fontsize=8, color=INK)
    ax.set_ylabel("Number of variants")
    ax.set_ylim(0, ymax * 1.24)
    ax.yaxis.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.set_title("B", loc="left", fontweight="bold", fontsize=11, x=-0.30)

    # (C) composition of false positives
    ax = fig.add_subplot(gs[1, 1])
    cats = [("SNV", "BCFtools"), ("SNV", "FreeBayes"), ("Indel", "BCFtools"), ("Indel", "FreeBayes")]
    parts = [("FP.gt", "Genotype mismatch (FP.gt)", "#c9c2f5"),
             ("FP.al", "Allele mismatch (FP.al)", "#8f82e6"),
             ("FP.other", "No matching truth variant", "#3f2fa3")]
    bottom = np.zeros(4)
    for ci, (col, name, colour) in enumerate(parts):
        vals = np.array([d.loc[c, col] / d.loc[c, "FP"] * 100 for c in cats])
        ax.bar(x, vals, 0.7, bottom=bottom, color=colour, label=name, edgecolor="white", linewidth=0.8)
        for xi, v, b in zip(x, vals, bottom):
            if v >= 8:
                ax.text(xi, b + v / 2, f"{v:.0f}%", ha="center", va="center", fontsize=6.5,
                        color=INK if ci == 0 else "white")
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([c for _, c in cats], fontsize=6.8)
    for xc, lab in ((0.5, "SNV"), (2.9, "Indel")):
        ax.text(xc, -0.20, lab, transform=ax.get_xaxis_transform(), ha="center", fontsize=8)
    ax.set_ylabel("Share of false positives (%)")
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=6.3, loc="upper center", bbox_to_anchor=(0.5, -0.27))
    ax.set_title("C", loc="left", fontweight="bold", fontsize=11, x=-0.30)

    os.makedirs(out_dir, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out_dir, f"Figure1.{ext}"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.join(out_dir, "Figure1.tiff"), dpi=300, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    d = load(a.results)
    table = d[["TP", "FN", "FP", "Recall", "Precision", "F1"]].copy()
    for c in ("Recall", "Precision", "F1"):
        table[c] = [str(Decimal(str(v)).quantize(Decimal("0.0001"), ROUND_HALF_UP)) for v in table[c]]
    table.to_csv(os.path.join(a.out if os.path.isdir(a.out) else ".", "Table1_performance.tsv"), sep="\t")
    print(table.to_string())
    plot(d, a.out)


if __name__ == "__main__":
    main()
