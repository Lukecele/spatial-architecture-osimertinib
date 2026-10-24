#!/usr/bin/env python3
"""Generates publication-quality figure for GSE301973 patient PFS spatial analysis.

Outputs:
  figures/figure_gse301973_patient_pfs.png
"""

import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data" / "public_validation" / "GSE301973_patient_pfs_analysis.json"
FIG_PATH = ROOT / "figures" / "figure_gse301973_patient_pfs.png"
BASE = ROOT / "data" / "public_validation" / "GSE301973" / "extracted"

COLOR_PALETTE = {
    "Tumor cell": "#e41a1c",
    "Endothelial cell": "#377eb8",
    "Stromal cell": "#4daf4a",
    "Plasmablast": "#984ea3",
    "Lymphocyte": "#ff7f00",
    "Myeloid cell": "#ffff33",
    "Alveolar epithelial cell": "#a65628",
    "Airway epithelial cell": "#f781bf",
}


def decode(v):
    return [x.decode() if isinstance(x, bytes) else x for x in v]


def main():
    with open(JSON_PATH) as f:
        data = json.load(f)

    patients = data["patients"]
    # Sort by PFS descending: s1 (69), s2 (22), s3 (15), s4 (12)
    patients = sorted(patients, key=lambda p: p["PFS_months"], reverse=True)

    fig = plt.figure(figsize=(16, 12), dpi=300)
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.28)

    # -------------------------------------------------------------
    # Panel A: Spatial Subsampling of the 4 Patient Regions
    # -------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.set_title("A. Visium HD Subcellular Maps (Subsampled 15,000 cells)", fontsize=13, fontweight="bold", loc="left")

    # Load and downsample cells for fast, clean rendering
    for slide, split_x, p_left, p_right in [("A", 12000, "s1_PFS_69", "s2_PFS_22"), ("C", 12500, "s3_PFS_15", "s4_PFS_12")]:
        h5_file = BASE / slide / f"{slide}_stardist_integrated.h5ad"
        with h5py.File(h5_file, "r") as h:
            coords = h["obsm"]["spatial"][:]
            cats = decode(h["obs"]["cell_type"]["categories"][:])
            codes = h["obs"]["cell_type"]["codes"][:]
            labels = np.array([cats[c] for c in codes])

        # Subsample ~3,500 cells per patient
        x = coords[:, 0]
        for mask, pid, off_x, off_y in [
            (x < split_x, p_left, 0 if slide == "A" else 25000, 0),
            (x >= split_x, p_right, 0 if slide == "A" else 25000, -25000 if slide == "A" else -25000),
        ]:
            idx = np.where(mask)[0]
            if len(idx) > 3000:
                np.random.seed(42)
                idx = np.random.choice(idx, size=3000, replace=False)
            
            sub_coords = coords[idx]
            sub_labels = labels[idx]

            for ct, col in COLOR_PALETTE.items():
                ct_mask = (sub_labels == ct)
                if np.sum(ct_mask) > 0:
                    ax_a.scatter(
                        sub_coords[ct_mask, 0] + off_x,
                        sub_coords[ct_mask, 1] + off_y,
                        c=col,
                        s=0.8,
                        alpha=0.6,
                        edgecolors="none",
                    )

    ax_a.text(6000, 24000, "Slide A: s1 (PFS 69 mo)", fontsize=10, fontweight="bold", color="#1f77b4")
    ax_a.text(6000, -2000, "Slide A: s2 (PFS 22 mo)", fontsize=10, fontweight="bold", color="#1f77b4")
    ax_a.text(31000, 24000, "Slide C: s3 (PFS 15 mo)", fontsize=10, fontweight="bold", color="#d62728")
    ax_a.text(31000, -2000, "Slide C: s4 (PFS 12 mo)", fontsize=10, fontweight="bold", color="#d62728")
    ax_a.set_xticks([])
    ax_a.set_yticks([])
    ax_a.axis("off")

    # Legend for Cell Types
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=col, markersize=7, label=ct)
        for ct, col in COLOR_PALETTE.items()
    ]
    ax_a.legend(handles=handles, loc="lower right", fontsize=7.5, framealpha=0.9)

    # -------------------------------------------------------------
    # Panel B: Homotypic Tumor Contact vs PFS Months
    # -------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.set_title("B. Tumor Confluence vs Progression-Free Survival (PFS)", fontsize=13, fontweight="bold", loc="left")

    pfs_vals = [p["PFS_months"] for p in patients]
    homo_vals = [p["key_tumor_contacts"]["homotypic_tumor"] * 100 for p in patients]
    colors = ["#2ca02c", "#2ca02c", "#d62728", "#d62728"]
    markers = ["o", "o", "s", "s"]

    for i, p in enumerate(patients):
        ax_b.scatter(
            pfs_vals[i],
            homo_vals[i],
            color=colors[i],
            s=180,
            zorder=4,
            label=f"{p['patient_id']} ({p['PFS_group']})",
        )
        ax_b.annotate(
            f" {p['patient_id']}\n ({homo_vals[i]:.1f}%)",
            (pfs_vals[i], homo_vals[i]),
            textcoords="offset points",
            xytext=(10, -5),
            fontsize=9,
            fontweight="semibold",
        )

    ax_b.axhline(75.9, color="#d62728", linestyle="--", alpha=0.5, label="Mean Short PFS (75.9%)")
    ax_b.axhline(64.3, color="#2ca02c", linestyle="--", alpha=0.5, label="Mean Long PFS (64.3%)")

    ax_b.set_xlabel("Progression-Free Survival (Months)", fontsize=11, fontweight="bold")
    ax_b.set_ylabel("Homotypic Tumor-Tumor Contact (%)", fontsize=11, fontweight="bold")
    ax_b.set_ylim(50, 95)
    ax_b.set_xlim(5, 78)
    ax_b.grid(True, linestyle=":", alpha=0.6)
    ax_b.legend(loc="upper right", fontsize=8)

    # -------------------------------------------------------------
    # Panel C: Microenvironment Neighborhood Contact Breakdown
    # -------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.set_title("C. Tumor Neighborhood Contact Breakdown by Cell Type", fontsize=13, fontweight="bold", loc="left")

    p_ids = [p["patient_id"] for p in patients]
    x_pos = np.arange(len(p_ids))
    width = 0.55

    bottom = np.zeros(len(patients))
    contact_order = [
        ("Tumor cell", "#e41a1c"),
        ("Endothelial cell", "#377eb8"),
        ("Stromal cell", "#4daf4a"),
        ("Plasmablast", "#984ea3"),
        ("Lymphocyte", "#ff7f00"),
        ("Myeloid cell", "#ffff33"),
        ("Alveolar epithelial cell", "#a65628"),
        ("Airway epithelial cell", "#f781bf"),
    ]

    for ct, col in contact_order:
        vals = np.array([p["tumor_neighborhood_contacts"].get(ct, 0.0) * 100 for p in patients])
        ax_c.bar(x_pos, vals, width, bottom=bottom, label=ct, color=col, edgecolor="black", linewidth=0.5)
        bottom += vals

    ax_c.set_xticks(x_pos)
    ax_c.set_xticklabels([f"{p['patient_id']}\n(PFS {p['PFS_months']:.0f}m)" for p in patients], fontsize=9.5, fontweight="bold")
    ax_c.set_ylabel("Neighborhood Contact Fraction (%)", fontsize=11, fontweight="bold")
    ax_c.set_ylim(0, 100)
    ax_c.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8, title="Neighbor Cell Type")

    # -------------------------------------------------------------
    # Panel D: Tumor Molecular Programs
    # -------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.set_title("D. Pre-Treatment Tumor Cell Expression Modules", fontsize=13, fontweight="bold", loc="left")

    mods = ["egfr_bypass", "emt", "stem_like", "spp1_cd44", "proliferation"]
    mod_labels = ["EGFR Bypass", "EMT", "Stem-Like", "SPP1-CD44", "Proliferation"]
    x_mod = np.arange(len(mods))
    bar_width = 0.18

    offsets = [-1.5 * bar_width, -0.5 * bar_width, 0.5 * bar_width, 1.5 * bar_width]
    colors_pat = ["#1f77b4", "#aec7e8", "#ff9896", "#d62728"]

    for i, p in enumerate(patients):
        vals = [p["tumor_module_means"][m] for m in mods]
        ax_d.bar(
            x_mod + offsets[i],
            vals,
            bar_width,
            label=f"{p['patient_id']} ({p['PFS_months']:.0f}m)",
            color=colors_pat[i],
            edgecolor="black",
            linewidth=0.5,
        )

    ax_d.set_xticks(x_mod)
    ax_d.set_xticklabels(mod_labels, fontsize=9.5, fontweight="bold")
    ax_d.set_ylabel("Mean Expression (Log-CPM Normalized)", fontsize=11, fontweight="bold")
    ax_d.grid(True, axis="y", linestyle=":", alpha=0.6)
    ax_d.legend(loc="upper right", fontsize=8.5)

    plt.suptitle("GSE301973: Pre-Treatment Visium HD Architecture in EGFR-Mutant NSCLC", fontsize=15, fontweight="bold", y=0.98)

    FIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(FIG_PATH, bbox_inches="tight")
    plt.close()
    print(f"Figure generated successfully: {FIG_PATH}")


if __name__ == "__main__":
    main()
