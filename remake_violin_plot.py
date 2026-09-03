#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


def require_columns(df, required, input_name):
    missing = set(required) - set(df.columns)
    if missing:
        raise SystemExit(
            f"{input_name} is missing required columns: {sorted(missing)}"
        )


def plot_top_bottom_violin(read_level, selected_summary, out_png):
    require_columns(
        read_level,
        {"dmr_id", "sample_group", "read_m_frac"},
        "Read-level TSV",
    )
    require_columns(
        selected_summary,
        {
            "dmr_id",
            "dmr_rank_label",
            "low_methylation_excess_patient_minus_control",
        },
        "Selected-summary TSV",
    )

    selected_summary = selected_summary.copy()
    plot_order = selected_summary["dmr_id"].tolist()
    n_dmrs = len(plot_order)

    if n_dmrs == 0:
        raise SystemExit("The selected-summary TSV contains no DMRs.")

    available_groups = set(read_level["sample_group"].dropna().unique())
    missing_groups = {"control", "patient"} - available_groups
    if missing_groups:
        raise SystemExit(
            "The read-level TSV must use the internal sample_group values "
            f"'control' and 'patient'. Missing: {sorted(missing_groups)}"
        )

    control_data = []
    urpl_data = []

    for dmr_id in plot_order:
        sub = read_level[read_level["dmr_id"] == dmr_id]

        control_vals = sub.loc[
            sub["sample_group"] == "control",
            "read_m_frac",
        ].dropna().to_numpy()

        urpl_vals = sub.loc[
            sub["sample_group"] == "patient",
            "read_m_frac",
        ].dropna().to_numpy()

        if len(control_vals) == 0 or len(urpl_vals) == 0:
            raise SystemExit(
                f"DMR {dmr_id!r} does not have read_m_frac values for both "
                "control and patient groups."
            )

        control_data.append(control_vals)
        urpl_data.append(urpl_vals)

    fig_width = max(12, n_dmrs * 0.55)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    positions = np.arange(1, n_dmrs + 1)

    vp_control = ax.violinplot(
        control_data,
        positions=positions,
        widths=0.75,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )

    vp_urpl = ax.violinplot(
        urpl_data,
        positions=positions,
        widths=0.75,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )

    control_color = "#1b9e77"
    urpl_color = "#d95f02"

    for body in vp_control["bodies"]:
        body.set_facecolor(control_color)
        body.set_edgecolor(control_color)
        body.set_alpha(0.35)

    for body in vp_urpl["bodies"]:
        body.set_facecolor(urpl_color)
        body.set_edgecolor(urpl_color)
        body.set_alpha(0.35)

    vp_control["cmedians"].set_color(control_color)
    vp_control["cmedians"].set_linewidth(1.2)
    vp_urpl["cmedians"].set_color(urpl_color)
    vp_urpl["cmedians"].set_linewidth(1.2)

    labels = []
    for _, row in selected_summary.iterrows():
        labels.append(
            f"{row['dmr_rank_label']}\n"
            f"{row['dmr_id']}\n"
            f"Δlow={row['low_methylation_excess_patient_minus_control']:.2f}"
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("Read-level m_frac within DMR")
    ax.set_xlabel("Selected DMRs")
    ax.set_title(
        "Overlaid read-level m_frac distributions for top and bottom DMRs\n"
        "ranked by uRPL excess of low-methylation reads"
    )
    ax.grid(axis="y", alpha=0.25)

    handles = [
        Patch(
            facecolor=control_color,
            edgecolor=control_color,
            alpha=0.35,
            label="Control",
        ),
        Patch(
            facecolor=urpl_color,
            edgecolor=urpl_color,
            alpha=0.35,
            label="uRPL",
        ),
    ]
    ax.legend(handles=handles, loc="upper right")

    n_bottom = selected_summary["dmr_rank_label"].str.startswith("bottom").sum()
    if 0 < n_bottom < n_dmrs:
        ax.axvline(n_bottom + 0.5, color="black", linewidth=1, alpha=0.4)

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Regenerate the top/bottom DMR violin plot from existing "
            "read-level and selected-summary TSV files."
        )
    )
    ap.add_argument(
        "--read-level",
        required=True,
        help="Existing top_bottom_low_methylation_excess_read_level.tsv",
    )
    ap.add_argument(
        "--selected-summary",
        required=True,
        help="Existing top_bottom_low_methylation_excess_dmrs.tsv",
    )
    ap.add_argument("--out-png", required=True)
    args = ap.parse_args()

    read_level = pd.read_csv(args.read_level, sep="\t")
    selected_summary = pd.read_csv(args.selected_summary, sep="\t")

    plot_top_bottom_violin(
        read_level=read_level,
        selected_summary=selected_summary,
        out_png=args.out_png,
    )

    print(f"Wrote: {args.out_png}")


if __name__ == "__main__":
    main()
