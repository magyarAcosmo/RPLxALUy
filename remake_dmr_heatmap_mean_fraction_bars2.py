#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def make_5pct_bin_labels(n_bins):
    labels = []

    if n_bins == 20:
        for i in range(n_bins):
            if i == 0:
                labels.append("0-5%")
            elif i == n_bins - 1:
                labels.append("96-100%")
            else:
                labels.append(f"{i * 5 + 1}-{(i + 1) * 5}%")
        return labels

    for i in range(n_bins):
        start = i / n_bins * 100
        end = (i + 1) / n_bins * 100
        labels.append(f"{start:.0f}-{end:.0f}%")

    return labels


def plot_three_panel_heatmap_from_bin_fractions(
    bin_fractions,
    out_png,
    n_bins=None,
    low_bin_id=0,
):
    required_cols = {"dmr_id", "sample_group", "bin_id", "read_count", "total_reads", "fraction"}
    missing = required_cols - set(bin_fractions.columns)
    if missing:
        raise ValueError(f"Input TSV is missing required columns: {sorted(missing)}")

    bin_fractions = bin_fractions.copy()
    bin_fractions["bin_id"] = bin_fractions["bin_id"].astype(int)

    if n_bins is None:
        n_bins = int(bin_fractions["bin_id"].max()) + 1

    bin_labels = make_5pct_bin_labels(n_bins)

    # Rank DMRs by the percentage of uRPL reads in the 0-5% methylation bin,
    # rather than by the raw number of reads in that bin.
    low = bin_fractions[bin_fractions["bin_id"] == low_bin_id].copy()

    low_wide = (
        low
        .pivot_table(
            index="dmr_id",
            columns="sample_group",
            values="fraction",
            aggfunc="first",
        )
        .fillna(0.0)
    )

    if "patient" not in low_wide.columns:
        low_wide["patient"] = 0.0
    if "control" not in low_wide.columns:
        low_wide["control"] = 0.0

    low_wide["urpl_low_read_percent"] = low_wide["patient"] * 100

    ordered_dmrs = (
        low_wide
        .sort_values(
            "urpl_low_read_percent",
            ascending=False,
        )
        .index
        .tolist()
    )

    full_index = pd.MultiIndex.from_product(
        [ordered_dmrs, ["control", "patient"], range(n_bins)],
        names=["dmr_id", "sample_group", "bin_id"],
    )

    bf = (
        bin_fractions
        .set_index(["dmr_id", "sample_group", "bin_id"])
        .reindex(full_index)
        .reset_index()
    )

    bf["fraction"] = bf["fraction"].fillna(0.0)
    bf["read_count"] = bf["read_count"].fillna(0)

    control = (
        bf[bf["sample_group"] == "control"]
        .pivot(index="dmr_id", columns="bin_id", values="fraction")
        .reindex(ordered_dmrs)
        .fillna(0.0)
    )

    patient = (
        bf[bf["sample_group"] == "patient"]
        .pivot(index="dmr_id", columns="bin_id", values="fraction")
        .reindex(ordered_dmrs)
        .fillna(0.0)
    )

    diff = patient - control

    # Mean per-DMR bin fractions for the top bar plots.
    # These bars summarize the heatmap rows directly.
    control_bar = (
        control
        .mean(axis=0)
        .reindex(range(n_bins))
        .fillna(0.0)
        .values
    )

    patient_bar = (
        patient
        .mean(axis=0)
        .reindex(range(n_bins))
        .fillna(0.0)
        .values
    )

    diff_bar = patient_bar - control_bar

    n_dmrs = len(ordered_dmrs)
    fig_height = max(8, min(18, n_dmrs * 0.035 + 4))

    fig = plt.figure(figsize=(15, fig_height))

    gs = fig.add_gridspec(
        2,
        3,
        height_ratios=[1.0, 6.0],
        width_ratios=[1, 1, 1.1],
        hspace=0.08,
        wspace=0.40,
    )

    bar_axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[0, 2]),
    ]

    heat_axes = [
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
        fig.add_subplot(gs[1, 2]),
    ]

    diff_abs_max = np.nanmax(np.abs(diff.values))
    if not np.isfinite(diff_abs_max) or diff_abs_max == 0:
        diff_abs_max = 1.0

    bar_abs_max = np.nanmax(np.abs(diff_bar))
    if not np.isfinite(bar_abs_max) or bar_abs_max == 0:
        bar_abs_max = 1.0

    fraction_bar_max = max(
        np.nanmax(control_bar),
        np.nanmax(patient_bar),
    )

    if not np.isfinite(fraction_bar_max) or fraction_bar_max == 0:
        fraction_bar_max = 1.0

    x = np.arange(n_bins)

    # Top bar plots: mean per-DMR bin fraction.
    bar_axes[0].bar(x, control_bar, width=0.85)
    bar_axes[0].set_title("Control fraction")
    bar_axes[0].set_ylabel("Mean\nfraction")
    bar_axes[0].set_ylim(0, fraction_bar_max * 1.1)
    bar_axes[0].grid(axis="y", alpha=0.25)

    bar_axes[1].bar(x, patient_bar, width=0.85)
    bar_axes[1].set_title("uRPL fraction")
    bar_axes[1].set_ylabel("Mean\nfraction")
    bar_axes[1].set_ylim(0, fraction_bar_max * 1.1)
    bar_axes[1].grid(axis="y", alpha=0.25)

    bar_axes[2].bar(x, diff_bar, width=0.85)
    bar_axes[2].axhline(0, color="black", linewidth=0.8, alpha=0.6)
    bar_axes[2].set_title("uRPL - Control")
    bar_axes[2].set_ylabel("Delta mean\nfraction", labelpad=2)
    bar_axes[2].set_ylim(-bar_abs_max * 1.1, bar_abs_max * 1.1)
    bar_axes[2].grid(axis="y", alpha=0.25)

    for ax in bar_axes:
        ax.set_xlim(-0.5, n_bins - 0.5)
        ax.set_xticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Heatmaps.
    im0 = heat_axes[0].imshow(
        control.values,
        aspect="auto",
        interpolation="nearest",
        cmap="viridis",
        vmin=0,
        vmax=1,
    )

    im1 = heat_axes[1].imshow(
        patient.values,
        aspect="auto",
        interpolation="nearest",
        cmap="viridis",
        vmin=0,
        vmax=1,
    )

    im2 = heat_axes[2].imshow(
        diff.values,
        aspect="auto",
        interpolation="nearest",
        cmap="RdBu_r",
        vmin=-diff_abs_max,
        vmax=diff_abs_max,
    )

    for ax in heat_axes:
        ax.set_xlabel("Read m_frac bin")
        ax.set_xticks(range(n_bins))
        ax.set_xticklabels(bin_labels, rotation=90, fontsize=8)

    heat_axes[0].set_ylabel("DMRs")

    if n_dmrs <= 60:
        heat_axes[0].set_yticks(range(n_dmrs))
        heat_axes[0].set_yticklabels(ordered_dmrs, fontsize=6)
    else:
        heat_axes[0].set_yticks([])

    heat_axes[1].tick_params(axis="y", left=False, labelleft=False)
    heat_axes[2].tick_params(axis="y", left=False, labelleft=False)

    cbar0 = fig.colorbar(im0, ax=heat_axes[0], fraction=0.046, pad=0.02)
    cbar0.set_label("Read fraction")

    cbar1 = fig.colorbar(im1, ax=heat_axes[1], fraction=0.046, pad=0.02)
    cbar1.set_label("Read fraction")

    cbar2 = fig.colorbar(im2, ax=heat_axes[2], fraction=0.046, pad=0.02)
    cbar2.set_label("Fraction difference")

    fig.suptitle(
        "Read-level methylation-bin fractions per DMR with mean per-DMR bin summaries",
        y=0.995,
        fontsize=14,
    )

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Remake the three-panel DMR bin-fraction heatmap from dmr_5pct_bin_fractions.tsv, "
            "using mean per-DMR bin fraction bars above each panel."
        )
    )
    ap.add_argument("--bin-fractions", required=True, help="TSV with dmr_id, sample_group, bin_id, read_count, total_reads, fraction")
    ap.add_argument("--out-png", required=True)
    ap.add_argument("--n-bins", type=int, default=None)
    ap.add_argument(
        "--low-bin-id",
        type=int,
        default=0,
        help="Bin used to sort DMRs by percent of uRPL reads. Default 0 = 0-5 percent bin.",
    )
    args = ap.parse_args()

    bin_fractions = pd.read_csv(args.bin_fractions, sep="\t")
    plot_three_panel_heatmap_from_bin_fractions(
        bin_fractions=bin_fractions,
        out_png=args.out_png,
        n_bins=args.n_bins,
        low_bin_id=args.low_bin_id,
    )

    print("Done.")
    print(f"Wrote: {args.out_png}")


if __name__ == "__main__":
    main()