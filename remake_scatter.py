#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


GROUP_ALIASES = {
    "control": "control",
    "patient": "patient",
    "urpl": "patient",
}


def load_variance_table(path, value_column, table_name):
    frame = pd.read_csv(path, sep="\t")
    required = {"dmr_id", "sample_group", value_column}
    missing = required - set(frame.columns)
    if missing:
        raise SystemExit(
            f"{table_name} is missing required columns: {sorted(missing)}"
        )

    normalized = frame["sample_group"].astype(str).str.strip().str.lower()
    unknown = sorted(set(normalized) - set(GROUP_ALIASES))
    if unknown:
        raise SystemExit(
            f"{table_name} has unrecognized sample_group values: {unknown}"
        )

    frame = frame.copy()
    frame["sample_group"] = normalized.map(GROUP_ALIASES)
    frame[value_column] = pd.to_numeric(frame[value_column], errors="coerce")
    frame = frame.dropna(subset=["dmr_id", "sample_group", value_column])
    return frame


def make_paired_table(frame, value_column, aggregation, table_name):
    paired = (
        frame.groupby(["dmr_id", "sample_group"], as_index=False)[value_column]
        .agg(aggregation)
        .pivot(index="dmr_id", columns="sample_group", values=value_column)
    )

    if not {"control", "patient"}.issubset(paired.columns):
        raise SystemExit(
            f"{table_name} does not contain paired Control and uRPL values."
        )

    paired = paired[["control", "patient"]].dropna().reset_index()
    paired = paired.rename(
        columns={"control": "Control", "patient": "uRPL"}
    )
    paired = paired[
        (paired["Control"] >= 0)
        & (paired["uRPL"] >= 0)
    ].copy()

    if paired.empty:
        raise SystemExit(f"No paired DMR values remain for {table_name}.")
    return paired


def plot_variance_scatter(
    paired,
    out_png,
    title,
    axis_metric,
    point_color,
    edge_color,
    point_alpha,
    point_size,
):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(
        paired["Control"],
        paired["uRPL"],
        s=point_size,
        alpha=point_alpha,
        color=point_color,
        edgecolors=edge_color,
        linewidth=0.3,
    )

    observed_max = max(paired["Control"].max(), paired["uRPL"].max())
    axis_max = observed_max * 1.05 if observed_max > 0 else 1.0
    ax.plot(
        [0, axis_max],
        [0, axis_max],
        linestyle="--",
        linewidth=1.25,
        color="#333333",
        alpha=0.75,
        label="uRPL = Control",
        zorder=0,
    )

    ax.set(xlim=(0, axis_max), ylim=(0, axis_max))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"Control {axis_metric}")
    ax.set_ylabel(f"uRPL {axis_metric}")
    ax.set_title(title)
    ax.grid(alpha=0.18, linewidth=0.7)
    ax.legend(frameon=False, loc="upper left")
    ax.text(
        0.98,
        0.02,
        f"n = {len(paired):,} DMRs",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
    )

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate only the inter-sample and intra-sample DMR variance "
            "scatter plots with darker, more opaque points."
        )
    )
    parser.add_argument("--inter-tsv", required=True)
    parser.add_argument("--intra-tsv", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--inter-value-column",
        default="inter_sample_var_mean_read_m_frac",
    )
    parser.add_argument(
        "--intra-value-column",
        default="intra_sample_var_read_m_frac",
    )
    parser.add_argument(
        "--intra-aggregation",
        choices=["mean", "median"],
        default="mean",
        help="How to combine per-sample intra-sample variances within each group.",
    )
    parser.add_argument(
        "--inter-color",
        default="#117a65",
        help="Point color for the inter-sample plot. Default: dark teal.",
    )
    parser.add_argument(
        "--intra-color",
        default="#c44e52",
        help="Point color for the intra-sample plot. Default: dark salmon.",
    )
    parser.add_argument(
        "--point-alpha",
        type=float,
        default=0.82,
        help="Point opacity from 0 to 1. Default: 0.82.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=30,
        help="Scatter-point area. Default: 30.",
    )
    args = parser.parse_args()

    if not 0 < args.point_alpha <= 1:
        raise SystemExit("--point-alpha must be in (0, 1].")
    if args.point_size <= 0:
        raise SystemExit("--point-size must be greater than zero.")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    inter = load_variance_table(
        args.inter_tsv,
        args.inter_value_column,
        "Inter-sample TSV",
    )
    intra = load_variance_table(
        args.intra_tsv,
        args.intra_value_column,
        "Intra-sample TSV",
    )

    inter_paired = make_paired_table(
        inter,
        value_column=args.inter_value_column,
        aggregation="mean",
        table_name="Inter-sample TSV",
    )
    intra_paired = make_paired_table(
        intra,
        value_column=args.intra_value_column,
        aggregation=args.intra_aggregation,
        table_name="Intra-sample TSV",
    )

    inter_paired.to_csv(
        out_dir / "inter_sample_variance_scatter_data.tsv",
        sep="\t",
        index=False,
    )
    intra_paired.to_csv(
        out_dir / "intra_sample_variance_scatter_data.tsv",
        sep="\t",
        index=False,
    )

    plot_variance_scatter(
        inter_paired,
        out_png=out_dir / "inter_sample_variance_scatter.png",
        title="Inter-sample variance by DMR",
        axis_metric="inter-sample variance",
        point_color=args.inter_color,
        edge_color="#0b4f42",
        point_alpha=args.point_alpha,
        point_size=args.point_size,
    )
    plot_variance_scatter(
        intra_paired,
        out_png=out_dir / "intra_sample_variance_scatter.png",
        title=(
            "Intra-sample variance by DMR "
            f"({args.intra_aggregation} across samples)"
        ),
        axis_metric="intra-sample variance",
        point_color=args.intra_color,
        edge_color="#7f3035",
        point_alpha=args.point_alpha,
        point_size=args.point_size,
    )

    print(f"Wrote: {out_dir / 'inter_sample_variance_scatter.png'}")
    print(f"Wrote: {out_dir / 'intra_sample_variance_scatter.png'}")
    print(f"Inter-sample paired DMRs: {len(inter_paired):,}")
    print(f"Intra-sample paired DMRs: {len(intra_paired):,}")


if __name__ == "__main__":
    main()