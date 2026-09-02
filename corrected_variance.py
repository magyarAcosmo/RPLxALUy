#!/usr/bin/env python3

import argparse
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter


GROUP_COLORS = {"control": "#1b9e77", "patient": "#d95f02"}
GROUP_LABELS = {"control": "Control", "patient": "uRPL"}


def sql_quote(value):
    return "'" + str(value).replace("'", "''") + "'"


def sql_list(paths):
    return "[" + ",".join(sql_quote(path) for path in paths) + "]"


def collect_ch3_files(folder, pattern):
    files = sorted(Path(folder).glob(pattern))
    if not files:
        raise SystemExit(f"No files matching {pattern!r} found in {folder}")
    return [str(path) for path in files]


def load_dmrs(path):
    dmrs = pd.read_csv(path, sep=r"\s+", header=None, comment="#", dtype={0: str})
    if dmrs.shape[1] < 3:
        raise SystemExit("DMR file must have at least 3 columns: chrom start end")

    dmrs = dmrs.iloc[:, :4].copy()
    if dmrs.shape[1] == 3:
        dmrs.columns = ["chrom", "dmr_start", "dmr_end"]
        dmrs["dmr_id"] = (
            dmrs["chrom"].astype(str)
            + ":"
            + dmrs["dmr_start"].astype(str)
            + "-"
            + dmrs["dmr_end"].astype(str)
        )
    else:
        dmrs.columns = ["chrom", "dmr_start", "dmr_end", "dmr_id"]
        generated = (
            dmrs["chrom"].astype(str)
            + ":"
            + dmrs["dmr_start"].astype(str)
            + "-"
            + dmrs["dmr_end"].astype(str)
        )
        dmrs["dmr_id"] = dmrs["dmr_id"].fillna(generated)

    dmrs["dmr_start"] = dmrs["dmr_start"].astype("int64")
    dmrs["dmr_end"] = dmrs["dmr_end"].astype("int64")
    return dmrs[["chrom", "dmr_start", "dmr_end", "dmr_id"]]


def create_input_tables(
    con,
    patient_files,
    control_files,
    dmrs,
    methylated_sql,
    min_cpgs_per_read,
):
    con.register("dmrs_df", dmrs)
    con.execute(
        f"""
        CREATE TEMP TABLE dmrs AS
        SELECT
            chrom::VARCHAR AS chrom,
            dmr_start::BIGINT AS dmr_start,
            dmr_end::BIGINT AS dmr_end,
            dmr_id::VARCHAR AS dmr_id
        FROM dmrs_df;

        CREATE TEMP VIEW ch3 AS
        SELECT
            'patient' AS sample_group,
            regexp_replace(
                regexp_extract(filename, '[^/]+$'),
                '\\.ch3$',
                ''
            ) AS sample_name,
            * EXCLUDE (filename, sample_name)
        FROM read_parquet(
            {sql_list(patient_files)},
            union_by_name=true,
            filename=true
        )

        UNION ALL

        SELECT
            'control' AS sample_group,
            regexp_replace(
                regexp_extract(filename, '[^/]+$'),
                '\\.ch3$',
                ''
            ) AS sample_name,
            * EXCLUDE (filename, sample_name)
        FROM read_parquet(
            {sql_list(control_files)},
            union_by_name=true,
            filename=true
        );

        CREATE TEMP TABLE read_level AS
        SELECT
            c.sample_group,
            c.sample_name,
            d.dmr_id,
            d.chrom,
            d.dmr_start,
            d.dmr_end,
            c.read_id,
            COUNT(*) AS n_cpg_calls,
            SUM(
                CASE WHEN c.call_code IN {methylated_sql} THEN 1 ELSE 0 END
            ) AS n_methylated_calls,
            SUM(
                CASE WHEN c.call_code IN {methylated_sql} THEN 1 ELSE 0 END
            )::DOUBLE / COUNT(*) AS read_m_frac
        FROM ch3 c
        JOIN dmrs d
          ON c.chrom = d.chrom
         AND c.start >= d.dmr_start
         AND c.start < d.dmr_end
        WHERE
            (LENGTH(c.query_kmer) = 2 AND c.query_kmer = 'CG')
            OR
            (LENGTH(c.query_kmer) = 5 AND SUBSTR(c.query_kmer, 3, 2) = 'CG')
        GROUP BY
            c.sample_group,
            c.sample_name,
            d.dmr_id,
            d.chrom,
            d.dmr_start,
            d.dmr_end,
            c.read_id
        HAVING COUNT(*) >= {min_cpgs_per_read};
        """
    )

    observed = dict(
        con.execute(
            """
            SELECT sample_group, COUNT(DISTINCT sample_name) AS n_samples
            FROM ch3
            GROUP BY sample_group
            """
        ).fetchall()
    )
    expected = {"patient": len(patient_files), "control": len(control_files)}
    if observed != expected:
        raise SystemExit(
            "Filename-derived sample count did not match the number of input files. "
            f"Expected {expected}; observed {observed}."
        )


def calculate_variance_tables(con, min_reads_per_sample_dmr):
    intra = con.execute(
        f"""
        SELECT
            sample_group,
            sample_name,
            dmr_id,
            chrom,
            dmr_start,
            dmr_end,
            COUNT(*) AS n_reads,
            SUM(n_cpg_calls) AS n_cpg_calls,
            AVG(read_m_frac) AS mean_read_m_frac,
            VAR_SAMP(read_m_frac) AS intra_sample_var_read_m_frac,
            STDDEV_SAMP(read_m_frac) AS intra_sample_sd_read_m_frac
        FROM read_level
        GROUP BY
            sample_group, sample_name, dmr_id, chrom, dmr_start, dmr_end
        HAVING COUNT(*) >= {min_reads_per_sample_dmr}
        ORDER BY sample_group, sample_name, chrom, dmr_start, dmr_end
        """
    ).df()

    con.execute(
        f"""
        CREATE TEMP TABLE sample_level AS
        SELECT
            sample_group,
            sample_name,
            dmr_id,
            chrom,
            dmr_start,
            dmr_end,
            COUNT(*) AS n_reads,
            SUM(n_cpg_calls) AS n_cpg_calls,
            SUM(n_methylated_calls) AS n_methylated_calls,
            SUM(n_methylated_calls)::DOUBLE / SUM(n_cpg_calls) AS pooled_m_frac,
            AVG(read_m_frac) AS mean_read_m_frac
        FROM read_level
        GROUP BY
            sample_group, sample_name, dmr_id, chrom, dmr_start, dmr_end
        HAVING COUNT(*) >= {min_reads_per_sample_dmr};
        """
    )

    inter = con.execute(
        """
        SELECT
            sample_group,
            dmr_id,
            chrom,
            dmr_start,
            dmr_end,
            COUNT(*) AS n_samples,
            AVG(pooled_m_frac) AS group_mean_pooled_m_frac,
            VAR_SAMP(pooled_m_frac) AS inter_sample_var_pooled_m_frac,
            STDDEV_SAMP(pooled_m_frac) AS inter_sample_sd_pooled_m_frac,
            AVG(mean_read_m_frac) AS group_mean_read_m_frac,
            VAR_SAMP(mean_read_m_frac) AS inter_sample_var_mean_read_m_frac,
            STDDEV_SAMP(mean_read_m_frac) AS inter_sample_sd_mean_read_m_frac
        FROM sample_level
        GROUP BY sample_group, dmr_id, chrom, dmr_start, dmr_end
        HAVING COUNT(*) >= 2
        ORDER BY sample_group, chrom, dmr_start, dmr_end
        """
    ).df()
    return intra, inter


def calculate_subpopulation_tables(
    con,
    n_bins,
    low_threshold,
    min_reads_per_dmr_group,
    top_n,
):
    bin_fractions = con.execute(
        f"""
        WITH binned AS (
            SELECT
                dmr_id,
                sample_group,
                CASE
                    WHEN read_m_frac >= 1.0 THEN {n_bins - 1}
                    WHEN read_m_frac <= 0.0 THEN 0
                    ELSE FLOOR(read_m_frac * {n_bins})::INTEGER
                END AS bin_id
            FROM read_level
        ),
        counts AS (
            SELECT dmr_id, sample_group, bin_id, COUNT(*) AS read_count
            FROM binned
            GROUP BY dmr_id, sample_group, bin_id
        ),
        totals AS (
            SELECT dmr_id, sample_group, SUM(read_count) AS total_reads
            FROM counts
            GROUP BY dmr_id, sample_group
        )
        SELECT
            c.dmr_id,
            c.sample_group,
            c.bin_id,
            (c.bin_id + 0.5) / {n_bins} AS bin_center,
            c.read_count,
            t.total_reads,
            c.read_count::DOUBLE / t.total_reads AS fraction
        FROM counts c
        JOIN totals t USING (dmr_id, sample_group)
        ORDER BY c.dmr_id, c.sample_group, c.bin_id
        """
    ).df()

    dmr_summary = con.execute(
        f"""
        WITH group_summary AS (
            SELECT
                dmr_id,
                sample_group,
                COUNT(*) AS n_reads,
                AVG(read_m_frac) AS mean_read_m_frac,
                MEDIAN(read_m_frac) AS median_read_m_frac,
                SUM(CASE WHEN read_m_frac <= {low_threshold} THEN 1 ELSE 0 END)
                    ::DOUBLE / COUNT(*) AS frac_low_methylation_reads,
                SUM(CASE WHEN read_m_frac = 0 THEN 1 ELSE 0 END)
                    ::DOUBLE / COUNT(*) AS frac_zero_methylation_reads,
                SUM(CASE WHEN read_m_frac >= 0.90 THEN 1 ELSE 0 END)
                    ::DOUBLE / COUNT(*) AS frac_high_methylation_reads
            FROM read_level
            GROUP BY dmr_id, sample_group
        ),
        wide AS (
            SELECT
                dmr_id,
                MAX(CASE WHEN sample_group = 'patient' THEN n_reads END)
                    AS patient_n_reads,
                MAX(CASE WHEN sample_group = 'control' THEN n_reads END)
                    AS control_n_reads,
                MAX(CASE WHEN sample_group = 'patient' THEN mean_read_m_frac END)
                    AS patient_mean_read_m_frac,
                MAX(CASE WHEN sample_group = 'control' THEN mean_read_m_frac END)
                    AS control_mean_read_m_frac,
                MAX(CASE WHEN sample_group = 'patient' THEN median_read_m_frac END)
                    AS patient_median_read_m_frac,
                MAX(CASE WHEN sample_group = 'control' THEN median_read_m_frac END)
                    AS control_median_read_m_frac,
                MAX(CASE WHEN sample_group = 'patient' THEN frac_low_methylation_reads END)
                    AS patient_frac_low_methylation_reads,
                MAX(CASE WHEN sample_group = 'control' THEN frac_low_methylation_reads END)
                    AS control_frac_low_methylation_reads,
                MAX(CASE WHEN sample_group = 'patient' THEN frac_zero_methylation_reads END)
                    AS patient_frac_zero_methylation_reads,
                MAX(CASE WHEN sample_group = 'control' THEN frac_zero_methylation_reads END)
                    AS control_frac_zero_methylation_reads,
                MAX(CASE WHEN sample_group = 'patient' THEN frac_high_methylation_reads END)
                    AS patient_frac_high_methylation_reads,
                MAX(CASE WHEN sample_group = 'control' THEN frac_high_methylation_reads END)
                    AS control_frac_high_methylation_reads
            FROM group_summary
            GROUP BY dmr_id
        )
        SELECT
            dmr_id,
            patient_n_reads,
            control_n_reads,
            patient_mean_read_m_frac,
            control_mean_read_m_frac,
            patient_mean_read_m_frac - control_mean_read_m_frac
                AS mean_read_m_frac_patient_minus_control,
            patient_median_read_m_frac,
            control_median_read_m_frac,
            patient_median_read_m_frac - control_median_read_m_frac
                AS median_read_m_frac_patient_minus_control,
            patient_frac_low_methylation_reads,
            control_frac_low_methylation_reads,
            patient_frac_low_methylation_reads - control_frac_low_methylation_reads
                AS low_methylation_excess_patient_minus_control,
            patient_frac_zero_methylation_reads,
            control_frac_zero_methylation_reads,
            patient_frac_zero_methylation_reads - control_frac_zero_methylation_reads
                AS zero_methylation_excess_patient_minus_control,
            patient_frac_high_methylation_reads,
            control_frac_high_methylation_reads,
            patient_frac_high_methylation_reads - control_frac_high_methylation_reads
                AS high_methylation_excess_patient_minus_control
        FROM wide
        WHERE patient_n_reads IS NOT NULL AND control_n_reads IS NOT NULL
        ORDER BY low_methylation_excess_patient_minus_control DESC
        """
    ).df()

    eligible = dmr_summary[
        (dmr_summary["patient_n_reads"] >= min_reads_per_dmr_group)
        & (dmr_summary["control_n_reads"] >= min_reads_per_dmr_group)
    ].copy()
    if eligible.empty:
        raise SystemExit(
            "No DMRs passed --min-reads-per-dmr-group in both groups."
        )

    top = eligible.nlargest(
        top_n, "low_methylation_excess_patient_minus_control"
    ).copy()
    top["dmr_rank_label"] = [f"top {i}" for i in range(1, len(top) + 1)]
    bottom = eligible.nsmallest(
        top_n, "low_methylation_excess_patient_minus_control"
    ).copy()
    bottom["dmr_rank_label"] = [f"bottom {i}" for i in range(1, len(bottom) + 1)]
    selected_summary = pd.concat([bottom, top], ignore_index=True)

    con.register(
        "selected_dmrs_df",
        selected_summary[["dmr_id"]].drop_duplicates(),
    )
    selected_read_level = con.execute(
        """
        SELECT
            r.sample_group,
            r.sample_name,
            r.dmr_id,
            r.read_id,
            r.n_cpg_calls,
            r.n_methylated_calls,
            r.read_m_frac
        FROM read_level r
        JOIN selected_dmrs_df s USING (dmr_id)
        ORDER BY r.dmr_id, r.sample_group, r.sample_name, r.read_m_frac
        """
    ).df()
    return bin_fractions, dmr_summary, selected_summary, selected_read_level


def make_5pct_bin_labels(n_bins):
    if n_bins == 20:
        return [
            "0-5%" if i == 0 else "96-100%" if i == 19
            else f"{i * 5 + 1}-{(i + 1) * 5}%"
            for i in range(n_bins)
        ]
    return [
        f"{i / n_bins * 100:.0f}-{(i + 1) / n_bins * 100:.0f}%"
        for i in range(n_bins)
    ]


def plot_selected_histograms(read_level, dmr_order, out_png, bins):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, dmr_id in zip(axes, dmr_order):
        sub = read_level[read_level["dmr_id"] == dmr_id]
        for group in ["control", "patient"]:
            values = sub.loc[
                sub["sample_group"] == group, "read_m_frac"
            ].dropna().to_numpy()
            if len(values):
                ax.hist(
                    values,
                    bins=bins,
                    range=(0, 1),
                    weights=np.full(len(values), 100.0 / len(values)),
                    histtype="step",
                    linewidth=2,
                    alpha=0.95,
                    label=GROUP_LABELS[group],
                    color=GROUP_COLORS[group],
                )
        ax.set_title(dmr_id)
        ax.set_xlabel("Read-level m_frac")
        ax.set_xlim(0, 1)
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
        ax.grid(axis="y", alpha=0.25)
        ax.legend(frameon=False)

    axes[0].set_ylabel("Percent of reads")
    fig.suptitle("Read-level m_frac distributions for selected DMRs", y=1.02, fontsize=14)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def make_paired_variance_table(df, value_column, aggregation):
    paired = (
        df.groupby(["dmr_id", "sample_group"], as_index=False)[value_column]
        .agg(aggregation)
        .pivot(index="dmr_id", columns="sample_group", values=value_column)
    )
    if not {"control", "patient"}.issubset(paired.columns):
        raise SystemExit(f"Cannot pair Control and uRPL values for {value_column}.")
    paired = paired[["control", "patient"]].dropna().reset_index()
    return paired.rename(columns={"control": "Control", "patient": "uRPL"})


def plot_variance_scatter(paired, out_png, title, axis_metric):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(
        paired["Control"], paired["uRPL"],
        s=24, alpha=0.55, color="#4c78a8", edgecolors="none",
    )
    observed_max = max(paired["Control"].max(), paired["uRPL"].max())
    axis_max = observed_max * 1.05 if observed_max > 0 else 1.0
    ax.plot(
        [0, axis_max], [0, axis_max], "--",
        linewidth=1.2, color="black", alpha=0.6, label="uRPL = Control",
    )
    ax.set(xlim=(0, axis_max), ylim=(0, axis_max))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"Control {axis_metric}")
    ax.set_ylabel(f"uRPL {axis_metric}")
    ax.set_title(title)
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, loc="upper left")
    ax.text(
        0.98, 0.02, f"n = {len(paired):,} DMRs",
        transform=ax.transAxes, ha="right", va="bottom",
    )
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def plot_dmr_heatmap(bin_fractions, out_png, n_bins, low_bin_id=0):
    bin_fractions = bin_fractions.copy()
    bin_fractions["bin_id"] = bin_fractions["bin_id"].astype(int)
    all_dmrs = pd.Index(bin_fractions["dmr_id"].drop_duplicates())
    low = bin_fractions[bin_fractions["bin_id"] == low_bin_id]
    low_wide = low.pivot_table(
        index="dmr_id", columns="sample_group", values="fraction", aggfunc="first"
    ).reindex(all_dmrs).fillna(0.0)
    for group in ["control", "patient"]:
        if group not in low_wide:
            low_wide[group] = 0.0
    ordered_dmrs = low_wide.sort_values("patient", ascending=False).index.tolist()

    full_index = pd.MultiIndex.from_product(
        [ordered_dmrs, ["control", "patient"], range(n_bins)],
        names=["dmr_id", "sample_group", "bin_id"],
    )
    bf = (
        bin_fractions.set_index(["dmr_id", "sample_group", "bin_id"])
        .reindex(full_index).reset_index()
    )
    bf["fraction"] = bf["fraction"].fillna(0.0)

    matrices = {}
    for group in ["control", "patient"]:
        matrices[group] = (
            bf[bf["sample_group"] == group]
            .pivot(index="dmr_id", columns="bin_id", values="fraction")
            .reindex(ordered_dmrs).fillna(0.0)
        )
    difference = matrices["patient"] - matrices["control"]
    bars = {
        group: matrices[group].mean(axis=0).reindex(range(n_bins)).fillna(0).to_numpy()
        for group in ["control", "patient"]
    }
    difference_bar = bars["patient"] - bars["control"]

    diff_max = np.nanmax(np.abs(difference.to_numpy()))
    diff_max = diff_max if np.isfinite(diff_max) and diff_max > 0 else 1.0
    bar_diff_max = np.nanmax(np.abs(difference_bar))
    bar_diff_max = bar_diff_max if np.isfinite(bar_diff_max) and bar_diff_max > 0 else 1.0
    fraction_max = max(np.nanmax(bars["control"]), np.nanmax(bars["patient"]))
    fraction_max = fraction_max if np.isfinite(fraction_max) and fraction_max > 0 else 1.0

    n_dmrs = len(ordered_dmrs)
    fig = plt.figure(figsize=(15, max(8, min(18, n_dmrs * 0.035 + 4))))
    grid = fig.add_gridspec(
        2, 3, height_ratios=[1, 6], width_ratios=[1, 1, 1.1],
        hspace=0.08, wspace=0.40,
    )
    bar_axes = [fig.add_subplot(grid[0, i]) for i in range(3)]
    heat_axes = [fig.add_subplot(grid[1, i]) for i in range(3)]
    x = np.arange(n_bins)

    bar_axes[0].bar(x, bars["control"], width=0.85)
    bar_axes[0].set_title("Control fraction")
    bar_axes[0].set_ylabel("Mean\nfraction")
    bar_axes[0].set_ylim(0, fraction_max * 1.1)
    bar_axes[1].bar(x, bars["patient"], width=0.85)
    bar_axes[1].set_title("uRPL fraction")
    bar_axes[1].set_ylabel("Mean\nfraction")
    bar_axes[1].set_ylim(0, fraction_max * 1.1)
    bar_axes[2].bar(x, difference_bar, width=0.85)
    bar_axes[2].axhline(0, color="black", linewidth=0.8, alpha=0.6)
    bar_axes[2].set_title("uRPL - Control")
    bar_axes[2].set_ylabel("Delta mean\nfraction", labelpad=2)
    bar_axes[2].set_ylim(-bar_diff_max * 1.1, bar_diff_max * 1.1)
    for ax in bar_axes:
        ax.set_xlim(-0.5, n_bins - 0.5)
        ax.set_xticks([])
        ax.grid(axis="y", alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    im0 = heat_axes[0].imshow(
        matrices["control"].to_numpy(), aspect="auto", interpolation="nearest",
        cmap="viridis", vmin=0, vmax=1,
    )
    im1 = heat_axes[1].imshow(
        matrices["patient"].to_numpy(), aspect="auto", interpolation="nearest",
        cmap="viridis", vmin=0, vmax=1,
    )
    im2 = heat_axes[2].imshow(
        difference.to_numpy(), aspect="auto", interpolation="nearest",
        cmap="RdBu_r", vmin=-diff_max, vmax=diff_max,
    )
    labels = make_5pct_bin_labels(n_bins)
    for ax in heat_axes:
        ax.set_xlabel("Read m_frac bin")
        ax.set_xticks(range(n_bins))
        ax.set_xticklabels(labels, rotation=90, fontsize=8)
    heat_axes[0].set_ylabel("DMRs")
    if n_dmrs <= 60:
        heat_axes[0].set_yticks(range(n_dmrs))
        heat_axes[0].set_yticklabels(ordered_dmrs, fontsize=6)
    else:
        heat_axes[0].set_yticks([])
    heat_axes[1].tick_params(axis="y", left=False, labelleft=False)
    heat_axes[2].tick_params(axis="y", left=False, labelleft=False)
    fig.colorbar(im0, ax=heat_axes[0], fraction=0.046, pad=0.02).set_label("Read fraction")
    fig.colorbar(im1, ax=heat_axes[1], fraction=0.046, pad=0.02).set_label("Read fraction")
    fig.colorbar(im2, ax=heat_axes[2], fraction=0.046, pad=0.02).set_label("Fraction difference")
    fig.suptitle(
        "Read-level methylation-bin fractions per DMR with mean per-DMR bin summaries",
        y=0.995, fontsize=14,
    )
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def plot_violin(read_level, selected_summary, out_png):
    plot_order = selected_summary["dmr_id"].tolist()
    positions = np.arange(1, len(plot_order) + 1)
    control_data, urpl_data = [], []
    for dmr_id in plot_order:
        sub = read_level[read_level["dmr_id"] == dmr_id]
        control_data.append(
            sub.loc[sub["sample_group"] == "control", "read_m_frac"].dropna().to_numpy()
        )
        urpl_data.append(
            sub.loc[sub["sample_group"] == "patient", "read_m_frac"].dropna().to_numpy()
        )

    fig, ax = plt.subplots(figsize=(max(12, len(plot_order) * 0.55), 6))
    control_violin = ax.violinplot(
        control_data, positions=positions, widths=0.75,
        showmeans=False, showmedians=True, showextrema=False,
    )
    urpl_violin = ax.violinplot(
        urpl_data, positions=positions, widths=0.75,
        showmeans=False, showmedians=True, showextrema=False,
    )
    for violin, group in [(control_violin, "control"), (urpl_violin, "patient")]:
        for body in violin["bodies"]:
            body.set_facecolor(GROUP_COLORS[group])
            body.set_edgecolor(GROUP_COLORS[group])
            body.set_alpha(0.35)
        violin["cmedians"].set_color(GROUP_COLORS[group])
        violin["cmedians"].set_linewidth(1.2)

    labels = [
        f"{row.dmr_rank_label}\n{row.dmr_id}\n"
        f"Δlow={row.low_methylation_excess_patient_minus_control:.2f}"
        for row in selected_summary.itertuples(index=False)
    ]
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
    ax.legend(
        handles=[
            Patch(
                facecolor=GROUP_COLORS[group], edgecolor=GROUP_COLORS[group],
                alpha=0.35, label=GROUP_LABELS[group],
            )
            for group in ["control", "patient"]
        ],
        loc="upper right",
    )
    n_bottom = selected_summary["dmr_rank_label"].str.startswith("bottom").sum()
    if 0 < n_bottom < len(plot_order):
        ax.axvline(n_bottom + 0.5, color="black", linewidth=1, alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate corrected DMR variance/subpopulation TSVs and all final "
            "histogram, scatter, heatmap, and violin figures in one run."
        )
    )
    parser.add_argument("--patients", required=True)
    parser.add_argument("--controls", required=True)
    parser.add_argument("--dmrs", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--pattern", default="*.ch3")
    parser.add_argument("--methylated-codes", default="m")
    parser.add_argument("--min-cpgs-per-read", type=int, default=1)
    parser.add_argument("--min-reads-per-sample-dmr", type=int, default=1)
    parser.add_argument("--min-reads-per-dmr-group", type=int, default=20)
    parser.add_argument("--heatmap-bins", type=int, default=20)
    parser.add_argument("--low-threshold", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--histogram-bins", type=int, default=20)
    parser.add_argument(
        "--dmr-ids", nargs=2, required=True,
        help="Two DMR IDs for the percent-of-reads histogram figure.",
    )
    parser.add_argument(
        "--inter-value-column", default="inter_sample_var_mean_read_m_frac"
    )
    parser.add_argument(
        "--intra-value-column", default="intra_sample_var_read_m_frac"
    )
    parser.add_argument(
        "--intra-aggregation", choices=["mean", "median"], default="mean"
    )
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    figure_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    patient_files = collect_ch3_files(args.patients, args.pattern)
    control_files = collect_ch3_files(args.controls, args.pattern)
    dmrs = load_dmrs(args.dmrs)

    missing_dmr_ids = set(args.dmr_ids) - set(dmrs["dmr_id"])
    if missing_dmr_ids:
        raise SystemExit(f"Histogram DMR IDs not found: {sorted(missing_dmr_ids)}")
    methylated_codes = [
        code.strip() for code in args.methylated_codes.split(",") if code.strip()
    ]
    if not methylated_codes:
        raise SystemExit("At least one methylated call code is required.")
    methylated_sql = "(" + ",".join(sql_quote(code) for code in methylated_codes) + ")"

    con = duckdb.connect()
    con.execute(f"PRAGMA threads = {args.threads}")
    create_input_tables(
        con, patient_files, control_files, dmrs,
        methylated_sql, args.min_cpgs_per_read,
    )

    intra, inter = calculate_variance_tables(con, args.min_reads_per_sample_dmr)
    bin_fractions, dmr_summary, selected_summary, selected_read_level = (
        calculate_subpopulation_tables(
            con,
            args.heatmap_bins,
            args.low_threshold,
            args.min_reads_per_dmr_group,
            args.top_n,
        )
    )
    histogram_read_level = con.execute(
        f"""
        SELECT sample_group, sample_name, dmr_id, read_id,
               n_cpg_calls, n_methylated_calls, read_m_frac
        FROM read_level
        WHERE dmr_id IN ({','.join(sql_quote(dmr_id) for dmr_id in args.dmr_ids)})
        ORDER BY dmr_id, sample_group, sample_name, read_m_frac
        """
    ).df()
    con.close()

    tables = {
        "intra_sample_variance.tsv": intra,
        "inter_sample_variance.tsv": inter,
        "dmr_5pct_bin_fractions.tsv": bin_fractions,
        "dmr_low_methylation_excess_summary.tsv": dmr_summary,
        "top_bottom_low_methylation_excess_dmrs.tsv": selected_summary,
        "top_bottom_low_methylation_excess_read_level.tsv": selected_read_level,
        "selected_dmr_read_level.tsv": histogram_read_level,
    }
    for filename, dataframe in tables.items():
        dataframe.to_csv(out_dir / filename, sep="\t", index=False)

    inter_paired = make_paired_variance_table(
        inter, args.inter_value_column, "mean"
    )
    intra_paired = make_paired_variance_table(
        intra, args.intra_value_column, args.intra_aggregation
    )
    inter_paired.to_csv(
        out_dir / "inter_sample_variance_scatter_data.tsv", sep="\t", index=False
    )
    intra_paired.to_csv(
        out_dir / "intra_sample_variance_scatter_data.tsv", sep="\t", index=False
    )

    plot_selected_histograms(
        histogram_read_level, args.dmr_ids,
        figure_dir / "selected_dmr_read_mfrac_histograms.png",
        args.histogram_bins,
    )
    plot_variance_scatter(
        inter_paired,
        figure_dir / "inter_sample_variance_scatter.png",
        "Inter-sample variance by DMR",
        "inter-sample variance",
    )
    plot_variance_scatter(
        intra_paired,
        figure_dir / "intra_sample_variance_scatter.png",
        f"Intra-sample variance by DMR ({args.intra_aggregation} across samples)",
        "intra-sample variance",
    )
    plot_dmr_heatmap(
        bin_fractions,
        figure_dir / "dmr_5pct_bin_heatmap_mean_fraction_bars.png",
        args.heatmap_bins,
    )
    plot_violin(
        selected_read_level,
        selected_summary,
        figure_dir / "top_bottom_low_methylation_excess_violin.png",
    )

    print(f"Wrote corrected TSVs to: {out_dir}")
    print(f"Wrote final figures to: {figure_dir}")
    print(f"Control samples (filename-derived): {len(control_files):,}")
    print(f"uRPL samples (filename-derived): {len(patient_files):,}")


if __name__ == "__main__":
    main()
