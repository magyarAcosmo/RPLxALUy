#!/usr/bin/env python3

import argparse
import gzip
import math
import re
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


GROUP_DISPLAY = {"control": "Control", "patient": "uRPL"}
GROUP_COLORS = {"control": "#2166ac", "patient": "#d73027"}


def sql_quote(value):
    return "'" + str(value).replace("'", "''") + "'"


def sql_list(values):
    return "[" + ",".join(sql_quote(value) for value in values) + "]"


def collect_files(folder, pattern):
    files = sorted(Path(folder).glob(pattern))
    if not files:
        raise SystemExit(f"No files matching {pattern!r} found in {folder}")
    return [str(path) for path in files]


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def parse_gtf_attributes(attributes):
    return {
        key: value
        for key, value in re.findall(r'(\S+)\s+"([^"]+)"', attributes)
    }


def load_gencode_gene_promoters(gtf_path, upstream, downstream):
    records = []

    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue

            chrom, _, _, start, end, _, strand, _, attributes = fields
            attrs = parse_gtf_attributes(attributes)
            gene_id = attrs.get("gene_id")
            if not gene_id or strand not in {"+", "-"}:
                continue

            # GTF is 1-based closed; ch3 start coordinates are treated as 0-based.
            gene_start0 = int(start) - 1
            gene_end0 = int(end)
            tss0 = gene_start0 if strand == "+" else gene_end0

            if strand == "+":
                promoter_start = max(0, tss0 - upstream)
                promoter_end = tss0 + downstream
            else:
                promoter_start = max(0, tss0 - downstream)
                promoter_end = tss0 + upstream

            records.append(
                {
                    "promoter_id": gene_id,
                    "chrom": chrom,
                    "promoter_start": promoter_start,
                    "promoter_end": promoter_end,
                }
            )

    promoters = pd.DataFrame.from_records(records)
    if promoters.empty:
        raise SystemExit(
            "No gene features were read from the GENCODE GTF. Confirm that the "
            "file is a GTF containing rows whose third column is 'gene'."
        )

    promoters = promoters.drop_duplicates("promoter_id")
    return promoters


def query_distribution(
    con,
    region_kind,
    n_bins,
    min_calls,
    methylated_sql,
    tile_size=None,
):
    if region_kind == "promoter":
        region_cte = f"""
            SELECT
                c.sample_group,
                c.sample_name,
                p.promoter_id AS region_id,
                COUNT(*) AS n_calls,
                SUM(
                    CASE WHEN c.call_code IN {methylated_sql} THEN 1 ELSE 0 END
                )::DOUBLE / COUNT(*) AS m_frac
            FROM ch3 c
            JOIN promoters p
              ON c.chrom = p.chrom
             AND c.start >= p.promoter_start
             AND c.start < p.promoter_end
            WHERE
                (LENGTH(query_kmer) = 2 AND query_kmer = 'CG')
                OR
                (LENGTH(query_kmer) = 5 AND SUBSTR(query_kmer, 3, 2) = 'CG')
            GROUP BY c.sample_group, c.sample_name, p.promoter_id
            HAVING COUNT(*) >= {min_calls}
        """
    elif region_kind == "tile":
        region_cte = f"""
            SELECT
                sample_group,
                sample_name,
                chrom || ':' ||
                    CAST(FLOOR(start / {tile_size}) * {tile_size} AS BIGINT)::VARCHAR
                    AS region_id,
                COUNT(*) AS n_calls,
                SUM(
                    CASE WHEN call_code IN {methylated_sql} THEN 1 ELSE 0 END
                )::DOUBLE / COUNT(*) AS m_frac
            FROM ch3
            WHERE
                (LENGTH(query_kmer) = 2 AND query_kmer = 'CG')
                OR
                (LENGTH(query_kmer) = 5 AND SUBSTR(query_kmer, 3, 2) = 'CG')
            GROUP BY
                sample_group,
                sample_name,
                chrom,
                CAST(FLOOR(start / {tile_size}) * {tile_size} AS BIGINT)
            HAVING COUNT(*) >= {min_calls}
        """
    else:
        raise ValueError(f"Unknown region kind: {region_kind}")

    result = con.execute(
        f"""
        WITH region_values AS (
            {region_cte}
        ),
        binned AS (
            SELECT
                sample_group,
                sample_name,
                CASE
                    WHEN m_frac >= 1.0 THEN {n_bins - 1}
                    WHEN m_frac <= 0.0 THEN 0
                    ELSE FLOOR(m_frac * {n_bins})::INTEGER
                END AS bin_id
            FROM region_values
        ),
        counts AS (
            SELECT
                sample_group,
                sample_name,
                bin_id,
                COUNT(*) AS region_count
            FROM binned
            GROUP BY sample_group, sample_name, bin_id
        )
        SELECT
            sample_group,
            sample_name,
            bin_id,
            region_count,
            region_count::DOUBLE
                / SUM(region_count) OVER (PARTITION BY sample_group, sample_name)
                AS fraction
        FROM counts
        ORDER BY sample_group, sample_name, bin_id
        """
    ).df()

    if result.empty:
        raise SystemExit(
            f"No {region_kind} values passed the minimum-call threshold."
        )
    return result


def complete_distribution_bins(distribution, n_bins):
    samples = distribution[["sample_group", "sample_name"]].drop_duplicates()
    full_index = pd.MultiIndex.from_tuples(
        [
            (row.sample_group, row.sample_name, bin_id)
            for row in samples.itertuples(index=False)
            for bin_id in range(n_bins)
        ],
        names=["sample_group", "sample_name", "bin_id"],
    )

    return (
        distribution.set_index(["sample_group", "sample_name", "bin_id"])
        .reindex(full_index)
        .reset_index()
        .assign(
            region_count=lambda frame: frame["region_count"].fillna(0),
            fraction=lambda frame: frame["fraction"].fillna(0.0),
        )
    )


def plot_distribution(distribution, out_png, title, y_label, panel_label, n_bins):
    distribution = complete_distribution_bins(distribution, n_bins)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    x = np.arange(n_bins) / n_bins

    for (group, sample_name), sub in distribution.groupby(
        ["sample_group", "sample_name"], sort=True
    ):
        sub = sub.sort_values("bin_id")
        ax.plot(
            x,
            sub["fraction"],
            marker="o",
            markersize=2.8,
            linewidth=1.1,
            color=GROUP_COLORS[group],
            alpha=0.38,
        )

    handles = [
        Line2D(
            [0],
            [0],
            color=GROUP_COLORS[group],
            linewidth=2,
            marker="o",
            markersize=4,
            label=GROUP_DISPLAY[group],
        )
        for group in ["control", "patient"]
    ]

    ax.set_xticks(x)
    ax.set_xticklabels([f"{value:.1f}" for value in x])
    ax.set_xlim(-0.035, x[-1] + 0.035)
    ax.set_xlabel("Methylation fraction bin")
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.legend(handles=handles, frameon=False)
    fig.text(0.015, 0.975, panel_label, fontsize=28, fontweight="bold", va="top")

    plt.tight_layout(rect=[0.03, 0, 1, 1])
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def query_10kb_values(con, tile_size, min_calls, methylated_sql):
    con.execute(
        f"""
        CREATE OR REPLACE TEMP TABLE tile10_values AS
        SELECT
            sample_group,
            sample_name,
            chrom,
            CAST(FLOOR(start / {tile_size}) * {tile_size} AS BIGINT) AS tile_start,
            COUNT(*) AS n_calls,
            SUM(
                CASE WHEN call_code IN {methylated_sql} THEN 1 ELSE 0 END
            )::DOUBLE / COUNT(*) AS m_frac
        FROM ch3
        WHERE
            (LENGTH(query_kmer) = 2 AND query_kmer = 'CG')
            OR
            (LENGTH(query_kmer) = 5 AND SUBSTR(query_kmer, 3, 2) = 'CG')
        GROUP BY sample_group, sample_name, chrom, tile_start
        HAVING COUNT(*) >= {min_calls}
        """
    )

    samples = con.execute(
        """
        SELECT DISTINCT sample_group, sample_name
        FROM ch3
        ORDER BY sample_group, sample_name
        """
    ).df()

    return samples


def build_pca_matrix(con, samples, min_sample_fraction):
    n_samples = len(samples)
    minimum_samples = max(2, math.ceil(n_samples * min_sample_fraction))

    values = con.execute(
        f"""
        WITH eligible_tiles AS (
            SELECT chrom, tile_start
            FROM tile10_values
            GROUP BY chrom, tile_start
            HAVING COUNT(*) >= {minimum_samples}
        )
        SELECT
            t.sample_group,
            t.sample_name,
            t.chrom,
            t.tile_start,
            t.m_frac
        FROM tile10_values t
        JOIN eligible_tiles e
          ON t.chrom = e.chrom
         AND t.tile_start = e.tile_start
        ORDER BY t.sample_group, t.sample_name, t.chrom, t.tile_start
        """
    ).df()

    if values.empty:
        raise SystemExit(
            "No 10 kb windows passed the requested cross-sample coverage threshold. "
            "Try lowering --min-window-sample-fraction or --min-calls-10kb."
        )

    matrix = values.pivot_table(
        index=["sample_group", "sample_name"],
        columns=["chrom", "tile_start"],
        values="m_frac",
        aggfunc="first",
    )

    matrix = matrix.dropna(axis=1, how="all")
    matrix = matrix.fillna(matrix.mean(axis=0))
    matrix = matrix.loc[:, matrix.var(axis=0, ddof=0) > 0]

    if matrix.shape[0] < 3 or matrix.shape[1] < 2:
        raise SystemExit(
            "The filtered 10 kb matrix needs at least three samples and two "
            "variable windows for PCA."
        )

    return matrix, minimum_samples


def pca_by_svd(matrix):
    x = matrix.to_numpy(dtype=float)
    centered = x - x.mean(axis=0, keepdims=True)
    u, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
    scores = u[:, :2] * singular_values[:2]
    eigenvalues = singular_values**2
    explained = eigenvalues / eigenvalues.sum()
    return centered, scores, explained[:2]


def pseudo_f_from_distances(distance_sq, labels):
    labels = np.asarray(labels)
    n = len(labels)
    groups = np.unique(labels)
    n_groups = len(groups)

    if n_groups < 2 or n <= n_groups:
        raise ValueError("PERMANOVA requires at least two groups and residual degrees of freedom.")

    total_ss = distance_sq.sum() / (2 * n)
    within_ss = 0.0
    for group in groups:
        idx = np.flatnonzero(labels == group)
        within_ss += distance_sq[np.ix_(idx, idx)].sum() / (2 * len(idx))

    between_ss = max(0.0, total_ss - within_ss)
    pseudo_f = (between_ss / (n_groups - 1)) / (within_ss / (n - n_groups))
    r_squared = between_ss / total_ss if total_ss > 0 else np.nan
    return pseudo_f, r_squared


def permanova(centered_matrix, labels, permutations, seed):
    gram = centered_matrix @ centered_matrix.T
    diagonal = np.diag(gram)
    distance_sq = diagonal[:, None] + diagonal[None, :] - 2 * gram
    distance_sq = np.maximum(distance_sq, 0.0)

    observed_f, r_squared = pseudo_f_from_distances(distance_sq, labels)
    rng = np.random.default_rng(seed)
    exceedances = 0

    for _ in range(permutations):
        permuted_f, _ = pseudo_f_from_distances(
            distance_sq,
            rng.permutation(labels),
        )
        if permuted_f >= observed_f - 1e-12:
            exceedances += 1

    p_value = (exceedances + 1) / (permutations + 1)
    return observed_f, r_squared, p_value


def plot_pca(
    matrix,
    scores,
    explained,
    pseudo_f,
    r_squared,
    p_value,
    permutations,
    out_png,
):
    metadata = matrix.index.to_frame(index=False)
    fig, ax = plt.subplots(figsize=(7.2, 6.2))

    for group in ["control", "patient"]:
        mask = metadata["sample_group"].to_numpy() == group
        ax.scatter(
            scores[mask, 0],
            scores[mask, 1],
            s=60,
            alpha=0.8,
            color=GROUP_COLORS[group],
            edgecolor="black",
            linewidth=0.4,
            label=GROUP_DISPLAY[group],
        )

    ax.axhline(0, color="0.8", linewidth=0.8)
    ax.axvline(0, color="0.8", linewidth=0.8)
    ax.set_xlabel(f"PC1 ({explained[0] * 100:.2f}%)")
    ax.set_ylabel(f"PC2 ({explained[1] * 100:.2f}%)")
    ax.set_title("Genome-wide methylation PCA (10 kb windows)")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)

    annotation = (
        "PERMANOVA (Euclidean)\n"
        f"pseudo-F = {pseudo_f:.3f}\n"
        f"R² = {r_squared:.3f}\n"
        f"p = {p_value:.4f} ({permutations:,} permutations)"
    )
    ax.text(
        0.98,
        0.02,
        annotation,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "alpha": 0.85},
    )
    fig.text(0.015, 0.975, "f", fontsize=28, fontweight="bold", va="top")

    plt.tight_layout(rect=[0.03, 0, 1, 1])
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Recreate promoter distributions (1c), genome-wide 1 kb "
            "distributions (1d), and a 10 kb-window PCA with PERMANOVA (1f) "
            "from Control and uRPL ch3 parquet files."
        )
    )
    ap.add_argument("--patients", required=True, help="Folder containing uRPL ch3 files")
    ap.add_argument("--controls", required=True, help="Folder containing Control ch3 files")
    ap.add_argument("--gencode-gtf", required=True, help="GENCODE gene annotation GTF or GTF.gz")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--pattern", default="*.ch3")
    ap.add_argument("--methylated-codes", default="m")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--distribution-bins", type=int, default=10)
    ap.add_argument("--promoter-upstream", type=int, default=2000)
    ap.add_argument("--promoter-downstream", type=int, default=500)
    ap.add_argument("--min-calls-promoter", type=int, default=10)
    ap.add_argument("--min-calls-1kb", type=int, default=5)
    ap.add_argument("--min-calls-10kb", type=int, default=20)
    ap.add_argument("--tile-size-1kb", type=int, default=1000)
    ap.add_argument("--tile-size-10kb", type=int, default=10000)
    ap.add_argument(
        "--min-window-sample-fraction",
        type=float,
        default=0.8,
        help="Minimum fraction of samples with data required for a PCA window.",
    )
    ap.add_argument("--permutations", type=int, default=999)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    if not 0 < args.min_window_sample_fraction <= 1:
        raise SystemExit("--min-window-sample-fraction must be in (0, 1].")
    if args.distribution_bins < 2:
        raise SystemExit("--distribution-bins must be at least 2.")
    if args.permutations < 1:
        raise SystemExit("--permutations must be at least 1.")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    patient_files = collect_files(args.patients, args.pattern)
    control_files = collect_files(args.controls, args.pattern)
    promoters = load_gencode_gene_promoters(
        args.gencode_gtf,
        upstream=args.promoter_upstream,
        downstream=args.promoter_downstream,
    )

    methylated_codes = [
        value.strip() for value in args.methylated_codes.split(",") if value.strip()
    ]
    if not methylated_codes:
        raise SystemExit("At least one --methylated-codes value is required.")
    methylated_sql = "(" + ",".join(sql_quote(code) for code in methylated_codes) + ")"

    con = duckdb.connect()
    con.execute(f"PRAGMA threads = {args.threads}")
    con.register("promoters_df", promoters)
    con.execute(
        f"""
        CREATE TEMP TABLE promoters AS
        SELECT
            promoter_id::VARCHAR AS promoter_id,
            chrom::VARCHAR AS chrom,
            promoter_start::BIGINT AS promoter_start,
            promoter_end::BIGINT AS promoter_end
        FROM promoters_df;

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
        """
    )

    group_counts = con.execute(
        """
        SELECT sample_group, COUNT(DISTINCT sample_name) AS n_samples
        FROM ch3
        GROUP BY sample_group
        ORDER BY sample_group
        """
    ).df()
    found_groups = set(group_counts["sample_group"])
    if found_groups != {"control", "patient"}:
        raise SystemExit("Both Control and uRPL folders must contain usable samples.")

    promoter_distribution = query_distribution(
        con,
        region_kind="promoter",
        n_bins=args.distribution_bins,
        min_calls=args.min_calls_promoter,
        methylated_sql=methylated_sql,
    )
    promoter_distribution.to_csv(
        out_dir / "figure_1c_promoter_distribution_data.tsv",
        sep="\t",
        index=False,
    )
    plot_distribution(
        promoter_distribution,
        out_png=out_dir / "figure_1c_promoters.png",
        title="Promoters",
        y_label="Fraction of promoters",
        panel_label="c",
        n_bins=args.distribution_bins,
    )

    tile1_distribution = query_distribution(
        con,
        region_kind="tile",
        n_bins=args.distribution_bins,
        min_calls=args.min_calls_1kb,
        methylated_sql=methylated_sql,
        tile_size=args.tile_size_1kb,
    )
    tile1_distribution.to_csv(
        out_dir / "figure_1d_1kb_distribution_data.tsv",
        sep="\t",
        index=False,
    )
    plot_distribution(
        tile1_distribution,
        out_png=out_dir / "figure_1d_genomewide_1kb.png",
        title="Genome-wide (1 kb bins)",
        y_label="Fraction of 1 kb windows",
        panel_label="d",
        n_bins=args.distribution_bins,
    )

    samples = query_10kb_values(
        con,
        tile_size=args.tile_size_10kb,
        min_calls=args.min_calls_10kb,
        methylated_sql=methylated_sql,
    )
    matrix, minimum_samples = build_pca_matrix(
        con,
        samples=samples,
        min_sample_fraction=args.min_window_sample_fraction,
    )
    centered, scores, explained = pca_by_svd(matrix)
    labels = matrix.index.get_level_values("sample_group").to_numpy()
    observed_f, r_squared, p_value = permanova(
        centered,
        labels=labels,
        permutations=args.permutations,
        seed=args.seed,
    )

    pca_scores = matrix.index.to_frame(index=False)
    pca_scores["PC1"] = scores[:, 0]
    pca_scores["PC2"] = scores[:, 1]
    pca_scores.to_csv(
        out_dir / "figure_1f_10kb_pca_scores.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        [
            {
                "distance": "Euclidean",
                "pseudo_F": observed_f,
                "R_squared": r_squared,
                "p_value": p_value,
                "permutations": args.permutations,
                "seed": args.seed,
                "n_samples": matrix.shape[0],
                "n_windows": matrix.shape[1],
                "minimum_samples_per_window": minimum_samples,
            }
        ]
    ).to_csv(
        out_dir / "figure_1f_permanova.tsv",
        sep="\t",
        index=False,
    )
    plot_pca(
        matrix,
        scores=scores,
        explained=explained,
        pseudo_f=observed_f,
        r_squared=r_squared,
        p_value=p_value,
        permutations=args.permutations,
        out_png=out_dir / "figure_1f_genomewide_10kb_pca.png",
    )

    con.close()
    print(f"Wrote outputs to: {out_dir}")
    print(f"Control files: {len(control_files):,}")
    print(f"uRPL files: {len(patient_files):,}")
    print(f"GENCODE gene promoters: {len(promoters):,}")
    print(f"PCA matrix: {matrix.shape[0]:,} samples x {matrix.shape[1]:,} windows")
    print(f"PERMANOVA: pseudo-F={observed_f:.4f}, R^2={r_squared:.4f}, p={p_value:.4g}")


if __name__ == "__main__":
    main()
