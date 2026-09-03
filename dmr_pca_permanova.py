#!/usr/bin/env python3

import argparse
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
    return [str(path.resolve()) for path in files]


def load_dmrs(path):
    dmrs = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        comment="#",
    )
    dmrs = dmrs.iloc[:, :4].copy()
    if dmrs.shape[1] < 3:
        raise SystemExit(
            "The DMR BED file must contain at least: chromosome, start, end."
        )

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
        missing_id = dmrs["dmr_id"].isna()
        dmrs.loc[missing_id, "dmr_id"] = (
            dmrs.loc[missing_id, "chrom"].astype(str)
            + ":"
            + dmrs.loc[missing_id, "dmr_start"].astype(str)
            + "-"
            + dmrs.loc[missing_id, "dmr_end"].astype(str)
        )

    dmrs["chrom"] = dmrs["chrom"].astype(str)
    dmrs["dmr_start"] = pd.to_numeric(dmrs["dmr_start"], errors="raise").astype(
        "int64"
    )
    dmrs["dmr_end"] = pd.to_numeric(dmrs["dmr_end"], errors="raise").astype(
        "int64"
    )
    dmrs["dmr_id"] = dmrs["dmr_id"].astype(str)

    invalid = (dmrs["dmr_start"] < 0) | (dmrs["dmr_end"] <= dmrs["dmr_start"])
    if invalid.any():
        examples = dmrs.loc[
            invalid, ["chrom", "dmr_start", "dmr_end", "dmr_id"]
        ].head()
        raise SystemExit(f"Invalid DMR coordinates found:\n{examples.to_string(index=False)}")

    duplicate_ids = dmrs.loc[dmrs["dmr_id"].duplicated(False), "dmr_id"].unique()
    if len(duplicate_ids):
        raise SystemExit(
            "DMR identifiers must be unique. Duplicates include: "
            + ", ".join(map(str, duplicate_ids[:10]))
        )

    return dmrs[["dmr_id", "chrom", "dmr_start", "dmr_end"]]


def make_sample_table(patient_files, control_files):
    rows = []
    for group, files in [
        ("patient", patient_files),
        ("control", control_files),
    ]:
        for filename in files:
            rows.append(
                {
                    "sample_group": group,
                    "sample_name": Path(filename).stem,
                    "filename": filename,
                }
            )

    samples = pd.DataFrame(rows)
    duplicate_keys = samples.duplicated(
        ["sample_group", "sample_name"], keep=False
    )
    if duplicate_keys.any():
        duplicates = samples.loc[
            duplicate_keys, ["sample_group", "sample_name", "filename"]
        ]
        raise SystemExit(
            "Duplicate filename-derived sample identifiers were found:\n"
            + duplicates.to_string(index=False)
        )
    return samples


def create_ch3_view(con, patient_files, control_files):
    con.execute(
        f"""
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


def query_dmr_values(con, methylated_sql):
    con.execute(
        f"""
        CREATE TEMP TABLE dmr_values AS
        SELECT
            c.sample_group,
            c.sample_name,
            d.dmr_id,
            d.chrom,
            d.dmr_start,
            d.dmr_end,
            COUNT(*) AS n_calls,
            SUM(
                CASE WHEN c.call_code IN {methylated_sql} THEN 1 ELSE 0 END
            )::DOUBLE / COUNT(*) AS m_frac
        FROM ch3 c
        JOIN dmrs d
          ON c.chrom = d.chrom
         AND c.start >= d.dmr_start
         AND c.start < d.dmr_end
        WHERE
            (LENGTH(c.query_kmer) = 2 AND c.query_kmer = 'CG')
            OR
            (
                LENGTH(c.query_kmer) = 5
                AND SUBSTR(c.query_kmer, 3, 2) = 'CG'
            )
        GROUP BY
            c.sample_group,
            c.sample_name,
            d.dmr_id,
            d.chrom,
            d.dmr_start,
            d.dmr_end;
        """
    )


def build_pca_matrix(con, samples):
    n_samples = len(samples)

    coverage = con.execute(
        """
        SELECT
            d.dmr_id,
            d.chrom,
            d.dmr_start,
            d.dmr_end,
            COUNT(DISTINCT v.sample_group || ':' || v.sample_name)
                AS n_samples_with_data,
            COUNT(DISTINCT CASE
                WHEN v.sample_group = 'control' THEN v.sample_name
            END) AS n_control_with_data,
            COUNT(DISTINCT CASE
                WHEN v.sample_group = 'patient' THEN v.sample_name
            END) AS n_urpl_with_data
        FROM dmrs d
        LEFT JOIN dmr_values v USING (dmr_id)
        GROUP BY d.dmr_id, d.chrom, d.dmr_start, d.dmr_end
        ORDER BY d.chrom, d.dmr_start, d.dmr_end
        """
    ).df()
    coverage["sample_fraction"] = coverage["n_samples_with_data"] / n_samples
    values = con.execute(
        """
        SELECT
            sample_group,
            sample_name,
            dmr_id,
            m_frac
        FROM dmr_values
        ORDER BY sample_group, sample_name, dmr_id
        """
    ).df()
    if values.empty:
        raise SystemExit(
            "No retained CG calls overlapped the supplied DMRs. Check that "
            "the BED and .ch3 files use the same chromosome naming and "
            "coordinate system."
        )

    matrix = values.pivot(
        index=["sample_group", "sample_name"],
        columns="dmr_id",
        values="m_frac",
    )
    matrix = matrix.reindex(columns=coverage["dmr_id"].tolist())
    expected_index = pd.MultiIndex.from_frame(
        samples[["sample_group", "sample_name"]]
    )
    matrix = matrix.reindex(expected_index)

    sample_coverage = pd.DataFrame(index=matrix.index)
    sample_coverage["n_dmrs_with_data"] = matrix.notna().sum(axis=1)
    sample_coverage["dmr_fraction_with_data"] = matrix.notna().mean(axis=1)
    sample_coverage = sample_coverage.reset_index()

    empty_samples = sample_coverage.loc[
        sample_coverage["n_dmrs_with_data"] == 0,
        ["sample_group", "sample_name"],
    ]
    if not empty_samples.empty:
        names = ", ".join(
            f"{row.sample_group}:{row.sample_name}"
            for row in empty_samples.itertuples(index=False)
        )
        raise SystemExit(
            "Some samples had no measurements at any supplied DMR and "
            f"cannot be placed meaningfully in the PCA: {names}"
        )

    matrix = matrix.fillna(matrix.mean(axis=0))
    variable = matrix.var(axis=0, ddof=0) > 0
    variable_ids = set(variable.index[variable])
    coverage["included_in_pca"] = coverage["dmr_id"].isin(variable_ids)
    matrix = matrix.loc[:, variable]

    if matrix.shape[0] < 3 or matrix.shape[1] < 2:
        raise SystemExit(
            "The filtered matrix needs at least three samples and two variable "
            "DMRs for PCA."
        )

    return matrix, coverage, sample_coverage


def pca_by_svd(matrix):
    x = matrix.to_numpy(dtype=float)
    centered = x - x.mean(axis=0, keepdims=True)
    _, singular_values, vt = np.linalg.svd(centered, full_matrices=False)
    scores = centered @ vt[:2].T
    eigenvalues = singular_values**2
    explained = eigenvalues / eigenvalues.sum()
    return centered, scores, explained[:2]


def pseudo_f_from_distances(distance_sq, labels):
    labels = np.asarray(labels)
    n = len(labels)
    groups = np.unique(labels)
    n_groups = len(groups)

    if n_groups < 2 or n <= n_groups:
        raise ValueError(
            "PERMANOVA requires at least two groups and residual degrees of freedom."
        )

    total_ss = distance_sq.sum() / (2 * n)
    within_ss = 0.0
    for group in groups:
        indices = np.flatnonzero(labels == group)
        within_ss += (
            distance_sq[np.ix_(indices, indices)].sum()
            / (2 * len(indices))
        )

    between_ss = max(0.0, total_ss - within_ss)
    pseudo_f = (
        (between_ss / (n_groups - 1))
        / (within_ss / (n - n_groups))
    )
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
    panel_label,
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
    ax.set_title("Methylation PCA at DMRs")
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
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": "white",
            "alpha": 0.85,
        },
    )

    if panel_label:
        fig.text(
            0.015,
            0.975,
            panel_label,
            fontsize=28,
            fontweight="bold",
            va="top",
        )

    plt.tight_layout(rect=[0.03, 0, 1, 1])
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create a Control-versus-uRPL methylation PCA and Euclidean "
            "PERMANOVA using DMR regions from a BED file."
        )
    )
    parser.add_argument(
        "--patients",
        required=True,
        help="Folder containing uRPL .ch3 files.",
    )
    parser.add_argument(
        "--controls",
        required=True,
        help="Folder containing Control .ch3 files.",
    )
    parser.add_argument(
        "--dmrs",
        required=True,
        help="BED-like DMR file: chrom start end [dmr_id].",
    )
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--pattern", default="*.ch3")
    parser.add_argument("--methylated-codes", default="m")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument(
        "--panel-label",
        default="f",
        help="Panel label drawn at upper left; pass an empty string to omit it.",
    )
    args = parser.parse_args()

    if args.permutations < 1:
        raise SystemExit("--permutations must be at least 1.")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    patient_files = collect_files(args.patients, args.pattern)
    control_files = collect_files(args.controls, args.pattern)
    samples = make_sample_table(patient_files, control_files)
    dmrs = load_dmrs(args.dmrs)

    methylated_codes = [
        value.strip()
        for value in args.methylated_codes.split(",")
        if value.strip()
    ]
    if not methylated_codes:
        raise SystemExit("At least one --methylated-codes value is required.")
    methylated_sql = (
        "(" + ",".join(sql_quote(code) for code in methylated_codes) + ")"
    )

    con = duckdb.connect()
    con.execute(f"PRAGMA threads = {args.threads}")
    con.register("dmrs_df", dmrs)
    con.execute(
        """
        CREATE TEMP TABLE dmrs AS
        SELECT
            dmr_id::VARCHAR AS dmr_id,
            chrom::VARCHAR AS chrom,
            dmr_start::BIGINT AS dmr_start,
            dmr_end::BIGINT AS dmr_end
        FROM dmrs_df
        """
    )
    create_ch3_view(con, patient_files, control_files)
    query_dmr_values(
        con,
        methylated_sql=methylated_sql,
    )

    dmr_values = con.execute(
        """
        SELECT *
        FROM dmr_values
        ORDER BY chrom, dmr_start, dmr_end, sample_group, sample_name
        """
    ).df()
    dmr_values.to_csv(
        out_dir / "dmr_sample_methylation_values.tsv",
        sep="\t",
        index=False,
    )

    matrix, coverage, sample_coverage = build_pca_matrix(
        con,
        samples=samples,
    )

    coverage.to_csv(
        out_dir / "dmr_coverage_summary.tsv",
        sep="\t",
        index=False,
    )
    sample_coverage.to_csv(
        out_dir / "sample_dmr_coverage_summary.tsv",
        sep="\t",
        index=False,
    )

    matrix_output = matrix.copy()
    matrix_output.columns.name = None
    matrix_output.reset_index().to_csv(
        out_dir / "dmr_pca_matrix_imputed.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
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
    pca_scores["group_label"] = pca_scores["sample_group"].map(GROUP_DISPLAY)
    pca_scores["PC1"] = scores[:, 0]
    pca_scores["PC2"] = scores[:, 1]
    pca_scores.to_csv(
        out_dir / "dmr_pca_scores.tsv",
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
                "n_control": int((labels == "control").sum()),
                "n_urpl": int((labels == "patient").sum()),
                "n_variable_dmrs": matrix.shape[1],
                "n_input_dmrs": len(dmrs),
                "n_dmrs_with_any_data": int(
                    (coverage["n_samples_with_data"] > 0).sum()
                ),
            }
        ]
    ).to_csv(
        out_dir / "dmr_permanova.tsv",
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
        panel_label=args.panel_label,
        out_png=out_dir / "dmr_pca_permanova.png",
    )

    con.close()
    print(f"Wrote outputs to: {out_dir}")
    print(f"Control files: {len(control_files):,}")
    print(f"uRPL files: {len(patient_files):,}")
    print(
        f"PCA matrix: {matrix.shape[0]:,} samples x "
        f"{matrix.shape[1]:,} variable DMRs"
    )
    print("No minimum call-count or cross-sample coverage filter was applied.")
    print(
        f"PERMANOVA: pseudo-F={observed_f:.4f}, "
        f"R^2={r_squared:.4f}, p={p_value:.4g}"
    )


if __name__ == "__main__":
    main()