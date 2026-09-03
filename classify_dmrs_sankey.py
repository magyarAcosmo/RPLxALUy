#!/usr/bin/env python3

import gzip
import urllib.request
from pathlib import Path
import pandas as pd

DMR_BED = "/Users/jordanmoore/Desktop/THEWAY/RPL/dmrs_2026-01-17_woAGE.bed"
OUTDIR = Path("/Users/jordanmoore/Desktop/THEWAY/RPL/dmr_sankey_hg38_refseq")
OUTDIR.mkdir(exist_ok=True)

PROM_UP = 2000
PROM_DOWN = 500

NCBI_REFSEQ_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz"
RMSK_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"

refseq_gz = OUTDIR / "ncbiRefSeq.txt.gz"
rmsk_gz = OUTDIR / "rmsk.txt.gz"

def download_if_needed(url, path):
    if not path.exists():
        print(f"Downloading {url}")
        urllib.request.urlretrieve(url, path)


def load_bed(path):
    rows = []
    with open(path) as f:
        for i, line in enumerate(f):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 else f"DMR_{i+1}"
            rows.append((chrom, start, end, name))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "dmr_id"])


def merge_intervals(intervals):
    """
    intervals: list of (start, end)
    returns merged list of (start, end)
    """
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def overlaps_any(chrom, start, end, interval_dict):
    """
    Any >=1 bp overlap.
    """
    for s, e in interval_dict.get(chrom, []):
        if e <= start:
            continue
        if s >= end:
            break
        return True
    return False


def add_interval(interval_dict, chrom, start, end):
    if end > start:
        interval_dict.setdefault(chrom, []).append((start, end))


def finalize_interval_dict(interval_dict):
    return {chrom: merge_intervals(vals) for chrom, vals in interval_dict.items()}


download_if_needed(NCBI_REFSEQ_URL, refseq_gz)
download_if_needed(RMSK_URL, rmsk_gz)

dmrs = load_bed(DMR_BED)

promoters = {}
exons = {}
genes = {}

print("Parsing ncbiRefSeq...")

with gzip.open(refseq_gz, "rt") as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")

        # UCSC genePred-like columns for ncbiRefSeq:
        # bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd,
        # exonCount, exonStarts, exonEnds, score, name2, ...
        chrom = fields[2]
        strand = fields[3]
        tx_start = int(fields[4])
        tx_end = int(fields[5])
        exon_starts = [int(x) for x in fields[9].rstrip(",").split(",") if x]
        exon_ends = [int(x) for x in fields[10].rstrip(",").split(",") if x]

        add_interval(genes, chrom, tx_start, tx_end)

        if strand == "+":
            p_start = max(0, tx_start - PROM_UP)
            p_end = tx_start + PROM_DOWN
        else:
            p_start = max(0, tx_end - PROM_DOWN)
            p_end = tx_end + PROM_UP

        add_interval(promoters, chrom, p_start, p_end)

        for s, e in zip(exon_starts, exon_ends):
            add_interval(exons, chrom, s, e)

promoters = finalize_interval_dict(promoters)
exons = finalize_interval_dict(exons)
genes = finalize_interval_dict(genes)

print("Building introns...")

introns = {}

for chrom, gene_intervals in genes.items():
    exon_intervals = exons.get(chrom, [])
    for g_start, g_end in gene_intervals:
        cursor = g_start
        for e_start, e_end in exon_intervals:
            if e_end <= g_start:
                continue
            if e_start >= g_end:
                break
            if e_start > cursor:
                add_interval(introns, chrom, cursor, min(e_start, g_end))
            cursor = max(cursor, e_end)
            if cursor >= g_end:
                break
        if cursor < g_end:
            add_interval(introns, chrom, cursor, g_end)

introns = finalize_interval_dict(introns)

print("Parsing RepeatMasker...")

repeats = {}
aluy = {}

with gzip.open(rmsk_gz, "rt") as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")

        # UCSC rmsk:
        # bin, swScore, milliDiv, milliDel, milliIns,
        # genoName, genoStart, genoEnd, genoLeft, strand,
        # repName, repClass, repFamily, ...
        chrom = fields[5]
        start = int(fields[6])
        end = int(fields[7])
        rep_name = fields[10]

        add_interval(repeats, chrom, start, end)

        # Includes exact AluY and AluY subfamilies such as AluYa5, AluYb8, etc.
        if rep_name.startswith("AluY"):
            add_interval(aluy, chrom, start, end)

repeats = finalize_interval_dict(repeats)
aluy = finalize_interval_dict(aluy)

print("Classifying DMRs...")

classified = []

for _, row in dmrs.iterrows():
    chrom = row["chrom"]
    start = row["start"]
    end = row["end"]
    dmr_id = row["dmr_id"]

    if overlaps_any(chrom, start, end, promoters):
        feature = "Promoter"
    elif overlaps_any(chrom, start, end, exons):
        feature = "Exon"
    elif overlaps_any(chrom, start, end, introns):
        feature = "Intron"
    else:
        feature = "Intergenic"

    is_repeat = overlaps_any(chrom, start, end, repeats)
    is_aluy = overlaps_any(chrom, start, end, aluy)

    if is_aluy:
        repeat_class = "AluY"
    elif is_repeat:
        repeat_class = "Other repeat"
    else:
        repeat_class = "Non-repeat"

    classified.append({
        "dmr_id": dmr_id,
        "chrom": chrom,
        "start": start,
        "end": end,
        "feature": feature,
        "is_repeat": is_repeat,
        "is_aluy": is_aluy,
        "repeat_class": repeat_class,
    })

classified_df = pd.DataFrame(classified)

breakdown = (
    classified_df
    .groupby(["feature", "repeat_class"])
    .size()
    .reset_index(name="n_dmrs")
)

feature_totals = (
    classified_df
    .groupby("feature")
    .size()
    .reset_index(name="total_dmrs")
)

wide = (
    breakdown
    .pivot(index="feature", columns="repeat_class", values="n_dmrs")
    .fillna(0)
    .astype(int)
    .reset_index()
)

for col in ["AluY", "Other repeat", "Non-repeat"]:
    if col not in wide.columns:
        wide[col] = 0

wide["repeat_dmrs"] = wide["AluY"] + wide["Other repeat"]
wide = wide.merge(feature_totals, on="feature", how="left")

wide = wide[[
    "feature",
    "total_dmrs",
    "repeat_dmrs",
    "Non-repeat",
    "AluY",
    "Other repeat",
]]

wide = wide.rename(columns={
    "Non-repeat": "non_repeat_dmrs",
    "AluY": "aluy_dmrs",
    "Other repeat": "other_repeat_dmrs",
})

feature_order = ["Promoter", "Exon", "Intron", "Intergenic"]
wide["feature"] = pd.Categorical(wide["feature"], categories=feature_order, ordered=True)
wide = wide.sort_values("feature")

total_row = pd.DataFrame([{
    "feature": "TOTAL",
    "total_dmrs": int(wide["total_dmrs"].sum()),
    "repeat_dmrs": int(wide["repeat_dmrs"].sum()),
    "non_repeat_dmrs": int(wide["non_repeat_dmrs"].sum()),
    "aluy_dmrs": int(wide["aluy_dmrs"].sum()),
    "other_repeat_dmrs": int(wide["other_repeat_dmrs"].sum()),
}])

wide = pd.concat([wide, total_row], ignore_index=True)

classified_df.to_csv(OUTDIR / "dmr_classification_per_region.tsv", sep="\t", index=False)
wide.to_csv(OUTDIR / "dmr_sankey_breakdown.tsv", sep="\t", index=False)

print("\nSankey breakdown:")
print(wide.to_string(index=False))

print(f"\nWrote:")
print(f"  {OUTDIR / 'dmr_classification_per_region.tsv'}")
print(f"  {OUTDIR / 'dmr_sankey_breakdown.tsv'}")