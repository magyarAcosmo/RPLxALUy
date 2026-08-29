# DMR EFFECT
## match up noAGE DMRs with AGE DMRs in dmrs_2026-01-17_effect.bed
library(GenomicRanges)
library(rtracklayer)

# read both files
woAGE <- read.delim("DMRs_woAGE/dmrs_2026-01-17_woAGE.bed", header = FALSE,
                    col.names = c("chr", "start", "end"))

effect <- read.delim("DMRs_wAGE/dmrs_2026-01-17_effect.bed", header = FALSE,
                     col.names = c("chr", "start", "end", "x_counter", "effect_value", "strand"))

# convert to GRanges
# bed files are 0-based half-open; GRanges is 1-based inclusive, so add 1 to start

gr_woAGE <- GRanges(seqnames = woAGE$chr,
                    ranges = IRanges(start = woAGE$start + 1, end = woAGE$end))
gr_effect <- GRanges(seqnames = effect$chr,
                     ranges = IRanges(start = effect$start + 1, end = effect$end),
                     effect_value = effect$effect_value)

# Find overlaps: for each woAGE DMR, find matching effect DMR
hits <- findOverlaps(gr_woAGE, gr_effect, type = "equal")   # exact coordinate match

# Build Output
df_out <- woAGE[queryHits(hits), ]
df_out$effect_value <- effect$effect_value[subjectHits(hits)]

# check if every woAGE DMR got a match
cat("woAGE DMRs:", nrow(woAGE), "\n")
cat("Matched DMRs:", nrow(df_out), "\n")

# write output file
write.table(df_out,
            file = "dmrs_2026-01-17_woAGE_effect.bed",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
######
#  Results: all 294 woAGE DMRs matched!
######

# AluY coordinates
## rmsk.txt.gz contains all AluY coordinates
## include all subfamilies that start with AluY (many subfamilies will be included)

# read rmsk.txt file
rmsk <- read.delim("rmsk.txt", header = FALSE)

# quick check on column 11 to make sure family is present -- looks good!
head(rmsk[, 11], 10)

# filter: repName (col 11) starts with "AluY"
aluY <- rmsk[grepl("^AluY", rmsk[, 11]), ]

# check what subfamilies were included
table(aluY[, 11])   # 28 subfamilies

cat("Total AluY elements found:", nrow(aluY), "\n")   #151627

# write output
write.table(aluY,
            file = "rmsk_AluY.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# clean AluY elements to relevant columns
aluY_clean <- data.frame(
  chr = aluY[, 6],        # genoName
  start = aluY[, 7],      # genoStart (0-based, BED-style)
  end = aluY[, 8],        # genoEnd
  strand = aluY[, 10],    # strand
  repName = aluY[, 11],   # aluY subfamily
  repStart = aluY[, 14],  # start position within AluY consensus
  repEnd = aluY[, 15],    # end position wihtin AluY consensus
  repLeft = aluY[, 16]    # bases remaining to end of consensus
)

head(aluY_clean, 5)
head(aluY_clean$repLeft, 5)

# Normalize consensus coordinates by strand
## + strand: true footprint = reStart -> repEnd
## - strand: true footprint = repLeft -> repEnd
## UCSC reference: https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema

aluY_clean$consensus_start <- ifelse(aluY_clean$strand == "+", aluY_clean$repStart, aluY_clean$repLeft)
aluY_clean$consensus_end <- aluY_clean$repEnd   # always the middle/true end coordinate

# check: consensus_start always <= consensus_end
sum(aluY_clean$consensus_start > aluY_clean$consensus_end)

### Establish a sensible minimum AluY element length to exlude other very short/low-quality fragments? (anything under 50bp?)

head(aluY_clean[, c("chr", "start", "end", "strand", "repStart", "repEnd", "repLeft",
                    "consensus_start", "consensus_end")], 5)

write.table(aluY_clean,
            file = "AluY_elements_annot.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)