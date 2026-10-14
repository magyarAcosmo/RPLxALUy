# RPLxALUy
Code used to generate results of manuscript.

Each FigureX.R script is self-contained: it includes the pre-figure data prep steps
and the plotting code for that figure. Run any single script top-to-bottom to reproduce
that figure from the input files.

### tools/
This folder holds third-party software used by this project. Each tool is included as a **git submodule**,
so this repo stores only a pointer to a specific commit of the tool's own repository rather than a copy of its code.

#### Contents
| Folder    | Upstream repo                                | Version | Purpose                                                                       |
|-----------|----------------------------------------------|---------|-------------------------------------------------------------------------------|
| `dnm/`    | https://github.com/smithlabcode/dnmtools.git | v1.5.1  | High-throughput analysis of DNA methylation data                              |
| `dorado/` | https://github.com/nanoporetech/dorado       | v1.1.1  | Oxford Nanopore (converts raw signal to reads with modification-base calling) |
| `modkit/` | https://github.com/nanoporetech/modkit       | v0.5.0  | Nanopore tool for summarizing and analyzing base modifications from BAM files |


### Script redundancy
Several figures use the same upstream prep steps (e.g. matching without age (woAGE) DMRs
to their effect sizes, or building the AluY consensus notation). That prep code is duplicated
across the relevant figure scripts rather than factored into a shared file. This is intentional:
it keeps each figure fully traceable and runnable on its own, at the cost of some duplicated code between 
files.

### Files

|            |                                                             |
|------------|-------------------------------------------------------------|
| **Script** | **Output**                                                  | 
| Figure1a.R | Fraction of sites covered (Control vs uRPL)                 |                     
| Figure1b.R | Read coverage per site (Control vs uRPL)                    |                     
| Figure1c.R | Fraction of promoters per methylation bin (Control vs uRPL) |                     
| Figure1d.R | Genome wide methylation faction bin (Control vs uRPL)       |                     
| Figure1e.R | Mean methylation per sample (Control vs uRPL)               |                     
| Figure1f.R | PCA (Control vs uRPL)                                       |                     
| Figure2a.R | Volcano plot (methylation difference vs -log10 p-value)     |
| Figure2b.R | Methylation PCA at DMRs                                     |                     
| Figure2c.R | DMR Classification                                          |                     
| Figure3a.R | Mean DNA methylation across CpG positions (Control vs uRPL) |                     
| Figure3b.R | DMR coverage across AluY consensus sequence                 |                     
| Figure3c.R | DMR and overlapping AluY length concordance                 |                     
| Figure3d.R |                                                             |                     
| Figure4a.R |                                                             |                     
| Figure4b.R | GO developmental genes                                      |                     
| Figure5a.R |                                                             |                     
| Figure5b.R |                                                             |                     
| Figure5c.R | Read methylation fraction within DMRs                       |                    
