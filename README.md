# DADA2_Pipeline

An automated DADA2 pipeline for 16S rRNA paired-end amplicon sequencing. Built for V3-V4 region data (Illumina paired-end), but works fine for V4 with minor parameter adjustments.

This was written as part of my master's thesis project on gut microbiome characterization. The goal was to have a single script that takes raw FASTQs all the way to an ASV table + taxonomy, without having to re-run things manually every time. Ended up being pretty reusable so posting it here.

---

## What it does

Two-step workflow:
1. **QC first** - generates quality profile plots so you can decide on trimming parameters
2. **Pipeline run** - filtering, denoising, merging, chimera removal, taxonomy

The pipeline is modular (each step is its own function), so you can also run individual steps if something goes wrong or you want to tweak one part.

---

## Requirements

- R >= 4.0
- `dada2` package (Bioconductor)
- `parallel` (base R, already included)

Install dada2 if you haven't:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")
```

You also need the Silva reference databases for taxonomy assignment. Download them from Zenodo:
- [Silva v138 training set](https://zenodo.org/record/4587955)

Place the `.fa.gz` files in your project root and update the paths in the config section at the top of the script if needed.

---

## Project structure

```
project/
├── data/               # your raw FASTQ files go here
├── output/
│   ├── qc/             # quality profile PDFs, error rate plots
│   ├── filtered/       # filtered reads (intermediate, can be deleted after)
│   └── results/        # final outputs
├── dada2_pipeline.R
├── silva_nr99_v138.2_toGenus_trainset.fa.gz
└── silva_species_assignment_v138.1.fa.gz
```

FASTQ files should follow standard Illumina naming: `samplename_R1.fastq` / `samplename_R2.fastq`. If your files use a different pattern (e.g. `_1.fastq`, `_2.fastq`), just change `pattern_fwd` and `pattern_rev` in `discover_fastq_files()`.

---

## Usage

### Step 1 - QC

```r
source("dada2_pipeline.R")

explore_and_suggest(
  data_dir = "data/",
  n_samples = 4           # how many samples to plot
)
```

This generates two PDFs in `output/qc/`:
- `quality_profiles_FORWARD.pdf`
- `quality_profiles_REVERSE.pdf`

Open them and look at where quality starts dropping. That's where you set your `truncLen` values for step 2. A good rule of thumb: truncate where median quality falls below Q30. Make sure the combined length still gives you enough overlap with your amplicon (default assumes V3-V4, ~460 bp).

### Step 2 - Run the pipeline

```r
results <- run_dada2_pipeline(
  truncLen_fwd = 240,    # set based on your quality profiles
  truncLen_rev = 160,    # set based on your quality profiles
  maxEE_fwd = 2,
  maxEE_rev = 2
)
```

That's it. When it's done, `output/results/` will contain:

| File | Description |
|------|-------------|
| `asv_table.csv` | ASV counts per sample |
| `asv_sequences.fasta` | FASTA with all ASV sequences |
| `taxonomy.csv` | taxonomy assignment for each ASV |
| `read_tracking.csv` | reads at each step (good for QC reporting) |
| `seqtab_nochim.rds` | R object for downstream analysis |
| `taxa.rds` | taxonomy as R object |
| `track.rds` | tracking table as R object |

---

## Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `truncLen_fwd` | required | truncation length for forward reads |
| `truncLen_rev` | required | truncation length for reverse reads |
| `maxEE_fwd` | 2 | max expected errors, forward |
| `maxEE_rev` | 2 | max expected errors, reverse |
| `trimLeft` | c(0, 0) | bp to trim from 5' end (use if primers weren't removed) |
| `pool` | "pseudo" | pooling strategy for ASV inference (FALSE / "pseudo" / TRUE) |
| `n_threads` | auto | defaults to all cores - 1 |

For V4 only (instead of V3-V4), set `amplicon_length = 253` inside `filter_and_trim_reads()` and `explore_and_suggest()`.

---

## Notes

- **Primer removal**: this pipeline assumes primers have already been removed (e.g. with cutadapt). If they haven't, use the `trimLeft` parameter to trim the appropriate number of bp from each end.
- **Pooling**: `pool = "pseudo"` is a good default. Use `pool = TRUE` if you care about detecting rare taxa, but it's noticeably slower on large datasets.
- **Low chimera retention warning**: if less than 50% of reads survive chimera removal, the most common culprit is leftover primers. Check with cutadapt first.

---

## Downstream

The `.rds` files are ready to load into phyloseq for community analysis:

```r
library(phyloseq)

seqtab <- readRDS("output/results/seqtab_nochim.rds")
taxa   <- readRDS("output/results/taxa.rds")

ps <- phyloseq(
  otu_table(seqtab, taxa_are_rows = FALSE),
  tax_table(taxa)
)
```

---

## References

Callahan et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581-583. https://doi.org/10.1038/nmeth.3869

Quast et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Research*, 41, D590-D596.
