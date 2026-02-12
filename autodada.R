# dada2_pipeline.R
# author: leom00
# 16S rRNA paired-end amplicon pipeline (V3-V4)
#
# two-step workflow:
#   1. explore_and_suggest()   --> check quality, pick truncLen
#   2. run_dada2_pipeline()    --> run everything

# ---- load libraries ----------------------------------------------------------

suppressPackageStartupMessages({
  library(dada2)
  library(parallel)
})

cat("pipeline loaded. run explore_and_suggest() first!\n\n")

# ---- config ------------------------------------------------------------------

# directories
DATA_DIR     <- "data/"
OUTPUT_DIR   <- "output/"
FILTERED_DIR <- file.path(OUTPUT_DIR, "filtered/")
RESULTS_DIR  <- file.path(OUTPUT_DIR, "results/")
QC_DIR       <- file.path(OUTPUT_DIR, "qc/")

# create dirs if missing
for (dir_path in c(DATA_DIR, OUTPUT_DIR, FILTERED_DIR, RESULTS_DIR, QC_DIR)) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    cat("Created directory:", dir_path, "\n")
  }
}

cat("\n")

# performance
N_THREADS        <- detectCores() - 1  
MEMORY_EFFICIENT <- TRUE 

# Silva v138 - update paths if needed
TAXONOMY_DB         <- "silva_nr99_v138.2_toGenus_trainset.fa.gz"
TAXONOMY_SPECIES_DB <- "silva_species_assignment_v138.1.fa.gz"

cat("config:\n")
cat("  - threads:", N_THREADS, "\n")
cat("  - memory-efficient mode:", MEMORY_EFFICIENT, "\n")
cat("  - data dir:", DATA_DIR, "\n")
cat("  - output dir:", OUTPUT_DIR, "\n\n")

# ---- helpers -----------------------------------------------------------------

`%R%` <- function(char, n) paste0(rep(char, n), collapse = "")

# ---- 1. file discovery -------------------------------------------------------

discover_fastq_files <- function(data_dir, 
                                 pattern_fwd = "_R1.fastq", 
                                 pattern_rev = "_R2.fastq") {
  # finds and pairs forward/reverse FASTQ files in data_dir
  # assumes standard Illumina naming (samplename_R1.fastq / _R2.fastq)
  # change pattern_fwd/rev if your files use a different convention
  
  cat("looking for FASTQ files...\n")
  
  fnFs <- sort(list.files(data_dir, pattern = pattern_fwd, full.names = TRUE))
  fnRs <- sort(list.files(data_dir, pattern = pattern_rev, full.names = TRUE))
  
  if (length(fnFs) == 0 || length(fnRs) == 0) {
    stop("ERROR: No FASTQ files found. Check data directory and patterns.")
  }
  
  if (length(fnFs) != length(fnRs)) {
    stop("ERROR: Unequal number of forward and reverse files!")
  }
  
  # takes everything before the first "_" as sample name
  sample_names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  sample_names<- make.unique(sample_names)
  
  cat("    Found", length(fnFs), "paired-end samples\n")
  cat("    Sample names:", paste(head(sample_names, 3), collapse = ", "), "...\n\n")
  
  return(list(
    forward = fnFs,
    reverse = fnRs,
    sample_names = sample_names
  ))
}

# ---- 2. QC - run this first! -------------------------------------------------

explore_and_suggest <- function(data_dir = "data/", 
                                output_dir = "qc/",
                                n_samples = 4,
                                min_quality = 30,
                                amplicon_length = 460) { # adjust for your region (253 for V4)
  # generates quality profile PDFs so you can decide on truncLen values
  # look at where quality drops below Q30 and truncate there
  # make sure fwd + rev - amplicon_length gives you enough overlap (>20 bp min)
  
  cat("\n")
  cat("=" %R% 80, "\n")
  cat("STEP 1: quality assessment\n")
  cat("=" %R% 80, "\n\n")
  
  files <- discover_fastq_files(data_dir)
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  n_plot <- min(n_samples, length(files$forward))
  
  cat("generating quality profiles for", n_plot, "samples...\n\n")
  
  # forward quality profiles
  pdf_fwd <- file.path(output_dir, "quality_profiles_FORWARD.pdf")
  tryCatch({
    pdf(pdf_fwd, width = 12, height = 7)
    print(plotQualityProfile(files$forward[1:n_plot]))
    dev.off()
    cat("  forward profiles saved:", pdf_fwd, "\n")
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("  could not generate forward profiles:", e$message, "\n")
  })
  
  # reverse quality profiles
  pdf_rev <- file.path(output_dir, "quality_profiles_REVERSE.pdf")
  tryCatch({
    pdf(pdf_rev, width = 12, height = 7)
    print(plotQualityProfile(files$reverse[1:n_plot]))
    dev.off()
    cat("  reverse profiles saved:", pdf_rev, "\n\n")
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("  could not generate reverse profiles:", e$message, "\n\n")
  })
  
  
  
  # return file paths for step 2
  invisible(list(
    files = files,
    qc_dir = output_dir,
    pdf_forward = pdf_fwd,
    pdf_reverse = pdf_rev
  ))
}

# ---- 3. filter & trim --------------------------------------------------------

filter_and_trim_reads <- function(fnFs, fnRs, sample_names, filtered_dir,
                                  truncLen_fwd, truncLen_rev,
                                  maxEE_fwd = 2, maxEE_rev = 2,
                                  trimLeft = c(0, 0), n_threads = 1) {
  # standard DADA2 filtering
  # truncLen_fwd and truncLen_rev are required - set them after checking QC plots
  # use trimLeft if primers weren't removed upstream (e.g. with cutadapt)
  
  cat("step 3: filtering and trimming reads...\n")
  
  # truncLen is required - can't guess this for you
  if (missing(truncLen_fwd) || missing(truncLen_rev)) {
    stop(paste0(
      "\ntrunLen_fwd and truncLen_rev are required!\n\n",
      "run explore_and_suggest() first to check quality profiles.\n"
    ))
  }
  
  # check overlap - important for merging later
  amplicon_length <- 460  # V3-V4 region default
  overlap <- truncLen_fwd + truncLen_rev - amplicon_length
  
  cat("    parameters:\n")
  cat("      - forward truncLen:", truncLen_fwd, "\n")
  cat("      - reverse truncLen:", truncLen_rev, "\n")
  cat("      - expected overlap:", overlap, "bp\n")
  cat("      - maxEE (F, R):", maxEE_fwd, ",", maxEE_rev, "\n")
  cat("      - trimLeft (F, R):", trimLeft[1], ",", trimLeft[2], "\n")
  
  if (overlap < 20) {
    warning(paste0(
      "\noverlap (", overlap, " bp) is below 20 bp - merging will likely fail!\n",
      "try increasing truncLen or check that amplicon_length is correct for your region\n"
    ))
  } else if (overlap < 50) {
    cat("      overlap is a bit tight but should be ok\n")
  } else {
    cat("      overlap looks good\n")
  }
  
  cat("\n")
  
  dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
  
  filtFs <- file.path(filtered_dir, paste0(sample_names, "_F_filt.fastq.gz"))
  filtRs <- file.path(filtered_dir, paste0(sample_names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample_names
  names(filtRs) <- sample_names
  
  out <- filterAndTrim(
    fnFs, filtFs, fnRs, filtRs,
    truncLen = c(truncLen_fwd, truncLen_rev),
    trimLeft = trimLeft,
    maxN = 0,                    # DADA2 requires no Ns
    maxEE = c(maxEE_fwd, maxEE_rev),
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = n_threads,
    verbose = TRUE
  )
  
  retention_rate <- mean(out[,2] / out[,1]) * 100
  
  cat("\n    filtering done for", nrow(out), "samples\n")
  cat("    mean retention rate:", round(retention_rate, 2), "%\n")
  
  if (retention_rate < 50) {
    warning(paste0(
      "\nretention rate is low (<50%)!\n",
      "possible causes:\n",
      "  1. truncLen too aggressive\n",
      "  2. maxEE too strict\n",
      "  3. poor sequencing quality\n",
      "go back and check quality profiles\n"
    ))
  } else if (retention_rate > 90) {
    cat("    excellent retention rate\n")
  } else {
    cat("    good retention rate\n")
  }
  
  cat("\n")
  
  # flag samples with zero reads after filtering
  failed_samples <- rownames(out)[out[,2] == 0]
  if (length(failed_samples) > 0) {
    warning("The following samples had zero reads after filtering: ",
            paste(failed_samples, collapse = ", "))
  }
  
  return(list(
    stats = out,
    filtered_forward = filtFs,
    filtered_reverse = filtRs
  ))
}

# ---- 4. error rate learning --------------------------------------------------

learn_error_rates <- function(filtFs, filtRs, n_threads = 1, 
                              plot_errors = TRUE, output_dir = NULL) {
  # learns the sequencing error model from the data
  # always check the plots - lines should follow the dots reasonably well
  
  cat("step 4: learning error rates (grab a coffee, this takes a while)...\n")
  
  # skip samples that got filtered out completely
  valid_idx    <- file.exists(filtFs) & file.exists(filtRs)
  filtFs_valid <- filtFs[valid_idx]
  filtRs_valid <- filtRs[valid_idx]
  
  if (length(filtFs_valid) == 0) {
    stop("ERROR: No valid filtered files found!")
  }
  
  errF <- learnErrors(filtFs_valid, multithread = n_threads, verbose = 1)
  cat("    forward error rates learned\n")
  
  errR <- learnErrors(filtRs_valid, multithread = n_threads, verbose = 1)
  cat("    reverse error rates learned\n\n")
  
  # save diagnostic plots
  if (plot_errors && !is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    tryCatch({
      pdf(file.path(output_dir, "error_rates_forward.pdf"), width = 10, height = 10)
      print(plotErrors(errF, nominalQ = TRUE))
      dev.off()
      
      pdf(file.path(output_dir, "error_rates_reverse.pdf"), width = 10, height = 10)
      print(plotErrors(errR, nominalQ = TRUE))
      dev.off()
      
      cat("    error rate plots saved to", output_dir, "\n\n")
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      warning("Could not generate error plots: ", e$message)
    })
  }
  
  return(list(forward = errF, reverse = errR))
}

# ---- 5. dereplication --------------------------------------------------------

dereplicate_sequences <- function(filtFs, filtRs, sample_names) {
  # collapse identical sequences - reduces memory usage significantly
  
  cat("step 5: dereplicating sequences...\n")
  
  # skip samples with no filtered reads
  valid_idx          <- file.exists(filtFs) & file.exists(filtRs)
  filtFs_valid       <- filtFs[valid_idx]
  filtRs_valid       <- filtRs[valid_idx]
  sample_names_valid <- sample_names[valid_idx]
  
  if (length(filtFs_valid) == 0) {
    stop("ERROR: No valid samples for dereplication!")
  }
  
  derepFs <- derepFastq(filtFs_valid, verbose = FALSE)
  derepRs <- derepFastq(filtRs_valid, verbose = FALSE)
  
  names(derepFs) <- sample_names_valid
  names(derepRs) <- sample_names_valid
  
  cat("    dereplication done for", length(derepFs), "samples\n\n")
  
  return(list(forward = derepFs, reverse = derepRs))
}

# ---- 6. ASV inference --------------------------------------------------------

infer_asvs <- function(derepFs, derepRs, errF, errR, n_threads = 1,
                       pool = "pseudo") {
  # core DADA2 step
  # pool = "pseudo" is a good default - faster than TRUE but catches more rare taxa than FALSE
  
  cat("step 6: inferring ASVs (pooling:", pool, ")...\n")
  cat("    this is the most computationally intensive step...\n")
  
  dadaFs <- dada(derepFs, err = errF, multithread = n_threads, 
                 pool = pool, verbose = 1)
  cat("    forward ASV inference done\n")
  
  dadaRs <- dada(derepRs, err = errR, multithread = n_threads,
                 pool = pool, verbose = 1)
  cat("    reverse ASV inference done\n\n")
  
  return(list(forward = dadaFs, reverse = dadaRs))
}

# ---- 7. merge paired reads ---------------------------------------------------

merge_paired_reads <- function(dadaFs, dadaRs, derepFs, derepRs,
                               minOverlap = 12, maxMismatch = 0) {
  # merge F and R reads into full amplicons
  
  cat("step 7: merging paired reads...\n")
  
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                        minOverlap = minOverlap,
                        maxMismatch = maxMismatch,
                        verbose = TRUE)
  
  cat("    merging done\n\n")
  
  return(mergers)
}

# ---- 8. sequence table -------------------------------------------------------

construct_sequence_table <- function(mergers) {
  # builds the ASV table (samples x ASVs)
  # check sequence length distribution - should cluster around expected amplicon size
  
  cat("step 8: constructing sequence table...\n")
  
  seqtab <- makeSequenceTable(mergers)
  
  cat("    sequence table dimensions:", dim(seqtab)[1], "samples x", 
      dim(seqtab)[2], "ASVs\n")
  
  seq_lengths <- nchar(getSequences(seqtab))
  cat("    sequence length distribution:\n")
  print(table(seq_lengths))
  cat("\n")
  
  return(seqtab)
}

# ---- 9. chimera removal ------------------------------------------------------

remove_chimeras <- function(seqtab, method = "consensus", n_threads = 1) {
  # removes bimeric sequences
  # if >50% of reads are removed here, check that primers were properly trimmed
  
  cat("step 9: removing chimeras...\n")
  
  seqtab_nochim <- removeBimeraDenovo(seqtab, method = method,
                                      multithread = n_threads,
                                      verbose = TRUE)
  
  pct_retained_seqs  <- (ncol(seqtab_nochim) / ncol(seqtab)) * 100
  pct_retained_reads <- (sum(seqtab_nochim) / sum(seqtab)) * 100
  
  cat("    chimera removal done\n")
  cat("    Retained", ncol(seqtab_nochim), "of", ncol(seqtab), 
      "ASVs (", round(pct_retained_seqs, 2), "%)\n")
  cat("    Retained", sum(seqtab_nochim), "of", sum(seqtab),
      "reads (", round(pct_retained_reads, 2), "%)\n")
  
  if (pct_retained_reads < 50) {
    warning(paste0(
      "\n<50% reads retained after chimera removal!\n",
      "possible causes:\n",
      "  1. primers not properly removed\n",
      "  2. too many PCR cycles (favors chimera formation)\n",
      "  3. unusual amplicon design\n"
    ))
  } else {
    cat("    chimera removal looks good\n")
  }
  
  cat("\n")
  
  return(seqtab_nochim)
}

# ---- 10. read tracking -------------------------------------------------------

track_reads <- function(filter_stats, dadaFs, dadaRs, mergers, seqtab_nochim,
                        sample_names) {
  # tracks reads through each step - useful for QC reporting
  # big drops at any step usually indicate a problem
  
  cat("step 10: tracking reads through pipeline...\n")
  
  getN <- function(x) sum(getUniques(x))
  
  track <- cbind(
    filter_stats,
    sapply(dadaFs, getN),
    sapply(dadaRs, getN),
    sapply(mergers, getN),
    rowSums(seqtab_nochim)
  )
  
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", 
                       "merged", "nonchim")
  rownames(track) <- sample_names
  
  cat("    read tracking table:\n")
  print(head(track))
  cat("\n")
  
  cat("    summary:\n")
  cat("      mean retention:", 
      round(mean(track[,"nonchim"] / track[,"input"]) * 100, 2), "%\n")
  cat("      median retention:",
      round(median(track[,"nonchim"] / track[,"input"]) * 100, 2), "%\n\n")
  
  return(track)
}

# ---- 11. taxonomy ------------------------------------------------------------

assign_taxonomy <- function(seqtab_nochim, taxonomy_db, 
                            species_db = NULL, n_threads = 1) {
  # naive Bayesian classifier against Silva
  # genus-level first, then species if db is available
  
  cat("step 11: assigning taxonomy...\n")
  
  if (!file.exists(taxonomy_db)) {
    warning("Taxonomy database not found at: ", taxonomy_db)
    cat("    Skipping taxonomic assignment.\n")
    cat("    Download Silva database from:\n")
    cat("    https://zenodo.org/record/4587955\n\n")
    return(NULL)
  }
  
  cat("    Using database:", taxonomy_db, "\n")
  cat("    This step may take 10-30 minutes depending on dataset size...\n")
  
  taxa <- assignTaxonomy(seqtab_nochim, taxonomy_db, 
                         multithread = n_threads, verbose = TRUE)
  
  cat("    Genus-level assignment completed\n")
  
  if (!is.null(species_db) && file.exists(species_db)) {
    cat("    Adding species-level annotation...\n")
    taxa <- addSpecies(taxa, species_db, verbose = TRUE)
    cat("    Species-level assignment completed\n")
  }
  
  cat("\n    Taxonomy summary:\n")
  print(head(taxa))
  cat("\n")
  
  return(taxa)
}

# ---- 12. save results --------------------------------------------------------

save_results <- function(seqtab_nochim, taxa, track, results_dir) {
  # saves everything in RDS (for R) and CSV (for everything else)
  
  cat("step 12: saving results...\n")
  
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  
  # RDS - preserves exact R objects
  saveRDS(seqtab_nochim, file.path(results_dir, "seqtab_nochim.rds"))
  cat("    Saved: seqtab_nochim.rds\n")
  
  if (!is.null(taxa)) {
    saveRDS(taxa, file.path(results_dir, "taxa.rds"))
    cat("    Saved: taxa.rds\n")
  }
  
  saveRDS(track, file.path(results_dir, "track.rds"))
  cat("    Saved: track.rds\n")
  
  # CSV - human-readable, works with Excel/Python
  
  # replace sequences with ASV IDs in the table
  asv_seqs <- colnames(seqtab_nochim)
  asv_ids  <- paste0("ASV", seq_along(asv_seqs))
  colnames(seqtab_nochim) <- asv_ids
  write.csv(seqtab_nochim, file.path(results_dir, "asv_table.csv"))
  cat("    Saved: asv_table.csv\n")
  
  # FASTA with ASV sequences
  fasta_file <- file.path(results_dir, "asv_sequences.fasta")
  fasta_conn <- file(fasta_file, "w")
  for (i in seq_along(asv_ids)) {
    writeLines(paste0(">", asv_ids[i]), fasta_conn)
    writeLines(asv_seqs[i], fasta_conn)
  }
  close(fasta_conn)
  cat("    Saved: asv_sequences.fasta\n")
  
  # taxonomy table
  if (!is.null(taxa)) {
    taxa_df          <- as.data.frame(taxa)
    taxa_df$ASV_ID   <- asv_ids
    taxa_df$Sequence <- asv_seqs
    write.csv(taxa_df, file.path(results_dir, "taxonomy.csv"), row.names = FALSE)
    cat("    Saved: taxonomy.csv\n")
  }
  
  # tracking table
  write.csv(track, file.path(results_dir, "read_tracking.csv"))
  cat("    Saved: read_tracking.csv\n")
  
  cat("\n    All results saved to:", results_dir, "\n\n")
}

# ---- main pipeline -----------------------------------------------------------

run_dada2_pipeline <- function(data_dir = "data/",
                               output_dir = "output/",
                               truncLen_fwd,         # required - set after QC
                               truncLen_rev,         # required - set after QC
                               maxEE_fwd = 2,
                               maxEE_rev = 2,
                               trimLeft = c(0, 0),
                               taxonomy_db = TAXONOMY_DB,
                               species_db = TAXONOMY_SPECIES_DB,
                               n_threads = NULL,
                               pool = "pseudo") {
  # runs the full pipeline end-to-end
  # truncLen values are required - run explore_and_suggest() first!
  
  if (missing(truncLen_fwd) || missing(truncLen_rev)) {
    stop(paste0(
      "\n", "=" %R% 80, "\n",
      "ERROR: truncLen_fwd and truncLen_rev are REQUIRED!\n",
      "=" %R% 80, "\n\n",
      "You must run explore_and_suggest() FIRST to determine these values.\n\n",
      "Workflow:\n",
      "  1. explore_and_suggest()              # Generates quality profiles\n",
      "  2. Review PDFs in qc/ directory\n",
      "  3. run_dada2_pipeline(truncLen_fwd = X, truncLen_rev = Y)\n\n"
    ))
  }
  
  start_time <- Sys.time()
  
  cat("\n")
  cat("=" %R% 80, "\n")
  cat("STEP 2: DADA2 PIPELINE EXECUTION\n")
  cat("=" %R% 80, "\n\n")
  
  if (is.null(n_threads)) {
    n_threads <- max(1, detectCores() - 1)
  }
  
  filtered_dir <- file.path(output_dir, "filtered/")
  results_dir  <- file.path(output_dir, "results/")
  qc_dir       <- file.path(output_dir, "qc/")
  
  files <- discover_fastq_files(data_dir)
  
  filtered <- filter_and_trim_reads(
    files$forward, files$reverse, files$sample_names,
    filtered_dir, truncLen_fwd, truncLen_rev,
    maxEE_fwd, maxEE_rev, trimLeft, n_threads
  )
  
  errors <- learn_error_rates(
    filtered$filtered_forward,
    filtered$filtered_reverse,
    n_threads, TRUE, qc_dir
  )
  
  derep <- dereplicate_sequences(
    filtered$filtered_forward,
    filtered$filtered_reverse,
    files$sample_names
  )
  
  dada_results <- infer_asvs(
    derep$forward, derep$reverse,
    errors$forward, errors$reverse,
    n_threads, pool
  )
  
  mergers <- merge_paired_reads(
    dada_results$forward, dada_results$reverse,
    derep$forward, derep$reverse
  )
  
  seqtab        <- construct_sequence_table(mergers)
  seqtab_nochim <- remove_chimeras(seqtab, "consensus", n_threads)
  
  track <- track_reads(
    filtered$stats, dada_results$forward,
    dada_results$reverse, mergers,
    seqtab_nochim, files$sample_names
  )
  
  taxa <- NULL
  if (!is.null(taxonomy_db) && file.exists(taxonomy_db)) {
    taxa <- assign_taxonomy(seqtab_nochim, taxonomy_db, species_db, n_threads)
  }
  
  save_results(seqtab_nochim, taxa, track, results_dir)
  
  end_time <- Sys.time()
  elapsed  <- difftime(end_time, start_time, units = "mins")
  
  cat("=" %R% 80, "\n")
  cat("PIPELINE COMPLETED SUCCESSFULLY!\n")
  cat("=" %R% 80, "\n")
  cat("Total runtime:", round(elapsed, 2), "minutes\n")
  cat("Final ASV count:", ncol(seqtab_nochim), "\n")
  cat("Results saved to:", results_dir, "\n")
  cat("=" %R% 80, "\n\n")
  
  return(list(
    seqtab       = seqtab_nochim,
    taxonomy     = taxa,
    tracking     = track,
    elapsed_time = elapsed
  ))
}

