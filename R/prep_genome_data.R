#!/usr/bin/env Rscript

## prep_genome_data.R
##
## Prepare core genome annotations for human (GRCh38 / hg38):
##  - Protein-coding genes as GRanges
##  - Protein-coding exons as GRanges
##  - Simple genome GRanges with chromosome lengths
##
## Output:
##  - data/gr_genes_pc.rds
##  - data/gr_exons_pc.rds
##  - data/gr_genome_hg38.rds

suppressPackageStartupMessages({
  library(EnsDb.Hsapiens.v86)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

## ---------- Config --------------------------------------------------------

# Where to save RDS files (relative to project root)
output_dir <- "R/data"

# Standard chromosomes to keep
std_chrs <- paste0("chr", c(1:22, "X", "Y", "MT"))

# Filenames
file_genes  <- file.path(output_dir, "gr_genes_pc.rds")
file_exons  <- file.path(output_dir, "gr_exons_pc.rds")
file_genome <- file.path(output_dir, "gr_genome_hg38.rds")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

edb <- EnsDb.Hsapiens.v86

## ---------- Protein-coding genes -----------------------------------------

if (!file.exists(file_genes)) {
  message("Preparing protein-coding gene GRanges...")
  
  filters <- list(
    GeneBiotypeFilter("protein_coding")
  )
  
  genes_pc <- genes(
    edb,
    filter  = filters,
    columns = c("gene_id", "gene_name", "gene_biotype")
  )
  
  # Add "chr" prefix to match typical UCSC-style naming
  seqlevels(genes_pc) <- paste0("chr", seqlevels(genes_pc))
  
  # Keep only standard chromosomes
  genes_pc <- keepSeqlevels(
    genes_pc,
    std_chrs,
    pruning.mode = "coarse"
  )
  
  # Sort for sanity
  genes_pc <- sort(genes_pc)
  
  saveRDS(genes_pc, file_genes)
  message("  -> Saved: ", file_genes)
} else {
  message("Skipping genes: file already exists at ", file_genes)
}

## ---------- Protein-coding exons -----------------------------------------

if (!file.exists(file_exons)) {
  message("Preparing protein-coding exon GRanges...")
  
  filters_exons <- list(
    TxBiotypeFilter("protein_coding"),
    GeneBiotypeFilter("protein_coding")
  )
  
  exons_pc <- exons(
    edb,
    filter  = filters_exons,
    columns = c("exon_id", "gene_id", "tx_id", "gene_name", "tx_biotype")
  )
  
  seqlevels(exons_pc) <- paste0("chr", seqlevels(exons_pc))
  
  exons_pc <- keepSeqlevels(
    exons_pc,
    std_chrs,
    pruning.mode = "coarse"
  )
  
  exons_pc <- sort(exons_pc)
  
  saveRDS(exons_pc, file_exons)
  message("  -> Saved: ", file_exons)
} else {
  message("Skipping exons: file already exists at ", file_exons)
}

## ---------- Genome GRanges (chr lengths) ---------------------------------

if (!file.exists(file_genome)) {
  message("Preparing genome GRanges with chromosome lengths...")
  
  # Use seqlengths from the EnsDb object (after UCSC-style renaming)
  seqinfo_edb <- seqinfo(edb)
  names(seqinfo_edb) <- paste0("chr", names(seqinfo_edb))
  
  # Keep just the standard chromosomes
  seqinfo_std <- seqinfo_edb[std_chrs]
  
  chr_lengths <- seqlengths(seqinfo_std)
  
  gr_genome <- GRanges(
    seqnames = names(chr_lengths),
    ranges   = IRanges(start = 1, end = chr_lengths)
  )
  
  seqinfo(gr_genome) <- seqinfo_std
  
  saveRDS(gr_genome, file_genome)
  message("  -> Saved: ", file_genome)
} else {
  message("Skipping genome: file already exists at ", file_genome)
}

## -------- Prepare membrane protein annotations from UniProt TSV export

# Input and output paths
uniprot_file <- "data/uniprotkb_taxonomy_id_2759_AND_model_or_2025_03_04.tsv"        # <-- change if needed
file_membrane <- "R/data/membrane_proteins.rds"

if (!file.exists(file_membrane)) {
  
  if (!file.exists(uniprot_file)) {
    warning("UniProt TSV not found at ", uniprot_file,
            ". Skipping membrane protein processing.")
  } else {
    
    message("Preparing membrane protein list from UniProt...")
    
    dt <- fread(uniprot_file)
    
    # Expected column names from UniProt export
    rename_map <- c("Gene Names"                 = "gene_names",
                    "Subcellular location [CC]" = "subcellular_location",
                    "Entry Name"                = "entry_name",
                    "Protein names"             = "protein_names")
    
    # Check required columns
    if (!all(names(rename_map) %in% colnames(dt))) {
      stop("ERROR: Required UniProt columns missing. ",
           "Expected: ", paste(names(rename_map), collapse = ", "))
    }
    
    # Rename to cleaner names
    setnames(dt, old = names(rename_map), new = rename_map)
    
    # Extract main gene symbol (first token)
    dt[, gene_name := sub(" .*", "", gene_names)]
    
    # Mark membrane-localized proteins (broad definition)
    dt[, membrane := grepl("membrane", subcellular_location, ignore.case = TRUE)]
    
    # Keep unique membrane genes
    dt_mem <- dt[membrane == TRUE & nchar(gene_name) > 0][!duplicated(gene_name)]
    dt_mem <- dt_mem[order(gene_name)]
    
    message("Membrane protein count detected: ", nrow(dt_mem))
    
    # Save RDS
    saveRDS(dt_mem, file_membrane)
    message("  -> Saved: ", file_membrane)
  }
  
} else {
  message("Skipping membrane proteins: file already exists at ", file_membrane)
}

## ======================================================================
## GTEx preprocessing: melt TPM matrix + annotate with tissue metadata
## ======================================================================

suppressPackageStartupMessages(library(data.table))

file_gtex_input  <- "data/dt_gtex.rds"            # existing GTEx wide-format input
file_gtex_output <- "R/data/gtex_melt.rds"          # output long-format annotated object

if (!file.exists(file_gtex_output)) {
  
  if (!file.exists(file_gtex_input)) {
    warning("GTEx RDS not found at ", file_gtex_input,
            ". Skipping GTEx processing.")
  } else {
    
    message("Preparing GTEx melted expression object...")
    
    dt.gtex <- readRDS(file_gtex_input)
    
    # Restrict to protein-coding genes present in genome GRanges
    gene_ids_pc <- genes_pc$gene_id
    dt.gtex <- dt.gtex[gene %in% gene_ids_pc]
    
    # Clean columns
    dt.gtex[, median_tpm := NULL]
    dt.gtex[, Description := as.factor(Description)]
    dt.gtex[, gene := as.factor(gene)]
    dt.gtex[, Name := NULL]
    
    # Melt (value column will be TPM)
    dt.gtex.melt <- melt(dt.gtex, id.vars = c("gene", "Description"))
    setnames(dt.gtex.melt, "variable", "SAMPID")
    
    # ------------------ Add metadata ---------------------
    metadata_url <- "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
    metadata <- fread(metadata_url)
    metadataFiltered <- metadata[, .(SAMPID, SMTS, SMTSD)]   # major tissue + specific tissue
    
    metadataFiltered[, SAMPID := as.factor(SAMPID)]
    metadataFiltered[, SMTS   := as.factor(SMTS)]
    metadataFiltered[, SMTSD  := as.factor(SMTSD)]
    
    # Ensure correspondence
    stopifnot(all(dt.gtex.melt$SAMPID %in% metadataFiltered$SAMPID))
    
    # Merge annotation
    dt.gtex.melt <- merge(dt.gtex.melt, metadataFiltered, by="SAMPID", all.x=TRUE)
    
    # ------------------ Per-gene summaries ------------------
    dt.gtex.melt[, gene_median := median(value), by = gene]
    dt.gtex.melt[, gene_smts_median := median(value), by = .(gene, SMTS)]
    
    # Save processed output
    saveRDS(dt.gtex.melt, file = file_gtex_output, compress = FALSE)
    message("  -> Saved GTEx melted expression table to: ", file_gtex_output)
  }
  
} else {
  message("Skipping GTEx processing: file already exists at ", file_gtex_output)
}
