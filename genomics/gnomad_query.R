#!/usr/bin/env Rscript

###
# this queries donwloaded gnomad files for gnomad hits for a given list of variantss

library(VariantAnnotation)
library(GenomicRanges)
library(data.table)

# 0) parse command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: gnomad_query.R <variants.csv> <output.csv>")
}
vcfin <- args[1]
csvout <- args[2]

# 1) Load your batch of variants (must have columns chr, pos, ref, alt)
variants <- fread(vcfin,colClasses = c(chr="character", pos="integer", ref="character",alt="character"))
variants[, chr := paste0("chr",chr)]

# 2) Build a GRanges of exactly the positions you care about
gr <- GRanges(variants$chr,IRanges(start=variants$pos, end=variants$pos))
mcols(gr)$ref <- variants$ref
mcols(gr)$alt <- variants$alt


# 3) Point to your local gnomAD VCFs  
#    (assumes files named like gnomad.genomes.r3.1.sites.chr1.vcf.bgz, chr2, ..., chrX)
chroms   <- unique(variants$chr)
vcf_paths <- file.path("/mnt/andy/gnomad_exome",sprintf("gnomad.exomes.v4.0.sites.%s.vcf.bgz", chroms))

vcf_list <- lapply(vcf_paths, function(vcf_file) {
  # 1. figure out which chromosome this file is
    contig <- sub(".*\\.sites\\.(chr[^.]+)\\.vcf\\.bgz$", "\\1", vcf_file)
    cat("...", contig, "\n")

  # 2. select just the SNPs on that chromosome
  subset_gr <- gr[seqnames(gr) == contig]
  if (length(subset_gr) == 0) return(NULL)

  # 3. read all VCF records at those positions
  param <- ScanVcfParam(which = subset_gr, info = "AF")
  rv    <- readVcf(vcf_file, "hg38", param)
  rr    <- rowRanges(rv)

# 3) build a table of all VCF rows: include the ID
  dt <- data.table(
    ID    = names(rr),                           # the variant name (rsID)
    CHROM = as.character(seqnames(rr)),
    POS   = start(rr),
    REF   = as.character(ref(rv)),
    ALT   = sapply(alt(rv), function(x) as.character(x)[1]),
    AF    = unlist(info(rv)$AF)
  )

  # 4) build the true SNP table to join against
  snp_dt <- data.table(
    CHROM = as.character(seqnames(subset_gr)),
    POS   = start(subset_gr),
    REF   = subset_gr$ref,
    ALT   = subset_gr$alt
  )

  # 5) join on all four keys, bringing along ID and AF
  setkey(dt,    CHROM, POS, REF, ALT)
  setkey(snp_dt, CHROM, POS, REF, ALT)
  dt[snp_dt, nomatch=0L]  # this will have columns ID, CHROM, POS, REF, ALT, AF
})

# 6) combine per-chromosome results
matched <- rbindlist(Filter(Negate(is.null), vcf_list))
                                        # 7) check we got one row per input SNP
# actually dont, because sometiems we are asking for a bad variant
#stopifnot(nrow(matched) == length(gr))

fwrite(matched, file=csvout, sep    = ",",
  quote  = FALSE,
  bom    = FALSE)
