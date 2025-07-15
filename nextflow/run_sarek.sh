nextflow run nf-core/sarek --input samplesheet.csv --genome  GATK.GRCh38 --outdir alignment -profile docker --save_mapped   --save_output_as_bam


nextflow run nf-core/sarek --input bam_samplesheet.csv --genome  GATK.GRCh38 --outdir results_strelka --step variant_calling --tools strelka --somatic -profile docker -c nextflow.config
