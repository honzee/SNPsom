# Wrapper functions

# Trimm Illumina Universal Adapters and run fastqc
trimmomatic_fastqc <- function(metadata, input_dir, output_dir, t = 2, min_length = 20) {

  # Perform Illumina Universal Adpater trimming
  dir.create(paste0(output_dir, "/trimmed_fastq/"))
  output_dir = paste0(output_dir, "/trimmed_fastq/")
  trimmomatic(metadata, input_dir, output_dir, t, min_length)

  # Perform fastqc
  metadata$fastq1 = sub(".fastq", "_trimmed.fastq", metadata$fastq1)
  metadata$fastq2 = sub(".fastq", "_trimmed.fastq", metadata$fastq2)
  dir.create(paste0(output_dir, "/fastqc/"))
  output_dir = paste0(output_dir, "/fastqc/")
  fastqc(metadata, input_dir, output_dir, t)
}

# Perform somatic variant analysis from pre-processed fastq files
SNPsom_wrapper <- function(metadata, input_dir, output_dir, reference, t, parallel, reference_version, annotation_dir) {

  # Perform alignment with bwa-mem2
  dir.create(paste0(output_dir, "/alignment/"))
  output_dir_alignment = paste0(output_dir, "/alignment/")
  bwa_mem2(metadata, input_dir, output_dir_alignment, reference, t)

  # Mark duplicates with MarkDuplicatesSpark
  dir.create(paste0(output_dir, "/MarkDuplicates/"))
  output_dir_MarkDuplicates = paste0(output_dir, "/MarkDuplicates/")
  MarkDuplicatesSpark(metadata, output_dir_alignment, output_dir_MarkDuplicates, t)

  # Call somatic SNV with Mutect2 and filter variants with FilterMutectcalls
  dir.create(paste0(output_dir, "/vcf/"))
  output_dir_vcf = paste0(output_dir, "/vcf/")
  Mutect2(metadata, output_dir_MarkDuplicates, output_dir_vcf, reference, parallel)
  FilterMutectCalls(metadata, output_dir_vcf, output_dir_vcf, reference)

  # Annotate variants - to be tested
  #Funcotator(metadata, output_dir_vcf, output_dir_vcf, reference, reference_version, annotation_dir)
}
