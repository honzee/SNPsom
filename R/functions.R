#' Perform fastqc on provided files
#'
#'@param metadata data frame;
#'@param input_dir character; path to the input directory
#'@param output_dir character; path to output directory
#'@param t numeric; number of threads to use; default=2

fastqc = function(metadata, input_dir, output_dir, t = 2)  {
  # create a vector of file names
  files = c(metadata$fastq1, metadata$fastq2)
  system(paste0("docker run --rm -i -v ", input_dir, ":", input_dir, " -v ", output_dir, ":", output_dir, " biocontainers/fastqc:v0.11.9_cv8 fastqc ", paste0(input_dir, files, collapse = " "), " -t ", t, " -o ", output_dir))
}

#' Perform adapter trimming
#'
#' Perform Illumina Universal adapter trimming and filter all reads under selected length (default 20).
#'
#'@param metadata data frame;
#'@param input_dir character; path to the input directory
#'@param output_dir character; path to output directory
#'@param t numeric; number of threads to use; default=2
#'@param min_length numeric; minimum read length to be kept; default=20

trimmomatic = function(metadata, input_dir, output_dir, t = 2, min_length = 20)  {
  for (i in 1:nrow(metadata)) {
    system(paste0("docker run --rm -i -v ", input_dir, ":", input_dir, " -v ", output_dir, ":", output_dir,
                  " honzik1/trimmomatic:0.39 TrimmomaticPE -threads ", t, " ", paste0(input_dir, metadata$fastq1[i]),
                  " ", paste0(input_dir, metadata$fastq2[i]), " ",
                  paste0(output_dir, sub(".fastq", "", c(metadata$fastq1[i], metadata$fastq2[i]))[c(1,1,2,2)],
                         rep(c("_trimmed.fastq", "_unpaired.fastq"), 2), collapse = " "),
                  " ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10:5:true MINLEN:", min_length))
  }
  system(paste0("rm ", output_dir, "/*_unpaired.fastq"))
}

#' Align the reads with bwa-mem2 and directly output bam file
#'
#' Bwa-mem2 is faster than bwa-mem and produces the same results. This function also assigns Read group ID and sample name to the bam header,
#' which is necessary for downstream GATK4 analysis. Checks for reference index file. If it does not exit, it is automatically created in the reference
#' directory.
#'
#'@param metadata data frame;
#'@param input_dir character; path to the input directory
#'@param output_dir character; path to output directory
#'@param reference character; path to the fasta genome reference file. Bwa-mem2 created index must be in the same directory.
#'@param t numeric; number of threads to use; default-2

bwa_mem2 = function(metadata, input_dir, output_dir, reference, t = 2) {

  if(!file.exists(sub(".fa$|.fasta$", ".bwt.2bit.64", reference))) {
    print("NO INDEX FILES FOR REFERENCE DETECTED, REFERENCE INDEX WILL BE AUTOMATICALLY GENERATED")
    system(paste0("docker run --rm -i -v ", dirname(reference), ":", dirname(reference),
                  " honzik1/bwa-mem2_samtools:2.2.1_1.13 bwa-mem2 index -p ", sub(".fa$|.fasta$", "", reference), " ", reference))
  }

  for (i in 1:nrow(metadata)) {
    system(paste0("docker run --rm -i -v ", input_dir, ":", input_dir, " -v ", output_dir, ":", output_dir, " -v ",
                  dirname(reference), ":", dirname(reference), " honzik1/bwa-mem2_samtools:2.2.1_1.13 /bin/bash -c 'bwa-mem2 mem -t ", t, " -M ",
                  "-R \"@RG\tID:", metadata$sample[i], "\tSM:", metadata$sample[i], "\tLB:lib1\tPL:illumina\" ",
                  sub(".fa|.fasta", "", reference), " ", input_dir, metadata$fastq1[i], " ", input_dir, metadata$fastq2[i], " | ",
                  "samtools sort -@ ", t, " -o ", output_dir, metadata$sample[i], ".bam -'"))
  }
}

#' Mark duplicates with MarkDuplicatesSpark
#'
#' MarkDuplicatesSpark is faster than standard MarkDuplicates from picard tools and produces same results
#'
#'@param metadata data frame;
#'@param input_dir character; path to the input directory
#'@param output_dir character; path to output directory
#'@param reference character; path to the fasta genome reference file. Bwa-mem2 created index must be in the same directory.
#'@param t numeric; number of threads to use; default-2

MarkDuplicatesSpark = function(metadata, input_dir, output_dir, t = 2) {
  for (i in 1:nrow(metadata)) {
    system(paste0("docker run --rm -i -v ", input_dir, ":", input_dir, " -v ", output_dir, ":", output_dir,
                  " broadinstitute/gatk:4.2.2.0 gatk MarkDuplicatesSpark -I ", input_dir, metadata$sample[i],
                  ".bam", " -O ", output_dir, metadata$sample[i], "_markdup.bam", " --conf 'spark.executor.cores=", t, "'"))
  }
}

#' Detect variants with Mutect2
#'
#' Run Mutect2 with the option to run multiple samples in parallel. Chekcs for
#' GATK sequence dictionary and samtools index, if the files do not exist, they
#' will be automatically created.
#'
#'@param metadata data frame;
#'@param input_dir character; path to the input directory
#'@param output_dir character; path to output directory
#'@param reference character; path to the fasta genome reference file. Bwa-mem2 created index must be in the same directory.
#'@param t numeric; number of threads to use; default-2
#'@param parallel; how many samples should be ran in parallel, each one will need around 15 GB of RAM

Mutect2 = function(metadata, input_dir, output_dir, reference, parallel) {

  # check for reference dictionary file and samtools index
  if(!file.exists(sub(".fa$|.fasta$", ".dict", reference))) {
    print("NO DICTIONARY FILE FOR REFERENCE DETECTED, REFERENCE DICTIONARY WILL BE AUTOMATICALLY GENERATED")
    system(paste0("docker run --rm -i -v ", dirname(reference), ":", dirname(reference),
                  " broadinstitute/gatk:4.2.2.0 gatk CreateSequenceDictionary -R ", reference))
  }
  if(!file.exists(paste0(reference, ".fai"))) {
    print("NO FASTA INDEX FOR REFERENCE DETECTED, FASTA INDEX WILL BE AUTOMATICALLY GENERATED")
    system(paste0("docker run --rm -i -v ", dirname(reference), ":", dirname(reference),
                  " broadinstitute/gatk:4.2.2.0 samtools faidx ", reference))
  }

  # check if the input bam files exist and change the metadata table accordingly
  files_exist = file.exists(paste0(input_dir, metadata$sample, "_markdup.bam"))
  if (any(!files_exist)) {
    print(paste0("Files: ", paste0(input_dir, metadata$sample[!files_exist], "_markdup.bam", collapse = ', '), " are missing. Analysis will be modified accordingly"))
    metadata = metadata[files_exist, ]
  }

  metadata_t = metadata[metadata$status == "t", ]
  metadata_c = metadata[metadata$status == "c", ]
  # check whether all tumour samples have thier control and if not filter them out
  metadata_t = metadata_t[metadata_t$patient %in% metadata_c$patient, ]

  for (i in 1:ceiling(nrow(metadata_t)/parallel)) {
    if (nrow(metadata_t) != 0) {
      sample_start = ((i-1)*parallel)+1
      if (i == ceiling(nrow(metadata_t)/parallel)) {
        sample_stop = nrow(metadata_t)
      } else {
        sample_stop  = i*parallel
      }

      for (ii in sample_start:sample_stop) {
        metadata_c_i <- metadata_c[which(metadata_c$patient == metadata_t$patient[ii]), ]
        system(paste0("docker run --rm -d --name ", metadata_t$sample[ii]," -v ", input_dir, ":", input_dir, " -v ", output_dir, ":", output_dir, " -v ",
                      dirname(reference), ":", dirname(reference), " broadinstitute/gatk:4.2.2.0 gatk Mutect2 -R ", reference,
                      " -I ", input_dir, metadata_t$sample[ii], "_markdup.bam ", "-I ", input_dir, metadata_c_i$sample, "_markdup.bam ",
                      "-normal ", metadata_c_i$sample, " -O ", output_dir, metadata_t$sample[ii], ".vcf"))
      }
      k = 1
      while(any(grepl(paste0(metadata_t$sample[sample_start:sample_stop], collapse = "|"), system("docker ps", intern = TRUE)[-1]))) {
        dockerps = system("docker ps", intern = TRUE)[-1]
        for (iii in sample_start:sample_stop) {
          if (any(grepl(metadata_t$sample[iii], dockerps))) {
            print(paste0("Container ", metadata_t$sample[iii], " still running"))
          }
        }
        Sys.sleep(300)
        print(paste0("Analysis running for: ", k*5, " minutes"))
        k = k+1
      }
    }
  }
}

#' Add filter information to somatic variants produced by Mutect2
#'
#' Adds information into the filter column of the vcf file produced by Mutect2
#'
#'@param metadata data frame;
#'@param input_dir character; path to the input directory
#'@param output_dir character; path to output directory
#'@param reference character; path to the fasta genome reference file. Bwa-mem2 created index must be in the same directory.

FilterMutectCalls <- function(metadata, input_dir, output_dir, reference) {
  metadata_t = metadata[metadata$status == "t", ]
  for (i in 1:nrow(metadata_t)) {
    system(paste0("docker run --rm -i --name ", metadata_t$sample[i]," -v ", input_dir, ":", input_dir, " -v ", output_dir, ":", output_dir, " -v ",
                  dirname(reference), ":", dirname(reference), " broadinstitute/gatk:4.2.2.0 gatk FilterMutectCalls -R ", reference,
                  " -V ", input_dir, metadata_t$sample[i], ".vcf", " -O ", output_dir, metadata_t$sample[i], "_filtered.vcf"))
  }
}

#' Annotate SNPs with Funcotator
#'
#' Funcotator annotates the SNPs with information stored in the annotation_dir. This directory has to be in specific format described more in detail here:
#' https://gatk.broadinstitute.org/hc/en-us/articles/4404604573851-Funcotator
#' The reference chromosome names have to be in format: chr1, chr2 ... in order for the Funcotator default database to work.
#'
#'@param metadata data frame;
#'@param input_dir character; path to the input directory
#'@param output_dir character; path to output directory
#'@param reference character; path to the fasta genome reference file. Bwa-mem2 created index must be in the same directory.
#'@param reference_version character; either "hg38" or "hg19"
#'@annotation_dir character; path to the directory with annotation data according to the Funcotator format

Funcotator <- function(metadata, input_dir, output_dir, reference, reference_version, annotation_dir) {
  metadata_t = metadata[metadata$status == "t", ]
  for (i in 1:nrow(metadata_t)) {
    system(paste0("docker run --rm -i --name ", metadata_t$sample[i]," -v ", input_dir, ":", input_dir, " -v ", output_dir, ":", output_dir, " -v ",
                  dirname(reference), ":", dirname(reference), " -v ", annotation_dir, ":", annotation_dir, " broadinstitute/gatk:4.2.2.0 gatk Funcotator -R ", reference,
                  " -V ", input_dir, metadata_t$sample[i], "_filtered.vcf", " -O ", output_dir, metadata_t$sample[i], "_filtered_an.vcf ",
                  "--output-file-format VCF ", "--data-sources-path ", annotation_dir, " --ref-version ", reference_version))
  }
}
