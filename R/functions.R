#' Perform fastqc on provided files
#'
#'@param metadata data frame;
#'@param input_dir character; path to the input directory
#'@param output_dir character; path to output directory
#'@param t numeric; number of threads to use; default=2

fastqc = function(metadata, input_dir, output_dir, t = 2)  {
  system(paste0("mkdir ", output_dir, "/fastqc/"))
  # create a vector of file names
  files = c(metadata$fastq1, metadata$fastq2)
  system(paste0("docker run --rm -i -v ", input_dir, ":", input_dir, " -v ", output_dir, ":", output_dir, " biocontainers/fastqc:v0.11.9_cv8 fastqc ", paste0(files, collapse = " "), " -t ", t, " -o ", output_dir, "/fastqc/"))
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
#' which is necessary for downstream GATK4 analysis.
#'
#'@param metadata data frame;
#'@param input_dir character; path to the input directory
#'@param output_dir character; path to output directory
#'@param reference character; path to the fasta genome reference file. Bwa-mem2 created index must be in the same directory.
#'@param t numeric; number of threads to use; default-2
bwa_mem2 = function(metadata, input_dir, output_dir, reference, t = 2) {
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
#'
MarkDuplicatesSpark = function(metadata, input_dir, output_dir, t = 2) {
  for (i in 1:nrow(metadata)) {
    system(paste0("docker run --rm -i -v ", input_dir, ":", input_dir, " -v ", output_dir, ":", output_dir,
                  " broadinstitute/gatk:4.2.2.0 gatk MarkDuplicatesSpark -I ", input_dir, metadata$sample[i],
                  ".bam", " -O ", output_dir, metadata$sample[i], "_markdup.bam", " --conf 'spark.executor.cores=", t, "'"))
  }
}

#' Detect variants with Mutect2
#'
#' Run Mutect2 with the option to run multiple samples in parallel.
#'
#'@param metadata data frame;
#'@param input_dir character; path to the input directory
#'@param output_dir character; path to output directory
#'@param reference character; path to the fasta genome reference file. Bwa-mem2 created index must be in the same directory.
#'@param t numeric; number of threads to use; default-2
#'@param parallel; how many samples should be ran in parallel, each one will need around 15 GB of RAM
#'
Mutect2 = function(metadata, input_dir, output_dir, reference, parallel) {
  metadata_t = metadata[metadata$status == "t", ]
  metadata_c = metadata[metadata$status == "c", ]

  for (i in 1:ceiling(nrow(metadata_t)/parallel)) {
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
                    " -I ", input_dir, metadata_t$sample[i], "_markdup.bam ", "-I ", input_dir, metadata_c_i$sample, "_markdup.bam ",
                    "-normal ", metadata_c_i$sample, " -O ", output_dir, metadata_t$sample[i], ".vcf"))
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
