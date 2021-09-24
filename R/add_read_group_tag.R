add_read_group_tag <- function(metadata, input_dir, output_dir) {
  #for (i in 1:nrow(metadata)) {
  for (i in 1:nrow(metadata)) {
    system(paste0("docker run --rm -i -v ", input_dir, ":", input_dir, " -v ", output_dir, ":", output_dir,
                  " broadinstitute/picard:2.26.2 ", " java -jar /usr/picard/picard.jar AddOrReplaceReadGroups I=",
                  input_dir, metadata$sample[i], ".bam ", "O=", output_dir, metadata$sample[i], "_tagged.bam ",
                  "-RGID ", metadata$sample[i], " ",
                  "-RGPL illumina -RGPU unit1 -RGLB lib1 -RGSM ", metadata$sample[i]))
    system(paste0("rm ", input_dir, metadata$sample[i], ".bam"))
  }
}
