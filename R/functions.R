fastqc = function(files, input_dir, output_dir, t)  {
system(paste0("docker run --rm -d -v ", input_dir, ":", input_dir, " -v ", output_dir, ":", output_dir, " biocontainers/fastqc:v0.11.9_cv8 fastqc ", paste0(input_dir, files, collapse = " "), " -t ", t, " -o ", output_dir))
}

system(paste0("docker run --rm -d -v ", input_dir, ":", input_dir, " biocontainers/fastqc:v0.11.9_cv8", " echo 'testtesxt' > ", input_dir, "/test.txt"))
