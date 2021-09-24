# SNPsom 

### DEVELOPMENT
This package is under development.

### An R package providing a containerized GATK-based pipeline for analyzing somatic variants
The package is based on GATK4 guidelines with docker containerized tools and shiny application interface. Containerization through docker provides high level of reproducibility and easy deployment without the need of installation of each tool and its requirements. In addition, implemented deployment of multiple docker containers enables parallel processing for functions without multi-threading options - Mutect2

### Table of contents
1.[Installation](#installation)<br><TD style="FONT-SIZE:13px; COLOR:#000000; LINE-HEIGHT:20px; FONT-FAMILY:Arial,Helvetica,sans-serif">
2.[Functionality](#functionality)<br><TD style="FONT-SIZE:13px; COLOR:#000000; LINE-HEIGHT:20px; FONT-FAMILY:Arial,Helvetica,sans-serif">
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.[Input](#input)<br><TD style="FONT-SIZE:13px; COLOR:#000000; LINE-HEIGHT:20px; FONT-FAMILY:Arial,Helvetica,sans-serif">

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.[Config](#config)<br><TD style="FONT-SIZE:13px; COLOR:#000000; LINE-HEIGHT:20px; FONT-FAMILY:Arial,Helvetica,sans-serif">

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.[Metadata](#metadata)<br><TD style="FONT-SIZE:13px; COLOR:#000000; LINE-HEIGHT:20px; FONT-FAMILY:Arial,Helvetica,sans-serif">

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.[SNPsom wrapper](#SNPsom_wrapper)<br><TD style="FONT-SIZE:13px; COLOR:#000000; LINE-HEIGHT:20px; FONT-FAMILY:Arial,Helvetica,sans-serif">

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.[Shiny](#shiny)<br><TD style="FONT-SIZE:13px; COLOR:#000000; LINE-HEIGHT:20px; FONT-FAMILY:Arial,Helvetica,sans-serif">


###  1. Installation <a name="installation"></a>
Users must have [R](https://www.r-project.org/) installed. [Rstudio](https://rstudio.com/products/rstudio/) is optional but recommended IDE.

To install SNPsom package from GitHub, it is convenient to use the package devtools, namely function: install_github.
```
# install devtools
install.packages("devtools")

# install RNAseqCNV package
devtools::install_github(repo = "honzee/SNPsom")
```
The tools are deployed in separate docker containers. The images are automatically downloaded from [dockerhub](https://hub.docker.com/). Docker must be installed and run in order for the package to work. More information on docker installation here: https://docs.docker.com/get-docker/

### 2. Functionality <a name="functionality"></a>
Each SNPsom tool can be run separately or via function wrappers. trimmomatic_fastqc wrapper provides options for adapter trimming and subsequent checking by fastqc. SNPsom_wrapper deploys tool for somatic variant analysis. For more information on functions and their parameters, please refer to the function documentation.

#### 2.1. Input <a name="input"></a>
The necessary for all of the functions are config and metadata files.

#### 2.1.1 Config <a name="config"></a>
Config file defines the input directory (input_dir), output directory (out_dir), reference fasta file (reference), reference version (reference_version) and Funcotator formated annotation directory (annotation_dir). Change the file directories accordingly but keep the key words identical as the example below:

```
input_dir = "/Path/to/input_dir"
output_dir = "/Path/to/output_dir"
reference = "/Path/to/reference.fasta"
reference_version = "hg19orhg38"
annotation_dir = "/Path/to/funcotator_formatted_dir"
```

#### 2.1.1 Metadata <a name="metadata"></a>
Metadata file contains information about the analyzed samples. The required columns are: sample - sample name; status - "t" or "c" for tumour of control sample, patient - patient id, fastq1 - file name for the forward fastq file, fastq2 file name for the reverse fastq file

| sample       | patient | status | fastq1                | fastq2                | 
|--------------|---------|--------|-----------------------|-----------------------| 
| PE01_tumour  | PE01    | t      | PE01_tumour_R1.fastq  | PE01_tumour_R2.fastq  | 
| PE01_control | PE01    | c      | PE01_control_R1.fastq | PE01_control_R2.fastq | 
| PE02_tumour  | PE02    | t      | PE02_tumour_1.fastq   | PE02_tumour_2.fastq   | 
| PE02_control | PE02    | c      | PE02_control_R1.fastq | PE02_control_2.fastq  | 


#### 2.1. SNPsom wrapper <a name="SNPsom_wrapper"></a>
SNPsom_wrapper allows deployment of alignment (bwa-mem2), marking duplicates (MarkDuplicatesSpark), variant calling (Mutect2), variant filtering (FilterMutectCalls) and variant annotation (Funcotator) in single command.

```
SNPsom_wrapper(metadata = metadata, input_dir = input_dir, output_dir = output_dir, reference = reference, t = t, parallel = parallel, reference_version = "hg19/hg38", annotation_dir = annotation_dir)
```

#### 2.3. Shiny <a name="shiny"></a>
Shiny app allows for deployment of SNPsom_wrapper by providing the config and metadata files through interactive interface.

```
launchApp()
```


