# DNBelab<sup>TM</sup> C Series Single-Cell RNA Analysis Software

# Introduction
- **Propose**
  - An open source and flexible pipeline to analyze DNBelab C Series<sup>TM</sup> single-cell RNA datasets. 
- **Language**
  - Workflow Description Language (WDL), Python3 and R scripts.
- **Hardware/Software requirements** 
  - x86-64 compatible processors
  - require at least 16GB of RAM, ideally 32GB. 
  - 64bit Linux
- **Workflow**
![](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software/blob/master/doc/fig/workflow.jpg)
# Directory contents
- **bin**        pre-compiled executables for Linux
- **config**     read structure configure files
- **pipelines**  WDL pipeline
- **scripts**    miscellaneous scripts

# Install

## Downlod the lastest binary [release](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software.git) and uncompress it
```
$ wget -c https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software.git
$ tar zxvf scRNA-pipe-demo_v0.1.0_Linux_static.tar.gz
```
or you can download with `git` 
```
$ git clone https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software.git
```

# Third-party software 

## Users need to install manually
* [java](www.java.com)
* [Cromwell](https://github.com/broadinstitute/cromwell/releases/download/35/cromwell-35.jar)
* [R](https://www.r-project.org/)  (3.5+) #  with following R packages installed
  * ggplot2
  * getopt
  * data.table
  * cowplot
* [python3](https://www.python.org/downloads/)  (3.6+) #  with following Python3 packages installed
  * numpy
  * pandas
  * scanpy(1.4.3+)
  * jinja2(2.10.3+)
* Tipe: You can use Conda to manage third-party software and create related system environment
* [Conda](https://anaconda.org/anaconda/conda) # with following anaconda installed
- Step 1: Create your environment name
```
conda create -n <env_name>
```
- Step 2: Activate your environment name and download third-party software
```
# Activate environment
source activate <env_name>

# Download R
conda install -c http://mirror/anaconda/pkgs/r/ r

# Download R packages
conda install -c https://repo.anaconda.com/pkgs/r/ ggplot2
conda install -c https://repo.anaconda.com/pkgs/r/ getopt
conda install -c https://repo.anaconda.com/pkgs/r/ data.table
conda install -c https://repo.anaconda.com/pkgs/r/ devtools

# Download cowplot package by R
R
library('devtools')
devtools::install_github('wilkelab/cowplot')

# Download Python3 packages by Conda
conda install numpy
conda install pandas
conda install scanpy
conda install jinja2
```
- Step 3: Load the Conda library instead of system library
```
echo 'export LD_LIBRARY_PATH="Your_Anaconda3_Path/envs/python37/lib":$LD_LIBRARY_PATH' >> ~/.bash_profile
source ~/.bash_profile
```
  
## Pre-compiled executables within binary releases
* [PISA](https://github.com/shiquan/LISA)
* [sambamba](https://lomereiter.github.io/sambamba/)
* [STAR](https://github.com/alexdobin/STAR)


# Database

## Download ready-made datasets
We provide human(GRCh38), mouse(mm10) and mixed dual species databases for download, including fasta, gtf, and STAR(V2.7.3a) index files.
- **Human(GRCh38)** [GRCh38](http://ftp.cngb.org/pub/CNSA/CNP0000906/)
- **Mouse(mm10)** [mm10](http://ftp.cngb.org/pub/CNSA/CNP0000906/)
- **Mixed Database(GRCh38 & mm10)** [GRCh38 & mm10](http://ftp.cngb.org/pub/CNSA/CNP0000906/)

## or you can build the database youself
Firstly, you need to prepare the fasta and gtf files of the reference database. And then build STAR index files. Please refer to the following command lines.
```
### Goto pipeline directory
$ ls
bin/  config/  doc/  example/  LICENSE  pipelines/  README.md  scripts/

### Create fold and prepare the reference files
$ cd example/database && mkdir star_index
$ gzip -d gtf/genes.gtf.gz

### build STAR index 
$ ../../bin/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles ./fasta/genome.fa --sjdbGTFfile ./gtf/genes.gtf
Dec 29 20:03:24 ..... started STAR run
*** logs ignored
Dec 29 20:04:33 ..... finished successfully
```
For a 3G reference sequence file, bulid index takes about 1 hour. It is worth noting that the STAR version corresponds to the STAR index, we default use the V2.7.3 STAR.

# Input JSON file
An input JSON file includes all input parameters and genome reference index directory for running pipelines. Always use absolute paths in an input JSON.

[Input JSON file specification](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software/blob/master/doc/input.md)

# Usage

## Example - single_Species

- Step 0: Build reference index

Please refer to [Database](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software#Database)

- Step 1: Prepare fastq
We provide mouse(mm10) pairs sequencing fastq for download[fastq](http://ftp.cngb.org/pub/CNSA/CNP0000906/CNS0196716/CNX0144388/CNR0177383/)

- Step 2: Setup configure file.
```
# Goto test directory
cd ./example/single_Species

# Check configure file
$ cat config.json
  {
    "main.fastq1": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/example/single_Species/read_1.fq.gz",
    "main.fastq2": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/example/single_Species/read_2.fq.gz",
    "main.root": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software",
    "main.gtf": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/databases/mm10/gtf/genes.gtf",
    "main.ID": "demo",
    "main.outdir": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/example/single_Species/result",
    "main.config": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/config/BGI_Droplet_scRNA_readstruct_v2.json",
    "main.Rscript":"/User/Pub/third_party/Rscript",
    "main.refdir": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/databases/mm10/star_index",
    "main.Python3": "/User/Pub/third_party/python3",
    "main.species":"mm10",
    "main.original":"cell line",
    "main.SampleTime":"2019-12-25",
    "main.ExperimentalTime":"2019-12-25"
  }
  
```

- Step 3:  Run this pipeline.

```
java -jar cromwell-35.jar run -i config.json ../../pipelines/Droplet_single.wdl
```

- Step 4: Check results.
```
# After all analysis processes ending, you will get these files below:
$ cd result && ls
outs/  report/  temp/  /symbol workflowtime.log

$ ls out
cell_barcodes.txt  cluster.h5ad  count_mtx.tsv.gz  final.bam

$ ls report
alignment_report.csv  annotated_report.csv  cell_report.csv  cluster.csv  cutoff.csv  iDrop_demo.html  marker.csv  RNA_counts.pdf  sample.csv  sequencing_report.csv  vln.csv

$ ls symbol
# In single_Species result,there are some follow files will be generated:

makedir_sigh.txt parseFastq_sigh.txt fastq2bam_sigh.txt sortBam_sigh.txt cellCount_sigh.txt cellCalling_sigh.txt countMatrix_sigh.txt report_sigh.txt
```

So the final html report is at `outdir Path`/report/iDrop_*.html

## Example - double_Species

- Step 0: Build reference index

Please refer to [Database](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software#Database)

- Step 1: Prepare fastq
We provide Mixed Sample(GRCh38 & mm10) pairs sequencing fastq for download[fastq](http://ftp.cngb.org/pub/CNSA/CNP0000906/CNS0196715/CNX0144387/CNR0177382/)

- Step 2: Setup configure file.
```
# Goto test directory
cd ./example/double_Species

# Check configure file
$ cat config.json
  {
    "main.fastq1": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/example/double_Species/read_1.fq.gz",
    "main.fastq2": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/example/double_Species/read_2.fq.gz",
    "main.root": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software",
    "main.gtf": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/databases/GRCh38_mm10/gtf/genes.gtf",
    "main.ID": "demo",
    "main.chrom":"/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/config/species_binding.txt",
    "main.outdir": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/example/double_Species/result",
    "main.config": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/config/BGI_Droplet_scRNA_readstruct_v2.json",
    "main.Rscript":"/User/Pub/third_party/Rscript",
    "main.refdir": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/databases/GRCh38_mm10/star_index",
    "main.Python3": "/User/Pub/third_party/python3",
    "main.species":"GRCh38_mm10",
    "main.original":"cell line",
    "main.SampleTime":"2019-12-25",
    "main.ExperimentalTime":"2019-12-25"
  }
```

- Step 3:  Run this pipeline.

```
java -jar cromwell-35.jar run -i config.json ../../pipelines/Droplet_double.wdl
```

- Step 4: Check results.
```
# After all analysis processes ending, you will get these files below:
$ cd result && ls
outs/  report/  temp/  /symbol workflowtime.log

$ ls out
anno_species.bam

$ ls symbol
# In double_Species result, there are some follow files willbe generated:

makedir_sigh.txt parseFastq_sigh.txt fastq2bam_sigh.txt sortBam_sigh.txt cellCount_sigh.txt cellCalling_sigh.txt report_sigh.txt

$ ls report
alignment_report.csv  cell_barcodes.txt cell_report.csv  mix_report.csv sample.csv iDrop_Demo.html annotated_report.csv  cell_count_summary.png  cutoff.csv mixture_cells.png  sequencing_report.csv  vln.csv
```
So the final html report is at `outdir Path`/report/iDrop_*.html

# License
[MIT](LICENSE) © Richard Littauer

# Frequently Asked Questions
1. Does this pipeline correct UMI errors?

   Yes, UMIs from same cell in same gene will be corrected using Hamming distance and frequency. This method retains only the UMI with the highest counts.

2. Does this pipeline correct Cell Barcode errors?

   Yes. Cell barcode not contained in white list will be corrected using levenshtein distance. Users can edit this distance parameter at the read structure configure file. If distance set to 0, pipeline will skip the correction.
   
3. Can you introduce how many QC steps performed at this pipeline? And for each step what parameters are used?

   For fastqs, this pipeline filter reads with low quality (< Q20) or if 2 bases < Q10 at first 15 bases. For bam, read with low mapping quality (< 20) will be filtered. For UMIs, check question 1.
   
4. Can I use this pipeline to analysis 10X Genomes single cell gene expression data ?

   Yes,  the read structure configure file can be found at `config/10X_3end_readstruct.json`. Other steps are same with the demonstration.
   
5. Can I continuse to run the workflow if some errors were happended in the process?
   Yes,  the result/symbol directory records the symbol for each step, you can delete the lastest symbol.txt file then keep the output path unchanged and run this pipeline after correct the error.

6. Why the inflection point is inaccurate on the total count curve?
   You can specify "main.umilow" in the configure file like "main.umilow": "1000". "main.umilow" is a numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets, default 1000.


