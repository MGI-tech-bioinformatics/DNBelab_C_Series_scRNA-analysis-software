# DNBelab C Series Single-Cell RNA Analysis Software

# Introduction
- **Propose**
  - An open source and flexible pipeline to analyze DNBelab C Series<sup>TM</sup> single-cell RNA datasets. 
- **Language**
  - Workflow Description Language (WDL), Python3 and R scripts.
- **Hardware/Software requirements** 
  - x86-64 compatible processors
  - require at least 36GB of RAM and 10 CPU. 
  - 64bit Linux
- **Workflow**
![](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software/blob/master/doc/fig/workflow.jpg)
# Directory contents
- **bin**        pre-compiled executables for Linux
- **config**     read structure configure files
- **pipelines**  WDL pipeline
- **scripts**    miscellaneous scripts

# With Docker
```
Setup
1. Install docker follow the official website
    https://www.docker.com/
2. Then do the following for the workflow:
    docker pull huangshunkai/dnbelab_c4:latest
    
    Notes:
    1. Please make sure that you run the docker container with at least 36GB memory and 10 CPU.
    2. The input is sample list and output directory which descripted below (Main progarm arguments).
    
Prepare

   cat config.json
   
   {
    "main.fastq1": "/DNBelab_C4/rawfq/demo_1.fq.gz",                 
    "main.fastq2": "/DNBelab_C4/rawfq/demo_2.fq.gz",                 
    "main.ID": "Demo_single",                                        
    "main.forceCell": "0",                                          
    "main.umilow": "1000",                                           
    "main.species":"GRCh38",                                    
    "main.original":"cell lines",                            
    "main.SampleTime":"2020-06-25",                         
    "main.ExperimentalTime":"2020-06-25"                    
   }
    
Running
1. Please set the following variables on your machine:
   (a) $DB_LOCAL: directory on your local machine that has the database files. Make sure that the directory must contains two subdirectories, "gtf" and "star_index". The gene annotation file named "genes.gtf" must be included under "gtf"; the genome index file for STAR under the "star_index". If you build the database youself, make sure the format of the directory path is correct.
   (b) $DATA_LOCAL: directory on your local machine that has the sequence data and "config.json" file.
      "config.json" must follow the format descripted bellow,
      and the *PATH* in "config.json" must be absolute dicrtory of $DATA_LOCAL.
   (c) $RESULT_LOCAL: directory for result.
  
2. Run the command:
   10x sequence data:
   docker run -d -P \
   --name $scRNANAME \
   -v $DB_LOCAL:/DNBelab_C4/database \
   -v $DATA_LOCAL:/DNBelab_C4/rawfq \
   -v $RESULT_LOCAL:/DNBelab_C4/result \
   huangshunkai/dnbelab_c4:latest \
   /bin/bash \
   /DNBelab_C4/bin/10xRun.sh
  
   mgi sequence data:
   docker run -d -P \
   --name $scRNANAME \
   -v $DB_LOCAL:/DNBelab_C4/database \
   -v $DATA_LOCAL:/DNBelab_C4/rawfq \
   -v $RESULT_LOCAL:/DNBelab_C4/result \
   huangshunkai/dnbelab_c4:latest \
   /bin/bash \
   /DNBelab_C4/bin/mgiRun.sh
  
3. After satisfactory result was generated:
   docker rm $scRNANAME
```

# Without Docker but run in local server
```
$ git clone https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software.git
```

# Third-party software 

## Users need to install manually
* [java](https://www.oracle.com/java/)
* [Cromwell](https://github.com/broadinstitute/cromwell/releases)
* [R](https://www.r-project.org/)  (3.5+) #  with following R packages installed
  * ggplot2
  * getopt
  * data.table
  * cowplot
  * DropletUtils(1.6.1+)
* [python3](https://www.python.org/downloads/)  (3.6+) #  with following Python3 packages installed
  * numpy
  * pandas
  * python-igraph
  * louvain
  * scanpy(1.4.3+)
  * jinja2(2.10.3+)

  
## Pre-compiled executables within binary releases
* [PISA](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software/blob/master/bin/PISA)
* [sambamba](https://lomereiter.github.io/sambamba/)
* [STAR](https://github.com/alexdobin/STAR)


# Database

## Download ready-made datasets
We provide the following databases for [download](http://ftp.cngb.org/pub/CNSA/data2/CNP0000906/), including fasta, gtf, and STAR(V2.7.3a) index files.
- **human(GRCh38)**
- **mouse(GRCm38)** 
- **Mixed Database(GRCh38 & GRCm38)** 

Note: Mixed dual species databases only for double species sample analysis.

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

- Step 1: Prepare fastq
We provide 100MB [Demo](https://bgitech-my.sharepoint.com/:f:/g/personal/huangshunkai_genomics_cn/EhbcIp5iKe1Ksq1Pf-7v0yIBkWNJvH08JX89rdlPOU1KVQ?e=dlHa4F) sequencing data for testing.
We also provide 36GB PBMCs sample fastq by pairs sequencing for download [fastq](http://ftp.cngb.org/pub/CNSA/data2/CNP0000906/CNS0232089/CNX0190688/CNR0248280/).

- Step 2: Setup configure file.
```
# Goto test directory
cd ./example/single_Species

# Check configure file
$ cat config.json
  {
    "main.fastq1": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/example/single_Species/fastq/Demo.human.fq.1.gz",
    "main.fastq2": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/example/single_Species/fastq/Demo.human.fq.1.gz",
    "main.root": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software",
    "main.gtf": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/databases/GRCh38/gtf/genes.gtf",
    "main.ID": "demo",
    "main.outdir": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/example/single_Species/result",
    "main.config": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/config/DNBelabC4_scRNA_readStructure.json",
    "main.Rscript":"/User/Pub/third_party/Rscript",
    "main.refdir": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/databases/GRCh38/star_index",
    "main.Python3": "/User/Pub/third_party/python3",
    "main.species":"human",
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
We provide 135MB [Demo](https://bgitech-my.sharepoint.com/:f:/g/personal/huangshunkai_genomics_cn/EneE2O8UOZxCj5S7FTSfdQ8BGhtk2KxuSbqClHJ7l-nKLw?e=mOJjiF) sequencing data for testing.
We also provide 52GB Mixed Sample(GRCh38 & mm10) pairs sequencing fastq for download[fastq](http://ftp.cngb.org/pub/CNSA/data2/CNP0000906/CNS0196715/CNX0144387/CNR0177382/)

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
    "main.config": "/User/pipeline/DNBelab_C_Series_scRNA-analysis-software/config/DNBelabC4_scRNA_readStructure.json",
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
[MIT](LICENSE)

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

   Yes,  the result/symbol directory records the symbol for each step, you can then keep the output path unchanged and run this pipeline after correct the error. If you use docker images, run the command `docker start $scRNANAME && docker exec -d $scRNANAME /bin/bash /DNBelab_C4/bin/10xRun.sh` or `docker start $scRNANAME && docker exec -d $scRNANAME /bin/bash /DNBelab_C4/bin/mgiRun.sh`.

6. Why the inflection point is inaccurate on the total count curve?

   You can specify "main.umilow" in the configure file like "main.umilow": "1000". "main.umilow" is a numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets, default 50.


