# Input JSON

An input JSON file includes all genomic data files, parameters and metadata for running pipelines. Our pipeline will use default values if they are not defined in an input JSON file.

We provide two example JSON files including [single-species](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software/tree/master/example/single_Species/config.json) and [mixed-species](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software/tree/master/example/double_Species/config.json).


## Pipeline parameters

| Parameter         | Type 		| Description                                                  |
| ------------------| --------- | ------------------------------------------------------------ |
| `main.ID` 		| String  	| MANDATORY. Sample id.                 |
| `main.fastq1` 	| Fastq   	| MANDATORY. Read 1 in fastq format. Can be gzipped. Fastqs from different lanes can be seperated with comma. For example, "L01_read_1.fq.gz,L02_read_1.fq.gz,..." |
| `main.fastq2` 	| Fastq   	| MANDATORY. Read 2 in fastq format. Can be gzipped. Fastqs from different lanes can be seperated with comma. For example, "L01_read_2.fq.gz,L02_read_2.fq.gz,..."|
| `main.outdir`	 	| Directory | MANDATORY. Output directory               |
| `main.root`   	| Directory | MANDATORY. Directory of this pipeline.               |
| `main.expectCell` | Integer 	| Optional, default: 1000. Expected cell number.               |
| `main.umilow` | Integer 	| Optional, default: 1000. Expected UMIs count per cell threshold.               |
| `main.config` 	| JSON file | MANDATORY. Read structure configure. UMI and Cell Barcode specified in this file. Can be found at `config/` directory.              |
| `main.Rscript` 	| Path 		| MANDATORY. Path to Rscript.               |
| `main.Python3` 	| Path 		| MANDATORY. Path to Python3.               |
| `main.refdir` 	| Directory | MANDATORY. STAR index directory of genome reference.               |
| `main.gtf` 		| File Path | MANDATORY. gtf file of genome reference.               |
| `main.chrom`		| File Path | Optional, default: `config/species_binding.txt` directory. Chromosome and species correspondence files in the mixed species process. And Single species process don't need this file.|
| `main.species` 	| String	| Optional, default: Null. Species.               |
| `main.original` 	| String 	| Optional, default: Null. original.               |
| `main.SampleTime` | String 	| Optional, default: Null. Sample time.               |
| `main.ExperimentalTime` | String | Optional, default: Null. Experimental time.               |


