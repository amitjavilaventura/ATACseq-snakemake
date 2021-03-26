# Pasini's lab ATAC-seq pipeline for Snakemake

**ONLY PAIRED-END**. To use an ATAC-seq pipeline compatible with single-end, go to the [se_and_pe branch](https://github.com/amitjavilaventura/ATACseq-snakemake/tree/se_and_pe) (*not fully updated*).

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.9.1-brightgreen.svg)](https://snakemake.bitbucket.io)


Developers: [amitjavilaventura](https://github.com/amitjavilaventura) & [dfernandezperez](https://github.com/dfernandezperez)

Snakemake-based ATAC-seq pipeline to be run in our PBS-based HPC using singularity containers. The singularity image that is used to run this pipeline is based in [this](https://github.com/amitjavilaventura/Dockerfiles/blob/main/atacseq_snakemake/Dockerfile) docker container.

Many features from this pipeline have been retrieved from or are based in [*deferernandezperez*'s ChIPseq pipeline for Snakemake](https://github.com/dfernandezperez/ChIPseq-snakemake/tree/development), so if you have some doubts that are not solved in this readme file or in the comments within the code, you may find the solutions in th

## Download

To download this repository run in the terminal:

`git clone https://github.com/amitjavilaventura/ATACseq-snakemake.git *empty_folder*`

It is recommended that you download the [docker container](https://github.com/amitjavilaventura/Dockerfiles/blob/main/atacseq_snakemake/Dockerfile) as a singularity image using `singularity pull`.


## Setup

The following files are located inside the folder `configuration/`:

* Raw data paths (`units.tsv`).
* Sample metadata (`samples.tsv`).
* Cluster configuration (`cluster.yaml`).
* Pipeline parameters -alignment, peak calling...- (`config.yaml`).

### Raw data

Paths to raw data are located in the file `units.tsv`. The file has the following structure:

| sample | lane | fq1 | fq2 |
|:------:|:----:|:---:|:----:|
| name_of_sample | name_of_lane_or_resequencing | path/to/forward.fastq | path/to/reverse.fastq |

* The first field correspond to the sample name. This field has to be the same as the sample name that is specified in the `samples.tsv` file (see below). It is recommended to NOT use underscores in the name of the samples, dashed are prefered. I still don't understand why sometimes I get errors if I use them so before fixing I strongly recommend to use dashes instead.

* The second field corresponds to `lane`. The idea of this field is to group fastq files corresponding to the same sample (or to samples that have to be merged). For example, if 1 sample arrived in 2 different lanes from a PE experiment, in total there will be 4 fastqs (2 forward and 2 reverse). In this case, one should enter the same sample 2 times, putting in the `lane` field the corresponding lanes (lane1 and lane2, for example). Actually one can write any word in this field, the idea is to group fastqs from the same sample. All the entries with the same name in the `sample` field with different `lane` will be merged in the same fastq. Here an example of how it would be with 1 sample that arrived in 2 lanes:

| sample | lane | fq1 | fq2 |
|:------:|:----:|:---:|:----:|
| foo | lane1 | path/to/forward_lane1.fastq | path/to/reverse_lane1.fastq |
| foo | lane2 | path/to/forward_lane2.fastq | path/to/reverse_lane2.fastq |

Here I am using lane1 and lane2 for consistency and making things more clear, but the following would also work:

| sample | lane | fq1 | fq2 |
|:------:|:----:|:---:|:----:|
| foo | potato | path/to/forward_lane1.fastq | path/to/reverse_lane1.fastq |
| foo | checazzo | path/to/forward_lane2.fastq | path/to/reverse_lane2.fastq |

* Finally the last 2 fields `fq1` and `fq2` correspond to the paths to the fastq files. `fq1` is the FORWARD read and  `fq2` the REVERSE. The order is very important because they will be sent in that order to the aligner. Both `fq1` and `fq2` must be provided, because this pipeline is intended **for PAIRED-END only**. If you want a pipeline compatible with single end, you can go to the [*se_and_pe branch*](https://github.com/amitjavilaventura/ATACseq-snakemake/tree/se_and_pe), although it may not be fully updated.


### Sample metadata

All metadata and information regarding every sample is located in `samples.tsv`. The file has the following structure:

| NAME | INPUT |  USER | GENOME | RUN | IS_INPUT | CONDITION | DOWNSAMPLE_GROUP |
|:----:|:-----:|:-----:|:--:    |:---:|:------:  |:---:      |:--------:        |
| name_of_sample | input_to_use | user | version of genome (i.e: mm10) | run of the sequencing | if the sample is an input | sample_condition | group whose bamfiles will be downsampled to the bamfile with the lower number of reads |

For example:

| NAME | INPUT |  USER | GENOME | RUN    | IS_INPUT | CONDITION | DOWNSAMPLE_GROUP |
|:----:|:-----:|:-----:|:--:    |:---:   |:------:  |:---------:|:----------------:|
| foo1 | False | user  | mm10   | 210303 | FALSE	| foo       | *not used*       |
| foo2 | False | user  | mm10   | 210303 | FALSE	| foo       | *not used*       |
| bar1 | False | user  | mm10   | 210303 | FALSE	| bar       | *not used*       |
| bar2 | False | user  | mm10   | 210303 | FALSE	| bar       | *not used*       |


* For every sample, the `NAME` field has to contain exactly the same name that was written in the `sample` column of the `units.tsv`.

* The `INPUT` field contains the name of the input corresponding to the given sample. It has to be the name of the input written in the fields `sample` and `NAME` from `units.tsv` and `samples.tsv`. At the moment this is set to FALSE because there is not an option to put the input.

* `GENOME`: Version of the genome used for the alignment. It will be used for peak annotation with ChIPseeker. Right now the accepted values are mm9, mm10, hg19 and hg38.

* `IS_INPUT`: The options are TRUE or FALSE. If the sample is an input set it to TRUE. Also, in case the sample is an input sequenced just to calculate the ratio sample/spike-in that won't be used to call peaks, set it to TRUE. At the moment this is set to FALSE because there is not an option to put the input.

* `CONDITION`: Condition in which the sample belongs. This feature has been added because there is an option to merge replicates inside the Genrich peak calling step.

* `DOWNSAMPLE_GROUP`: Group of different samples whose bamfiles will be downsampled to the file with lower number of reads. This function has not been added yet. 


### Configuration of pipeline parameters

To configure the pipeline parameters go to `configuration/config.yaml`. 
This file contains the configuration of the software and parameters used in the pipeline. Modify them as you wish. 
Check always that you are using the correct genome files corresponding to the version that you want to use.

Structure of the `configuration/config.yaml`:

* `units` has the path to the `configuration/units.tsv`. It must not change. 
* `samples` has the path to the `configuration/samples.tsv`. It must not change. 
* `cluster` has the path to the `configuration/cluster.yaml`. It must not change. 
* `tmp` contains the path to a folder where fastq files will be saved. This folder can be any folder you want.
* `params` contains the parameters of some software used:
	
	+ `bowtie2`
	+ `samblaster`
	+ `samtools`
	+ `fastp`
	+ `genrich` is the peak caller used in this pipeline, look at [it's repository](https://github.com/jsh58/Genrich) for more information. Parameters in Genrich:

		+ `atacmode`: whether to use ATAC mode (-j, default) or "normal" mode. One of "-j" or "".
		+ `pe_or_se`: wheter samples are PE or SE. This is required only when `options: genrich_merge: True`, otherwise is not used.
		+ `se`: parameters for peak calling from single-end reads. The default is "-y -d 150". "-y" is compulsory.
		+ `p_or_q`: whether to filter by p-value or q-value. One of "p" (default) or "q". 
		+ `pqval`: p- or q-value threshold in the peak calling step. Default is 5e-4
		+ `filt_peaks_pqval`: p- or q- value threshold to filter the peaks in a second filtering step. In log10 format. Default is "5"
		+ `chrM`: whether to remove reads aligning to chrM (or others) in the peak calling step. One of "-e chrM" (default) or "".
		+ `rm_pcr_dups`: whether to remove the pcr duplicates. One of "-r" or "". The default is "" because PCR duplicates are removed by samblaster.
	
	+ `bam2bigwig`

* `promoter`: how many bases after and before the TSS define the promoter.
* `ref`: reference files (genome indexes...) used.
* `options`: options for the pipeline. The ones available are:

	+ `genrich_pool`: if you are working with replicates, setting this to True will make Genrich to call the peaks using the bam files of both replicates. The output will be narrowPeak file for the condition specified in `samples.tsv`. If False, a .narrowPeak file for each replicate will be generated.
	+ `rm_chrM`: remove the reads mapping to the mitochondrial chromosome in the BAM files. *Not available yet*



### Cluster configuration

`configuration/cluster.yaml` contains the per rule cluster parameters (ncpus, ram, walltime...). It can be modified as desired. In the future I want to remove this file in favour of the new [snakemake profiles](https://github.com/Snakemake-Profiles) system (see below), but I still need to understand a little bit better how it works and how to properly do the migration.


### Snakemake profiles

In Snakemake 4.1 [snakemake profiles](https://github.com/Snakemake-Profiles) were introduced. They are supposed to substitute the classic cluster.json file and make the execution of snakemake more simple. 
The parameters that will be passed to the `snakemake`command (i.e: `--cluster`, `--use-singularity`...) now are inside a yaml file (`config.yaml`) inside the profile folder (in the case of this repository is `snakemake_profile/`). 
The `snakemake_profile/config.yaml` file contains the parameters passed to snakemake. So if you were executing snakemake as `snakemake --cluster qsub --use-singularity` the new `config.yaml` would be like this:

```yaml
cluster: qsub
use-singularity: true
```

## Execution of the pipeline

Once you have all the configuration files as desired, it's time to execute the pipeline. For that you have to execute the `execute_pipeline.sh` script, followed by the name of the rule that you want to execute. If any rule is given it will automatically execute the rule `all` (which would execute the standard pipeline). Examples:

```bash
./execute_pipeline.sh all
```

is equivalent to 

```bash
./execute_pipeline.sh
```

At the end of the `configuration/Snakefile` you will find all the possible target rules and their corresponding output files:

* `all`: will do the alignment (BAMs), peak calling and summits (NarrowPeak and BEDs), peak annotation and quality control.
* `get_peaks`: will generate the NarrowPeak files and the corresponding p-value filtered BEDs.
* `get_peakanno`: will generate the annotated peaks from the p-value filtered BEDs.
* `get_summits`: will generate the summits of the peaks and the filtered peaks.
* `get_bams`: will generate the BAM files
* `get_bw`: will generate the BIGWIG files from the BAMs.
* `get_multiqc`: will do all the QC analyses and will generate a MultiQC report.
* `get_fastqc`: will do a fastQC analysis.
* `get_fingerprint`: will perform a QC with `plotFingerPrint` from `deepTools`.


## To Do's

* Add an error message if forward and reverse fq files are not found
* Migrate 100% to snakemake profiles and stop using the `cluster.yaml` configuration.
* Add options to downsample the *bamfiles* by desired groups.
* Add options to remove chrM reads from bamfiles.
* Add options to split bam files and call peaks from that files.
* Adding the control/input option in Genrich in case it is provided.
