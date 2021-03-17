# Python functions to be used for the snakemake pipeline.

# ----- Function get_lanes() ----- #
# It will read the variable "units" which has the names of the sample, the lane and the paths to the fastq files.
# It will make a list with the paths to fq1 and fq2 for each sample and each lane. 
def get_lanes(wildcards):
    return units.loc[(wildcards.sample, wildcards.lane), ["fq1", "fq2"]].dropna()

# ----- Function is_single_end() ----- #
# It will say if the sample we are analyzing is single end.
def is_single_end(sample):
    return pd.isnull(units.loc[(sample), "fq2"][0])

# ----- Function get_fq() ----- #
# Get raw or trimmed reads based on trimming configuration.
# It has to be used as {input} of the rule align.
# If the function "is_single_end()" is not true, it will be a paired-end sample and get_fq() will return a list where there are two strings: 
#   the first string will be the forward read and the second string will be the reverse read. 
# In the case of a single-end sample, the list will contain only one string.
# To make this work, the {params.pe} has to change to {params.reads}, which will be the output of the function "set_reads()": reads = set_reads.
# When using this, the arguments "-1 {input.fw} -2 {input.rv}" in the shell command have to be deleted, because they are included in the function set_reads().
# When using this, the {input.fw} and {input.rv} have to be substituted by "get_fq"
def get_fq(wildcards):
    if config["options"]["trimming"]:
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("{tmp}/fastq/trimmed/{sample}.{group}.fastq.gz", group=[1, 2], **wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/trimmed/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)
    else:
        # no trimming, use raw reads
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("{tmp}/fastq/{sample}.{group}.fastq.gz", group=[1, 2], **wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)

# ----- Function set_reads() ----- #
# To use as one of the params in align.
# It will return a list where there are the names of the reads to aling that will serve as input.
# It will look at how long is the input of the rule.
# If it is paired-end, it will return the two reads, if it's singl-end it will return only one.
# The {params.pe} in the {params} of the rule align will have to change to "reads = set_reads", but now the pipeline is to do paired-end only. 
def set_reads(wildcards, input):
        n = len(input)
        if n == 1:
            reads = "{}".format(*input)
            return reads
        else:
            reads = config["params"]["bowtie2"]["pe"] + " -1 {} -2 {}".format(*input)
            return reads



### FUNCTIONS FOR RULE BAM2BW NO SUBTRACT ###
### ===================================== ###
# ----- Function get_bam() ----- #
# To get bams for bam2bw
def get_bam(wildcards):
    return { "case": "results/02_bam/{sample}.bam".format(sample=wildcards.sample) }


# ----- Function set_reads_extension() ----- #
def set_read_extension(wildcards):
    if is_single_end(wildcards.sample):
        return "--extendReads " + str(config['bam2bigwig']['read_extension'])
    return "--extendReads"

# ----- Function set_reads_spike2() ----- #
def set_reads_spike2(wildcards, input):
        n = len(input)
        assert n == 1 or n == 2, "input->sample must have 1 (sample) or 2 (sample + spike) elements"
        if n == 1:
            reads = "scripts/bam2bigwig_noSubtract.py"
            return reads
        if n == 2:
            reads = "scripts/bam2bigwig_spike_noSubtract.py --spike {}".format(input.spike)
            return reads  