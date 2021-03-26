# ----- Aligning with BOWTIE2 and processing with SAMTOOLS ----- #
# Aligning with bowtie
# Reomving PCR duplicates with samblaster
# Sorting and indexing with samtools

rule align:
    input:
        get_fq,
    output:
        bam   = "results/02_bam/{sample}.bam",
        index = "results/02_bam/{sample}.bam.bai"
    threads:
        CLUSTER["align"]["cpu"]
    params:
        index  	     = config["ref"]["index"],
        bowtie2 	 = config["params"]["bowtie2"]["global"],
        samblaster   = config["params"]["samblaster"],
        reads  	     = set_reads,
        rm_chrM      = config["options"]["rm_chrM"],
        samtools_mem = config["params"]["samtools"]["memory"]
    message:
        "Aligning {input} with parameters {params.bowtie2}"
    log:
       align   = "results/00_log/align/{sample}.log",
       rm_dups = "results/00_log/align/rm_dup/{sample}.log",
    benchmark:
        "results/.benchmarks/{sample}.align.benchmark.txt"
    shell:
        """
        bowtie2 -p {threads} {params.bowtie2} -x {params.index} {params.reads} 2> {log.align} \
        | samblaster {params.samblaster} 2> {log.rm_dups} \
        | samtools view -q2 -Sb -F 4 - \
        | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
        samtools index {output.bam}
        """


## add an option to remove chromosome M reads

## add a rule for bam splitting.
