# ----- Aligning with BOWTIE2 and processing with SAMTOOLS ----- #
# Aligning with bowtie
# Reomving PCR duplicates with samblaster
# Sorting and indexing with samtool

rule align:
    input:
        get_fq
    output:
         bam   = temp("results/02_bam/{sample}.bam"),
         index = temp("results/02_bam/{sample}.bam.bai")
    threads:
        CLUSTER["align"]["cpu"]
    params:
        index  	     = config["ref"]["index"],
        bowtie2 	 = config["params"]["bowtie2"]["global"],
        samblaster   = config["params"]["samblaster"],
        reads  	     = set_reads,
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


rule genrich_sort:
    input:
       "results/02_bam/{sample}.bam",
    output:
        bam = temp("results/02_bam/{sample}.bam.tmp"),
    threads:
        CLUSTER["align"]["cpu"]
    params:
        samtools_mem = config["params"]["samtools"]["memory"],
    log:
        sort   = "results/00_log/genrich_sort/{sample}_bam_sort_genrich.log",
    shell:
        """
        samtools sort -n -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} {input} 2>> {log.sort}
        """

