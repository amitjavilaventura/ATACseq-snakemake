# ----- Aligning with BOWTIE and processing with SAMTOOLS ----- #
# Aligning with bowtie
# Reomving PCR duplicates with samblaster
# Sorting and indexing with samtool
rule align:
    input:
        get_fq,
    output:
        bam = "results/02_bam/{sample}.bam",
        bai = "results/02_bam/{sample}.bam.bai",
    threads:
        CLUSTER["align"]["cpu"]
    params:
        index        = config["ref"]["index"],
        bowtie       = config["params"]["bowtie"]["global"],
        reads        = set_reads,
        samtools_mem = config["params"]["samtools"]["memory"],
    log:
        align   = "results/00_log/align/{sample}_align.log",
        rm_dups = "results/00_log/align/{sample}_rm-dup.log",
        index   = "results/00_log/align/{sample}_index.log",
    shell:
        """
        bowtie -p {threads} {params.bowtie} {params.index} {params.reads} 2> {log.align} \
        | samblaster --removeDups 2> {log.rm_dups} \
        | samtools view -Sb -F 4 - \
        | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
        samtools index {output.bam} 2> {log.index}
        """

rule genrich_sort:
    input:
       "results/02_bam/{sample}.bam",
    output:
        bam = "results/02_bam/{sample}.sorted.bam",
    threads:
        CLUSTER["align"]["cpu"]
    params:
        index        = config["ref"]["index"],
        bowtie       = config["params"]["bowtie"]["global"],
        reads        = set_reads,
        samtools_mem = config["params"]["samtools"]["memory"],
    log:
        sort   = "results/00_log/genrich_sort/{sample}_bam_sort_genrich.log",
    shell:
        """
        samtools sort -n -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} {input} 2>> {log.sort}
        """

