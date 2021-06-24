# ----- Aligning with BOWTIE2 and processing with SAMTOOLS ----- #
# Aligning with bowtie
# Reomving PCR duplicates with samblaster
# Sorting and indexing with samtools

rule align:
    input:
        get_fq,
    output:
        bam   = "results/02_bam/{sample}.bam",
        index = "results/02_bam/{sample}.bam.bai",
    threads:
        CLUSTER["align"]["cpu"]
    params:
        index  	     = config["ref"]["index"],
        bowtie2 	 = config["params"]["bowtie2"]["global"],
        samblaster   = config["params"]["samblaster"],
        reads  	     = set_reads,
        rm_chrM      = config["options"]["rm_chrM"],
        chrM_name    = config["ref"]["chrM_name"],
        samtools_mem = config["params"]["samtools"]["memory"],
    message:
        "Aligning {input} with parameters {params.bowtie2} and {params.chrM_name} removing == {params.rm_chrM}"
    log:
       align   = "results/00_log/align/{sample}.log",
       rm_dups = "results/00_log/align/rm_dup/{sample}.log",
    benchmark:
        "results/.benchmarks/{sample}.align.benchmark.txt"
    shell:
        """
        # -q2 removes multimapping, -F 4 removes unmapped, -f 2 removes unpaired reads.

        bowtie2 -p {threads} {params.bowtie2} -x {params.index} {params.reads} 2> {log.align} \
            | samtools view -h - \
            | samblaster {params.samblaster} 2> {log.rm_dups} \
            | samtools view -q2 -Sb -F 4 - \
            | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}

        samtools index {output.bam}          

        """

### Remove mitochondrial chromosome and change name of the BAMS in order to be 
###   able to have the filtered BAMS with the same name as the original.
if config["options"]["rm_chrM"]:
    rule remove_chrM:
        input:
            bam = rules.align.output.bam
        output:
            bam_filt = temp("results/02_bam/{sample}_noChrM.bam"),
            bam_indx = temp("results/02_bam/{sample}_noChrM.bam.bai"),
        log:
            log = "results/00_log/remove_mt/{sample}"
        params:
            chrM_name = config["ref"]["chrM_name"]
        shell:
            """
            # remove reads mapping to chromosome M 
            samtools view -h {input} | fgrep -w -v {params.chrM_name} | samtools view -b > {output.bam_filt} 2> {log}
            samtools index {output.bam_filt}
            """