# ----- Aligning with BOWTIE2 and processing with SAMTOOLS ----- #
# Aligning with bowtie
# Reomving PCR duplicates with samblaster
# Sorting and indexing with samtool

rule align:
    input:
        get_fq
    output:
         bam   = "results/02_bam/{sample}.bam",
         idx   = "results/02_bam/{sample}.bam.bai",
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
    shell:
        """
        # -q2 removes multimapping, -F 4 removes unmapped, -f 2 removes unpaired reads.

        bowtie2 -p {threads} {params.bowtie2} -x {params.index} {params.reads} 2> {log.align} \
        | samblaster --ignoreUnmated {params.samblaster} 2> {log.rm_dups} \
        | samtools view -q2 -Sb -F 4 -f 2 - \
        | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
        samtools index {output.bam}
        """


rule remove_chrM:
    input:
        bam = "results/02_bam/{sample}.bam"
    output:
        bam_filt = "results/02_bam/{sample}_noMT.bam",
        bam_indx = "results/02_bam/{sample}_noMT.bam.bai"
    log:
        log = "results/00_log/remove_mt/{sample}"
    params:
        chrM_name = config["ref"]["chrM_name"]
    shell:
        """
        samtools view -h {input} | fgrep -w -v {params.chrM_name} | samtools view -b > {output.bam_filt} 2> {log}
        samtools index {output.bam_filt}
        """


# Rule to split bam into nucleosome free regions (nfr) and nucleosomes. 
# The bam is divided into fragments of 150bp or less (nfr) or 151bp or more (nucleosomes).
rule split_bam:
    input:
        bam      = "results/02_bam/{sample}_noMT.bam",
        bam_indx = "results/02_bam/{sample}_noMT.bam.bai"
    output:
        nfr        = "results/02_bam/{sample}_nfr.bam",
        nucleo     = "results/02_bam/{sample}_nucleosomes.bam",
        nfr_idx    = "results/02_bam/{sample}_nfr.bam.bai",
        nucleo_idx = "results/02_bam/{sample}_nucleosomes.bam.bai",
        metrics    = "results/01_QCs/split_bam/{sample}_metrics.txt"
    params:
        samtools_mem = config["params"]["samtools"]["memory"],
    threads:
        CLUSTER["split_bam"]["cpu"]
    log:
        log = "results/00_log/split_bam/{sample}.log"
    shadow:
        "minimal"
    shell:
        """
        alignmentSieve --bam {input.bam} \
        --outFile {output.nfr}.tmp -p {threads} \
        --ATACshift \
        --filterMetrics {output.metrics} \
        --maxFragmentLength 150 \
        --filteredOutReads {output.nucleo}.tmp \
        2> {log}

        # Index the output bam files because it will be used to create the bigwig
        samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.nucleo}_tmp -o {output.nucleo} {output.nucleo}.tmp 2>> {log}
        samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.nfr}_tmp -o {output.nfr} {output.nfr}.tmp 2>> {log}

        samtools index {output.nucleo}
        samtools index {output.nfr}
        """


rule genrich_sort:
    input:
        nfr     = "results/02_bam/{sample}_nfr.bam",
        nucleo  = "results/02_bam/{sample}_nucleosomes.bam",
    output:
        nfr_sort    = "results/02_bam/{sample}_nfr.bam.sorted",
        nucleo_sort = "results/02_bam/{sample}_nucleosomes.bam.sorted"
    threads:
        CLUSTER["genrich_sort"]["cpu"]
    params:
        samtools_mem = config["params"]["samtools"]["memory"],
    log:
        sort   = "results/00_log/genrich_sort/{sample}_bam_sort_genrich.log",
    shell:
        """
        samtools sort -n -m {params.samtools_mem}G -@ {threads} -T {output.nfr_sort}.tmp -o {output.nfr_sort} {input.nfr} 2> {log.sort}
        samtools sort -n -m {params.samtools_mem}G -@ {threads} -T {output.nucleo_sort}.tmp -o {output.nucleo_sort} {input.nucleo} 2>> {log.sort}
        """

