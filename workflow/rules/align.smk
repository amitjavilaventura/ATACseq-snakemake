# ----- Aligning with BOWTIE and processing with SAMTOOLS ----- #
# Aligning with bowtie
# Reomving PCR duplicates with samblaster
# Sorting and indexing with samtools
rule align:
    input:
        get_trimmed,
    output:
        bam = "results/02_bam/{sample}.bam",
        bai = "results/02_bam/{sample}.bam.bai",
    threads: 10
    params:
        index  = config["ref"]["index"],
        bowtie = config["params"]["bowtie"]["global"],
        reads  = set_reads,
    log:
        align   = "results/00_log/align/{sample}_align.log",
        rm_dups = "results/00_log/align/{sample}_rm-dup.log",
        index   = "results/00_log/align/{sample}_index.log",
    shell:
        """
        bowtie -p {threads} {params.bowtie} {params.index} {params.reads} 2> {log.align} \
        | samblaster --removeDups 2> {log.rm_dups} \
        | samtools view -Sb -F 4 - \
        | samtools sort -m 5G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
        samtools index {output.bam} 2> {log.index}
        """

# # ----- Get rid of reads that align to chrM ----- #
# # Try it before applying it to the pipeline. 
#
# rule align_noChrM:
#    input:
#        get_trimmed,
#    output:
#        bam = "results/02_bam/{sample}/{sample}_noChrM.bam",
#        bai = "results/02_bam/{sample}/{sample}_noChrM.bam.bai",
#     threads: 10
#     params:
#         index  = config["ref"]["index"],
#         bowtie = config["params"]["bowtie"]["global"],
#         reads  = set_reads,
#     log:
#         align   = "results/00_log/alig_noChrM/{sample}_align.log",
#         rm_dups = "results/00_log/align_noChrM/{sample}_rm-dup.log",
#         index   = "results/00_log/align_noChrM/{sample}_index.log",
#     shell:
#         """
#         bowtie -p {threads} {params.bowtie} {params.index} {params.reads} 2> {log.align} \
#         | samblaster --removeDups 2> {log.rm_dups} \
#         | samtools view -Sb -F 4 - \
#         | samtools sort -m 5G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
#         samtools view -h {output.bam} \
#         | awk '{if($3 != "chrM" && $3 != "chrUn"){print $0}}' \
#         | samtools view -Shb - > {output.bam} 2> {log.align}
#         samtools index {output.bam} 2> {log.index}
#         """

# # ----- Shift reads in (+) strand by +4 bp and reads in (-) strand by -5 bp ----- #
# # This is done because Tn5 binds as a dimer and inserts to adaptors separated by 9bp.
# #   If we want to know the exact Tn5 insertion site, we want to shift the reads consequently.
# #   From Buenrostro et al. https://www.nature.com/articles/nmeth.2688#ref11
# # Rule needs to be deffined.
