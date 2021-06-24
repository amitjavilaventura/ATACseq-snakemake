# SPLIT BAMS AND CALL PEAKS TO FIND NFR AND NUCLEOSOME REGIONS

# --- Split bam files --- #
# Rule to split bam into nucleosome free regions (nfr) and nucleosomes. 
# The bam is divided into fragments of 150bp or less (nfr) or 151bp or more (nucleosomes).
if config["options"]["splitBam"]:
    rule split_bam:
        input:
            bam   = "results/02_bam/{sample}_noChrM.bam" if config["options"]["rm_chrM"] else "results/02_bam/{sample}.bam",
            index = "results/02_bam/{sample}_noChrM.bam.bai" if config["options"]["rm_chrM"] else "results/02_bam/{sample}.bam.bai",
        output:
            nfr     = "results/02_bam/splitBam/nfr/{sample}_NFR.bam",
            nuc     = "results/02_bam/splitBam/nuc/{sample}_NUC.bam",
            nfr_idx = "results/02_bam/splitBam/nfr/{sample}_NFR.bam.bai",
            nuc_idx = "results/02_bam/splitBam/nuc/{sample}_NUC.bam.bai",
            metrics = "results/01_QCs/splitBam/{sample}_metrics.txt",
        params:
            samtools_mem = config["params"]["samtools"]["memory"],
            nfr_length   = config["params"]["splitbam"]["nfr_length"]
        threads:
            CLUSTER["split_bam"]["cpu"]
        log:
            log = "results/00_log/split_bam/{sample}.log"
        shadow:
            "minimal"
        shell:
            """
            # Correct Tn5 insertion bias 
            # Filter nucleosome free regions and nucleosome regions using fragment lenght = 150
            alignmentSieve --bam {input.bam} \
            --outFile {output.nfr}.tmp -p {threads} \
            --ATACshift \
            --filterMetrics {output.metrics} \
            --maxFragmentLength {params.nfr_length} \
            --filteredOutReads {output.nuc}.tmp \
            2> {log.log}

            # Sort the BAM files to sort with Genrich
            samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.nuc}_tmp -o {output.nuc} {output.nuc}.tmp 2>> {log.log}
            samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.nfr}_tmp -o {output.nfr} {output.nfr}.tmp 2>> {log.log}
            samtools index {output.nuc}
            samtools index {output.nfr}
            """

    rule splitBam_genrich_sort:
        input:
            nfr     = "results/02_bam/splitBam/nfr/{sample}_NFR.bam",
            nuc     = "results/02_bam/splitBam/nuc/{sample}_NUC.bam",
        output:
            nfr_sort = temp("results/02_bam/splitBam/nfr/{sample}_NFR.bam.tmp"),
            nuc_sort = temp("results/02_bam/splitBam/nuc/{sample}_NUC.bam.tmp"),
        threads:
            CLUSTER["genrich_sort"]["cpu"]
        params:
            samtools_mem = config["params"]["samtools"]["memory"],
        log:
            sort = "results/00_log/genrich_sort/{sample}_splitBam_sort_genrich.log",
        shell:
            """
            samtools sort -n -m {params.samtools_mem}G -@ {threads} -T {output.nfr_sort}.tmp -o {output.nfr_sort} {input.nfr} 2> {log.sort}
            samtools sort -n -m {params.samtools_mem}G -@ {threads} -T {output.nuc_sort}.tmp -o {output.nuc_sort} {input.nuc} 2>> {log.sort}
            """



    rule genrich_splitBam:
        input:
            bams_nfr = lambda w: expand("results/02_bam/splitBam/nfr/{condition}_NFR.bam.tmp", condition = SAMPLES.loc[SAMPLES["CONDITION"] == w.condition].NAME) if config["options"]["genrich_pool"] else "results/02_bam/splitBam/{sample}_NFR.bam.tmp",
            bams_nuc = lambda w: expand("results/02_bam/splitBam/nuc/{condition}_NUC.bam.tmp", condition = SAMPLES.loc[SAMPLES["CONDITION"] == w.condition].NAME) if config["options"]["genrich_pool"] else "results/02_bam/splitBam/{sample}_NUC.bam.tmp",
        output:
            peak_nfr = "results/03_genrich/" + "splitBam/pooled/{condition}/{condition}_NFR_peaks.narrowPeak" if config["options"]["genrich_pool"] else "splitBam/{sample}/{sample}_NFR_peaks.narrowPeak",
            peak_nuc = "results/03_genrich/" + "splitBam/pooled/{condition}/{condition}_NUC_peaks.narrowPeak" if config["options"]["genrich_pool"] else "splitBam/{sample}/{sample}_NUC_peaks.narrowPeak",
        params:
            # path to genrich since it's in the singularity image it's not needed.
            #genrich  = config["params"]["genrich"]["path"],
            # pvalue/qvalue threshold
            p_or_q   = config["params"]["genrich"]["p_or_q"],
            pqval    = config["params"]["genrich"]["pqval"],
            # rm chr M reads
            chrM     = config["params"]["genrich"]["chrM"],
            # rm pcr duplicates
            PRCdups = config["params"]["genrich"]["PRCdups"],
            # readme file
            readme   = "results/03_genrich/" + "splitBam/pooled/README.md" if config["options"]["genrich_pool"] else "splitBam/README.md"
        log:
            "results/00_log/genrich_pool/{condition}_splitBam_peakcalling.log" if config["options"]["genrich_pool"] else "results/00_log/genrich/{sample}_splitBam_peakcalling.log",
        shell:
            """
            ## CALL PEAKS WITHOUT ATAC MODE
            ## DO NOT USE UNMATCHED REAEDS

            ## PEAKS FROM NUC-FREE REGIONS
            Genrich \
            -t '{input.bams_nfr}' \
            -o {output.peak_nfr} \
            -{params.p_or_q} {params.pqval} \
            {params.chrM} {params.PRCdups} -v 2> {log}

            ## PEAKS FROM NUCLEOSOME REGIONS
            Genrich \
            -t '{input.bams_nuc}' \
            -o {output.peak_nuc} \
            -{params.p_or_q} {params.pqval} \
            {params.chrM} {params.PRCdups} -v 2> {log}

            echo '''

            # Genrich NFR/NUC peak calling

            This folder contains the peaks called from splitted BAM files into nucleosome-free and nucleosome regions

            ## Inputs

            * **NFR**: {input.bams_nfr}

            * **NUC**: {input.bams_nuc}

            ## Outputs

            * **NFR**: {output.peak_nfr}

            * **NUC**: {output.peak_nuc}

            ## Params

            * ATAC mode: NO

            * Use unmatched reads: NO

            * p- or q-value: {params.p_or_q} -> {params.pqval}

            * Remove PRC dups: {params.PRCdups}

            * Remove mitocondrial reads: {params.chrM}
            ''' > {params.readme}
            """