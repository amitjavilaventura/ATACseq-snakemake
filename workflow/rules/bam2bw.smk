rule bam2bigwig:
    input: 
        nfr  = "results/02_bam/{sample}_nfr.bam",
        nucl = "results/02_bam/{sample}_nucleosomes.bam",
    output:  
        nfr  = "results/05_bigwig/{sample}_nfr.bw",
        nucl = "results/05_bigwig/{sample}_nucleosomes.bw",
    params: 
        params     = config["params"]["bam2bigwig"]["other"]
    log: 
        "results/00_log/bam2bw/{sample}_bigwig.bam2bw"
    threads: 
        CLUSTER["bam2bigwig"]["cpu"]
    message: 
        "making input subtracted bigwig for sample {wildcards.sample}"
    shell:
        """
        python workflow/scripts/bam2bigwig.py \
        --case {input.nfr} \
        --bigwig {output.nfr} \
        --threads {threads} \
        --otherParams --extendReads {params.params} &> {log}

        python workflow/scripts/bam2bigwig.py \
        --case {input.nucl} \
        --bigwig {output.nucl} \
        --threads {threads} \
        --otherParams --extendReads {params.params} &> {log}
        """


# Create bigwig for pooled replicates
rule bam2bigwig_merge:
    input: 
        nfr  = lambda w: expand("results/02_bam/{sample}_nfr.bam", sample = SAMPLES.loc[SAMPLES["CONDITION"] == w.condition].NAME),
        nucl = lambda w: expand("results/02_bam/{sample}_nucleosomes.bam", sample = SAMPLES.loc[SAMPLES["CONDITION"] == w.condition].NAME),
    output:  
        nfr  = "results/05_bigwig/merged/{condition}_nfr_merged.bw",
        nucl = "results/05_bigwig/merged/{condition}_nucleosomes_merged.bw",
    params: 
        params     = config["params"]["bam2bigwig"]["other"],
        samtools_mem = config["params"]["samtools"]["memory"],
    log: 
        "results/00_log/bam2bw/{condition}_bigwig.bam2bw"
    threads: 
        CLUSTER["bam2bigwig"]["cpu"]
    shell:
        """
        # Merge, sort rand index replicates
        samtools merge -f {wildcards.condition}_nfr_merged_tmp.bam {input.nfr} 
        samtools merge -f {wildcards.condition}_nucleosomes_merged_tmp.bam {input.nucl}
        samtools index {wildcards.condition}_nfr_merged_tmp.bam
        samtools index {wildcards.condition}_nucleosomes_merged_tmp.bam

        # bam2bw
        python workflow/scripts/bam2bigwig.py \
        --case {wildcards.condition}_nfr_merged_tmp.bam \
        --bigwig {output.nfr} \
        --threads {threads} \
        --otherParams --extendReads {params.params} &> {log}
        rm {wildcards.condition}_nfr_merged_tmp.bam* # Remove tmp bam

        python workflow/scripts/bam2bigwig.py \
        --case {wildcards.condition}_nucleosomes_merged_tmp.bam \
        --bigwig {output.nucl} \
        --threads {threads} \
        --otherParams --extendReads {params.params} &> {log}
        rm {wildcards.condition}_nucleosomes_merged_tmp.bam*
        """


rule bigwig2server:
    input: 
        bw         = "results/05_bigwig/{sample}.bw",
        samblaster = "results/00_log/align/rm_dup/{sample}.log",
        bowtie2    = "results/00_log/align/{sample}.log"
    output:
        temp("results/temp_file_{sample}.txt")
    params:
        user     = lambda wildcards : SAMPLES.USER[wildcards.sample],
        genome   = lambda wildcards : config["ref"]["genome"],
        run      = lambda wildcards : SAMPLES.RUN[wildcards.sample],
        chip     = lambda wildcards : "ATACseq"
    run:
        # Get number of removed reported reads by bowtie
        with open(input.bowtie,"r") as fi:
            for ln in fi:
                lines = str.split(ln)
                if "exactly" in lines:
                    nreads = int(lines[0])
                if ">1" in lines:
                    nreads += int(lines[0])
        # Get number of removed reads
        with open(input.samblaster,"r") as fi:
            for ln in fi:
                if ln.startswith("samblaster: Marked "):
                    removed_reads = int( str.split(ln)[2] )

        # Total number of final reads is reported by bowtie minus duplicated removed
        total_reads = nreads-removed_reads
        shell(
            "cp {input} \
            /hpcnfs/data/DP/UCSC_tracks/Data/bigWig/{sample}_{control}_{user}_{nreads}_{chip}_{antibody}_{genome}_{run}.bigWig".format(
            input    = input.bw,
            sample   = wildcards.sample,
            control  = "X",
            user     = params.user,
            nreads   = total_reads,
            chip     = params.chip,
            antibody = "X",
            genome   = params.genome,
            run      = params.run)
            )
        shell("touch {output}".format(output = output))
