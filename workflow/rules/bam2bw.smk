## BAM2BW
## ======

# ----- .bam to .bigwig with custom Python script (no input subtraction) ----- #
rule bam2bigwig_noSubstract:
    input: 
        unpack(get_bam)
    output:  
        "results/02_bigwig/noSubtract/{sample}.bw"
    params: 
        read_exten = set_read_extension,
        reads      = set_reads_spike2,
        params     = config["params"]["bam2bigwig"]["other"]
    log: 
        "results/00_log/bam2bw/{sample}_bigwig.bam2bw"
    threads: 
        CLUSTER["bam2bigwig"]["cpu"]
    message: 
        "making input subtracted bigwig for sample {wildcards.sample}"
    shell:
        """
        python {params.reads} \
        --case {input.case} \
        --bigwig {output} \
        --threads {threads} \
        --otherParams {params.read_exten} {params.params} &> {log}
        """

# ================
# BigWig to server
# ================

# ----- Put .bw files to the folder where they are read to be visualized in UCSC ----- # 
rule bigwig2server:
    input: 
        bw         = "results/02_bigwig/noSubtract/{sample}.bw",
        samblaster = "results/00_log/alignments/rm_dup/{sample}.log",
        bowtie     = "results/00_log/alignments/{sample}.log"
    output:
        temp("results/temp_file_{sample}_{control}.txt")
    params:
        user     = lambda wildcards : SAMPLES.USER[wildcards.sample],
        antibody = lambda wildcards : SAMPLES.AB[wildcards.sample],
        genome   = lambda wildcards : SAMPLES.GENOME[wildcards.sample],
        run      = lambda wildcards : SAMPLES.RUN[wildcards.sample],
        chip     = lambda wildcards : str("ChIPseq") if SAMPLES.SPIKE[wildcards.sample] == False else str("ChIPseqSpike")
    run:
        # Get number of removed reported reads by bowtie
        with open(input.bowtie,"r") as fi:
            for ln in fi:
                if ln.startswith("# reads with at least one reported alignment:"):
                    nreads = int( str.split(ln)[8] )

        # Get number of removed reads
        with open(input.samblaster,"r") as fi:
            for ln in fi:
                if ln.startswith("samblaster: Removed "):
                    removed_reads = int( str.split(ln)[2] )

        # Total number of final reads is reported by bowtie minus duplicated removed
        total_reads = nreads-removed_reads

        shell(
            "cp {input} \
            /hpcnfs/data/DP/UCSC_tracks/Data/bigWig/{sample}_{control}_{user}_{nreads}_{chip}_{antibody}_{genome}_{run}.bigWig".format(
            input    = input.bw,
            sample   = wildcards.sample,
            control  = wildcards.control,
            user     = params.user,
            nreads   = total_reads,
            chip     = params.chip,
            antibody = params.antibody,
            genome   = params.genome,
            run      = params.run)
            )
        shell("touch {output}".format(output = output))
