################################
# ----- QUALITY CONTROLS ----- #
################################

# Here I want to create a rule to plot a bar graph with the number of reads aligned in each chromosome.
# Has to use a function of amitjavilaventura called chromReads() in the chromReads.R script.

## Most things are copied from dfernadezperez chipseq pipeline
# ----- Function get_fq_forward(). ----- #
# This function will return the forward read, in the case of Paired-End, in order to analize it with fastqc. 
# To use as input for the "rule fastqc" when updating the pipeline to differentiate single-end and paired-end. Now it is only for PE.
# Get raw or trimmed reads based on trimming configuration. Used for fastqc
def get_fq_forward(wildcards):
    if config["options"]["trimming"]:
        if not is_single_end(**wildcards):
            # paired-end sample
            return "{tmp}/fastq/trimmed{sample}.1.fastq.gz".format(**wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/trimmed{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)
    else:
        # no trimming, use raw reads
        if not is_single_end(**wildcards):
            # paired-end sample
            return "{tmp}/fastq/{sample}.1.fastq.gz".format(**wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)

# ------- FASTQC ------- #
rule fastqc:
    """This rule is a little bit tricky to fit PE fastq in multiqc.
    Basically the probelm is that fastqc needs to run R1 and R2 separatelly,
    which means 2 fastqc_zip files with different names. This will be recognized
    by multiqc as different samples, so the report will be a mess.
    My workaround has been to use just the forward fastq to create the report.
    For this I need to change the fastq file name (because it has .1.) to fit
    what multiqc is expecting as name. If multiqc reads A.1.fastq it won't know
    that that file must match A.bam in the report, and they will be in different
    rows. I know it's a pain of workaround but it's the only solution I found.
    For this I need to create a symlink to the fastq to be able to change it's name
    (without duplicating the file which would be less efficient). To make sure that 
    the symlink will work it needs to be created from the folder where it's going to be,
    that's why the cd command of the rule it's imporant. Since the fastq folder can change
    this step needs to work always, it's the only solution I cam up with.
    """
    input:  
        get_fq_forward
    output: 
        "results/01_QCs/fastQC/{sample}_fastqc.zip"
    log:    
        "results/00_log/fastQC/{sample}.log"
    params:
        folder_name = "results/01_QCs/fastQC/",
        tmp = "{sample}.fastq.gz"
    threads: 
        CLUSTER["fastqc"]["cpu"]
    message: 
        "Running fastqc for {input}"
    shadow: 
        "minimal"
    shell:
        """
        cd {params.folder_name} # Move to folder where symlink is going to be created
        ln -s {input} {params.tmp} # Create symlink to fastq file. Imporant to set the desired file name.
        cd - # Go back to workdir
        fastqc -o {params.folder_name} -f fastq -t {threads} --noextract {params.folder_name}/{params.tmp} 2> {log}
        """

# # ------- InsertSize calculation ------- #
rule insert_size:
    input:
        "results/02_bam/{sample}.bam"
    output:
        txt="results/01_QCs/insert_size/{sample}.isize.txt",
        pdf="results/01_QCs/insert_size/{sample}.isize.pdf"
    log:
        "results/00_log/picard/insert_size/{sample}.log"
    params:
        # optional parameters (e.g. relax checks as below)
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE"
    shell:
        """
        # Create the outfiles to handle
        touch {output}
        picard CollectInsertSizeMetrics {params} \
        INPUT={input} OUTPUT={output.txt} \
        HISTOGRAM_FILE={output.pdf} > {log}
        """

# # # ----- Function set_read_extension() ----- #
# # # This function will be used as parameter for the rule plotFingerprint.
# # # The params.read_exten have to change from "--extendReads" to "set_read_extension()"
# # # It will return --extendReads and, if the sample is single-end, it will return more parameters. 
# # # To use when the snakemake pipeline is update to be used with both single-end and paired-end. Not now.
# # def set_read_extension(wildcards):
# #     if is_single_end(wildcards.sample):
# #         return "--extendReads " + str(config['bam2bigwig']['read_extension'])
# #     return "--extendReads"

# # ------- Deeptools quality control ------- #
rule plotFingerprint:
    input: 
        case      = "results/02_bam/{sample}.bam", 
    output: 
        qualMetrics = "results/01_QCs/fingerPrint/{sample}.qualityMetrics.tsv",
        raw_counts  = "results/01_QCs/fingerPrint/{sample}.rawcounts.tsv",
        plot        = "results/01_QCs/fingerPrint/{sample}.plot.pdf",
    log:
        "results/00_log/plotFingerprint/{sample}.log"
    params:
        read_exten = "--extendReads",
    threads:
        CLUSTER["phantom_peak_qual"]["cpu"]
    shell:
        """
        plotFingerprint -b {input} \
        -p {threads} \
        --outQualityMetrics {output.qualMetrics} \
        --outRawCounts {output.raw_counts} \
        --plotFile {output.plot}
        """

rule GC_bias:
    input: 
        bam = "results/02_bam/{sample}.bam",
        bed = rules.filter_peaks.output.macs2_filt
    output: 
        pdf      = "results/01_QCs/GCbias/{sample}_GCbias.pdf",
        freq_txt = "results/01_QCs/GCbias/{sample}_GCbias.txt"
    log:
        "results/00_log/GCbias/{sample}_GCbias.log"
    params:
        repeatMasker = config["ref"]['rep_masker'],
        tempBed      = "results/01_QCs/GCbias/{sample}_Repeatmasker.bed.tmp",
        bit_file     = config["ref"]["2bit"],
        egenome_size = config["ref"]["egenome_size"]
    threads: 5
    message:
        "Computing GC bias for sample {wildcards.sample}"
    shell:
        """
        bedops -u {input.bed} {params.repeatMasker} > {params.tempBed}
        bp_peaks=$(bedops --merge {input.bed} | bedmap --bases - | awk "{{sum+=\$1}}END{{print sum}}")
        total_eGsize=$(({params.egenome_size}-$bp_peaks))

        computeGCBias -b {input.bam} \
            -p {threads} \
            --effectiveGenomeSize $total_eGsize \
            -g {params.bit_file} \
            -l 200 \
            -bl {params.tempBed} \
            --biasPlot {output.pdf} \
            --GCbiasFrequenciesFile {output.freq_txt} 2> {log}
        rm -f {params.tempBed}
        """

# ---------------- MultiQC report ----------------- #
# if config["options"]["genrich_merge"] == False:

rule multiQC_inputs:
    input:
        expand("results/00_log/align/{sample}.log", sample = ALL_SAMPLES),
        expand("results/01_QCs/fastQC/{sample}_fastqc.zip", sample = ALL_SAMPLES),
        expand("results/01_QCs/insert_size/{sample}.isize.txt", sample = ALL_SAMPLES),
        expand("results/00_log/align/rm_dup/{sample}.log", sample = ALL_SAMPLES),
        expand("results/01_QCs/fingerPrint/{sample}.qualityMetrics.tsv", zip, sample = ALL_SAMPLES),
        expand("results/01_QCs/fingerPrint/{sample}.rawcounts.tsv", zip, sample = ALL_SAMPLES),
        #expand("results/03_macs2/{sample}/{sample}_peaks.xls", zip, sample = ALL_SAMPLES) #this is done with genrich
    output: 
        file = "results/01_QCs/multiQC/multiQC_inputs.txt"
    message:
        "create file containing all multiqc input files"
    run:
        with open(output.file, 'w') as outfile:
            for fname in input:
                outfile.write(fname + "\n")

rule multiQC:
    input:
        "results/01_QCs/multiQC/multiQC_inputs.txt"
    output: 
        "results/01_QCs/multiQC/multiQC_report.html"
    params:
        log_name = "multiQC_report",
        folder   = "results/01_QCs/multiQC"
    log:
        "results/00_log/multiQC/multiQC.log"
    message:
        "multiQC for all logs"
    shell:
        """
        multiqc -o {params.folder} -l {input} -f -v -n {params.log_name} 2> {log}
        """

# elif config["options"]["genrich_merge"] == True:

#     rule multiQC_inputs:
#         input:
#             expand("results/00_log/align/{sample}_align.log", sample = ALL_SAMPLES),
#             expand("results/01_QCs/fastQC/{sample}_fastqc.zip", sample = ALL_SAMPLES),
#             expand("results/01_QCs/insert_size/{sample}.isize.txt", sample = ALL_SAMPLES),
#             expand("results/01_QCs/phantom_peak_qual/{sample}.spp.out", sample = ALL_SAMPLES),
#             expand("results/00_log/align/{sample}_rm-dup.log", sample = ALL_SAMPLES),
#             expand("results/01_QCs/fingerPrint/{sample}.qualityMetrics.tsv", zip, sample = ALL_SAMPLES),
#             expand("results/01_QCs/fingerPrint/{sample}.rawcounts.tsv", zip, sample = ALL_SAMPLES),
#             #expand("results/03_macs2/{sample}/{sample}_peaks.xls", zip, sample = ALL_SAMPLES) #this is done with genrich
#         output: 
#             file = "results/01_QCs/multiQC/multiQC_inputs_genrich_merge.txt"
#         message:
#             "create file containing all multiqc input files"
#         run:
#             with open(output.file, 'w') as outfile:
#                 for fname in input:
#                         outfile.write(fname + "\n")


#     rule multiQC:
#         input:
#             "results/01_QCs/multiQC/multiQC_inputs_genrich_merge.txt"
#         output: 
#             "results/01_QCs/multiQC/multiQC_report_genrich_merge.html"
#         params:
#             log_name = "multiQC_report",
#             folder   = "results/01_QCs/multiQC"
#         log:
#             "results/00_log/multiQC/multiQC.log"
#         message:
#             "multiQC for all logs"
#         shell:
#             """
#             multiqc -o {params.folder} -l {input} -f -v -n {params.log_name} 2> {log}
#             """