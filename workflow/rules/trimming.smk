# ----- Copy the fastq files to a temporary directory ----- #
# This rule will take the output of the function get_fastq() and will do a link (ln) to a temporary directory
rule get_fastq_pe:
    input:
        get_fastq,
    output:
        fastq1 = temp("fastq/{sample}-{lane}.r1.fq"),
        fastq2 = temp("fastq/{sample}-{lane}.r2.fq"),
    message:   
        "Copying .fastq files for {wildcards.sample}-{wildcards.lane}"
    shell:
        """
        ln -s {input[0]} {output.fastq1}
        ln -s {input[1]} {output.fastq2}
        """
rule get_fastq_se:
    input:
        get_fastq,
    output:
        temp("results/fastq/{sample}-{lane}.fastq.gz"),
    message:
        "Copying fastq files {input}"
    shell:
        """
        ln -s {input} {output}
        """

# ----- Merging the different lanes from the .fastq files and trimming with fastp ----- #
rule fastp_trim_pe:
    input:
        fw = lambda wildcards: expand("fastq/{row.sample}-{row.lane}.r1.fq", row=units.loc[wildcards.sample].itertuples()),
        rv = lambda wildcards: expand("fastq/{row.sample}-{row.lane}.r2.fq", row=units.loc[wildcards.sample].itertuples())
    output:
	    fastq1 = temp("results/fastq/{sample}.r1.fastq"),
	    fastq2 = temp("results/fastq/{sample}.r2.fastq")
    log:
        "results/00_log/fastp/{sample}_trim.log",
    threads: 3
    params:
        fastp_params = config["params"]["fastp"]["pe"],
        tmp_fw       = "results/fastq/{sample}.1.fastq.tmp.gz",
        tmp_rv       = "results/fastq/{sample}.2.fastq.tmp.gz"
    message:
        "Processing fastq files from {input}"
    shell:
        """
		cat {input.fw} > {params.tmp_fw}
		cat {input.rv} > {params.tmp_rv}
		fastp -i {params.tmp_fw} \
		-I {params.tmp_rv} \
		-o {output.fastq1} \
		-O {output.fastq2} \
		-w {threads} \
		{params.fastp_params} 2> {log}
        """
rule fastp_trim_se:
	input:
		lambda w: expand("results/fastq/{row.sample}-{row.lane}.fastq.gz", row=units.loc[w.sample].itertuples()),
	output:
		temp("results/fastq/{sample}.se.fastq")
	log:
		"results/00_log/fastp/{sample}_trim.log"
	threads: 5
	params:
		fastp_params = config["params"]["fastp"]["se"],
	message:
		"Processing fastq files from {input}"
	shell:
		"""
		zcat {input} | \
		fastp -o {output} \
		-w {threads} {params.fastp_params} 2> {log}
		"""
