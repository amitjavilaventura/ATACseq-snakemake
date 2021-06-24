#########################################
# ----- Peak calling with Genrich ----- #
#########################################

### Using GENRICH ###
### ============= ###

# https://github.com/jsh58/Genrich

### To do: 
###    allow to input a control --> make a function for input and a function for params

rule sortbam_genrich:
    input:
       rules.align.output.bam,
    output:
        bam = temp("results/02_bam/{sample}.bam.tmp"),
    threads:
        CLUSTER["genrich_sort"]["cpu"]
    params:
        samtools_mem = config["params"]["samtools"]["memory"],
    log:
        sort   = "results/00_log/genrich_sort/{sample}_bam_sort_genrich.log",
    shell:
        """
        samtools sort -n -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} {input} 2>> {log.sort}
        """


rule genrich:
	input:
		# if config["options"]["genrich_pool"] is true, it calls the all the bamfiles in each condition, if it's false, it calls each bamfile individually
		# if config["options"]["rm_chrM"] is true, calls the bam file without chromosome reads, else it will call the original bam (have in mind that only one is creeated in the alignment step)
		t = lambda w: expand("results/02_bam/{condition}.bam.tmp", condition = SAMPLES.loc[SAMPLES["CONDITION"] == w.condition].NAME) if config["options"]["genrich_pool"] else "results/02_bam/{sample}.bam.tmp"
	output:
		"results/03_genrich/" + "pooled/{condition}/{condition}_peaks.narrowPeak" if config["options"]["genrich_pool"] else "{sample}/{sample}_peaks.narrowPeak",
	params:
		# path to genrich since it's in the singularity image it's not needed.
		#genrich  = config["params"]["genrich"]["path"],
		atacmode = config["params"]["genrich"]["atac_mode"],
		unmatch  = lambda w: config["params"]["genrich"]["unmatch"],
		# pvalue/qvalue threshold
		p_or_q   = config["params"]["genrich"]["p_or_q"],
		pqval    = config["params"]["genrich"]["pqval"],
		# rm chr M reads
		chrM     = config["params"]["genrich"]["chrM"],
		# rm pcr duplicates
		PRCdups = config["params"]["genrich"]["PRCdups"],
		# readme file
		readme   = "results/03_genrich/" + "pooled/readme.md" if config["options"]["genrich_pool"] else "readme.md"
	log:
		"results/00_log/genrich_pool/{condition}_peakcalling.log" if config["options"]["genrich_pool"] else "results/00_log/genrich/{sample}_peakcalling.log",
	shell:
		"""
		Genrich {params.atacmode} \
		-t '{input.t}' \
		-o {output} \
		-{params.p_or_q} {params.pqval} \
		{params.unmatch} {params.chrM} \
		{params.PRCdups} -v 2> {log}

		echo '''
		
		# Genrich peak calling 
		
		This folder contains the peak called from the complete BAM files. 

		## Inputs

		{input.t}

		## Outputs

		{output}

		## Params

		* ATAC mode: {params.atacmode}. 
		
			+ If the ATAC mode is on (-j), it will mean that the fragments in the BAM correspond to the Tn5-cutting sites, which means that all the peaks called correspond to open chromatin regions. 
			+ If you want the NFR and NUC regions separately, try setting config["options"]["splitBam"] to True. The peaks from the complete BAM will be still called. 

		* Use unmatched reads: {params.unmatch}

			+ Unmatched reads (-y) should be used in case that ATAC mode (-j) is used. Otherways, feel free to use it or not.

		* p- or q-value: {params.p_or_q} -> {params.pqval}

		* Remove PRC dups: {params.PRCdups}

		* Remove mitocondrial reads: {params.chrM}
		''' > {params.readme}
 		"""


pval = config["params"]["genrich"]["filt"]

rule filter_genrich:
		input:
			narrow = rules.genrich.output,
		output:
			bed_filt = "results/03_genrich/" + "pooled/{{condition}}/{{condition}}_peaks_p{pvalue}.bed".format(pvalue = pval) if config["options"]["genrich_pool"] else "{{sample}}/{{sample}}_peaks_p{pvalue}.bed".format(pvalue = pval), 
		params:
			column_pqval   = lambda w: "8" if config["params"]["genrich"]["p_or_q"]=="p" else "9",
			pq_filt        = config["params"]["genrich"]["filt"]
		log:
			"results/00_log/genrich_pool/{condition}_filt_peaks.log" if config["options"]["genrich_pool"] else "results/00_log/genrich/{sample}_filt_peaks.log",
		shell:
			"""
			awk "\${params.column_pqval} >= {params.pq_filt}" {input.narrow} | cut -f1-4,8,9,10 > {output.bed_filt} 2> {log}
			"""

rule peakAnnot_genrich:
		input:
			rules.filter_genrich.output.bed_filt,
		output:
			annot             = "results/04_peakAnno/" + "pooled/{{condition}}/{{condition}}-peaks_logp{pvalue}.annot".format(pvalue = pval) if config["options"]["genrich_pool"] else "{{sample}}/{{sample}}-peaks_logp{pvalue}.annot".format(pvalue = pval),
			promo_bed_targets = "results/04_peakAnno/" + "pooled/{{condition}}/{{condition}}-peaks_logp{pvalue}_promoTargets.bed".format(pvalue = pval) if config["options"]["genrich_pool"] else  "{{sample}}/{{sample}}-peaks_logp{pvalue}_promoTargets.bed".format(pvalue = pval),
			promoTargets      = "results/04_peakAnno/" + "pooled/{{condition}}/{{condition}}-peaks_logp{pvalue}_promoTargets.txt".format(pvalue = pval) if config["options"]["genrich_pool"] else  "{{sample}}/{{sample}}-peaks_logp{pvalue}_promoTargets.txt".format(pvalue = pval),
			promoBed          = "results/04_peakAnno/" + "pooled/{{condition}}/{{condition}}-peaks_logp{pvalue}_promoPeaks.bed".format(pvalue = pval) if config["options"]["genrich_pool"] else  "{{sample}}/{{sample}}-peaks_logp{pvalue}_promoPeaks.bed".format(pvalue = pval),
			distalBed         = "results/04_peakAnno/" + "pooled/{{condition}}/{{condition}}-peaks_logp{pvalue}_distalPeaks.bed".format(pvalue = pval) if config["options"]["genrich_pool"] else  "{{sample}}/{{sample}}-peaks_logp{pvalue}_distalPeaks.bed".format(pvalue = pval),
		params:
			before = config["ref"]["promoter"]["bTSS"],
			after  = config["ref"]["promoter"]["aTSS"],
			genome = config["ref"]["genome"]
		log: 
			"results/00_log/peakAnnot/" + "pooled/{condition}_peakanot_genrich.log" if config["options"]["genrich_pool"] else "{sample}_peakanot_genrich.log",
		message:
			"Annotating peaks for {wildcards.condition}"  if config["options"]["genrich_pool"] else  "Annotating peaks for {wildcards.sample}",
		shell:
			"""
			Rscript --vanilla workflow/scripts/pipeline/peakAnno.R {input} {params.before} {params.after}   \
				{output.annot} {output.promo_bed_targets} {output.promoTargets} {output.promoBed} \
				{output.distalBed} {params.genome} 2> {log}
			"""

rule summits_genrich:
	input:
		narrow   = rules.genrich.output,
		bed_filt = rules.filter_genrich.output,
	output:
		summit   = "results/03_genrich/" + "pooled/{condition}/{condition}_summits.bed" if config["options"]["genrich_pool"] else "{sample}/{sample}_summits.bed",
		sum_filt = "results/03_genrich/" + "pooled/{{condition}}/{{condition}}_summits_p{pvalue}.bed".format(pvalue = config["params"]["genrich"]["filt"]) if config["options"]["genrich_pool"] else "{{sample}}/{{sample}}_summits_p{pvalue}.bed".format(pvalue = config["params"]["genrich"]["filt"]), 
	log:
		"results/00_log/genrich_pool/{condition}_peakcalling.log" if config["options"]["genrich_pool"] else "results/00_log/genrich/{sample}_peakcalling.log",
	shell:
		"""
		awk  'BEGIN {{OFS="\t"}}; {{ print $1, $2+$10, $2+$10+1, $4, $8, $9}}' {input.narrow} > {output.summit} 2> {log}
		awk 'BEGIN {{OFS="\t"}}; {{print $1, $2+$7, $2+$7+1, $4, $5, $6 }}' {input.bed_filt} > {output.sum_filt} 2>> {log}
		"""