#########################################
# ----- Peak calling with Genrich ----- #
#########################################

### Using GENRICH ###
### ============= ###

# https://github.com/jsh58/Genrich

### To do: 
###    allow to input a control --> make a function for input and a function for params
###    allow to input more than one replicate --> function in the input // different rule??

if config["options"]["genrich_pool"] == False:
	
	rule genrich:
		input:
			"results/02_bam/{sample}.bam.tmp",
		output:
			"results/03_genrich/{sample}/{sample}_peaks.narrowPeak",
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
			readme   = "results/03_genrich/readme.txt",
		threads: 1
		log:
			"results/00_log/genrich/{sample}_peakcalling.log",
		shell:
			"""
			Genrich {params.atacmode} \
			-t {input} \
			-o {output} \
			-{params.p_or_q} {params.pqval} \
			{params.unmatch} {params.chrM} \
			{params.PRCdups} -v 2>> {log}

			echo 'This peaks have been called with genrich using {{params.unmatch}} \
			and a {params.p_or_q}value threshold of {params.pqval}' > {params.readme}

			#rm {input}  
			"""

	rule summits_genrich:
		input:
			narrow = "results/03_genrich/{sample}/{sample}_peaks.narrowPeak",
		output:
			summit = "results/03_genrich/{sample}/{sample}_summits.bed",
		params:
			log_pqval_filt = config["params"]["genrich"]["filt"],
			column_pqval   = lambda w: "8" if config["params"]["genrich"]["p_or_q"]=="p" else "9",
			readme         = "results/03_genrich/readme.txt",
			p_or_q         = config["params"]["genrich"]["p_or_q"],
			pq_filt        = config["params"]["genrich"]["filt"]
		log:
			"results/00_log/genrich/{condition}_get_summits.log"
		shell:
			"""
			awk  'BEGIN {{OFS="\t"}}; {{ print $1, $2+$10, $2+$10+1, $4, $8, $9}}' {input.narrow} > {output.summit} 2>> {log}
			"""


	pq    = config["params"]["genrich"]["p_or_q"]
	pqval = config["params"]["genrich"]["filt"]

	# filter peaks by p or q value
	rule filter_peaks_genrich:
		input:
			narrow = "results/03_genrich/{sample}/{sample}_peaks.narrowPeak",
		output:
			bed_filt = "results/03_genrich/{{sample}}/{{sample}}_peaks_{p_or_q}{pqvalue}.bed".format(pqvalue = config["params"]["genrich"]["filt"], p_or_q = config["params"]["genrich"]["p_or_q"]),
			sum_filt = "results/03_genrich/{{sample}}/{{sample}}_peaks_{p_or_q}{pqvalue}.bed".format(pqvalue = config["params"]["genrich"]["filt"], p_or_q = config["params"]["genrich"]["p_or_q"]),
		params:
			log_pqval_filt = config["params"]["genrich"]["filt"],
			column_pqval   = lambda w: "8" if config["params"]["genrich"]["p_or_q"]=="p" else "9",
			readme         = "results/03_genrich/readme.txt",
			p_or_q         = config["params"]["genrich"]["p_or_q"],
			pq_filt        = config["params"]["genrich"]["filt"]
		log:
			"results/00_log/genrich/{sample}_genrich_filt_peaks.log"
		shell:
			"""
			awk "\${params.column_pqval} >= {params.log_pqval_filt}" {input.narrow} | cut -f1-4,8,9,10 > {output.bed_filt} 2> {log}

			awk 'BEGIN {{OFS="\t"}}; {{print$1, $2+$7, $2+$7+1, $5, $6 }}' {output.bed_filt} > {output.sum_filt}

			echo '
			Genrich narrowpeaks were filtered using a -log10({params.p_or_q}value) threshold of {params.pq_filt}.
			The' \
			> {params.readme}
			"""


		# annotate peaks as promoter and distal
	rule peakAnnot_genrich:
		input:
			rules.filter_peaks_genrich.output.bed_filt,
		output:
			annot             = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_log{p_or_q}{pqvalue}.annot".format(pqvalue = pqval, p_or_q = pq),
			promo_bed_targets = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_log{p_or_q}{pqvalue}_promoTargets.bed".format(pqvalue = pqval, p_or_q = pq),
			promoTargets      = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_log{p_or_q}{pqvalue}_promoTargets.txt".format(pqvalue = pqval, p_or_q = pq),
			promoBed          = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_log{p_or_q}{pqvalue}_promoPeaks.bed".format(pqvalue = pqval, p_or_q = pq),
			distalBed         = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_log{p_or_q}{pqvalue}_distalPeaks.bed".format(pqvalue = pqval, p_or_q = pq),
		params:
			before = config["promoter"]["bTSS"],
			after  = config["promoter"]["aTSS"],
			genome = lambda wildcards: SAMPLES.GENOME[wildcards.sample],
		log: 
			"results/00_log/peakAnnot/{sample}_peakanot_genrich.log",
		message:
			"Annotating peaks for {wildcards.sample}",
		shell:
			"""
			Rscript --vanilla workflow/scripts/pipeline/peakAnno.R {input} {params.before} {params.after}   \
				{output.annot} {output.promo_bed_targets} {output.promoTargets} {output.promoBed} \
				{output.distalBed} {params.genome} 2> {log}
			"""
	


## Merging more than one replicate in the same peak calling step
elif config["options"]["genrich_pool"] == True:

	rule genrich_pool:
		input:
			t = lambda w: expand("results/02_bam/{condition}.bam.tmp", condition = SAMPLES.loc[SAMPLES["CONDITION"] == w.condition].NAME),
		output:
			"results/03_genrich/pooled/{condition}/{condition}_peaks.narrowPeak",
		params:
			# path to genrich
			atacmode = config["params"]["genrich"]["atac_mode"],
			# use unmatched reads
			unmatch   = config["params"]["genrich"]["unmatch"],
			# pvalue/qvalue threshold
			p_or_q  = config["params"]["genrich"]["p_or_q"],
			pqval   = config["params"]["genrich"]["pqval"],
			# rm chr M reads
			chrM     = config["params"]["genrich"]["chrM"],
			# rm pcr duplicates
			PRCdups = config["params"]["genrich"]["PRCdups"],
			# readme file
			readme  = "results/03_genrich/pooled/{condition}/readme.txt",
		threads: 5
		log:
			"results/00_log/genrich_pool/{condition}_peakcalling.log",
		shell:
			"""
			Genrich {params.atacmode} \
			-t '{input.t}' \
			-o {output} \
			-{params.p_or_q} {params.pqval} \
			{params.unmatch} {params.chrM} \
			{params.PRCdups} -v 2>> {log}

			echo 'This peaks have been called using more than 1 replicates ({input}) for each condition, the parameters {params.unmatch} and a {params.p_or_q}value threshold of {params.pqval}.' > {params.readme}
			"""

	rule summits_genrich_pool:
		input:
			narrow = "results/03_genrich/pooled/{condition}/{condition}_peaks.narrowPeak",
		output:
			summit = "results/03_genrich/pooled/{condition}/{condition}_summits.bed",
		log: 
			"results/00_log/genrich_pool/{condition}_get_summits.log",
		shell:
			"""
			awk  'BEGIN {{OFS="\t"}}; {{ print $1, $2+$10, $2+$10+1, $4, $8, $9 }}'  {input.narrow} > {output.summit} 2>> {log}
			"""

	rule filter_peaks_genrich_pool:
		input:
			narrow = "results/03_genrich/pooled/{condition}/{condition}_peaks.narrowPeak",
		output:
			bed_filt = "results/03_genrich/pooled/{{condition}}/{{condition}}_peaks_p{pvalue}.bed".format(pvalue = config["params"]["genrich"]["filt"]),
			sum_filt = "results/03_genrich/pooled/{{condition}}/{{condition}}_summits_p{pvalue}.bed".format(pvalue = config["params"]["genrich"]["filt"]),
		params:
			log_pqval_filt = config["params"]["genrich"]["filt"],
			column_pqval   = lambda w: "8" if config["params"]["genrich"]["p_or_q"]=="p" else "9",
			readme         = "results/03_genrich/readme.txt",
			p_or_q         = config["params"]["genrich"]["p_or_q"],
			pq_filt        = config["params"]["genrich"]["filt"]
		log:
			"results/00_log/genrich_pool/{condition}_genrich_filt_peaks.log"
		shell:
			"""
			awk "\${params.column_pqval} >= {params.log_pqval_filt}" {input.narrow} | cut -f1-4,8,9,10 > {output.bed_filt} 2> {log}
			awk 'BEGIN {{OFS="\t"}}; {{print $1, $2+$7, $2+$7+1, $4, $5, $6 }}' {output.bed_filt} > {output.sum_filt}
			"""


	# ----- Annotate peaks with peakAnno from ChIPseeker in R ----- #
	pq    = config["params"]["genrich"]["p_or_q"]
	pqval = config["params"]["genrich"]["filt"]

	rule peakAnnot_genrich_pool:
		input:
			rules.filter_peaks_genrich_pool.output.bed_filt,
		output:
			annot             = "results/04_peakAnno/pooled/{{condition}}/{{condition}}-peaks_log{p_or_q}{pqvalue}.annot".format(pqvalue = pqval, p_or_q = pq),
			promo_bed_targets = "results/04_peakAnno/pooled/{{condition}}/{{condition}}-peaks_log{p_or_q}{pqvalue}_promoTargets.bed".format(pqvalue = pqval, p_or_q = pq),
			promoTargets      = "results/04_peakAnno/pooled/{{condition}}/{{condition}}-peaks_log{p_or_q}{pqvalue}_promoTargets.txt".format(pqvalue = pqval, p_or_q = pq),
			promoBed          = "results/04_peakAnno/pooled/{{condition}}/{{condition}}-peaks_log{p_or_q}{pqvalue}_promoPeaks.bed".format(pqvalue = pqval, p_or_q = pq),
			distalBed         = "results/04_peakAnno/pooled/{{condition}}/{{condition}}-peaks_log{p_or_q}{pqvalue}_distalPeaks.bed".format(pqvalue = pqval, p_or_q = pq),
		params:
			before = config["promoter"]["bTSS"],
			after  = config["promoter"]["aTSS"],
			genome = "mm10"
		log: 
			"results/00_log/peakAnnot/{condition}_peakanot.log"
		message:
			"Annotating peaks for {wildcards.condition}"
		shell:
			"""
			Rscript --vanilla workflow/scripts/pipeline/peakAnno.R {input} {params.before} {params.after}   \
				{output.annot} {output.promo_bed_targets} {output.promoTargets} {output.promoBed} \
				{output.distalBed} {params.genome} 2> {log}
			"""
