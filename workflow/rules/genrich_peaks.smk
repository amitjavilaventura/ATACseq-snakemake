#########################################
# ----- Peak calling with Genrich ----- #
#########################################

# https://github.com/jsh58/Genrich

### To do: 
###    allow to input a control --> make a function for input and a function for params
###    allow to input more than one replicate --> function in the input // different rule??

rule genrich:
	input:
		rules.align.output.bam,
	output:
		"results/03_genrich/{sample}/{sample}_peaks.narrowPeak",
	params:
		# path to genrich
		genrich = config["params"]["genrich"]["path"],
		# pe/se parameters
		pe      = lambda w: config["params"]["genrich"]["se"] if is_single_end(w.sample) else "",
		# pvalue/qvalue threshold
		p_or_q  = config["params"]["genrich"]["p_or_q"],
		pqval   = config["params"]["genrich"]["pqval"],
		# readme file
		readme  = "results/03_genrich/readme.txt",
	threads: 5
	log:
		"results/00_log/genrich/{sample}_peakcalling.log",
	shell:
		"""
		{params.genrich} -j \
		-t {input} \
		-o {output} \
		{params.p_or_q} {params.pqval} \
		{params.pe} 2>> {log}

		echo 'This peaks have been called using {params.pe}, and a pvalue threshold of {params.pvalue} and a qvalue threshold of {params.qval}' > {params.readme}
		"""


	
