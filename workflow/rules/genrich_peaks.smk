#########################################
# ----- Peak calling with Genrich ----- #
#########################################

### To do: 
###    allow to input a control --> make a function for input and a function for params
###    allow to input more than one replicate --> function in the input // different rule??

rule genrich:
	input:
		rules.align.output.bam,
	output:
		"results/03_genrich/{sample}/{sample}_peaks.narrowPeak",
	params: 
		genrich = config["params"]["genrich"]["path"],
		pe      = lambda w: config["params"]["genrich"]["se"] if is_single_end(w.sample) else "",
		qval    = config["params"]["genrich"]["qval"],
	threads: 5
	log:
		"results/00_log/genrich/{sample}_peakcalling.log",
	shell:
		"""
		{params.genrich} -j \
		-t {input} \
		-o {output} \
		-q {params.qval} \
		{params.pe} 2>> {log}
		"""

	
