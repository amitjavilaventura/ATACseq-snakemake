# ANNOTATE PEAKS
# ==============

# ----- Annotate peaks with peakAnno from ChIPseeker in R ----- #
rule peakAnnot:
    input:
        rules.filter_peaks.output.bed_filt,
    output:
        annot             = "results/04_peakAnno/{{sample}}/{{sample}}_peaks_logp{pvalue}.annot".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        promo_bed_targets = "results/04_peakAnno/{{sample}}/{{sample}}_peaks_logp{pvalue}_promoTargets.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        promoTargets      = "results/04_peakAnno/{{sample}}/{{sample}}_peaks_logp{pvalue}_promoTargets.txt".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        promoBed          = "results/04_peakAnno/{{sample}}/{{sample}}_peaks_logp{pvalue}_promoPeaks.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        distalBed         = "results/04_peakAnno/{{sample}}/{{sample}}_peaks_logp{pvalue}_distalPeaks.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
    params:
        before = config["promoter"]["bTSS"],
        after  = config["promoter"]["aTSS"],
        genome = lambda wildcards: SAMPLES.GENOME[wildcards.sample]
    log: 
        "results/00_log/peakAnnot/{sample}_peakanot.log"
    message:
        "Annotating peaks for {wildcards.sample}"
    shell:
        """
        Rscript --vanilla workflow/scripts/pipeline/peakAnno.R {input} {params.before} {params.after}   \
            {output.annot} {output.promo_bed_targets} {output.promoTargets} {output.promoBed} \
            {output.distalBed} {params.genome} 2> {log}
        """

