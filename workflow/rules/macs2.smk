############################
# ----- Peak calling ----- #
############################

### Using MACS2 ###
### =========== ###

### TO DO: allow input file as control.

if config["options"]["peakcaller"] == "macs":

    # Call peaks
    rule macs2_callpeak:
        input:
            "results/02_bam/{sample}.bam",
        output:
            narrow = "results/03_macs2/{sample}/{sample}_peaks.narrowPeak",
            xls    = "results/03_macs2/{sample}/{sample}_peaks.xls"
        params:
            callpeak = config["params"]["macs2"]["callpeak"],
            pvalue   = config["params"]["macs2"]["pvalue"],
            gsize    = config["params"]["macs2"]["gsize"],
            pe       = lambda w: "--format BAM --nomodel --extsize 200 --shift -100" if is_single_end(w.sample) else "--format BAMPE",
            outdir   = "results/03_macs2/{sample}/",
        threads:
            CLUSTER["macs2_callpeak"]["cpu"]
        log:
            "results/00_log/macs2/{sample}_macs2.log"
        shell:
            """
            macs2 callpeak {params.pe} \
                --treatment {input} \
                --gsize {params.gsize} \
                --outdir {params.outdir} \
                --name {wildcards.sample} \
                --pvalue {params.pvalue} \
                {params.callpeak} 2> {log} 
            """


    # ----- Filter peaks from .narrowPeak files ----- #
    rule filter_peaks:
            input:
                narrow = "results/03_macs2/{sample}/{sample}_peaks.narrowPeak",
            output:
                bed_filt = "results/03_macs2/{{sample}}/{{sample}}_peaks_p{pvalue}.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"])
            params:
                log_pval_filt = config["params"]["macs2"]["filt_peaks_pval"],
            log:
                "results/00_log/macs2/{sample}_macs2_filt_peaks.log"
            shell:
                """
                awk "\$8 >= {params.log_pval_filt}" {input.narrow} | cut -f1-4,8,9 > {output.bed_filt} 2> {log}
                """
    
    # ----- Annotate peaks with peakAnno from ChIPseeker in R ----- #
    rule peakAnnot:
        input:
            rules.filter_peaks.output.bed_filt,
        output:
            annot             = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_logp{pvalue}.annot".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
            promo_bed_targets = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_logp{pvalue}_promoTargets.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
            promoTargets      = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_logp{pvalue}_promoTargets.txt".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
            promoBed          = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_logp{pvalue}_promoPeaks.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
            distalBed         = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_logp{pvalue}_distalPeaks.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
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
