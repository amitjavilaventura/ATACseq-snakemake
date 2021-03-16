############################
# ----- Peak calling ----- #
############################

### Using MACS2 ###
### =========== ###

### TO DO: allow input file as control.

if config["options"]["peakcaller"] == "macs2":

    # Call peaks
    rule macs2_callpeak:
        input:
            "results/02_bam/{sample}.bam",
        output:
            narrow = "results/03_macs2/{sample}/{sample}_peaks.narrowPeak",
            xls    = "results/03_macs2/{sample}/{sample}_peaks.xls"
        params:
            callpeak = config["params"]["macs2"]["callpeak"],
            p_or_q   = config["params"]["macs2"]["p_or_q"],
            pqvalue  = config["params"]["macs2"]["pqval"],
            gsize    = config["params"]["macs2"]["gsize"],
            pe       = lambda w: "--format BAM --nomodel --extsize 200 --shift -100" if is_single_end(w.sample) else "--format BAMPE",
            outdir   = "results/03_macs2/{sample}/",
            readme   = "results/03_macs2/readme.txt",
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
                --{params.p_or_q}value {params.pqvalue} \
                {params.callpeak} 2> {log}
            

			echo 'This peaks have been called with genrich using {params.pe} \
			and a {params.p_or_q}value threshold of {params.pqvalue}' > {params.readme}
            """


    # ----- Filter peaks from .narrowPeak files ----- #
    rule filter_peaks_macs:
            input:
                narrow = "results/03_macs2/{sample}/{sample}_peaks.narrowPeak",
            output:
                bed_filt = "results/03_macs2/{{sample}}/{{sample}}_peaks_{p_or_q}{pqvalue}.bed".format(pqvalue = config["params"]["macs2"]["filt_peaks_pqval"], p_or_q = config["params"]["macs2"]["p_or_q"])
            params:
                log_pqval_filt = config["params"]["macs2"]["filt_peaks_pqval"],
                column_pqval   = lambda w: "8" if config["params"]["macs2"]["p_or_q"]=="p" else "9",
                readme         = "results/03_macs2/readme.txt",
                p_or_q         = config["params"]["macs2"]["p_or_q"],
                pq_filt        = config["params"]["macs2"]["filt_peaks_pqval"]
            log:
                "results/00_log/macs2/{sample}_macs2_filt_peaks.log"
            shell:
                """
                awk "\${params.column_pqval} >= {params.log_pqval_filt}" {input.narrow} | cut -f1-4,8,9 > {output.bed_filt} 2> {log}
                
                echo '
				Macs2 narrowpeaks were filtered using a -log10({params.p_or_q}value) threshold of {params.pq_filt}' \
				>> {params.readme}

                """
    
    pq   = config["params"]["macs2"]["p_or_q"]
    pqval = config["params"]["macs2"]["filt_peaks_pqval"]

    # ----- Annotate peaks with peakAnno from ChIPseeker in R ----- #
    rule peakAnnot_macs:
        input:
            rules.filter_peaks_macs.output.bed_filt,
        output:
            annot             = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_log{p_or_q}{pqvalue}.annot".format(pqvalue = pqval, p_or_q = pq),
            promo_bed_targets = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_log{p_or_q}{pqvalue}_promoTargets.bed".format(pqvalue = pqval, p_or_q = pq),
            promoTargets      = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_log{p_or_q}{pqvalue}_promoTargets.txt".format(pqvalue = pqval, p_or_q = pq),
            promoBed          = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_log{p_or_q}{pqvalue}_promoPeaks.bed".format(pqvalue = pqval, p_or_q = pq),
            distalBed         = "results/04_peakAnno/{{sample}}/{{sample}}-peaks_log{p_or_q}{pqvalue}_distalPeaks.bed".format(pqvalue = pqval, p_or_q = pq),
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
