rule genrich:
    input:
        nfr  = lambda w: expand("results/02_bam/{sample}_nfr.bam.sorted", sample = SAMPLES.loc[SAMPLES["CONDITION"] == w.condition].NAME),
        nucl = lambda w: expand("results/02_bam/{sample}_nucleosomes.bam.sorted", sample = SAMPLES.loc[SAMPLES["CONDITION"] == w.condition].NAME),
    output:
        narrowPeak_nfr  = "results/03_genrich/{condition}/{condition}_peaks_nfr.narrowPeak",
        narrowPeak_nucl = "results/03_genrich/{condition}/{condition}_peaks_nucleosomes.narrowPeak",
        summit_nfr      = "results/03_genrich/{condition}/{condition}_summits_nfr.bed",
        summit_nucl     = "results/03_genrich/{condition}/{condition}_summits_nucleosomes.bed",
    params:
        p_or_q  = config["params"]["genrich"]["p_or_q"],
        pqval   = config["params"]["genrich"]["pqval"],
        readme  = "results/03_genrich/{condition}/readme.txt",
    log:
        "results/00_log/genrich/{condition}_peakcalling.log",
    shell:
        """
        #------------------------------ Call peaks for NFR and generate summit file ------------------------------#
        Genrich \
        -t '{input.nfr}' \
        -o {output.narrowPeak_nfr} \
        -{params.p_or_q} {params.pqval} \
        -v 2> {log}

        awk  'BEGIN {{OFS="\t"}}; {{ print $1, $2+$10, $2+$10+1, $4, $8, $9 }}'  {output.narrowPeak_nfr} > {output.summit_nfr}

        #------------------------------ Call peaks for nucleosomes and generate summit file ------------------------------#
        Genrich \
        -t '{input.nucl}' \
        -o {output.narrowPeak_nucl} \
        -{params.p_or_q} {params.pqval} \
        -v 2> {log}

        awk  'BEGIN {{OFS="\t"}}; {{ print $1, $2+$10, $2+$10+1, $4, $8, $9 }}'  {output.narrowPeak_nucl} > {output.summit_nucl}
        """


rule filter_peaks_genrich:
    input:
        narrow_nfr  = "results/03_genrich/{condition}/{condition}_peaks_nfr.narrowPeak",
        narrow_nucl = "results/03_genrich/{condition}/{condition}_peaks_nucleosomes.narrowPeak",
    output:
        bed_filt_nfr         = "results/03_genrich/{condition}/{condition}_peaks_nfr_{threshold}{pqvalue}.bed",
        sum_filt_nfr         = "results/03_genrich/{condition}/{condition}_summits_nfr_{threshold}{pqvalue}.bed",
        bed_filt_nucleosomes = "results/03_genrich/{condition}/{condition}_peaks_nucleosomes_{threshold}{pqvalue}.bed",
        sum_filt_nucleosomes = "results/03_genrich/{condition}/{condition}_summits_nucleosomes_{threshold}{pqvalue}.bed"
    params:
        log_pqval_filt = config["params"]["genrich"]["filt_peaks_pqval"],
        column_pqval   = lambda w: "8" if config["params"]["genrich"]["p_or_q"]== "p" else "9",
        readme         = "results/03_genrich/readme.txt",
    log:
        "results/00_log/genrich/{condition}_genrich_filt_peaks_{threshold}{pqvalue}.log"
    shell:
        """
        #------------------------------ NFR ------------------------------#
        awk "\${params.column_pqval} >= {params.log_pqval_filt}" {input.narrow_nfr} | cut -f1-4,8,9,10 > {output.bed_filt_nfr} 2> {log}
        awk 'BEGIN {{OFS="\t"}}; {{print $1, $2+$7, $2+$7+1, $4, $5, $6 }}' {output.bed_filt_nfr} > {output.sum_filt_nfr}

        #------------------------------ Nucleosomes ------------------------------#
        awk "\${params.column_pqval} >= {params.log_pqval_filt}" {input.narrow_nucl} | cut -f1-4,8,9,10 > {output.bed_filt_nucleosomes} 2> {log}
        awk 'BEGIN {{OFS="\t"}}; {{print $1, $2+$7, $2+$7+1, $4, $5, $6 }}' {output.bed_filt_nucleosomes} > {output.sum_filt_nucleosomes}
        """


rule peakAnnot_genrich:
    input:
        rules.filter_peaks_genrich.output.bed_filt_nfr,
    output:
        annot             = "results/04_peakAnno/{condition}/{condition}_peaks_nfr_{threshold}{pqvalue}.annot",
        promo_bed_targets = "results/04_peakAnno/{condition}/{condition}_peaks_nfr_{threshold}{pqvalue}_promoTargets.bed",
        promoTargets      = "results/04_peakAnno/{condition}/{condition}_peaks_nfr_{threshold}{pqvalue}_promoTargets.txt",
        promoBed          = "results/04_peakAnno/{condition}/{condition}_peaks_nfr_{threshold}{pqvalue}_promoPeaks.bed",
        distalBed         = "results/04_peakAnno/{condition}/{condition}_peaks_nfr_{threshold}{pqvalue}_distalPeaks.bed",
    params:
        before = config["promoter"]["bTSS"],
        after  = config["promoter"]["aTSS"],
        genome = config["ref"]["genome"]
    log: 
        "results/00_log/peakAnnot/{condition}_{threshold}{pqvalue}_peakanot.log"
    message:
        "Annotating peaks for {wildcards.condition}"
    shell:
        """
        Rscript --vanilla workflow/scripts/peakAnno.R {input} {params.before} {params.after}   \
            {output.annot} {output.promo_bed_targets} {output.promoTargets} {output.promoBed} \
            {output.distalBed} {params.genome} 2> {log}
        """
