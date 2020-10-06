# ----- Peak calling with MACS2 ----- #

# I have to determine the best parameters for ATAC-seq. 

# In some cases, in Paired-end samples, each read can be on each side of one nucleosome (long fragments: aprox >200bp o >150bp). 
#   This will center the peak in the nucleosome, which is not accurate for ATAC-seq, since Tn5-insertion sites are in accessible regions.

# Some people use --format BAM even when treating with Paired-end samples (along with --nomodel --extsize 200/150 and --shift -100/-75, 
#   depending on the assumption of how long is a nucleosome-free region with one nucleosome removed, even some people use --extsize 73 --shift -37).
#   Using --format BAM (--extsize x --shift -x/2) will took the 5' end of the first read (forward), 
#       it will extend a x bp-window with --extsize and then will center the peak with --shift. Using --format BAM does not take into account the reverse read.
#   Other people use --format BAMPE.  Using --format BAMPE will take into account the left mate (forward read) 
#       and the insertion site (+- lenght of the fragment between left and right mates)

# Also, don't forget to look at HMMRATAC (https://github.com/LiuLabUB/HMMRATAC) 
# and Genrich (https://github.com/jsh58/Genrich), which are peak callers specifically designed for ATAC-seq.

rule macs2_callpeak:
    input:
        case = "results/02_bam/{sample}/{sample}.bam",
    output:
        narrow = "results/03_macs2/{sample}/{sample}_peaks.narrowPeak",
        xls    = "results/03_macs2/{sample}/{sample}_peaks.xls"
    params:
        callpeak = config["params"]["macs2"]["callpeak"],
        gsize    = config["params"]["macs2"]["gsize"],
        pe       = lambda w: "--format BAM --nomodel --extsize 200 --shift -100" if is_single_end(w.sample) else "--format BAMPE",
        outdir   = "results/03_macs2/{sample}/",
        name     = "{sample}",
    threads: 5
    log:
        "results/00_log/macs2/{sample}_macs2.log"
    shell:
        """
        macs2 callpeak \
        --treatment {input.case} \
        --outdir {params.outdir} \
        {params.callpeak} \
        {params.pe} --extsize 200 --shift -100 \
        --gsize {params.gsize} \
        --name {params.name} \
        -p 1e-5 2> {log} 
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
        awk "\$8 >= {params.log_pval_filt}" {input.narrow} | cut -f1-4,8 > {output.bed_filt} 2> {log}
        """


