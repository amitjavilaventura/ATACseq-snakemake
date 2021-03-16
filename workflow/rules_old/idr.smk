# IDR peak merging
# ================

if config["options"]["idr_merge"] == True:
    rule idr:
        input:
            lambda w: expand("results/03_macs2/{sample}/{sample}_peaks.narrowPeak", sample = SAMPLES.loc[SAMPLES["CONDITION"] == w.condition].NAME)
        output:
            "results/03_macs2/IDR/{condition}.txt"
        threads:
            CLUSTER["idr"]["cpu"]
        params:
        log:
            "results/00_log/idr/{condition}.log",
        shell:
            """
            echo "Using the samples {input} from condition {wildcards.condition} as input for IDR" > {output}
            """