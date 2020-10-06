#PBS -l select=1:ncpus=1:mem=5G

rep1=$1 
rep2=$2
condition=$3

# Usage: call from (qsub) call_scripts.sh --> bash other_scripts/idr.sh rep1 rep2 condition


# ¡¡ This script needs to have idr installed !! --> in the conda environment 'idr', idr is installed

mkdir -p \
05_IDR/ \
05_IDR/${condition}/


# ----- IRREPRODUCIBLE DISCOVERY RATE ----- #
# https://github.com/nboley/idr#irreproducible-discovery-rate-idr
# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html
# The output file format has the following columns.
#	#01: chromosome:	the name of the chromosome where the "common" feature is.
#	#02: chromStart: 	the starting position of the "common" feature in the chromosome.
#	#03: chromEnd: 		the end position of the "common" feature in the chromosome.
#	#04: name:			the name of the region. Can be '.' if no name is assigned.
#	#05: score:			the scaled idr value. min(int(log2(-125IDR), 1000) --> peaks with a idr of 0 = score of 1000, idr of 0.05 = score of 540; idr of 1 = score of 0.
#	#06: strandValue:	the strand (+ or -). '.' if no strand is assigned.
#	#07: signal:		the measurement of enrichment of the region. If peak list is provided, value assigned from the peak list.
#	#08: p-value:		the p-value of merged peaks.  If peak list is provided, value assigned from the peak list.
#	#09: q-value:		the q-value of merged peaks.  If peak list is provided, value assigned from the peak list.
#	#10: summit:		the summit of merged peaks. 
#	#11: local IDR: 	the -log10(local_IDR_value).
#	#12: global IDR:	the -log10(global_IDR_value).
#	#13: rep1_chrStart:	the starting position of feature in replicate 1.
#	#14: rep1_chrEnd:	the ending position of the feature in replicate 1. 
#	#15: rep1_signal:	the signal measure from replicate 1. It is determined by the "--rank" option. If "--rank p.value" --> column 15 = column 8 from narrowPeak,
#						if "--rank signal.value" --> column 15 in the output = column 8 in the narrowPeak from replicate 1. 
#	#16: rep1_summit:	the summit of the peak in replicate 1.
#	#17: rep2_chrStart:	the same as #13 but for replicate 2. 
#	#18: rep2_chrEnd:	the same as #14 but for replicate 2. 
#	#19: rep2_signal:	the same as #15 but for replicate 2. 
#	#20: rep2_summit:	the same as #16 but for replicate 2. 
 

# sort  narrowPeak files
sort -k8,8nr 03_macs2/${rep1}/${rep1}_peaks.narrowPeak > 03_macs2/${rep1}/${rep1}_peaks.sort.narrowPeak
sort -k8,8nr 03_macs2/${rep2}/${rep2}_peaks.narrowPeak > 03_macs2/${rep2}/${rep2}_peaks.sort.narrowPeak

# filter narrowPeak files
awk "\$8 >= 10" 03_macs2/${rep1}/${rep1}_peaks.sort.narrowPeak | cut -f1-10 > 03_macs2/${rep1}/${rep1}_peaks.filt.sort.narrowPeak
awk "\$8 >= 10" 03_macs2/${rep2}/${rep2}_peaks.sort.narrowPeak | cut -f1-10 > 03_macs2/${rep2}/${rep2}_peaks.filt.sort.narrowPeak

#################
#### run IDR ####
#################
idr --samples 03_macs2/${rep1}/${rep1}_peaks.filt.sort.narrowPeak 03_macs2/${rep2}/${rep2}_peaks.filt.sort.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--soft-idr-threshold 0.1 \
-o 05_IDR/${condition}/${condition}.idr \
--log-output-file 05_IDR/${condition}/log_${condition}.idr.log \
--plot

# filter merged-peak files by log10(pval) > 10 and score > 415 (= IDRvalue > 0.10)
awk '{if($8>=10){print $0}}' 05_IDR/${condition}/${condition}.idr |
awk '{if($5>=415){print $0}}' > 05_IDR/${condition}/${condition}_logp10_score415.filt.idr

# take columns 1 (chrom), 2 (peakStart), 3 (peakEnd), 5 (score), 8 (p-value) and 10 (summit).
cut -f 1-3,5,8,10 05_IDR/${condition}/${condition}_logp10_score415.filt.idr > 05_IDR/${condition}/${condition}_logp_10_score415.filt.coord.idr 
