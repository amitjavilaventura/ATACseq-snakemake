#PBS -l select=1:ncpus=1:mem=5G


##### -------------------------------------------------------------- #####
##### Script merge_all_peaks.sh have to be laucnhed before this one! #####
##### -------------------------------------------------------------- #####

# qsub has to be done from the main directory of the project.
# to do the qsub we will use a script called "call_scripts.sh" from where we will call the other scritps with their respective arguments.

# This script will be used to get a featureCounts output to analyse with DESeq2.
# It can only be used when we have more than one replicate per each condition, otherways, DESeq2 does not work. 


sample=$1

mkdir -p \
09_fCounts/ \
09_fCounts/${sample}/


# # First we run the "awk" command to pick the geneID column and the coordinates of the peak annotation bed files (promoter and distal targets, distal target is all things different to promoter):
# #	Column1 = peak id
# #	Column2 = chromosome number (column 1 of bed file)
# #	Column3 = start position of the peak (column 2 of the bed)
# #	Column4 = end position of the peak (column 3 of the bed)

awk 'OFS="\t" {print $1"."$2+1"."$3, $1, $2+1, $3, "."}' 03_macs2/all_filt_peaks.sort.bed > 09_fCounts/all_filt_peaks.sort.saf

# ----- Run featureCounts ----- #
featureCounts -a 9_fCounts/all_filt_peaks.sort.saf -F SAF -o 09_fCounts/${sample}/${sample}_allPeaks_countMatrix.featureCounts 02_bamfiles/${sample}/${sample}.bam
cut -f 1,7 09_fCounts/${sample}/${sample}_allPeaks_countMatrix.featureCounts | tail -n +2 > 09_fCounts/${sample}/${sample}_allPeaks.counts