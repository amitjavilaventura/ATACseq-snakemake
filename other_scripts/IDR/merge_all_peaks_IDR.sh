#PBS -l select=1:ncpus=1:mem=5G

# Need bedtools installed --> in my case conda environment: deeptools 

mkdir -p 05_IDR/ \
05_IDR/all_peaks/

# First pick all the .bed files obtained from the filtering of macs2 .narrowPeak files and put all them in one unique file
cat 05_IDR/*/*_logp_10.filt.coord.idr > 05_IDR/all_peaks/all_filt_peaks.bed

# Sort by chromosome and peakStart.
sort -k1,1 -k2,2n 05_IDR/all_peaks/all_filt_peaks.bed > 05_IDR/all_peaks/all_filt_peaks.sort.bed

# Run bedtools merge. 
bedtools merge -i 05_IDR/all_peaks/all_filt_peaks.sort.bed > 05_IDR/all_peaks/all_filt_peaks.sort.merge.bed
