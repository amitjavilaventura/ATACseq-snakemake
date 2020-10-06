#PBS -l select=1:ncpus=1:mem=5G
cd $PBS_O_WORKDIR

# qsub has to be done from the main directory of the project.
# to do the qsub we will use a script called "call_scripts.sh" from where we will call the other scritps with their respective arguments.

#bedtools need to be installed (can create an environement)
sample1=$1
sample2=$2

# ----- Create directories ----- #
mkdir -p 05_overlaps/ 05_overlaps/bedtools/ 05_overlaps/bedtools/${sample1}_vs_${sample2}/

######################
# bedtools intersect #
######################

# overlaps between sample1 and sample2
bedtools intersect -a 03_macs2/${sample1}/${sample1}_coord_filt_peaks.bed \
-b 03_macs2/${sample2}/${sample2}_coord_filt_peaks.bed \
> 05_overlaps/bedtools/${sample1}_vs_${sample2}/${sample1}_vs_${sample2}_overlap.bed

# no overlap. unique peaks from sample1
bedtools intersect -a 03_macs2/${sample1}/${sample1}_coord_filt_peaks.bed \
-b 03_macs2/${sample2}/${sample2}_coord_filt_peaks.bed -v \
> 05_overlaps/bedtools/${sample1}_vs_${sample2}/${sample1}_vs_${sample2}_no_overlap.bed

# no overlap. unique peaks from sample2
bedtools intersect -a 03_macs2/${sample2}/${sample2}_coord_filt_peaks.bed \
-b 03_macs2/${sample1}/${sample1}_coord_filt_peaks.bed -v \
> 05_overlaps/bedtools/${sample1}_vs_${sample2}/${sample2}_vs_${sample1}_no_overlap.bed
