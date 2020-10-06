#PBS -l select=1:ncpus=1:mem=5G

condition1=$1
condition2=$2
condition3=$3

mkdir -p \
05_IDR/ \
05_IDR/heatmaps/ \
05_IDR/heatmaps/all/ \
05_IDR/heatmaps/all/matrix/

source activate deeptools

# the regions file have the coordinates of the all mm10 genes. 

computeMatrix reference-point \
--referencePoint TSS \
--beforeRegionStartLength 2500 \
--afterRegionStartLength 2500 \
--scoreFileName  05_IDR/bigwig/${condition1}_merged.bw 05_IDR/bigwig/${condition2}_merged.bw 05_IDR/bigwig/${condition3}_merged.bw \
--regionsFileName /hpcnfs/scratch/DP/amitjavila/other/anno_files/mm10.annot.sorted.bed \
--outFileName 05_IDR/heatmaps/all/matrix/all_matrix_reads-in-peaks_TSS.gz

# the regions file is all the peaks called by macs2 in all samples, replicate-merged with IDR, filtered and merged in one file. 

computeMatrix reference-point \
--referencePoint Summit \
--beforeRegionStartLength 2500 \
--afterRegionStartLength 2500 \
--scoreFileName  05_IDR/bigwig/${condition1}_merged.bw 05_IDR/bigwig/${condition2}_merged.bw 05_IDR/bigwig/${condition3}_merged.bw \
--regionsFileName 05_IDR/all_peaks/all_filt_peaks.sort.merge.bed\
--outFileName 05_IDR/heatmaps/all/matrix/all_matrix_reads-in-peaks.gz

plotHeatmap \
--matrixFile 05_IDR/heatmaps/all/matrix/all_matrix_reads-in-peaks_TSS.gz \
--outFileName 05_IDR/heatmaps/all/all_heatmap_TSS_reads-in-peaks.png \
--xAxisLabel "Relative distance to TSS (bp)" \
--heatmapWidth 4 \
--heatmapHeight 20 \
--plotTitle  "Reads in TSS peaks (all)" \
--samplesLabel "${condition1}" "${condition2}" "${condition3}" \
--colorList "white,blue"

plotHeatmap \
--matrixFile 05_IDR/heatmaps/all/matrix/all_matrix_reads-in-peaks.gz \
--outFileName 05_IDR/heatmaps/all/all_heatmap_reads-in-peaks.png \
--xAxisLabel "Relative distance to peak center (bp)" \
--heatmapWidth 4 \
--heatmapHeight 20 \
--plotTitle  "Reads in peaks (all)" \
--samplesLabel "${condition1}" "${condition2}" "${condition3}" \
--colorList "white,blue"