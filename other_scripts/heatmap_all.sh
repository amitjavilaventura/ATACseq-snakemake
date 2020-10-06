#PBS -l select=1:ncpus=1:mem=5G

mkdir -p \
06_plotHeatmap/ \
06_plotHeatmap/all/

source activate deeptools

# the regions file have the coordinates of the genes that have a peak in the promoter

computeMatrix reference-point \
--referencePoint TSS \
--beforeRegionStartLength 2500 \
--afterRegionStartLength 2500 \
--scoreFileName 02_bigwigfiles/APC_KO_org_rep1/APC_KO_org_rep1.bw 02_bigwigfiles/APC_KO_org_rep2/APC_KO_org_rep2.bw \
02_bigwigfiles/bCAT_noEx3_org_rep1b/bCAT_noEx3_org_rep1b.bw 02_bigwigfiles/bCAT_noEx3_org_rep1c/bCAT_noEx3_org_rep1c.bw \
02_bigwigfiles/WT_WNT_org_rep1/WT_WNT_org_rep1.bw 02_bigwigfiles/WT_WNT_org_rep2/WT_WNT_org_rep2.bw \
--regionsFileName /hpcnfs/scratch/DP/amitjavila/other/anno_files/mm10.annot.sorted.bed \
--outFileName 06_plotHeatmap/all/all_matrix_reads-in-peaks_TSS.gz

computeMatrix reference-point \
--referencePoint TSS \
--beforeRegionStartLength 2500 \
--afterRegionStartLength 2500 \
--scoreFileName 02_bigwigfiles/APC_KO_org_rep1/APC_KO_org_rep1.bw 02_bigwigfiles/APC_KO_org_rep2/APC_KO_org_rep2.bw \
02_bigwigfiles/bCAT_noEx3_org_rep1b/bCAT_noEx3_org_rep1b.bw 02_bigwigfiles/bCAT_noEx3_org_rep1c/bCAT_noEx3_org_rep1c.bw \
02_bigwigfiles/WT_WNT_org_rep1/WT_WNT_org_rep1.bw 02_bigwigfiles/WT_WNT_org_rep2/WT_WNT_org_rep2.bw \
--regionsFileName 03_macs2/all_filt_peaks.sort.bed \
--outFileName 06_plotHeatmap/all/all_matrix_reads-in-peaks.gz

plotHeatmap \
--matrixFile 06_plotHeatmap/all/all_matrix_reads-in-peaks_TSS.gz \
--outFileName 06_plotHeatmap/all/all_heatmap_TSS_reads-in-peaks.png \
--xAxisLabel "Relative distance to TSS (bp)" \
--heatmapWidth 4 \
--heatmapHeight 20 \
--plotTitle  "Reads in TSS peaks (all)" \
--samplesLabel "APC KO rep1" "APC KO rep1b" "bCAT noEx3 1b" "bCAT noEx3 1c" "WT with WNT 1" "WT with WNT 1b" \
--colorList "white,blue"

plotHeatmap \
--matrixFile 06_plotHeatmap/all/all_matrix_reads-in-peaks.gz \
--outFileName 06_plotHeatmap/all/all_heatmap_TSS_reads-in-peaks.png \
--xAxisLabel "Relative distance to peak center (bp)" \
--heatmapWidth 4 \
--heatmapHeight 20 \
--plotTitle  "Reads in peaks (all)" \
--samplesLabel "APC KO rep1" "APC KO rep1b" "bCAT noEx3 1b" "bCAT noEx3 1c" "WT with WNT 1" "WT with WNT 1b" \
--colorList "white,blue"