#PBS -l select=1:ncpus=1:mem=5G

# qsub has to be done from the main directory of the project.
# to do the qsub we will use a script called "call_scripts.sh" from where we will call the other scritps with their respective arguments.

#deeptools needs to be installed (create a conda environment and do conda activate environmentname).

sample=$1

mkdir -p \
06_plotHeatmap/ \
06_plotHeatmap/${sample}/

# the regions file is all the peaks called by macs2 in all samples, filtered and merged in one file. 

computeMatrix reference-point \
--referencePoint TSS \
--beforeRegionStartLength 2500 \
--afterRegionStartLength 2500 \
--scoreFileName 02_bigwigfiles/${sample}/${sample}.bw \
--regionsFileName 03_macs2/all_filt_peaks.sort.bed \
--outFileName 06_plotHeatmap/${sample}/${sample}_matrix_reads-in-peaks.gz

plotHeatmap \
--matrixFile 06_plotHeatmap/${sample}/${sample}_matrix_reads-in-peaks.gz \
--outFileName 06_plotHeatmap/${sample}/${sample}_heatmap_TSS_reads-in-peaks.png \
--xAxisLabel "Relative distance to peak center (bp)" \
--heatmapWidth 4 \
--heatmapHeight 20 \
--plotTitle  "Reads in TSS peaks (${sample})" \
--samplesLabel "${sample}"