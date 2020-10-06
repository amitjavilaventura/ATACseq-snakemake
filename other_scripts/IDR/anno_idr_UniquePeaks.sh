#PBS -l select=1:ncpus=1:mem=5G
#cd $PBS_O_WORKDIR

# Usage: qsub of the call_scripts.sh
# qsub has to be done from the main directory of the project.
# to do the qsub we will use a script called "call_scripts.sh" from where we will call the other scritps with their respective arguments.

condition1=$1 #condition (ex. apc KO)
condition2=$2

mkdir -p  05_IDR/ \
05_IDR/annotation/ \
05_IDR/annotation/overlaps/ \
05_IDR/annotation/overlaps/${condition1}-vs-${condition2}/

# after doing the overlaps with bedtools intersect, this scripts calls the Rscript annoUniquePeaks.R that does the annotation with R annotatePeak function (package ChIPseeker).
# 	then, it filters the different outputs in order to have:
#		-peakAnnot.bed = The whole annotation file (with only one promoter region, from 0 to 3kb, and not separated in different distances from TSS)		
#		-promoTargets.bed = The genes with peaks in their promoter, and the coordinates of the genes. 
#		-promoTargets.txt = The names of the genes with peaks in their promoter
#		-promoPeaks = The coordinates of the peaks that fall in gene promoters, the peakID, the p-value, the gene name in ENSEMBL and SYMBOL formats.
#		-distalPeaks = The coordinates of the peaks that fall in other regions than promoters, the peakID, the p-value and the name of the nearest gene in ENSEMBL and SYMBOL formats. 
#		-distalTargets.bed = The nearest genes of distal peaks and the coordinates of the genes
#		-distalTargets.txt = The names of nearest genes from distal peaks

Rscript --vanilla other_scripts/annoUniquePeaks.R \
05_IDR/overlaps/${condition1}-vs-${condition2}/${condition1}_unique.bed \
05_IDR/annotation/overlaps/${condition1}-vs-${condition2}/${condition1}_unique_peakAnnot.bed \
05_IDR/annotation/overlaps/${condition1}-vs-${condition2}/${condition1}_unique_peakAnnot_promoTargets.bed \
05_IDR/annotation/overlaps/${condition1}-vs-${condition2}/${condition1}_unique_peakAnnot_promoTargets.txt \
05_IDR/annotation/overlaps/${condition1}-vs-${condition2}/${condition1}_unique_peakAnnot_promoPeaks.bed \
05_IDR/annotation/overlaps/${condition1}-vs-${condition2}/${condition1}_unique_peakAnnot_distalPeaks.bed \
05_IDR/annotation/overlaps/${condition1}-vs-${condition2}/${condition1}_unique_peakAnnot_distalTargets.bed \
05_IDR/annotation/overlaps/${condition1}-vs-${condition2}/${condition1}_unique_peakAnnot_distalTargets.txt
