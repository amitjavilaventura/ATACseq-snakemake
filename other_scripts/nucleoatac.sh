#PBS -l select=1:ncpus=10:mem=40G
cd $PBS_O_WORKDIR


##### IT DOES NOT WORK YET!!! ###

sample=$1

mkdir -p \
07_nucleoatac/ \
07_nucleoatac/genome/ \
07_nucleoatac/bedfiles/ \
07_nucleoatac/bedfiles/${sample} \
07_nucleoatac/bamfiles/ \
07_nucleoatac/bamfiles/${sample} \
07_nucleoatac/${sample}/ \
07_nucleoatac/${sample}/textfiles/

# To determine nucleosome positioning and occupancy with nucleoATAC
# NucleoATAC has to be installed --> https://nucleoatac.readthedocs.io/en/latest/
# Only if samples are paired-end

# ----- NECESSARY FILES ----- #

# The .bed file can be obtained by the macs2 callpeak function with the --broad flag set on (.broadPeak) it can also be obtained from a .narrowPeak file.

# The .bam sample is the one we generate after the alignment.
#   Bam file with aligned reads. These generally will be filtered for reads not mapping to mitochondria & reads with high mapping quality.
# The .fasta file is the indexed reference genome. It must be indexed with "samtools faidx".


# # ----- Filter chrM peaks from bam files ----- #

# # This is optional. .bam files without mitochondrial peaks can be used or we can use the original .bam files.

cp \
02_bamfiles/${sample}/${sample}.bam \
07_nucleoatac/bamfiles/${sample}/${sample}.bam

# samtools view -h 02_bamfiles/${sample}/${sample}.bam \
# | awk '{if($3 != "chrM" && $3 != "chrUn"){print $0}}' \
# | samtools view -Sb - > 07_nucleoatac/bamfiles/${sample}/${sample}_noChrM.bam


# ----- Index reference genome ----- #

# First make a link to the reference genome .fasta file and then generate the index.

cp \
/hpcnfs/techunits/bioinformatics/refdata/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa \
07_nucleoatac/genome/mm10_ref_genome.fa

samtools faidx 07_nucleoatac/genome/mm10_ref_genome.fa #it will create a file named ref_genome.fasta.fai


# ----- EXTEND PEAK REGIONS ----- #

#Create a copy of the filtered macs2 narrowPeak file (.bed format) to the nucleoatac directory.
cp 03_macs2/${sample}/${sample}_coord_filt_peaks.bed 07_nucleoatac/bedfiles/${sample}/${sample}_filt_peaks_ext.bed

# The peak regions should be extended a bit (ex. 100-200bp), specially if they come from .narrowPeak files,
#   in order to detect the nucleosomes in the case that the peak falls inside a NFR flanked by nucleosomes (the flanking nucleosomes would not be able to be detected)
bedtools slop -b 200 \
-i 07_nucleoatac/bedfiles/${sample}/${sample}_filt_peaks_ext.bed \
-g 07_nucleoatac/genome/mm10_ref_genome.fa

# The extended peak regions should not overlap, so we have to use bedtools merge in order to merge the overlaps
bedtools merge -i 07_nucleoatac/bedfiles/${sample}/${sample}_filt_peaks_ext.bed > 07_nucleoatac/bedfiles/${sample}/${sample}_filt_peaks_ext_merged.bed

##################
# RUN NUCLEOATAC #
##################

# Article for references https://genome.cshlp.org/content/25/11/1757

nucleoatac run \
--bed 07_nucleoatac/bedfiles/${sample}/${sample}_filt_peaks_ext_merged.bed \
--bam 02_bamfiles/${sample}/${sample}.bam \
--fasta 07_nucleoatac/genome/mm10_ref_genome.fa \
--out 07_nucleoatac/${sample}/${sample}_nucleoatac --cores 3

# we can also use this without nucleosomes: --bam 07_nucleoatac/bamfiles/${sample}/${sample}_noChrM.bam \

# ----- Remove files ----- #

# # To remove (from the working directory) the reference genome files (which are already in the cluster), the filtered .bam files (no chrM) and the extended and merged .bed files.
# rm -r 07_nucleoatac/genome 07_nucleoatac/bedfiles/ 07_nucleoatac/bamfiles/
