###################################
# Merge BAM files #### Bam2Bigwig #
###################################

condition=$1
replicate1=$2
replicate2=$3

mkdir -p 05_IDR/ \
05_IDR/bamfiles/ \
05_IDR/bigwig/

samtools merge 05_IDR/bamfiles/${condition}_merged.bam 02_bamfiles/${replicate1}/${replicate1}.bam 02_bamfiles/${replicate2}/${replicate2}.bam

samtools sort -m 5G -@ 5 -o 05_IDR/bamfiles/${condition}_merged.sort.bam 05_IDR/bamfiles/${condition}_merged.bam

samtools index 05_IDR/bamfiles/${condition}_merged.sort.bam

bamCoverage -b 05_IDR/bamfiles/${condition}_merged.sort.bam -o 05_IDR/bigwig/${condition}_merged.bw --normalizeUsing CPM

rm -r 05_IDR/bamfiles/
