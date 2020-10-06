##################################
### HOMER - FIND MOTIFS GENOME ###
##################################

# this script runs findMotifsGenome.pl from HOMER --> http://homer.ucsd.edu/homer/ngs/peakMotifs.html

sample=$1 #ex: APC_KO_org_rep1

mkdir -p 10_motifs/ \
10_motifs/${sample}/

# this script must be used with the summits.bed output from MACS2. 
# first argument = file with summits
# second argument = genome
# third argument = output directory
# -size = fragment size for motif finding. integer (ex. 50 -> 25 before and 25 after the summit) or 2 integers (ex. -50,50 -> 50 before and 50 after the summit)

findMotifsGenome.pl \
03_macs2/${sample}/${sample}_summits.bed \
/hpcnfs/scratch/DP/amitjavila/other/genomes/mm10/mm10_genome.fa \
10_motifs/${sample}/ \
-size 50

