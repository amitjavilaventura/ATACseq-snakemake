#PBS -l select=1:ncpus=1:mem=5G

# qsub has to be done from the main directory of the project.
# to do the qsub we will use a script called "call_scripts.sh" from where we will call the other scritps with their respective arguments.

# the script overlaps_reps.sh has to be run before this one. (if we have more than three reps we will have to do the second ant third, etc, intersects manually in the overlaps_reps.sh)

#bedtools need to be installed (can create an environement)
condition1=$1
condition2=$2

# ----- Create directories ----- #
mkdir -p 05_IDR/ 05_IDR/overlaps/ 05_IDR/overlaps/${condition1}-vs-${condition2}/

######################
# bedtools intersect #
######################

# overlaps between condition1 and condition2
bedtools intersect -a 05_IDR/${condition1}/${condition1}_logp_10.filt.coord.idr  \
-b 05_IDR/${condition2}/${condition2}_logp_10.filt.coord.idr  \
> 05_IDR/overlaps/${condition1}-vs-${condition2}/${condition1}-vs-${condition2}_overlap.bed

# no overlap. unique peaks from condition1
bedtools intersect -a 05_IDR/${condition1}/${condition1}_logp_10.filt.coord.idr  \
-b 05_IDR/${condition2}/${condition2}_logp_10.filt.coord.idr  -v \
> 05_IDR/overlaps/${condition1}-vs-${condition2}/${condition1}_unique.bed

# no overlap. unique peaks from condition2
bedtools intersect -a 05_IDR/${condition2}/${condition2}_logp_10.filt.coord.idr  \
-b 05_IDR/${condition1}/${condition1}_logp_10.filt.coord.idr  -v \
> 05_IDR/overlaps/${condition1}-vs-${condition2}/${condition2}_unique.bed
