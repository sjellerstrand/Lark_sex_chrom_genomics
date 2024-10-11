#!/bin/bash -l

#SBATCH -A snic2022-5-71
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 15:00:00
#SBATCH -J call_variants
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Rasolark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
REF=$WORKDIR/A4_make_new_reference/$PROJECT\_consensus_reference.fasta;
METADATA=$WORKDIR/metadata;
NONDIPLOIDS=No;

### Load modules
module load bioinfo-tools freebayes/1.3.2 \
bcftools/1.14 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/A6_call_variants;
OUTDIR=$OUTDIR/A6_call_variants;

### Setup file info
fasta_generate_regions.py $REF.fai 100000 > $OUTDIR/vcf_regions.txt;
regions=$OUTDIR/vcf_regions.txt;
INDS=$(cat $METADATA/sample_info.txt | tail -n+2 | cut -f1);
BAMS=$(for IND in ${INDS[@]}; do \
echo $WORKDIR/A5_align_reads/$IND\_qsorted_merged_nodup_filtered_sorted.bam; done);

## Set maximum coverage allowed
MAX_COV=$(echo "$(echo $INDS | tr ' ' '\n' | wc -l) * 200" | bc);
echo "A maximum coverage of $MAX_COV is allowed per site for alignments to be processed";

## Set ploidies
if [ $NONDIPLOIDS == "Yes" ]; then
PLOIDY=$(echo --cnv-map $METADATA/ploidy.txt);
GENOMES=$(echo "$(cat $METADATA/ploidy.txt | \
awk -F'\t' 'BEGIN {sum=0;} {sum+=$2;} END {print sum}')");
echo "Some or all samples are non-diploid. Individual ploidies are provided by a file.";
else
GENOMES=$(echo "$(echo $INDS | tr ' ' '\n' | wc -l) * 2 " | bc);
echo "All samples assumed to be diploid.";
fi;
echo "Calling variants from $GENOMES genomes";

## Set maximum number of alleles
if [ $(echo "$(echo $GENOMES) / 4" | bc) -lt 10 ]; then
MAX_ALLELES=10;
else
MAX_ALLELES=$(echo "$(echo $GENOMES) / 4" | bc);
fi;
echo "A maximum number of $MAX_ALLELES alleles will be processed per site";

### Call variants in parallel with freebayes
freebayes-parallel $regions 20 -f $REF -b $BAMS --standard-filters \
-n $MAX_ALLELES -g  $MAX_COV --min-coverage 3 $PLOIDY \
> $OUTDIR/$PROJECT\_raw.vcf;
bgzip $OUTDIR/$PROJECT\_raw.vcf;
tabix -p vcf $OUTDIR/$PROJECT\_raw.vcf.gz;
